#include <arbor/domain_decomposition.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/recipe.hpp>
#include <arbor/symmetric_recipe.hpp>
#include <arbor/context.hpp>

#include "cell_group_factory.hpp"
#include "execution_context.hpp"
#include "gpu_context.hpp"
#include "util/maputil.hpp"
#include "util/partition.hpp"
#include "util/span.hpp"

namespace arb {

domain_decomposition partition_load_balance(
    const recipe& rec,
    const context& ctx,
    partition_hint_map hint_map)
{
    const bool gpu_avail = ctx->gpu->has_gpu();

    struct partition_gid_domain {
        partition_gid_domain(gathered_vector<cell_gid_type> divs, unsigned domains):
            gid_divisions(std::move(divs)), num_domains(domains)
        {}

        int operator()(cell_gid_type gid) const {
            auto gid_part = gid_divisions.partition();
            for (unsigned i = 0; i < num_domains; i++) {
                auto rank_gids = util::subrange_view(gid_divisions.values(), gid_part[i], gid_part[i+1]);
                if(std::binary_search(rank_gids.begin(), rank_gids.end(), gid)) {
                    return i;
                }
            }
            return -1;
        }

        const gathered_vector<cell_gid_type> gid_divisions;
        unsigned num_domains;
    };

    struct cell_identifier {
        cell_gid_type id;
        bool is_super_cell;
    };

    using util::make_span;

    unsigned num_domains = ctx->distributed->size();
    unsigned domain_id = ctx->distributed->id();
    auto num_global_cells = rec.num_cells();

    auto dom_size = [&](unsigned dom) -> cell_gid_type {
        const cell_gid_type B = num_global_cells/num_domains;
        const cell_gid_type R = num_global_cells - num_domains*B;
        return B + (dom<R);
    };

    // Global load balance

    std::vector<cell_gid_type> gid_divisions;
    auto gid_part = make_partition(
        gid_divisions, transform_view(make_span(num_domains), dom_size));

    // Local load balance

    std::vector<std::vector<cell_gid_type>> super_cells; //cells connected by gj
    std::vector<cell_gid_type> reg_cells; //independent cells

    // Map to track visited cells (cells that already belong to a group)
    std::unordered_map<cell_gid_type, bool> visited;

    // Connected components algorithm using BFS
    std::queue<cell_gid_type> q;
    for (auto gid: make_span(gid_part[domain_id])) {
        if (!rec.gap_junctions_on(gid).empty()) {
            // If cell hasn't been visited yet, must belong to new super_cell
            // Perform BFS starting from that cell
            if (visited.find(gid) == visited.end()) {
                std::vector<cell_gid_type> cg;
                q.push(gid);
                visited[gid] = true;
                while (!q.empty()) {
                    auto element = q.front();
                    q.pop();
                    cg.push_back(element);
                    // Adjacency list
                    auto conns = rec.gap_junctions_on(element);
                    for (auto c: conns) {
                        if (visited.find(c.location.gid) == visited.end()) {
                            q.push(c.location.gid);
                            visited[c.location.gid] = true;
                        }
                    }
                }
                super_cells.push_back(cg);
            }
        }
        else {
            // If cell has no gap_junctions, put in separate group of independent cells
            reg_cells.push_back(gid);
        }
    }

    // Sort super_cell groups and only keep those where the first element in the group belongs to domain
    super_cells.erase(std::remove_if(super_cells.begin(), super_cells.end(),
            [gid_part, domain_id](std::vector<cell_gid_type>& cg)
            {
                std::sort(cg.begin(), cg.end());
                return cg.front() < gid_part[domain_id].first;
            }), super_cells.end());

    // Collect local gids that belong to this rank, and sort gids into kind lists
    // kind_lists maps a cell_kind to a vector of either:
    // 1. gids of regular cells (in reg_cells)
    // 2. indices of supercells (in super_cells)

    std::vector<cell_gid_type> local_gids;
    std::unordered_map<cell_kind, std::vector<cell_identifier>> kind_lists;
    for (auto gid: reg_cells) {
        local_gids.push_back(gid);
        kind_lists[rec.get_cell_kind(gid)].push_back({gid, false});
    }

    for (unsigned i = 0; i < super_cells.size(); i++) {
        auto kind = rec.get_cell_kind(super_cells[i].front());
        for (auto gid: super_cells[i]) {
            if (rec.get_cell_kind(gid) != kind) {
                throw arbor_internal_error("Cells of different kinds connected by gap_junctions: not allowed");
            }
            local_gids.push_back(gid);
        }
        kind_lists[kind].push_back({i, true});
    }


    // Create a flat vector of the cell kinds present on this node,
    // partitioned such that kinds for which GPU implementation are
    // listed before the others. This is a very primitive attempt at
    // scheduling; the cell_groups that run on the GPU will be executed
    // before other cell_groups, which is likely to be more efficient.
    //
    // TODO: This creates an dependency between the load balancer and
    // the threading internals. We need support for setting the priority
    // of cell group updates according to rules such as the back end on
    // which the cell group is running.

    auto has_gpu_backend = [&ctx](cell_kind c) {
        return cell_kind_supported(c, backend_kind::gpu, *ctx);
    };

    std::vector<cell_kind> kinds;
    for (auto l: kind_lists) {
        kinds.push_back(cell_kind(l.first));
    }
    std::partition(kinds.begin(), kinds.end(), has_gpu_backend);

    std::vector<group_description> groups;
    for (auto k: kinds) {
        partition_hint hint;
        if (auto opt_hint = util::value_by_key(hint_map, k)) {
            hint = opt_hint.value();
        }

        backend_kind backend = backend_kind::multicore;
        std::size_t group_size = hint.cpu_group_size;

        if (hint.prefer_gpu && gpu_avail && has_gpu_backend(k)) {
            backend = backend_kind::gpu;
            group_size = hint.gpu_group_size;
        }

        std::vector<cell_gid_type> group_elements;
        // group_deps provides information about super_cells
        // group_elements are sorted such that the gids of all members of a super_cell are consecutive.
        // group_deps identifies which cells belong to the same super_cell
        std::vector<int> group_deps;
        for (auto cell: kind_lists[k]) {
            if (cell.is_super_cell == false) {
                group_elements.push_back(cell.id);
                group_deps.push_back(0);
            } else {
                if (group_elements.size() + super_cells[cell.id].size() > group_size && !group_elements.empty()) {
                    groups.push_back({k, std::move(group_elements), backend});
                    group_elements.clear();
                    group_deps.clear();
                }
                for (auto gid: super_cells[cell.id]) {
                    group_elements.push_back(gid);
                    if (gid == super_cells[cell.id].front()) {
                       group_deps.push_back(super_cells[cell.id].size());
                    }
                    else {
                        group_deps.push_back(0);
                    }
                }
            }
            if (group_elements.size()>=group_size) {
                groups.push_back({k, std::move(group_elements), backend});
                group_elements.clear();
                group_deps.clear();
            }
        }
        if (!group_elements.empty()) {
            groups.push_back({k, std::move(group_elements), backend});
        }
    }

    cell_size_type num_local_cells = local_gids.size();

    // Exchange gid list with all other nodes

    util::sort(local_gids);

    // global all-to-all to gather a local copy of the global gid list on each node.
    auto global_gids = ctx->distributed->gather_gids(local_gids);

    domain_decomposition d;
    d.num_domains = num_domains;
    d.domain_id = domain_id;
    d.num_local_cells = num_local_cells;
    d.num_global_cells = num_global_cells;
    d.groups = std::move(groups);
    d.gid_domain = partition_gid_domain(std::move(global_gids), num_domains);

    return d;
}

} // namespace arb

