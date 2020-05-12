#include <utility>
#include <vector>

#include <arbor/assert.hpp>
#include <arbor/common_types.hpp>
#include <arbor/domain_decomposition.hpp>
#include <arbor/recipe.hpp>
#include <arbor/spike.hpp>

#include "algorithms.hpp"
#include "communication/gathered_vector.hpp"
#include "connection.hpp"
#include "distributed_context.hpp"
#include "execution_context.hpp"
#include "profile/profiler_macro.hpp"
#include "threading/threading.hpp"
#include "util/partition.hpp"
#include "util/rangeutil.hpp"
#include "util/span.hpp"

#include "communication/communicator.hpp"

namespace arb {

// For caching information about each cell
struct gid_info {
    using connection_list = decltype(std::declval<recipe>().connections_on(0));
    gid_info() = default;           // so we can in a std::vector
    gid_info(cell_gid_type g, cell_size_type di, connection_list c):
        gid(g), index_on_domain(di), conns(std::move(c)) {}
    
    cell_gid_type gid;              // global identifier of cell
    cell_size_type index_on_domain; // index of cell in this domain
    connection_list conns;          // list of connections terminating at this cell
};

struct chunk_info_type {
    using cell_list = communicator::cell_list;
    
    cell_list conns_part;  // partition local cells globally, and not by chunks (for connections_ext_)
    cell_list chunk_part;  // partition local cells into chunks [0, n-chunk-1, n-chunk-1+2 .... num-cells]
    cell_size_type num_chunks = 0;

    chunk_info_type(std::vector<gid_info>& gid_infos, int threads, cell_size_type n_cons)
    {
        using util::make_span;
        using util::make_partition;

        cell_list chunk_conns; // number of conns by chunk
        cell_list chunk_ncells; // number of ids by chunk

        const auto conns_per_thread = n_cons / threads;
        chunk_conns.reserve(threads+1); // at most, threads+1
        chunk_ncells.reserve(threads+1);
        
        // split cells into approximately equal numbers of connection chunks with approximately thread-number chunks
        cell_size_type conns_used = 0; // number of conns in current chunk    
        cell_size_type ncells = 0; // offset of next gid in chunk
        for (auto&& id: make_span(gid_infos.size())) {
            ++ncells;
            
            auto&& cell = gid_infos[id];
            conns_used += cell.conns.size();
            if (conns_used >= conns_per_thread) {
                chunk_conns.push_back(conns_used);
                chunk_ncells.push_back(ncells);
                
                conns_used = 0;
                ncells = 0;
                ++num_chunks;
            }
        }
        // any conns left over?
        if (ncells) {
            chunk_conns.push_back(conns_used);
            chunk_ncells.push_back(ncells);
            ++num_chunks;
        }

        make_partition(conns_part, chunk_conns);
        make_partition(chunk_part, chunk_ncells);
    }
};

communicator::communicator(const recipe& rec,
                           const domain_decomposition& dom_dec,
                           execution_context& ctx)
{
    using util::subrange_view;

    distributed_ = ctx.distributed;
    thread_pool_ = ctx.thread_pool;

    num_domains_ = distributed_->size();
    num_local_groups_ = dom_dec.groups.size();
    num_local_cells_ = dom_dec.num_local_cells;

    // Make a list of local gid with their group index and connections
    //   -> gid_infos
    // Count the number of local connections (i.e. connections terminating on this domain)
    //   -> n_cons: scalar
    // Calculate and store domain id of the presynaptic cell on each local connection
    //   -> src_domains: array with one entry for every local connection
    // Also the count of presynaptic sources from each domain
    //   -> src_counts: array with one entry for each domain

    // Record all the gid in a flat vector.
    // These are used to map from local index to gid in the parallel loop
    // that populates gid_infos.
    std::vector<cell_gid_type> gids;
    gids.reserve(num_local_cells_);
    for (auto&& g: dom_dec.groups) {
        util::append(gids, g.gids);
    }

    // Build cell partition by group for passing events to cell groups
    index_part_ = util::make_partition( // thus we have a partition of local cells by group id
        index_divisions_,
        util::transform_view(// index on dom_dec.groups -> number of cells in group
            dom_dec.groups,
            [](const group_description& g){return g.gids.size();}));
    
    // Build the connection information for local cells in parallel.
    std::vector<gid_info> gid_infos;
    gid_infos.resize(num_local_cells_);
    threading::parallel_for::apply(0, gids.size(), thread_pool_.get(),
        [&](cell_size_type i) {
            auto gid = gids[i];
            gid_infos[i] = gid_info{gid, i, rec.connections_on(gid)};
        });

    cell_local_size_type n_cons =
        util::sum_by(gid_infos, [](const gid_info& g){ return g.conns.size(); });
    const auto threads = thread_pool_->get_num_threads();

    chunk_info_type chunk_info{gid_infos, threads, n_cons};
    num_chunks_ = chunk_info.num_chunks;

    // now partitions applied by chunk
    connections_.resize(num_chunks_);
    connection_part_.resize(num_chunks_);
    connections_ext_.resize(n_cons);
    threading::parallel_for::apply(0, num_chunks_, thread_pool_.get(),
        [&](cell_size_type chunk) {
            auto chunk_gid_infos = subrange_view(gid_infos, chunk_info.chunk_part[chunk], chunk_info.chunk_part[chunk+1]);
            auto conn_off = chunk_info.conns_part[chunk];
            const auto chunk_n_cons = chunk_info.conns_part[chunk+1] - conn_off;

            auto& chunk_connections = connections_[chunk];
    
            std::vector<unsigned> src_domains;
            src_domains.reserve(chunk_n_cons);
            std::vector<cell_size_type> src_counts(num_domains_);
    
            for (auto&& cell: chunk_gid_infos) {
                for (auto&& c: cell.conns) {
                    auto src = dom_dec.gid_domain(c.source.gid);
                    src_domains.push_back(src);
                    src_counts[src]++;
                    connections_ext_[conn_off++] = {c.source, c.dest, c.weight, c.delay, cell.index_on_domain};
                }
            }

            // Construct the connections.
            // The loop above gave the information required to construct in place
            // the connections as partitioned by the domain of their source gid.
            chunk_connections.resize(chunk_n_cons);
            connection_part_[chunk] = algorithms::make_index(src_counts);
            auto& chunk_connection_part = connection_part_[chunk];
            
            auto offsets = chunk_connection_part;
            std::size_t pos = 0;
            for (auto&& cell: chunk_gid_infos) {
                for (auto&& c: cell.conns) {
                    const auto i = offsets[src_domains[pos]]++;
                    chunk_connections[i] = {c.source, c.dest, c.weight, c.delay, cell.index_on_domain};
                    ++pos;
                }
            }

            // Sort the connections for each domain.
            // This is num_domains_ independent sorts, so it can be parallelized trivially.
            const auto& cp = chunk_connection_part;
            threading::parallel_for::apply(0, num_domains_, thread_pool_.get(),
                [&](cell_size_type i) {
                    util::sort(util::subrange_view(chunk_connections, cp[i], cp[i+1]));
                });
        });
}

std::pair<cell_size_type, cell_size_type> communicator::group_queue_range(cell_size_type i) {
    arb_assert(i<num_local_groups_);
    return index_part_[i];
}

time_type communicator::min_delay() {
    auto local_min = std::numeric_limits<time_type>::max();
    for (auto&& chunk_connections : connections_) {
        for (auto&& con: chunk_connections) {
            local_min = std::min(local_min, con.delay());
        }
    }

    return distributed_->min(local_min);
}

gathered_vector<spike> communicator::exchange(std::vector<spike> local_spikes) {
    PE(communication_exchange_sort);
    // sort the spikes in ascending order of source gid
    util::sort_by(local_spikes, [](spike s){return s.source;});
    PL();

    PE(communication_exchange_gather);
    // global all-to-all to gather a local copy of the global spike list on each node.
    auto global_spikes = distributed_->gather_spikes(local_spikes);
    num_spikes_ += global_spikes.size();
    PL();

    return global_spikes;
}

static void make_queues_by_conns(
    std::vector<pse_vector>& queues,
    std::vector<connection>::iterator cn,
    const std::vector<connection>::iterator cend,
    std::vector<spike>::const_iterator sp,
    const std::vector<spike>::const_iterator send)
{
    using util::make_range;
    
    struct spike_pred {
        bool operator()(const spike& spk, const cell_member_type& src)
            {return spk.source<src;}
        bool operator()(const cell_member_type& src, const spike& spk)
            {return src<spk.source;}
    };

    for (auto&& c: make_range(cn, cend)) {
        if (sp == send) {break;};
        
        auto spikes = std::equal_range(sp, send, c.source(), spike_pred());
        if (spikes.first != spikes.second) {
            auto& q = queues[c.index_on_domain()];
            for (auto&& s: make_range(spikes)) {
                q.push_back(c.make_event(s));
            }
        }
        
        sp = spikes.first; // should be first, range of connections may have same source
    }
}

static void make_queues_by_spikes(
    std::vector<pse_vector>& queues,
    std::vector<connection>::iterator cn,
    const std::vector<connection>::iterator cend,
    std::vector<spike>::const_iterator sp,
    const std::vector<spike>::const_iterator send)
{
    using util::make_range;

    for (auto&& s: make_range(sp, send)) {
        if (cn == cend) {break;};
        
        auto targets = std::equal_range(cn, cend, s.source);
        for (auto&& c: make_range(targets)) {
            auto& q = queues[c.index_on_domain()];
            q.push_back(c.make_event(s));
        }

        cn = targets.second; // range of connections with this source handled
    }
}

void communicator::make_event_queues(
        const gathered_vector<spike>& global_spikes,
        std::vector<pse_vector>& queues)
{
    arb_assert(queues.size()==num_local_cells_);

    using util::subrange_view;
    using util::make_span;

    const auto& sp = global_spikes.partition();
    threading::parallel_for::apply(0, num_chunks_, thread_pool_.get(),
        [&](cell_size_type chunk) {
            const auto& cp = connection_part_[chunk];
            for (auto dom: make_span(num_domains_)) {
                auto cons = subrange_view(connections_[chunk], cp[dom], cp[dom+1]);
                auto spks = subrange_view(global_spikes.values(), sp[dom], sp[dom+1]);
                
                // We have a choice of whether to walk spikes or connections:
                // i.e., we can iterate over the spikes, and for each spike search
                // the for connections that have the same source; or alternatively
                // for each connection, we can search the list of spikes for spikes
                // with the same source.
                //
                // We iterate over whichever set is the smallest, which has
                // complexity of order max(S log(C), C log(S)), where S is the
                // number of spikes, and C is the number of connections.
                if (cons.size()<spks.size()) {
                    make_queues_by_conns(queues, cons.begin(), cons.end(), spks.begin(), spks.end()); 
                }
                else {
                    make_queues_by_spikes(queues, cons.begin(), cons.end(), spks.begin(), spks.end());
                }
            }
        });
}

std::uint64_t communicator::num_spikes() const {
    return num_spikes_;
}

cell_size_type communicator::num_local_cells() const {
    return num_local_cells_;
}

const std::vector<connection>& communicator::connections() const {
    return connections_ext_;
}

void communicator::reset() {
    num_spikes_ = 0;
}

} // namespace arb

