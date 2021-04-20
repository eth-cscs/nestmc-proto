#include <algorithm>
#include <string>
#include <vector>

#include <arbor/spike.hpp>

#include "distributed_context.hpp"
#include "label_resolver.hpp"
#include "threading/threading.hpp"

namespace arb {

struct dry_run_context_impl {

    explicit dry_run_context_impl(unsigned num_ranks, unsigned num_cells_per_tile):
        num_ranks_(num_ranks), num_cells_per_tile_(num_cells_per_tile) {};

    gathered_vector<arb::spike>
    gather_spikes(const std::vector<arb::spike>& local_spikes) const {
        using count_type = typename gathered_vector<arb::spike>::count_type;

        count_type local_size = local_spikes.size();

        std::vector<arb::spike> gathered_spikes;
        gathered_spikes.reserve(local_size*num_ranks_);

        for (count_type i = 0; i < num_ranks_; i++) {
            gathered_spikes.insert(gathered_spikes.end(), local_spikes.begin(), local_spikes.end());
        }

        for (count_type i = 0; i < num_ranks_; i++) {
            for (count_type j = i*local_size; j < (i+1)*local_size; j++){
                gathered_spikes[j].source.gid += num_cells_per_tile_*i;
            }
        }

        std::vector<count_type> partition;
        for (count_type i = 0; i <= num_ranks_; i++) {
            partition.push_back(static_cast<count_type>(i*local_size));
        }

        return gathered_vector<arb::spike>(std::move(gathered_spikes), std::move(partition));
    }

    gathered_vector<cell_gid_type>
    gather_gids(const std::vector<cell_gid_type>& local_gids) const {
        using count_type = typename gathered_vector<cell_gid_type>::count_type;

        count_type local_size = local_gids.size();

        std::vector<cell_gid_type> gathered_gids;
        gathered_gids.reserve(local_size*num_ranks_);

        for (count_type i = 0; i < num_ranks_; i++) {
            gathered_gids.insert(gathered_gids.end(), local_gids.begin(), local_gids.end());
        }

        for (count_type i = 0; i < num_ranks_; i++) {
            for (count_type j = i*local_size; j < (i+1)*local_size; j++){
                gathered_gids[j] += num_cells_per_tile_*i;
            }
        }

        std::vector<count_type> partition;
        for (count_type i = 0; i <= num_ranks_; i++) {
            partition.push_back(static_cast<count_type>(i*local_size));
        }

        return gathered_vector<cell_gid_type>(std::move(gathered_gids), std::move(partition));
    }

    int id() const { return 0; }

    int size() const { return num_ranks_; }

    template <typename T>
    T min(T value) const { return value; }

    template <typename T>
    T max(T value) const { return value; }

    template <typename T>
    T sum(T value) const { return value * num_ranks_; }

    template <typename T>
    std::vector<T> gather(T value, int) const {
        return std::vector<T>(num_ranks_, value);
    }

    cell_labeled_ranges gather_labeled_range(const cell_labeled_ranges& local_ranges) const {
        std::vector<cell_gid_type> gids;
        std::vector<cell_tag_type> labels;
        std::vector<lid_range> ranges;
        std::vector<std::size_t> sorted_partitions = {0};

        gids.reserve(local_ranges.gids.size()*num_ranks_);
        labels.reserve(local_ranges.labels.size()*num_ranks_);
        ranges.reserve(local_ranges.ranges.size()*num_ranks_);

        for (unsigned i = 0; i < num_ranks_; i++) {
            arb_assert(local_ranges.is_one_partition());
            std::transform(local_ranges.gids.begin(), local_ranges.gids.end(), gids.begin(), [&](cell_gid_type gid){return gid+num_cells_per_tile_*i;});
            labels.insert(labels.end(), local_ranges.labels.begin(), local_ranges.labels.end());
            ranges.insert(ranges.end(), local_ranges.ranges.begin(), local_ranges.ranges.end());
            sorted_partitions.push_back(sorted_partitions.back() + local_ranges.gids.size());
        }
        return cell_labeled_ranges(gids, labels, ranges, sorted_partitions);
    }

    void barrier() const {}

    std::string name() const { return "dryrun"; }

    unsigned num_ranks_;
    unsigned num_cells_per_tile_;
};

std::shared_ptr<distributed_context> make_dry_run_context(unsigned num_ranks, unsigned num_cells_per_tile) {
    return std::make_shared<distributed_context>(dry_run_context_impl(num_ranks, num_cells_per_tile));
}

} // namespace arb
