#pragma once

#include <vector>

#include <arbor/common_types.hpp>
#include <arbor/domain_decomposition.hpp>
#include <arbor/recipe.hpp>
#include <arbor/spike.hpp>
#include <arbor/version.hpp>

#include "communication/gathered_vector.hpp"
#include "connection.hpp"
#include "execution_context.hpp"
#include "util/partition.hpp"

namespace arb {

// When the communicator is constructed the number of target groups and targets
// is specified, along with a mapping between local cell id and local
// target id.
//
// The user can add connections to an existing communicator object, where
// each connection is between any global cell and any local target.
//
// Once all connections have been specified, the construct() method can be used
// to build the data structures required for efficient spike communication and
// event generation.

class communicator {
public:
    explicit communicator(const recipe& rec,
                          const domain_decomposition& dom_dec,
                          execution_context& ctx);

    /// The range of event queues that belong to cells in group i.
    std::pair<cell_size_type, cell_size_type> group_queue_range(cell_size_type i);

    /// The minimum delay of all connections in the global network.
    time_type min_delay();

    /// Perform exchange of spikes.
    ///
    /// Takes as input the list of local_spikes that were generated on the calling domain.
    /// Returns the full global set of vectors, along with meta data about their partition
    gathered_vector<spike> exchange(std::vector<spike> local_spikes);

    /// Check each global spike in turn to see it generates local events.
    /// If so, make the events and insert them into the appropriate event list.
    ///
    /// Takes reference to a vector of event lists as an argument, with one list
    /// for each local cell group. On completion, the events in each list are
    /// all events that must be delivered to targets in that cell group as a
    /// result of the global spike exchange, plus any events that were already
    /// in the list.
    void make_event_queues(
            const gathered_vector<spike>& global_spikes,
            std::vector<pse_vector>& queues);

    /// Returns the total number of global spikes over the duration of the simulation
    std::uint64_t num_spikes() const;

    cell_size_type num_local_cells() const;

    const std::vector<connection>& connections() const; // do we need this?

    void reset();

    using connection_list = std::vector<connection>;
    using cell_list = std::vector<cell_size_type>;
    using group_partition = util::partition_view_type<cell_list>; 
    
private:
    cell_size_type num_local_cells_;
    cell_size_type num_local_groups_;
    cell_size_type num_domains_;
    cell_size_type num_chunks_;

    std::vector<connection_list> connections_; // connections broken into chunks, ordered by domains
    connection_list connections_ext_; // same list, but flattened for external use, ordered by cell id (is this needed?)

    std::vector<cell_list> connection_part_; // partition connections_ by cells, partitioned by chunks

    group_partition index_part_; // partition local cells by group index in domain_decomposition
    cell_list index_divisions_; // underlying division for index_part_

    distributed_context_handle distributed_;
    task_system_handle thread_pool_;
    std::uint64_t num_spikes_ = 0u;
};

} // namespace arb
