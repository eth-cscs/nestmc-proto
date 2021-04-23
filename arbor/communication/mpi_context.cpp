// Attempting to acquire an MPI context without MPI enabled will produce
// a link error.

#ifndef ARB_HAVE_MPI
#error "build only if MPI is enabled"
#endif

#include <string>
#include <vector>

#include <mpi.h>

#include <arbor/spike.hpp>

#include "communication/mpi.hpp"
#include "distributed_context.hpp"
#include "label_resolver.hpp"
#include "util/partition.hpp"

namespace arb {

// Throws arb::mpi::mpi_error if MPI calls fail.
struct mpi_context_impl {
    int size_;
    int rank_;
    MPI_Comm comm_;

    explicit mpi_context_impl(MPI_Comm comm): comm_(comm) {
        size_ = mpi::size(comm_);
        rank_ = mpi::rank(comm_);
    }

    gathered_vector<arb::spike>
    gather_spikes(const std::vector<arb::spike>& local_spikes) const {
        return mpi::gather_all_with_partition(local_spikes, comm_);
    }

    gathered_vector<cell_gid_type>
    gather_gids(const std::vector<cell_gid_type>& local_gids) const {
        return mpi::gather_all_with_partition(local_gids, comm_);
    }

    cell_labeled_ranges gather_cell_labeled_ranges(const cell_labeled_ranges& local_ranges) const {
        cell_labeled_ranges global_ranges;
        global_ranges.gids   = mpi::gather_all(local_ranges.gids, comm_);
        global_ranges.sizes  = mpi::gather_all(local_ranges.sizes, comm_);
        global_ranges.labels = mpi::gather_all(local_ranges.labels, comm_);
        global_ranges.ranges = mpi::gather_all(local_ranges.ranges, comm_);
        return global_ranges;
    }

    template <typename T>
    std::vector<T> gather(T value, int root) const {
        return mpi::gather(value, root, comm_);
    }

    std::string name() const { return "MPI"; }
    int id() const { return rank_; }
    int size() const { return size_; }

    template <typename T>
    T min(T value) const {
        return mpi::reduce(value, MPI_MIN, comm_);
    }

    template <typename T>
    T max(T value) const {
        return mpi::reduce(value, MPI_MAX, comm_);
    }

    template <typename T>
    T sum(T value) const {
        return mpi::reduce(value, MPI_SUM, comm_);
    }

    void barrier() const {
        mpi::barrier(comm_);
    }
};

template <>
std::shared_ptr<distributed_context> make_mpi_context(MPI_Comm comm) {
    return std::make_shared<distributed_context>(mpi_context_impl(comm));
}

} // namespace arb
