#include <arbor/fvm_types.hpp>
#include <arbor/gpu/gpu_common.hpp>

#include "matrix_common.hpp"

namespace arb {
namespace gpu {

namespace kernels {

/// GPU implementation of Hines Matrix solver.
/// Flat format
/// p: parent index for each variable. Needed for backward and forward sweep
template <typename T, typename I>
__global__
void solve_matrix_flat(
    T* rhs, T* d, const T* u, const I* p, const I* cell_cv_divs, int num_mtx)
{
    auto tid = threadIdx.x + blockDim.x*blockIdx.x;

    if (tid<num_mtx) {
        // get range of this thread's cell matrix
        const auto first = cell_cv_divs[tid];
        const auto last  = cell_cv_divs[tid+1];

        // Zero diagonal term implies dt==0; just leave rhs (for whole matrix)
        // alone in that case.
        if (d[last-1]==0) return;

        // backward sweep
        for (auto i=last-1; i>first; --i) {
            auto factor = u[i] / d[i];
            d[p[i]]   -= factor * u[i];
            rhs[p[i]] -= factor * rhs[i];
        }
        rhs[first] /= d[first];

        // forward sweep
        for(auto i=first+1; i<last; ++i) {
            rhs[i] -= u[i] * rhs[p[i]];
            rhs[i] /= d[i];
        }
    }
}

/// GPU implementation of Hines Matrix solver.
/// Block-interleaved format
template <typename T, typename I, int BlockWidth>
__global__
void solve_matrix_interleaved(
    T* rhs, T* d, const T* u, const I* p, const I* sizes, int padded_size, int num_mtx)
{
    auto tid = threadIdx.x + blockDim.x*blockIdx.x;

    if (tid<num_mtx) {
        const auto block       = tid/BlockWidth;
        const auto block_start = block*BlockWidth;
        const auto block_lane  = tid - block_start;

        // get range of this thread's cell matrix
        const auto first    = block_start*padded_size + block_lane;
        const auto last     = first + BlockWidth*(sizes[tid]-1);
        const auto last_max = first + BlockWidth*(sizes[block_start]-1);

        // Zero diagonal term implies dt==0; just leave rhs (for whole matrix)
        // alone in that case.
        if (d[last]==0) return;

        // backward sweep
        for(auto i=last_max; i>first; i-=BlockWidth) {
            if (i<=last) {
                auto factor = u[i] / d[i];
                d[p[i]]   -= factor * u[i];
                rhs[p[i]] -= factor * rhs[i];
            }
        }
        rhs[first] /= d[first];

        // forward sweep
        for(auto i=first+BlockWidth; i<=last; i+=BlockWidth) {
            rhs[i] -= u[i] * rhs[p[i]];
            rhs[i] /= d[i];
        }
    }
}

} // namespace kernels

void solve_matrix_flat(
    fvm_value_type* rhs,
    fvm_value_type* d,
    const fvm_value_type* u,
    const fvm_index_type* p,
    const fvm_index_type* cell_cv_divs,
    int num_mtx)
{
    constexpr unsigned block_dim = 128;
    const unsigned grid_dim = impl::block_count(num_mtx, block_dim);
    kernels::solve_matrix_flat
        <fvm_value_type, fvm_index_type>
        <<<grid_dim, block_dim>>>
        (rhs, d, u, p, cell_cv_divs, num_mtx);
}

void solve_matrix_interleaved(
    fvm_value_type* rhs,
    fvm_value_type* d,
    const fvm_value_type* u,
    const fvm_index_type* p,
    const fvm_index_type* sizes,
    int padded_size,
    int num_mtx)
{
    constexpr unsigned block_dim = impl::matrices_per_block();
    const unsigned grid_dim = impl::block_count(num_mtx, block_dim);
    kernels::solve_matrix_interleaved<fvm_value_type, fvm_index_type, block_dim>
        <<<grid_dim, block_dim>>>
        (rhs, d, u, p, sizes, padded_size, num_mtx);
}

} // namespace gpu
} // namespace arb
