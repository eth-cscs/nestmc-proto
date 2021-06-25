#include <arbor/fvm_types.hpp>

#include <arbor/gpu/gpu_common.hpp>
#include "matrix_common.hpp"

namespace arb {
namespace gpu {

namespace kernels {
/// GPU implementatin of Hines matrix assembly
/// Flat layout
/// For a given time step size dt
///     - use the precomputed alpha and alpha_d values to construct the diagonal
///       and off diagonal of the symmetric Hines matrix.
///     - compute the RHS of the linear system to solve
template <typename T, typename I>
__global__
void assemble_matrix_flat(
        T* d,
        T* rhs,
        const T* invariant_d,
        const T* voltage,
        const T* current,
        const T* conductivity,
        const T* cv_capacitance,
        const T* cv_area,
        const I* cv_to_cell,
        const T* dt_intdom,
        const I* cell_to_intdom,
        unsigned n)
{
    const unsigned tid = threadIdx.x + blockDim.x*blockIdx.x;

    if (tid<n) {
        auto cid = cv_to_cell[tid];
        auto dt = dt_intdom[cell_to_intdom[cid]];

        // Note: dt==0 case is expected only at the end of a mindelay/2
        // integration period, and consequently divergence is unlikely
        // to be a peformance problem.

        if (dt>0) {
            // The 1e-3 is a constant of proportionality required to ensure that the
            // conductance (gi) values have units μS (micro-Siemens).
            // See the model documentation in docs/model for more information.
            T oodt_factor = 1e-3/dt; // [1/μs]
            T area_factor = 1e-3*cv_area[tid]; // [1e-9·m²]

            auto gi = oodt_factor * cv_capacitance[tid] + area_factor*conductivity[tid]; // [μS]
            d[tid] = gi + invariant_d[tid];
            rhs[tid] = gi*voltage[tid] - area_factor*current[tid];
        }
        else {
            d[tid] = 0;
            rhs[tid] = voltage[tid];
        }
    }
}

/// GPU implementatin of Hines matrix assembly
/// Interleaved layout
/// For a given time step size dt
///     - use the precomputed alpha and alpha_d values to construct the diagonal
///       and off diagonal of the symmetric Hines matrix.
///     - compute the RHS of the linear system to solve
template <typename T, typename I, unsigned BlockWidth, unsigned LoadWidth, unsigned Threads>
__global__
void assemble_matrix_interleaved(
        T* d,
        T* rhs,
        const T* invariant_d,
        const T* voltage,
        const T* current,
        const T* conductivity,
        const T* cv_capacitance,
        const T* area,
        const I* sizes,
        const I* starts,
        const I* matrix_to_cell,
        const T* dt_intdom,
        const I* cell_to_intdom,
        unsigned padded_size, unsigned num_mtx)
{
    static_assert(BlockWidth*LoadWidth==Threads,
        "number of threads must equal number of values to process per block");
    __shared__ T buffer_v[Threads];
    __shared__ T buffer_i[Threads];
    __shared__ T buffer_g[Threads];

    const unsigned tid = threadIdx.x + blockIdx.x*blockDim.x;
    const unsigned lid = threadIdx.x;

    const unsigned mtx_id   = tid/LoadWidth;
    const unsigned mtx_lane = tid - mtx_id*LoadWidth;

    const unsigned blk_id   = tid/(BlockWidth*LoadWidth);
    const unsigned blk_row  = lid/BlockWidth;
    const unsigned blk_lane = lid - blk_row*BlockWidth;

    const unsigned blk_pos  = LoadWidth*blk_lane + blk_row;

    const bool do_load  = mtx_id<num_mtx;

    unsigned load_pos  = do_load? starts[mtx_id] + mtx_lane     : 0;
    const unsigned end = do_load? starts[mtx_id] + sizes[mtx_id]: 0;
    unsigned store_pos = blk_id*BlockWidth*padded_size + (blk_row*BlockWidth + blk_lane);

    const unsigned max_size = sizes[0];

    T oodt_factor = 0;
    T dt = 0;
    const unsigned permuted_cid = blk_id*BlockWidth + blk_lane;

    if (permuted_cid<num_mtx) {
        auto cid = matrix_to_cell[permuted_cid];
        dt = dt_intdom[cell_to_intdom[cid]];

        // The 1e-3 is a constant of proportionality required to ensure that the
        // conductance (gi) values have units μS (micro-Siemens).
        // See the model documentation in docs/model for more information.

        oodt_factor = dt>0? T(1e-3)/dt: 0;
    }

    for (unsigned j=0u; j<max_size; j+=LoadWidth) {
        if (do_load && load_pos<end) {
            buffer_v[lid] = voltage[load_pos];
            buffer_i[lid] = current[load_pos];
            buffer_g[lid] = conductivity[load_pos];
        }

        __syncthreads();

        if (j+blk_row<padded_size) {
            T area_factor = T(1e-3)*area[store_pos];
            const auto gi = oodt_factor*cv_capacitance[store_pos] + area_factor*buffer_g[blk_pos];

            if (dt>0) {
                d[store_pos]   = (gi + invariant_d[store_pos]);
                rhs[store_pos] = (gi*buffer_v[blk_pos] - area_factor*buffer_i[blk_pos]);
            }
            else {
                d[store_pos]   = 0;
                rhs[store_pos] = buffer_v[blk_pos];
            }
        }

        __syncthreads();

        store_pos += LoadWidth*BlockWidth;
        load_pos  += LoadWidth;
    }
}

} // namespace kernels

void assemble_matrix_flat(
        fvm_value_type* d,
        fvm_value_type* rhs,
        const fvm_value_type* invariant_d,
        const fvm_value_type* voltage,
        const fvm_value_type* current,
        const fvm_value_type* conductivity,
        const fvm_value_type* cv_capacitance,
        const fvm_value_type* area,
        const fvm_index_type* cv_to_cell,
        const fvm_value_type* dt_intdom,
        const fvm_index_type* cell_to_intdom,
        unsigned n)
{
    constexpr unsigned block_dim = 128;
    const unsigned grid_dim = impl::block_count(n, block_dim);

    kernels::assemble_matrix_flat
        <fvm_value_type, fvm_index_type>
        <<<grid_dim, block_dim>>>
        (d, rhs, invariant_d, voltage, current, conductivity, cv_capacitance,
         area, cv_to_cell, dt_intdom, cell_to_intdom, n);
}

//template <typename T, typename I, unsigned BlockWidth, unsigned LoadWidth, unsigned Threads>
void assemble_matrix_interleaved(
    fvm_value_type* d,
    fvm_value_type* rhs,
    const fvm_value_type* invariant_d,
    const fvm_value_type* voltage,
    const fvm_value_type* current,
    const fvm_value_type* conductivity,
    const fvm_value_type* cv_capacitance,
    const fvm_value_type* area,
    const fvm_index_type* sizes,
    const fvm_index_type* starts,
    const fvm_index_type* matrix_to_cell,
    const fvm_value_type* dt_intdom,
    const fvm_index_type* cell_to_intdom,
    unsigned padded_size, unsigned num_mtx)
{
    constexpr unsigned bd = impl::matrices_per_block();
    constexpr unsigned lw = impl::load_width();
    constexpr unsigned block_dim = bd*lw;

    // The number of threads is threads_per_matrix*num_mtx
    const unsigned grid_dim = impl::block_count(num_mtx*lw, block_dim);

    kernels::assemble_matrix_interleaved
        <fvm_value_type, fvm_index_type, bd, lw, block_dim>
        <<<grid_dim, block_dim>>>
        (d, rhs, invariant_d, voltage, current, conductivity, cv_capacitance, area,
         sizes, starts, matrix_to_cell,
         dt_intdom, cell_to_intdom, padded_size, num_mtx);
}

} // namespace gpu
} // namespace arb
