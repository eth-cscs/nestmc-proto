#include <arbor/fvm_types.hpp>
#include <arbor/gpu/gpu_common.hpp>

namespace arb {
namespace gpu {
namespace kernel {
__global__
void multiply_in_place_(fvm_value_type* s, const fvm_index_type* p, int n) {
    int tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid_<n) {
        s[tid_] *= p[tid_];
    }
}
}

void multiply_in_place_(fvm_value_type* s, const fvm_index_type* p, int n) {
        unsigned block_dim = 128;
        unsigned grid_dim = gpu::impl::block_count(n, block_dim);

        kernel::multiply_in_place_<<<grid_dim, block_dim>>>(s, p, n);
}
} // namespace gpu
} // namespace arb
