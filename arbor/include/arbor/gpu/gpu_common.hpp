#pragma once

#include "gpu_api.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#   define HOST_DEVICE_IF_GPU __host__ __device__
#else
#   define HOST_DEVICE_IF_GPU
#endif

namespace arb {
namespace gpu {

namespace impl {
// Number of threads per warp
// This has always been 32, however it may change in future NVIDIA gpus
HOST_DEVICE_IF_GPU
constexpr inline unsigned threads_per_warp() {
#ifdef ARB_HIP
    return 64u;
#else
    return 32u;
#endif
}

// The minimum number of bins required to store n values where the bins have
// dimension of block_size.
HOST_DEVICE_IF_GPU
constexpr inline unsigned block_count(unsigned n, unsigned block_size) {
    return (n+block_size-1)/block_size;
}

} // namespace impl

} // namespace gpu
} // namespace arb
