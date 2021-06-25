#include "memory.hpp"

#ifdef __linux__
extern "C" {
    #include <malloc.h>
}
#endif

#ifdef ARB_HAVE_GPU
    #include <arbor/gpu/gpu_api.hpp>
#endif

namespace arb {
namespace hw {

#if defined(__linux__)
memory_size_type allocated_memory() {
    auto m = mallinfo();
    return m.hblkhd + m.uordblks;
}
#else
memory_size_type allocated_memory() {
    return -1;
}
#endif

#ifdef ARB_HAVE_GPU
memory_size_type gpu_allocated_memory() {
    std::size_t free;
    std::size_t total;
    auto success = gpu::device_mem_get_info(&free, &total);

    return success? total-free: -1;
}
#else
memory_size_type gpu_allocated_memory() {
    return -1;
}
#endif

} // namespace hw
} // namespace arb
