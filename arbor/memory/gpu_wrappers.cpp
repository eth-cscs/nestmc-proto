#include <cstdlib>
#include <string>

#include <arbor/arbexcept.hpp>

#include "util.hpp"

#ifdef ARB_HAVE_GPU

#include <arbor/gpu/gpu_api.hpp>

#define HANDLE_GPU_ERROR(error, msg)\
throw arbor_exception("GPU memory:: "+std::string(__func__)+" "+std::string((msg))+": "+error.description());

namespace arb {
namespace memory {

using std::to_string;
using namespace gpu;
 
void gpu_memcpy_d2d(void* dest, const void* src, std::size_t n) {
    auto status = device_memcpy(dest, src, n, gpuMemcpyDeviceToDevice);
    if (!status) {
        HANDLE_GPU_ERROR(status, "n="+to_string(n));
    }
}

void gpu_memcpy_d2h(void* dest, const void* src, std::size_t n) {
    auto status = device_memcpy(dest, src, n, gpuMemcpyDeviceToHost);
    if (!status) {
        HANDLE_GPU_ERROR(status, "n="+to_string(n));
    }
}

void gpu_memcpy_h2d(void* dest, const void* src, std::size_t n) {
    auto status = device_memcpy(dest, src, n, gpuMemcpyHostToDevice);
    if (!status) {
        HANDLE_GPU_ERROR(status, "n="+to_string(n));
    }
}

void* gpu_host_register(void* ptr, std::size_t size) {
    auto status = host_register(ptr, size, gpuHostRegisterPortable);
    if (!status) {
        HANDLE_GPU_ERROR(status, "unable to register host memory");
    }
    return ptr;
}

void gpu_host_unregister(void* ptr) {
    host_unregister(ptr);
}

void* gpu_malloc(std::size_t n) {
    void* ptr;

    auto status = device_malloc(&ptr, n);
    if (!status) {
        HANDLE_GPU_ERROR(status, "unable to allocate "+to_string(n)+" bytes");
    }
    return ptr;
}

void gpu_free(void* ptr) {
    auto status = device_free(ptr);
    if (!status) {
        HANDLE_GPU_ERROR(status, "");
    }
}

} // namespace memory
} // namespace arb

#else

#define NOGPU \
LOG_ERROR("memory:: "+std::string(__func__)+"(): no GPU support")

namespace arb {
namespace memory {

void gpu_memcpy_d2d(void* dest, const void* src, std::size_t n) {
    NOGPU;
}

void gpu_memcpy_d2h(void* dest, const void* src, std::size_t n) {
    NOGPU;
}

void gpu_memcpy_h2d(void* dest, const void* src, std::size_t n) {
    NOGPU;
}

void* gpu_host_register(void* ptr, std::size_t size) {
    NOGPU;
    return 0;
}

void gpu_host_unregister(void* ptr) {
    NOGPU;
}

void* gpu_malloc(std::size_t n) {
    NOGPU;
    return 0;
}

void gpu_free(void* ptr) {
    NOGPU;
}

} // namespace memory
} // namespace arb

#endif // def ARB_HAVE_GPU

