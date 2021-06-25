#pragma once

#include <vector>

#include <arbor/gpu/gpu_api.hpp>

/*
 * Helpers for using GPU memory in unit tests.
 *
 * The memory helpers can't be used in .cu files, because we don't let nvcc
 * compile most of our headers to avoid compiler bugs and c++ version issues.
 */

template <typename T>
struct gpu_ref_proxy {
    T* ptr;

    gpu_ref_proxy(T* p): ptr(p) {}

    gpu_ref_proxy& operator=(const T& value) {
        arb::gpu::device_memcpy(ptr, &value, sizeof(T), arb::gpu::gpuMemcpyHostToDevice);
        return *this;
    }

    operator T() const {
        T tmp;
        arb::gpu::device_memcpy(&tmp, ptr, sizeof(T), arb::gpu::gpuMemcpyDeviceToHost);
        return tmp;
    }
};

template <typename T>
class gpu_vector {
    using value_type = T;
    using size_type = std::size_t;

public:
    gpu_vector() = default;

    gpu_vector(size_type n) {
        allocate(n);
    }

    gpu_vector(const std::vector<T>& other) {
        allocate(other.size());
        to_device(other.data());
    }

    ~gpu_vector() {
        if (data_) arb::gpu::device_free(data_);
    }

    std::vector<T> host_vector() const {
        std::vector<T> v(size());
        to_host(v.data());
        return v;
    }

    value_type* data() {
        return data_;
    }

    const value_type* data() const {
        return data_;
    }

    size_type size() const {
        return size_;
    }

    value_type operator[](size_type i) const {
        value_type tmp;
        arb::gpu::device_memcpy(&tmp, data_+i, sizeof(value_type), arb::gpu::gpuMemcpyDeviceToHost);
        return tmp;
    }

    gpu_ref_proxy<value_type> operator[](size_type i) {
        return gpu_ref_proxy<value_type>(data_+i);
    }

private:

    void allocate(size_type n) {
        size_ = n;
        arb::gpu::device_malloc(&data_, n*sizeof(T));
    }

    void to_device(const value_type* other) {
        arb::gpu::device_memcpy(data_, other, size_in_bytes(), arb::gpu::gpuMemcpyHostToDevice);
    }

    void to_host(value_type* other) const {
        arb::gpu::device_memcpy(other, data_, size_in_bytes(), arb::gpu::gpuMemcpyDeviceToHost);
    }

    size_type size_in_bytes() const {
        return size_*sizeof(value_type);
    }

    void free() {
        arb::gpu::device_free(data_);
        size_ = 0;
        data_ = nullptr;
    }

    size_type size_;
    value_type* data_;
};
