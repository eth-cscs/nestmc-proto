#include "../gtest.h"

#include <algorithm>
#include <vector>

#include <arbor/gpu/reduce_by_key.hpp>

#include "gpu_vector.hpp"

using namespace arb;

template <typename T, typename I>
__global__
void reduce_kernel(const T* src, T* dst, const I* index, int n) {
    unsigned tid = threadIdx.x + blockIdx.x*blockDim.x;

    unsigned mask = gpu::ballot(0xffffffff, tid<n);
    if (tid<n) {
        gpu::reduce_by_key(src[tid], dst, index[tid], mask);
    }
}

template <typename T>
std::vector<T> reduce(const std::vector<T>& in, size_t n_out, const std::vector<int>& index, unsigned block_dim=128) {
    EXPECT_EQ(in.size(), index.size());
    EXPECT_TRUE(std::is_sorted(index.begin(), index.end()));

    using array  = gpu_vector<T>;
    using iarray = gpu_vector<int>;

    int n = in.size();

    array  src(in);
    iarray idx(index);
    array  dst(std::vector<T>(n_out, 0));

    unsigned grid_dim = (n-1)/block_dim + 1;
    reduce_kernel<<<grid_dim, block_dim>>>(src.data(), dst.data(), idx.data(), n);

    return dst.host_vector();
}

TEST(reduce_by_key, no_repetitions)
{
    int n = 64;
    std::vector<int> index(n);
    for (int i=0; i<n; ++i) index[i] = i;

    {
        std::vector<float> in(n, 1);

        auto out = reduce(in, n, index);
        for (auto o: out) EXPECT_EQ(o, 1.0f);
    }
    {
        std::vector<double> in(n, 1);

        auto out = reduce(in, n, index);
        for (auto o: out) EXPECT_EQ(o, 1.0);
    }
}

TEST(reduce_by_key, single_repeated_index)
{
    // Perform reduction of a sequence of 1s of length n
    // The expected result is n
    for (auto n: {1, 2, 7, 31, 32, 33, 63, 64, 65, 128}) {
        std::vector<double> in(n, 1);
        std::vector<int> index(n, 0);

        auto out = reduce(in, 1, index, 32);
        EXPECT_EQ(double(n), out[0]);
    }
    // Perform reduction of an ascending sequence of {1,2,3,...,n}
    // The expected result is n*(n+1)/2
    for (auto n: {1, 2, 7, 31, 32, 33, 63, 64, 65, 128}) {
        std::vector<double> in(n);
        for (int i=0; i<n; ++i) in[i] = i+1;
        std::vector<int> index(n, 0);

        auto out = reduce(in, 1, index);
        EXPECT_EQ(out[0], double((n+1)*n/2));
    }
}

TEST(reduce_by_key, scatter)
{
    // A monotonic sequence of keys with repetitions and gaps, for a reduction
    // onto an array of length 12.
    std::size_t n = 12;
    std::vector<int> index = {0,0,0,1,2,2,2,2,3,3,7,7,7,7,7,11};
    std::vector<double> in(index.size(), 1);
    std::vector<double> expected = {3., 1., 4., 2., 0., 0., 0., 5., 0., 0., 0., 1.};

    EXPECT_EQ(n, expected.size());

    auto out = reduce(in, n, index);
    EXPECT_EQ(expected, out);

    // rerun with 7 threads per thread block, to test
    //  * using more than one thread block
    //  * thread blocks that are not a multiple of 32
    //  * thread blocks that are less than 32

    out = reduce(in, n, index, 7);
    EXPECT_EQ(expected, out);
}

// Test kernels that perform more than one reduction in a single invokation.
// Used to reproduce and test for synchronization issues on V100 GPUs.

template <typename T, typename I>
__global__
void reduce_twice_kernel(const T* src, T* dst, const I* index, int n) {
    unsigned tid = threadIdx.x + blockIdx.x*blockDim.x;

    unsigned mask = gpu::ballot(0xffffffff, tid<n);
    if (tid<n) {
        gpu::reduce_by_key(src[tid], dst, index[tid], mask);
        gpu::reduce_by_key(src[tid], dst, index[tid], mask);
    }
}

template <typename T>
std::vector<T> reduce_twice(const std::vector<T>& in, size_t n_out, const std::vector<int>& index, unsigned block_dim=128) {
    EXPECT_EQ(in.size(), index.size());
    EXPECT_TRUE(std::is_sorted(index.begin(), index.end()));

    using array  = gpu_vector<T>;
    using iarray = gpu_vector<int>;

    int n = in.size();

    array  src(in);
    iarray idx(index);
    array  dst(std::vector<T>(n_out, 0));

    unsigned grid_dim = (n-1)/block_dim + 1;
    reduce_twice_kernel<<<grid_dim, block_dim>>>(src.data(), dst.data(), idx.data(), n);

    return dst.host_vector();
}

TEST(reduce_by_key, scatter_twice)
{
    // A monotonic sequence of keys with repetitions and gaps, for a reduction
    // onto an array of length 12.
    std::size_t n = 12;
    std::vector<int> index = {0,0,0,1,2,2,3,7,7,7,11};
    std::vector<double> in(index.size(), 1);
    std::vector<double> expected = {6., 2., 4., 2., 0., 0., 0., 6., 0., 0., 0., 2.};

    EXPECT_EQ(n, expected.size());

    auto out = reduce_twice(in, n, index);
    EXPECT_EQ(expected, out);

    // rerun with 7 threads per thread block, to test
    //  * using more than one thread block
    //  * thread blocks that are not a multiple of 32
    //  * thread blocks that are less than 32

    out = reduce_twice(in, n, index, 7);
    EXPECT_EQ(expected, out);
}
