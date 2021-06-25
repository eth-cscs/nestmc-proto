#include <numeric>
#include <random>
#include <vector>

#ifdef ARB_HIP
#include <hip/hip_runtime.h>
#endif

#ifdef ARB_CUDA
#include <cuda.h>
#endif

#include <arbor/math.hpp>
#include <arbor/gpu/gpu_common.hpp>

#include "matrix.hpp"
#include "memory/memory.hpp"
#include "util/span.hpp"

#include "backends/gpu/matrix_state_flat.hpp"
#include "backends/gpu/matrix_state_fine.hpp"

#include "../gtest.h"
#include "common.hpp"


using namespace arb;

using util::make_span;
using util::assign_from;
using memory::on_gpu;
using memory::on_host;

using testing::seq_almost_eq;

using std::begin;
using std::end;

// test that the flat and interleaved storage back ends produce identical results
TEST(matrix, backends)
{
    using T = fvm_value_type;
    using I = fvm_index_type;

    using state_flat = gpu::matrix_state_flat<T, I>;
    using state_fine = gpu::matrix_state_fine<T, I>;

    using gpu_array  = memory::device_vector<T>;

    // There are two matrix structures:
    //
    // p_1: 3 branches, 6 compartments
    //
    //           3
    //          /.
    // 0 - 1 - 2
    //          \.
    //           4
    //            \.
    //             5
    //
    // p_2: 5 branches, 8 compartments
    //
    //             4
    //            /.
    //           3
    //          / \.
    // 0 - 1 - 2   5
    //          \.
    //           6
    //            \.
    //             7

    // The parent indexes that define the two matrix structures
    std::vector<std::vector<I>>
        p_base = { {0,0,1,2,2,4}, {0,0,1,2,3,3,2,6} };

    // Make a set of matrices based on repeating this pattern.
    // We assign the patterns round-robin, i.e. so that the input
    // matrices will have alternating sizes of 6 and 8, which will
    // test the solver with variable matrix size, and exercise
    // solvers that reorder matrices according to size.
    const int num_mtx = 200;

    std::vector<I> p;
    std::vector<I> cell_cv_divs;
    std::vector<I> cell_to_intdom;
    for (auto m=0; m<num_mtx; ++m) {
        auto &p_ref = p_base[m%2];
        auto first = p.size();
        for (auto i: p_ref) {
            p.push_back(i + first);
        }
        cell_cv_divs.push_back(first);
        cell_to_intdom.push_back(m);
    }
    cell_cv_divs.push_back(p.size());

    auto group_size = cell_cv_divs.back();

    // Build the capacitance, (axial) conductance, voltage, current density,
    // and membrane conductance vectors. Populate them with nonzero random values.
    auto gen  = std::mt19937();
    gen.seed(100);
    auto dist = std::uniform_real_distribution<T>(1, 200);

    std::vector<T> Cm(group_size);
    std::vector<T> g(group_size);
    std::vector<T> v(group_size);
    std::vector<T> i(group_size);
    std::vector<T> mg(group_size);
    std::vector<T> area(group_size, 1e3);

    std::generate(Cm.begin(), Cm.end(), [&](){return dist(gen);});
    std::generate(g.begin(), g.end(), [&](){return dist(gen);});
    std::generate(v.begin(), v.end(), [&](){return dist(gen);});
    std::generate(i.begin(), i.end(), [&](){return dist(gen);});
    std::generate(mg.begin(), mg.end(), [&](){return dist(gen);});

    // Make the reference matrix and the gpu matrix
    auto flat = state_flat(p, cell_cv_divs, Cm, g, area, cell_to_intdom); // flat
    auto fine = state_fine(p, cell_cv_divs, Cm, g, area, cell_to_intdom); // interleaved

    // Set the integration times for the cells to be between 0.01 and 0.02 ms.
    std::vector<T> dt(num_mtx, 0);

    auto dt_dist = std::uniform_real_distribution<T>(0.01, 0.02);
    std::generate(dt.begin(), dt.end(), [&](){return dt_dist(gen);});

    // Voltage, current, and membrane conductance values.
    auto gpu_dt = on_gpu(dt);
    auto gpu_v = on_gpu(v);
    auto gpu_i = on_gpu(i);
    auto gpu_mg = on_gpu(mg);

    auto x_flat_d = gpu_array(group_size);
    auto x_fine_d = gpu_array(group_size);

    flat.assemble(gpu_dt, gpu_v, gpu_i, gpu_mg);
    fine.assemble(gpu_dt, gpu_v, gpu_i, gpu_mg);

    flat.solve(x_flat_d);
    fine.solve(x_fine_d);

    // Compare the results.
    // We expect exact equality for the two gpu matrix implementations because both
    // perform the same operations in the same order on the same inputs.
    auto x_flat = on_host(x_flat_d);
    // as the fine algorithm contains atomics the solution might be slightly
    // different from flat and interleaved
    auto x_fine = on_host(x_fine_d);

    auto max_diff_fine =
        util::max_value(
            util::transform_view(
                util::count_along(x_flat),
                [&](unsigned i) {return std::abs(x_flat[i] - x_fine[i]);}));

    EXPECT_LE(max_diff_fine, 1e-12);
}
