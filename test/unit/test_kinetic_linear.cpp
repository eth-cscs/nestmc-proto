#include <vector>

#include <arbor/mechanism.hpp>
#include <arbor/version.hpp>

#include "backends/multicore/fvm.hpp"

#ifdef ARB_GPU_ENABLED
#include "backends/gpu/fvm.hpp"
#endif

#include "common.hpp"
#include "mech_private_field_access.hpp"
#include "fvm_lowered_cell.hpp"
#include "fvm_lowered_cell_impl.hpp"
#include "simple_recipes.hpp"
#include "unit_test_catalogue.hpp"

using namespace arb;

using backend = arb::multicore::backend;
using fvm_cell = arb::fvm_lowered_cell_impl<backend>;

using shared_state = backend::shared_state;
ACCESS_BIND(std::unique_ptr<shared_state> fvm_cell::*, private_state_ptr, &fvm_cell::state_)

template <typename backend>
void run_test(std::string mech_name,
        std::vector<std::string> state_variables,
        std::unordered_map<std::string, fvm_value_type> assigned_variables,
        std::vector<fvm_value_type> t0_values,
        std::vector<fvm_value_type> t1_values,
        fvm_value_type dt) {

    auto cat = make_unit_test_catalogue();

    fvm_size_type ncell = 1;
    fvm_size_type ncv = 1;
    std::vector<fvm_index_type> cv_to_intdom(ncv, 0);

    std::vector<fvm_gap_junction> gj = {};
    auto instance = cat.instance(backend::kind, mech_name);
    auto& test = instance.mech;

    std::vector<fvm_value_type> temp(ncv, 300.);
    std::vector<fvm_value_type> diam(ncv, 1.);
    std::vector<fvm_value_type> vinit(ncv, -65);
    std::vector<fvm_index_type> src_to_spike = {};

    auto shared_state = std::make_unique<typename backend::shared_state>(
            ncell, ncell, 0, cv_to_intdom, cv_to_intdom, gj, vinit, temp, diam, src_to_spike, test->data_alignment());

    mechanism_layout layout;
    mechanism_overrides overrides;

    layout.weight.assign(ncv, 1.);
    for (fvm_size_type i = 0; i<ncv; ++i) {
        layout.cv.push_back(i);
    }

    shared_state->instantiate(*test, 0, overrides, layout);

    for (auto a: assigned_variables) {
        shared_state->set_parameter(*test, a.first, std::vector<fvm_value_type>(ncv,a.second));
    }

    shared_state->reset();

    test->initialize();

    if (!t0_values.empty()) {
        for (unsigned i = 0; i < state_variables.size(); i++) {
            for (unsigned j = 0; j < ncv; j++) {
                EXPECT_NEAR(t0_values[i], mechanism_field(test.get(), state_variables[i]).at(j), 1e-6);
            }
        }
    }

    shared_state->update_time_to(dt, dt);
    shared_state->set_dt();

    test->update_state();

    if (!t1_values.empty()) {
        for (unsigned i = 0; i < state_variables.size(); i++) {
            for (unsigned j = 0; j < ncv; j++) {
                EXPECT_NEAR(t1_values[i], mechanism_field(test.get(), state_variables[i]).at(j), 1e-6);
            }
        }
    }
}

TEST(mech_kinetic, kinetic_linear_scaled) {
    std::vector<std::string> state_variables = {"s", "h", "d"};
    std::vector<fvm_value_type> t0_values = {0.5, 0.2, 0.3};
    std::vector<fvm_value_type> t1_0_values = {0.373297, 0.591621, 0.0350817};
    std::vector<fvm_value_type> t1_1_values = {0.329897, 0.537371, 0.132732};

    run_test<multicore::backend>("test0_kin_compartment", state_variables, {}, t0_values, t1_0_values, 0.5);
    run_test<multicore::backend>("test1_kin_compartment", state_variables, {}, t0_values, t1_1_values, 0.5);

}

TEST(mech_kinetic, kinetic_linear_1_conserve) {
    std::vector<std::string> state_variables = {"s", "h", "d"};
    std::vector<fvm_value_type> t0_values = {0.5, 0.2, 0.3};
    std::vector<fvm_value_type> t1_0_values = {0.380338, 0.446414, 0.173247};
    std::vector<fvm_value_type> t1_1_values = {0.218978, 0.729927, 0.0510949};

    run_test<multicore::backend>("test0_kin_diff", state_variables, {}, t0_values, t1_0_values, 0.5);
    run_test<multicore::backend>("test0_kin_conserve", state_variables, {}, t0_values, t1_0_values, 0.5);
    run_test<multicore::backend>("test0_kin_steadystate", state_variables, {}, t0_values, t1_1_values, 0.5);
}

TEST(mech_kinetic, kinetic_linear_2_conserve) {
    std::vector<std::string> state_variables = {"a", "b", "x", "y"};
    std::vector<fvm_value_type> t0_values = {0.2, 0.8, 0.6, 0.4};
    std::vector<fvm_value_type> t1_0_values = {0.217391304, 0.782608696, 0.33333333, 0.66666666};
    std::vector<fvm_value_type> t1_1_values = {0.230769, 0.769231, 0.189189, 0.810811};

    run_test<multicore::backend>("test1_kin_diff", state_variables, {}, t0_values, t1_0_values, 0.5);
    run_test<multicore::backend>("test1_kin_conserve", state_variables, {}, t0_values, t1_0_values, 0.5);
    run_test<multicore::backend>("test1_kin_steadystate", state_variables, {}, t0_values, t1_1_values, 0.5);
}

TEST(mech_kinetic, kinetic_nonlinear) {
    std::vector<std::string> state_variables = {"a", "b", "c"};
    std::vector<fvm_value_type> t0_values = {0.2, 0.3, 0.5};
    std::vector<fvm_value_type> t1_0_values = {0.222881, 0.31144, 0.48856};
    std::vector<fvm_value_type> t1_1_values = {0.2078873133, 0.34222075, 0.45777925};

    run_test<multicore::backend>("test2_kin_diff", state_variables, {}, t0_values, t1_0_values, 0.025);
    run_test<multicore::backend>("test3_kin_diff", state_variables, {}, t0_values, t1_1_values, 0.025);

}

TEST(mech_kinetic, kinetic_nonlinear_scaled) {
    std::vector<std::string> state_variables = {"A", "B", "C", "d", "e"};
    std::vector<fvm_value_type> t0_values = {4.5, 6.6, 0.28, 2, 0};
    std::vector<fvm_value_type> t1_values = {4.087281958014442,
                                             6.224088678118931,
                                             0.6559113218810689,
                                             1.8315624742412617,
                                             0.16843752575873824};

    run_test<multicore::backend>("test4_kin_compartment", state_variables, {}, t0_values, t1_values, 0.1);
}

TEST(mech_linear, linear_system) {
    std::vector<std::string> state_variables = {"h", "s", "d"};
    std::vector<fvm_value_type> values = {0.5, 0.2, 0.3};
    std::unordered_map<std::string, fvm_value_type> assigned_variables = {{"a0", 2.5}, {"a1",0.5}, {"a2",3}, {"a3",2.3}};

    run_test<multicore::backend>("test_linear_state", state_variables, assigned_variables, {}, values, 0.5);
    run_test<multicore::backend>("test_linear_init", state_variables, assigned_variables, values, {}, 0.5);
    run_test<multicore::backend>("test_linear_init_shuffle", state_variables, assigned_variables, values, {}, 0.5);
}

#ifdef ARB_GPU_ENABLED
TEST(mech_kinetic_gpu, kinetic_linear_scaled) {
    std::vector<std::string> state_variables = {"s", "h", "d"};
    std::vector<fvm_value_type> t0_values = {0.5, 0.2, 0.3};
    std::vector<fvm_value_type> t1_0_values = {0.373297, 0.591621, 0.0350817};
    std::vector<fvm_value_type> t1_1_values = {0.329897, 0.537371, 0.132732};

    run_test<gpu::backend>("test0_kin_compartment", state_variables, {}, t0_values, t1_0_values, 0.5);
    run_test<gpu::backend>("test1_kin_compartment", state_variables, {}, t0_values, t1_1_values, 0.5);
}

TEST(mech_kinetic_gpu, kinetic_linear_1_conserve) {
    std::vector<std::string> state_variables = {"s", "h", "d"};
    std::vector<fvm_value_type> t0_values = {0.5, 0.2, 0.3};
    std::vector<fvm_value_type> t1_0_values = {0.380338, 0.446414, 0.173247};
    std::vector<fvm_value_type> t1_1_values = {0.218978, 0.729927, 0.0510949};

    run_test<gpu::backend>("test0_kin_diff", state_variables, {}, t0_values, t1_0_values, 0.5);
    run_test<gpu::backend>("test0_kin_conserve", state_variables, {}, t0_values, t1_0_values, 0.5);
    run_test<gpu::backend>("test0_kin_steadystate", state_variables, {}, t0_values, t1_1_values, 0.5);
}

TEST(mech_kinetic_gpu, kinetic_linear_2_conserve) {
    std::vector<std::string> state_variables = {"a", "b", "x", "y"};
    std::vector<fvm_value_type> t0_values = {0.2, 0.8, 0.6, 0.4};
    std::vector<fvm_value_type> t1_0_values = {0.217391304, 0.782608696, 0.33333333, 0.66666666};
    std::vector<fvm_value_type> t1_1_values = {0.230769, 0.769231, 0.189189, 0.810811};

    run_test<gpu::backend>("test1_kin_diff", state_variables, {}, t0_values, t1_0_values, 0.5);
    run_test<gpu::backend>("test1_kin_conserve", state_variables, {}, t0_values, t1_0_values, 0.5);
    run_test<gpu::backend>("test1_kin_steadystate", state_variables, {}, t0_values, t1_1_values, 0.5);
}

TEST(mech_kinetic_gpu, kinetic_nonlinear) {
    std::vector<std::string> state_variables = {"a", "b", "c"};
    std::vector<fvm_value_type> t0_values = {0.2, 0.3, 0.5};
    std::vector<fvm_value_type> t1_0_values = {0.222881, 0.31144, 0.48856};
    std::vector<fvm_value_type> t1_1_values = {0.2078873133, 0.34222075, 0.45777925};

    run_test<gpu::backend>("test2_kin_diff", state_variables, {}, t0_values, t1_0_values, 0.025);
    run_test<gpu::backend>("test3_kin_diff", state_variables, {}, t0_values, t1_1_values, 0.025);
}

TEST(mech_kinetic_gpu, kinetic_nonlinear_scaled) {
    std::vector<std::string> state_variables = {"A", "B", "C", "d", "e"};
    std::vector<fvm_value_type> t0_values = {4.5, 6.6, 0.28, 2, 0};
    std::vector<fvm_value_type> t1_values = {4.087281958014442,
                                             6.224088678118931,
                                             0.6559113218810689,
                                             1.8315624742412617,
                                             0.16843752575873824};

    run_test<gpu::backend>("test4_kin_compartment", state_variables, {}, t0_values, t1_values, 0.1);
}

TEST(mech_linear_gpu, linear_system) {
    std::vector<std::string> state_variables = {"h", "s", "d"};
    std::vector<fvm_value_type> values = {0.5, 0.2, 0.3};
    std::unordered_map<std::string, fvm_value_type> assigned_variables = {{"a0", 2.5},{"a1",0.5},{"a2",3},{"a3",2.3}};

    run_test<gpu::backend>("test_linear_state", state_variables, assigned_variables, {}, values, 0.5);
    run_test<gpu::backend>("test_linear_init", state_variables, assigned_variables, values, {}, 0.5);
    run_test<gpu::backend>("test_linear_init_shuffle", state_variables, assigned_variables, values, {}, 0.5);
}

#endif
