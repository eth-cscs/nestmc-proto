#include <vector>

#include <arbor/mechanism.hpp>
#include <arbor/version.hpp>

#include "backends/multicore/fvm.hpp"
#ifdef ARB_GPU_ENABLED
#include "backends/gpu/fvm.hpp"
#endif

#include "common.hpp"
#include "mech_private_field_access.hpp"
#include "unit_test_catalogue.hpp"

using namespace arb;

template <typename backend>
void run_celsius_test() {
    auto cat = make_unit_test_catalogue();

    // one cell, three CVs:

    fvm_size_type ncell = 1;
    fvm_size_type ncv = 3;
    std::vector<fvm_index_type> cv_to_intdom(ncv, 0);

    std::vector<fvm_gap_junction> gj = {};
    auto instance = cat.instance(backend::kind, "celsius_test");
    auto& celsius_test = instance.mech;

    double temperature_K = 300.;
    double temperature_C = temperature_K-273.15;

    std::vector<fvm_value_type> temp(ncv, temperature_K);
    std::vector<fvm_value_type> diam(ncv, 1.);
    std::vector<fvm_value_type> vinit(ncv, -65);
    std::vector<fvm_index_type> src_to_spike = {};

    auto shared_state = std::make_unique<typename backend::shared_state>(
        ncell, ncell, 0, cv_to_intdom, cv_to_intdom, gj, vinit, temp, diam, src_to_spike, celsius_test->data_alignment());

    mechanism_layout layout;
    mechanism_overrides overrides;

    layout.weight.assign(ncv, 1.);
    for (fvm_size_type i = 0; i<ncv; ++i) {
        layout.cv.push_back(i);
    }

    shared_state->instantiate(*celsius_test, 0, overrides, layout);
    shared_state->reset();

    // expect 0 value in state 'c' after init:

    celsius_test->initialize();
    std::vector<fvm_value_type> expected_c_values(ncv, 0.);

    EXPECT_EQ(expected_c_values, mechanism_field(celsius_test.get(), "c"));

    // expect temperature_C value in state 'c' after state update:

    celsius_test->update_state();
    expected_c_values.assign(ncv, temperature_C);

    EXPECT_EQ(expected_c_values, mechanism_field(celsius_test.get(), "c"));
}

template <typename backend>
void run_diam_test() {
    auto cat = make_unit_test_catalogue();

    // one cell, three CVs:

    fvm_size_type ncell = 1;
    fvm_size_type ncv = 3;
    std::vector<fvm_index_type> cv_to_intdom(ncv, 0);

    std::vector<fvm_gap_junction> gj = {};
    auto instance = cat.instance(backend::kind, "diam_test");
    auto& celsius_test = instance.mech;

    std::vector<fvm_value_type> temp(ncv, 300.);
    std::vector<fvm_value_type> vinit(ncv, -65);
    std::vector<fvm_value_type> diam(ncv);
    std::vector<fvm_index_type> src_to_spike = {};

    mechanism_layout layout;
    mechanism_overrides overrides;

    layout.weight.assign(ncv, 1.);

    for (fvm_size_type i = 0; i < ncv; ++i) {
        diam[i] = i*2 + 0.1;
        layout.cv.push_back(i);
    }

    auto shared_state = std::make_unique<typename backend::shared_state>(
            ncell, ncell, 0, cv_to_intdom, cv_to_intdom, gj, vinit, temp, diam, src_to_spike, celsius_test->data_alignment());


    shared_state->instantiate(*celsius_test, 0, overrides, layout);
    shared_state->reset();

    // expect 0 value in state 'd' after init:

    celsius_test->initialize();
    std::vector<fvm_value_type> expected_d_values(ncv, 0.);

    EXPECT_EQ(expected_d_values, mechanism_field(celsius_test.get(), "d"));

    // expect original diam values in state 'd' after state update:

    celsius_test->update_state();
    expected_d_values = diam;

    EXPECT_EQ(expected_d_values, mechanism_field(celsius_test.get(), "d"));
}

TEST(mech_temperature, celsius) {
    run_celsius_test<multicore::backend>();
    run_diam_test<multicore::backend>();
}

#ifdef ARB_GPU_ENABLED
TEST(mech_temperature_gpu, celsius) {
    run_celsius_test<gpu::backend>();
    run_diam_test<gpu::backend>();
}
#endif
