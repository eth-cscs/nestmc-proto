#include <vector>
#include <string>

#include "../test/gtest.h"

#include <arbor/mechanism_abi.h>
#include <arbor/mechanism.hpp>

#include "backends/multicore/shared_state.hpp"

using namespace std::string_literals;

TEST(abi, multicore_initialisation) {
    std::vector<arb_field_info> globals = {{ "G0", "kg",  123.0,     0.0, 2000.0},
                                           { "G1", "lb",  456.0,     0.0, 2000.0},
                                           { "G2", "gr",  789.0,     0.0, 2000.0}};
    std::vector<arb_field_info> states  = {{ "S0", "nA",      0.123, 0.0, 2000.0},
                                           { "S1", "mV",      0.456, 0.0, 2000.0}};
    std::vector<arb_field_info> params  = {{ "P0", "lm", -123.0,     0.0, 2000.0}};

    arb_mechanism_type type{};
    type.globals    = globals.data(); type.n_globals    = globals.size();
    type.parameters = params.data();  type.n_parameters = params.size();
    type.state_vars = states.data();  type.n_state_vars = states.size();

    arb_mechanism_interface iface { arb_backend_kind_cpu,
                                    1,
                                    1,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    nullptr };

    auto mech = arb::mechanism(type, iface);

    arb_size_type ncell = 1;
    arb_size_type ncv = 1;
    std::vector<arb_index_type> cv_to_intdom(ncv, 0);
    std::vector<arb_value_type> temp(ncv, 23);
    std::vector<arb_value_type> diam(ncv, 1.);
    std::vector<arb_value_type> vinit(ncv, -65);
    std::vector<arb::fvm_gap_junction> gj = {};
    std::vector<arb_index_type> src_to_spike = {};

    arb::multicore::shared_state shared_state(ncell, ncell, 0,
                                              cv_to_intdom, cv_to_intdom,
                                              gj, vinit, temp, diam, src_to_spike,
                                              mech.data_alignment());

    arb::mechanism_layout layout;
    layout.weight.assign(ncv, 1.);
    for (arb_size_type i = 0; i<ncv; ++i) layout.cv.push_back(i);

    arb::mechanism_overrides overrides;

    shared_state.instantiate(mech, 42, overrides, layout);

    {
        auto tb = mech.global_table();
        EXPECT_EQ(tb.size(), globals.size());
        for (auto idx = 0ul; idx < globals.size(); ++idx) {
            const auto [k, v] = tb[idx];
            EXPECT_EQ(k, globals[idx].name);
            EXPECT_EQ(v, globals[idx].default_value);
            EXPECT_EQ(mech.field_data(globals[idx].name), nullptr);
        }
    }

    {
        auto tb = mech.state_table();
        EXPECT_EQ(tb.size(), states.size());
        for (auto idx = 0ul; idx < states.size(); ++idx) {
            const auto& [k, v]  = tb[idx];
            const auto& [vs, d] = v;
            EXPECT_EQ(k, states[idx].name);
            EXPECT_EQ(d, states[idx].default_value);
            for (auto cv = 0ul; cv < ncv; ++cv) EXPECT_EQ(d, vs[cv]);
            EXPECT_EQ(mech.field_data(states[idx].name), vs);
        }
    }

    {
        for (auto idx = 0ul; idx < params.size(); ++idx) {
            const auto& vs = mech.field_data(params[idx].name);
            for (auto cv = 0ul; cv < ncv; ++cv) EXPECT_EQ(vs[cv], params[idx].default_value);
        }
    }
}


#ifdef ARB_GPU_ENABLED
TEST(abi, gpu_initialisation) {
    std::vector<arb_field_info> globals = {{ "G0", "kg",  123.0,     0.0, 2000.0},
                                           { "G1", "lb",  456.0,     0.0, 2000.0},
                                           { "G2", "gr",  789.0,     0.0, 2000.0}};
    std::vector<arb_field_info> states  = {{ "S0", "nA",      0.123, 0.0, 2000.0},
                                           { "S1", "mV",      0.456, 0.0, 2000.0}};
    std::vector<arb_field_info> params  = {{ "P0", "lm", -123.0,     0.0, 2000.0}};

    arb_mechanism_type type{};
    type.globals    = globals.data(); type.n_globals    = globals.size();
    type.parameters = params.data();  type.n_parameters = params.size();
    type.state_vars = states.data();  type.n_state_vars = states.size();

    arb_mechanism_interface iface { arb_backend_kind_gpu,
                                    1,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    nullptr };

    auto mech = arb::mechanism(type, iface);

    arb_size_type ncell = 1;
    arb_size_type ncv = 1;
    std::vector<arb_index_type> cv_to_intdom(ncv, 0);
    std::vector<arb_value_type> temp(ncv, 23);
    std::vector<arb_value_type> diam(ncv, 1.);
    std::vector<arb_value_type> vinit(ncv, -65);
    std::vector<arb::fvm_gap_junction> gj = {};
    std::vector<arb_index_type> src_to_spike = {};

    arb::multicore::shared_state shared_state(ncell, ncell, 0,
                                              cv_to_intdom, cv_to_intdom,
                                              gj, vinit, temp, diam, src_to_spike,
                                              mech.data_alignment());

    arb::mechanism_layout layout;
    layout.weight.assign(ncv, 1.);
    for (arb_size_type i = 0; i<ncv; ++i) layout.cv.push_back(i);

    arb::mechanism_overrides overrides;

    mech.instantiate(42, shared_state, overrides, layout);

    {
        auto tb = mech.global_table();
        EXPECT_EQ(tb.size(), globals.size());
        for (auto idx = 0ul; idx < globals.size(); ++idx) {
            const auto [k, v] = tb[idx];
            EXPECT_EQ(k, globals[idx].name);
            EXPECT_EQ(v, globals[idx].default_value);
            EXPECT_EQ(mech.field_data(globals[idx].name), nullptr);
        }
    }

    {
        auto tb = mech.state_table();
        EXPECT_EQ(tb.size(), states.size());
        for (auto idx = 0ul; idx < states.size(); ++idx) {
            const auto& [k, v]  = tb[idx];
            const auto& [vs, d] = v;
            EXPECT_EQ(k, states[idx].name);
            EXPECT_EQ(d, states[idx].default_value);
            for (auto cv = 0; cv < ncv; ++cv) EXPECT_EQ(d, vs[cv]);
            EXPECT_EQ(mech.field_data(states[idx].name), vs);
        }
    }

    {
        for (auto idx = 0ul; idx < params.size(); ++idx) {
            const auto& vs = mech.field_data(params[idx].name);
            for (auto cv = 0ul; cv < ncv; ++cv) EXPECT_EQ(vs[cv], params[idx].default_value);
        }
    }
}
#endif
