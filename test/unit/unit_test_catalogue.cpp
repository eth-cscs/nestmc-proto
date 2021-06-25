#include <arbor/mechcat.hpp>
#include <arbor/version.hpp>

#ifdef ARB_GPU_ENABLED
#include "backends/gpu/fvm.hpp"
#endif
#include "backends/multicore/fvm.hpp"

#include "unit_test_catalogue.hpp"
#include "mechanisms/ca_linear.hpp"
#include "mechanisms/celsius_test.hpp"
#include "mechanisms/diam_test.hpp"
#include "mechanisms/non_linear.hpp"
#include "mechanisms/param_as_state.hpp"
#include "mechanisms/post_events_syn.hpp"
#include "mechanisms/test0_kin_diff.hpp"
#include "mechanisms/test_linear_state.hpp"
#include "mechanisms/test_linear_init.hpp"
#include "mechanisms/test_linear_init_shuffle.hpp"
#include "mechanisms/test0_kin_conserve.hpp"
#include "mechanisms/test0_kin_steadystate.hpp"
#include "mechanisms/test0_kin_compartment.hpp"
#include "mechanisms/test1_kin_compartment.hpp"
#include "mechanisms/test1_kin_diff.hpp"
#include "mechanisms/test1_kin_conserve.hpp"
#include "mechanisms/test2_kin_diff.hpp"
#include "mechanisms/test3_kin_diff.hpp"
#include "mechanisms/test4_kin_compartment.hpp"
#include "mechanisms/test1_kin_steadystate.hpp"
#include "mechanisms/fixed_ica_current.hpp"
#include "mechanisms/point_ica_current.hpp"
#include "mechanisms/linear_ca_conc.hpp"
#include "mechanisms/test_cl_valence.hpp"
#include "mechanisms/test_ca_read_valence.hpp"
#include "mechanisms/read_eX.hpp"
#include "mechanisms/write_Xi_Xo.hpp"
#include "mechanisms/write_multiple_eX.hpp"
#include "mechanisms/write_eX.hpp"
#include "mechanisms/read_cai_init.hpp"
#include "mechanisms/write_cai_breakpoint.hpp"
#include "mechanisms/test_ca.hpp"
#include "mechanisms/test_kin1.hpp"
#include "mechanisms/test_kinlva.hpp"

#include "../gtest.h"

#ifndef ARB_GPU_ENABLED
#define ADD_MECH(c, x)\
c.add(#x, make_testing_##x());   \
c.register_implementation(#x, std::make_unique<arb::mechanism>(make_testing_##x(), *make_testing_##x##_interface_multicore()));
#else
#define ADD_MECH(c, x)\
c.add(#x, make_testing_##x());                                      \
c.register_implementation(#x, std::make_unique<arb::mechanism>(make_testing_##x(), *make_testing_##x##_interface_multicore())); \
c.register_implementation(#x, std::make_unique<arb::mechanism>(make_testing_##x(), *make_testing_##x##_interface_gpu()));
#endif

using namespace arb;

mechanism_catalogue make_unit_test_catalogue(const mechanism_catalogue& from) {
    mechanism_catalogue cat(from);

    ADD_MECH(cat, test_ca)
    ADD_MECH(cat, test_kin1)
    ADD_MECH(cat, test_kinlva)
    ADD_MECH(cat, ca_linear)
    ADD_MECH(cat, celsius_test)
    ADD_MECH(cat, diam_test)
    ADD_MECH(cat, param_as_state)
    ADD_MECH(cat, post_events_syn)
    ADD_MECH(cat, test_linear_state)
    ADD_MECH(cat, test_linear_init)
    ADD_MECH(cat, test_linear_init_shuffle)
    ADD_MECH(cat, test0_kin_diff)
    ADD_MECH(cat, test0_kin_conserve)
    ADD_MECH(cat, test0_kin_steadystate)
    ADD_MECH(cat, test0_kin_compartment)
    ADD_MECH(cat, test1_kin_diff)
    ADD_MECH(cat, test1_kin_conserve)
    ADD_MECH(cat, test2_kin_diff)
    ADD_MECH(cat, test3_kin_diff)
    ADD_MECH(cat, test1_kin_steadystate)
    ADD_MECH(cat, test1_kin_compartment)
    ADD_MECH(cat, test4_kin_compartment)
    ADD_MECH(cat, fixed_ica_current)
    ADD_MECH(cat, non_linear)
    ADD_MECH(cat, point_ica_current)
    ADD_MECH(cat, linear_ca_conc)
    ADD_MECH(cat, test_cl_valence)
    ADD_MECH(cat, test_ca_read_valence)
    ADD_MECH(cat, read_eX)
    ADD_MECH(cat, write_Xi_Xo)
    ADD_MECH(cat, write_multiple_eX)
    ADD_MECH(cat, write_eX)
    ADD_MECH(cat, read_cai_init)
    ADD_MECH(cat, write_cai_breakpoint)

    return cat;
}

