#include "gtest.h"
#include "../test_util.hpp"

#include <cell.hpp>
#include <mechanisms/expsyn.hpp>
#include <mechanisms/exp2syn.hpp>
#include <backends/memory_traits.hpp>

// compares results with those generated by nrn/ball_and_stick.py
TEST(synapses, add_to_cell)
{
    using namespace nest::mc;

    nest::mc::cell cell;

    // setup global state for the mechanisms
    // nest::mc::mechanisms::setup_mechanism_helpers();

    // Soma with diameter 12.6157 um and HH channel
    auto soma = cell.add_soma(12.6157/2.0);
    soma->add_mechanism(hh_parameters());

    parameter_list exp_default("expsyn");
    parameter_list exp2_default("exp2syn");

    cell.add_synapse({0, 0.1}, exp_default);
    cell.add_synapse({1, 0.2}, exp2_default);
    cell.add_synapse({0, 0.3}, exp_default);

    EXPECT_EQ(3u, cell.synapses().size());
    const auto& syns = cell.synapses();

    EXPECT_EQ(syns[0].location.segment, 0u);
    EXPECT_EQ(syns[0].location.position, 0.1);
    EXPECT_EQ(syns[0].mechanism.name(), "expsyn");

    EXPECT_EQ(syns[1].location.segment, 1u);
    EXPECT_EQ(syns[1].location.position, 0.2);
    EXPECT_EQ(syns[1].mechanism.name(), "exp2syn");

    EXPECT_EQ(syns[2].location.segment, 0u);
    EXPECT_EQ(syns[2].location.position, 0.3);
    EXPECT_EQ(syns[2].mechanism.name(), "expsyn");
}

// compares results with those generated by nrn/ball_and_stick.py
TEST(synapses, expsyn_basic_state)
{
    using namespace nest::mc;

    using synapse_type = mechanisms::expsyn::mechanism_expsyn<multicore::memory_traits>;
    auto num_syn = 4;

    synapse_type::index_type indexes(num_syn);
    synapse_type::vector_type voltage(num_syn, -65.0);
    synapse_type::vector_type current(num_syn,   1.0);
    auto mech = mechanisms::make_mechanism<synapse_type>( voltage, current, indexes );

    auto ptr = dynamic_cast<synapse_type*>(mech.get());

    auto n = ptr->size();
    using view = synapse_type::view;

    // parameters initialized to default values
    for(auto e : view(ptr->e, n)) {
        EXPECT_EQ(e, 0.);
    }
    for(auto tau : view(ptr->tau, n)) {
        EXPECT_EQ(tau, 2.0);
    }

    // current and voltage vectors correctly hooked up
    for(auto v : view(ptr->vec_v_, n)) {
        EXPECT_EQ(v, -65.);
    }
    for(auto i : view(ptr->vec_i_, n)) {
        EXPECT_EQ(i, 1.0);
    }

    // should be initialized to NaN
    for(auto g : view(ptr->g, n)) {
        EXPECT_NE(g, g);
    }

    // initialize state then check g has been set to zero
    ptr->nrn_init();
    for(auto g : view(ptr->g, n)) {
        EXPECT_EQ(g, 0.);
    }

    // call net_receive on two of the synapses
    ptr->net_receive(1, 3.14);
    ptr->net_receive(3, 1.04);
    EXPECT_EQ(ptr->g[1], 3.14);
    EXPECT_EQ(ptr->g[3], 1.04);
}

TEST(synapses, exp2syn_basic_state)
{
    using namespace nest::mc;

    using synapse_type = mechanisms::exp2syn::mechanism_exp2syn<multicore::memory_traits>;
    auto num_syn = 4;

    synapse_type::index_type indexes(num_syn);
    synapse_type::vector_type voltage(num_syn, -65.0);
    synapse_type::vector_type current(num_syn,   1.0);
    auto mech = mechanisms::make_mechanism<synapse_type>( voltage, current, indexes );

    auto ptr = dynamic_cast<synapse_type*>(mech.get());

    auto n = ptr->size();
    using view = synapse_type::view;

    // parameters initialized to default values
    for(auto e : view(ptr->e, n)) {
        EXPECT_EQ(e, 0.);
    }
    for(auto tau1: view(ptr->tau1, n)) {
        EXPECT_EQ(tau1, 0.5);
    }
    for(auto tau2: view(ptr->tau2, n)) {
        EXPECT_EQ(tau2, 2.0);
    }

    // should be initialized to NaN
    for(auto factor: view(ptr->factor, n)) {
        EXPECT_NE(factor, factor);
    }

    // initialize state then check factor has sane (positive) value
    // and A and B are zero
    ptr->nrn_init();
    for(auto factor: view(ptr->factor, n)) {
        EXPECT_GT(factor, 0.);
    }
    for(auto A: view(ptr->A, n)) {
        EXPECT_EQ(A, 0.);
    }
    for(auto B: view(ptr->B, n)) {
        EXPECT_EQ(B, 0.);
    }

    // call net_receive on two of the synapses
    ptr->net_receive(1, 3.14);
    ptr->net_receive(3, 1.04);

    EXPECT_NEAR(ptr->A[1], ptr->factor[1]*3.14, 1e-6);
    EXPECT_NEAR(ptr->B[3], ptr->factor[3]*1.04, 1e-6);
}

