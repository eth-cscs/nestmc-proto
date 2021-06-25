#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "../gtest.h"

#include <arbor/common_types.hpp>
#include <arbor/domain_decomposition.hpp>
#include <arbor/fvm_types.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/math.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/recipe.hpp>
#include <arbor/sampling.hpp>
#include <arbor/simulation.hpp>
#include <arbor/schedule.hpp>
#include <arbor/mechanism.hpp>
#include <arbor/string_literals.hpp>
#include <arbor/util/any_ptr.hpp>

#include <arborenv/concurrency.hpp>

#include "backends/multicore/fvm.hpp"
#include "execution_context.hpp"
#include "fvm_lowered_cell.hpp"
#include "fvm_lowered_cell_impl.hpp"
#include "mech_private_field_access.hpp"
#include "util/meta.hpp"
#include "util/maputil.hpp"
#include "util/rangeutil.hpp"
#include "util/transform.hpp"

#include "common.hpp"
#include "unit_test_catalogue.hpp"
#include "../common_cells.hpp"
#include "../simple_recipes.hpp"

using namespace std::string_literals;

using backend = arb::multicore::backend;
using fvm_cell = arb::fvm_lowered_cell_impl<backend>;

// Access to fvm_cell private data:

using shared_state = backend::shared_state;
ACCESS_BIND(std::unique_ptr<shared_state> fvm_cell::*, private_state_ptr, &fvm_cell::state_)

using matrix = arb::matrix<arb::multicore::backend>;
ACCESS_BIND(matrix fvm_cell::*, private_matrix_ptr, &fvm_cell::matrix_)

ACCESS_BIND(std::vector<arb::mechanism_ptr> fvm_cell::*, private_mechanisms_ptr, &fvm_cell::mechanisms_)

arb::mechanism* find_mechanism(fvm_cell& fvcell, const std::string& name) {
    for (auto& mech: fvcell.*private_mechanisms_ptr) {
        if (mech->internal_name()==name) {
            return mech.get();
        }
    }
    return nullptr;
}

arb::mechanism* find_mechanism(fvm_cell& fvcell, int index) {
    auto& mechs = fvcell.*private_mechanisms_ptr;
    return index<(int)mechs.size()? mechs[index].get(): nullptr;
}

// Access to mechanism-internal data:
ACCESS_BIND(\
    arb::mechanism_global_table (arb::mechanism::*)(), \
    private_global_table_ptr,\
    &arb::mechanism::global_table)

ACCESS_BIND(\
    arb::mechanism_field_table (arb::mechanism::*)(), \
    private_field_table_ptr,\
    &arb::mechanism::field_table)

ACCESS_BIND(\
    arb::mechanism_ion_table (arb::mechanism::*)(), \
    private_ion_table_ptr,\
    &arb::mechanism::ion_table)

using namespace arb;

class gap_recipe_0: public recipe {
public:
    gap_recipe_0() {}

    cell_size_type num_cells() const override {
        return size_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        auto c = soma_cell_builder(20).make_cell();
        c.decorations.place(mlocation{0, 1}, gap_junction_site{}, "gj");
        return {cable_cell{c}};
    }

    cell_kind get_cell_kind(cell_gid_type gid) const override {
        return cell_kind::cable;
    }
    std::vector<gap_junction_connection> gap_junctions_on(cell_gid_type gid) const override {
        switch (gid) {
            case 0 :
                return {gap_junction_connection({5, "gj"}, {"gj"}, 0.1)};
            case 2 :
                return {
                        gap_junction_connection({3, "gj"}, {"gj"}, 0.1),
                };
            case 3 :
                return {
                        gap_junction_connection({7, "gj"}, {"gj"}, 0.1),
                        gap_junction_connection({2, "gj"}, {"gj"}, 0.1)
                };
            case 5 :
                return {gap_junction_connection({0, "gj"}, {"gj"}, 0.1)};
            case 7 :
                return {
                        gap_junction_connection({3, "gj"}, {"gj"}, 0.1),
                };
            default :
                return {};
        }
    }

private:
    cell_size_type size_ = 12;
};

class gap_recipe_1: public recipe {
public:
    gap_recipe_1() {}

    cell_size_type num_cells() const override {
        return size_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type) const override {
        return soma_cell_builder(20).make_cell();
    }

    cell_kind get_cell_kind(cell_gid_type gid) const override {
        return cell_kind::cable;
    }

private:
    cell_size_type size_ = 12;
};

class gap_recipe_2: public recipe {
public:
    gap_recipe_2() {}

    cell_size_type num_cells() const override {
        return size_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type) const override {
        auto c = soma_cell_builder(20).make_cell();
        c.decorations.place(mlocation{0,1}, gap_junction_site{}, "gj");
        return {cable_cell{c}};
    }

    cell_kind get_cell_kind(cell_gid_type gid) const override {
        return cell_kind::cable;
    }
    std::vector<gap_junction_connection> gap_junctions_on(cell_gid_type gid) const override {
        switch (gid) {
            case 0 :
                return {
                        gap_junction_connection({2, "gj"}, {"gj"}, 0.1),
                        gap_junction_connection({3, "gj"}, {"gj"}, 0.1),
                        gap_junction_connection({5, "gj"}, {"gj"}, 0.1)
                };
            case 2 :
                return {
                        gap_junction_connection({0, "gj"}, {"gj"}, 0.1),
                        gap_junction_connection({3, "gj"}, {"gj"}, 0.1),
                        gap_junction_connection({5, "gj"}, {"gj"}, 0.1)
                };
            case 3 :
                return {
                        gap_junction_connection({0, "gj"}, {"gj"}, 0.1),
                        gap_junction_connection({2, "gj"}, {"gj"}, 0.1),
                        gap_junction_connection({5, "gj"}, {"gj"}, 0.1)
                };
            case 5 :
                return {
                        gap_junction_connection({2, "gj"}, {"gj"}, 0.1),
                        gap_junction_connection({3, "gj"}, {"gj"}, 0.1),
                        gap_junction_connection({0, "gj"}, {"gj"}, 0.1)
                };
            default :
                return {};
        }
    }

private:
    cell_size_type size_ = 12;
};


TEST(fvm_lowered, matrix_init)
{
    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    auto isnan = [](auto v) { return std::isnan(v); };
    auto ispos = [](auto v) { return v>0; };
    auto isneg = [](auto v) { return v<0; };

    soma_cell_builder builder(12.6157/2.0);
    builder.add_branch(0, 200, 1.0/2, 1.0/2, 10, "dend"); // 10 compartments
    cable_cell cell = builder.make_cell();

    fvm_cell fvcell(context);
    fvcell.initialize({0}, cable1d_recipe(cell));

    auto& J = fvcell.*private_matrix_ptr;
    auto& S = fvcell.*private_state_ptr;
    EXPECT_EQ(J.size(), 12u);

    // Test that the matrix is initialized with sensible values

    fvcell.integrate(0.01, 0.01, {}, {});

    auto n = J.size();
    auto& mat = J.state_;

    EXPECT_FALSE(arb::util::any_of(util::subrange_view(mat.u, 1, n), isnan));
    EXPECT_FALSE(arb::util::any_of(mat.d, isnan));
    EXPECT_FALSE(arb::util::any_of(S->voltage, isnan));

    EXPECT_FALSE(arb::util::any_of(util::subrange_view(mat.u, 1, n), ispos));
    EXPECT_FALSE(arb::util::any_of(mat.d, isneg));
}

TEST(fvm_lowered, target_handles) {
    using namespace arb;

    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    cable_cell_description descriptions[] = {
        make_cell_ball_and_stick(),
        make_cell_ball_and_3stick()
    };

    // (in increasing target order)
    descriptions[0].decorations.place(mlocation{0, 0.7}, "expsyn", "syn0");
    descriptions[0].decorations.place(mlocation{0, 0.3}, "expsyn", "syn1");
    descriptions[1].decorations.place(mlocation{2, 0.2}, "exp2syn", "syn2");
    descriptions[1].decorations.place(mlocation{2, 0.8}, "expsyn", "syn3");

    descriptions[1].decorations.place(mlocation{0, 0}, threshold_detector{3.3}, "detector");

    cable_cell cells[] = {descriptions[0], descriptions[1]};

    EXPECT_EQ(cells[0].morphology().num_branches(), 1u);
    EXPECT_EQ(cells[1].morphology().num_branches(), 3u);

    auto test_target_handles = [&](fvm_cell& cell, const std::vector<target_handle>& targets) {
        mechanism* expsyn = find_mechanism(cell, "expsyn");
        ASSERT_TRUE(expsyn);
        mechanism* exp2syn = find_mechanism(cell, "exp2syn");
        ASSERT_TRUE(exp2syn);

        unsigned expsyn_id = expsyn->mechanism_id();
        unsigned exp2syn_id = exp2syn->mechanism_id();

        EXPECT_EQ(4u, targets.size());

        EXPECT_EQ(expsyn_id, targets[0].mech_id);
        EXPECT_EQ(1u, targets[0].mech_index);
        EXPECT_EQ(0u, targets[0].intdom_index);

        EXPECT_EQ(expsyn_id, targets[1].mech_id);
        EXPECT_EQ(0u, targets[1].mech_index);
        EXPECT_EQ(0u, targets[1].intdom_index);

        EXPECT_EQ(exp2syn_id, targets[2].mech_id);
        EXPECT_EQ(0u, targets[2].mech_index);
        EXPECT_EQ(1u, targets[2].intdom_index);

        EXPECT_EQ(expsyn_id, targets[3].mech_id);
        EXPECT_EQ(2u, targets[3].mech_index);
        EXPECT_EQ(1u, targets[3].intdom_index);
    };

    fvm_cell fvcell0(context);
    auto fvm_info0 = fvcell0.initialize({0, 1}, cable1d_recipe(cells, true));
    test_target_handles(fvcell0, fvm_info0.target_handles);

    fvm_cell fvcell1(context);
    auto fvm_info1 = fvcell1.initialize({0, 1}, cable1d_recipe(cells, false));
    test_target_handles(fvcell1, fvm_info1.target_handles);

}


TEST(fvm_lowered, stimulus) {
    // Ball-and-stick with two stimuli:
    //
    //           |stim0 |stim1
    // -----------------------
    // delay     |   5  |    1
    // duration  |  80  |    2
    // amplitude | 0.3  |  0.1
    // CV        |   5  |    0

    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    auto desc = make_cell_ball_and_stick(false);

    // At end of stick
    desc.decorations.place(mlocation{0,1},   i_clamp::box(5., 80., 0.3), "clamp0");
    // On the soma CV, which is over the approximate interval: (cable 0 0 0.1)
    desc.decorations.place(mlocation{0,0.05}, i_clamp::box(1., 2.,  0.1), "clamp1");

    std::vector<cable_cell> cells{desc};

    const fvm_size_type soma_cv = 0u;
    const fvm_size_type tip_cv = 5u;

    // The implementation of the stimulus is tested by creating a lowered cell, then
    // testing that the correct currents are injected at the correct control volumes
    // as during the stimulus windows.
    cable_cell_global_properties gprop;
    gprop.default_parameters = neuron_parameter_defaults;

    fvm_cv_discretization D = fvm_cv_discretize(cells, gprop.default_parameters, context);
    const auto& A = D.cv_area;

    fvm_cell fvcell(context);
    fvcell.initialize({0}, cable1d_recipe(cells));

    auto& state = *(fvcell.*private_state_ptr).get();
    auto& J = state.current_density;
    auto& T = state.time;

    // Test that no current is injected at t=0.
    memory::fill(J, 0.);
    memory::fill(T, 0.);
    state.add_stimulus_current();

    for (auto j: J) {
        EXPECT_EQ(j, 0.);
    }

    constexpr double reltol = 1e-10;
    using testing::near_relative;

    // Test that 0.1 nA current is injected at soma at t=1.
    memory::fill(J, 0.);
    memory::fill(T, 1.);
    state.add_stimulus_current();
    constexpr double unit_factor = 1e-3; // scale A/m²·µm² to nA
    EXPECT_TRUE(near_relative(-0.1, J[soma_cv]*A[soma_cv]*unit_factor, reltol));

    // Test that 0.1 nA is again injected at t=1.5, for a total of 0.2 nA.
    memory::fill(T, 1.);
    state.add_stimulus_current();
    EXPECT_TRUE(near_relative(-0.2, J[soma_cv]*A[soma_cv]*unit_factor, reltol));

    // Test that at t=10, no more current is injected at soma, and that
    // that 0.3 nA is injected at dendrite tip.
    memory::fill(T, 10.);
    state.add_stimulus_current();
    EXPECT_TRUE(near_relative(-0.2, J[soma_cv]*A[soma_cv]*unit_factor, reltol));
    EXPECT_TRUE(near_relative(-0.3, J[tip_cv]*A[tip_cv]*unit_factor, reltol));
}

TEST(fvm_lowered, ac_stimulus) {
    // Simple cell (one CV) with oscillating stimulus.

    arb::execution_context context; // Just use default context for this one!

    decor dec;
    segment_tree tree;
    tree.append(mnpos, {0., 0., 0., 1.}, {100., 0., 0., 1.}, 1);

    const double freq = 0.02; // (kHz)
    const double phase = 0.1; // (radian)
    const double max_amplitude = 30; // (nA)
    const double max_time = 8; // (ms)

    // Envelope is linear ramp from 0 to max_time.
    dec.place(mlocation{0, 0}, i_clamp({{0, 0}, {max_time, max_amplitude}, {max_time, 0}}, freq, phase), "clamp");
    std::vector<cable_cell> cells = {cable_cell(tree, {}, dec)};

    cable_cell_global_properties gprop;
    gprop.default_parameters = neuron_parameter_defaults;

    fvm_cv_discretization D = fvm_cv_discretize(cells, gprop.default_parameters, context);
    const auto& A = D.cv_area;

    fvm_cell fvcell(context);
    fvcell.initialize({0}, cable1d_recipe(cells));

    auto& state = *(fvcell.*private_state_ptr).get();
    auto& J = state.current_density;
    auto& T = state.time;

    // Current from t = 0 to max_time should be given by
    //     I = max_amplitude * t/max_time * sin(2π f t)
    // where t is in ms and f = freq/1000 is the frequency in kHz.

    constexpr double reltol = 1e-10;
    using testing::near_relative;

    for (double t: {0., 0.1*max_time, 0.7*max_time, 1.1*max_time}) {
        constexpr double unit_factor = 1e-3; // scale A/m²·µm² to nA

        memory::fill(J, 0.);
        memory::fill(T, t);
        state.add_stimulus_current();

        double expected_I = t<=max_time? max_amplitude*t/max_time*std::sin(2*math::pi<double>*t*freq+phase): 0;
        EXPECT_TRUE(near_relative(-expected_I, J[0]*A[0]*unit_factor, reltol));
    }
}

// Test derived mechanism behaviour.

TEST(fvm_lowered, derived_mechs) {
    // Create ball and stick cells with the 'test_kin1' mechanism, which produces
    // a voltage-independent current density of the form a + exp(-t/tau) as a function
    // of time t.
    //
    // 1. Default 'test_kin1': tau = 10 [ms].
    //
    // 2. Specialized version 'custom_kin1' with tau = 20 [ms].
    //
    // 3. Cell with both test_kin1 and custom_kin1.

    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }

    std::vector<cable_cell> cells;
    cells.reserve(3);
    soma_cell_builder builder(6);
    builder.add_branch(0, 100, 0.5, 0.5, 4, "dend");
    for (int i = 0; i<3; ++i) {
        auto cell = builder.make_cell();

        switch (i) {
            case 0:
                cell.decorations.paint(reg::all(), "test_kin1");
                break;
            case 1:
                cell.decorations.paint(reg::all(), "custom_kin1");
                break;
            case 2:
                cell.decorations.paint(reg::all(), "test_kin1");
                cell.decorations.paint(reg::all(), "custom_kin1");
                break;
        }
        cells.push_back(cell);
    }

    cable1d_recipe rec(cells);
    rec.catalogue() = make_unit_test_catalogue();
    rec.catalogue().derive("custom_kin1", "test_kin1", {{"tau", 20.0}});

    cable_probe_total_ion_current_density where{builder.location({1, 0.3})};
    rec.add_probe(0, 0, where);
    rec.add_probe(1, 0, where);
    rec.add_probe(2, 0, where);

    {
        // Test initialization and global parameter values.

        arb::execution_context context(resources);
        fvm_cell fvcell(context);
        fvcell.initialize({0, 1, 2}, rec);

        // Both mechanisms will have the same internal name, "test_kin1".

        using fvec = std::vector<fvm_value_type>;
        fvec tau_values;
        for (auto& mech: fvcell.*private_mechanisms_ptr) {
            EXPECT_EQ("test_kin1"s, mech->internal_name());

            ASSERT_TRUE(mech);

            auto opt_tau_ptr = util::value_by_key((mech.get()->*private_global_table_ptr)(), "tau"s);
            ASSERT_TRUE(opt_tau_ptr);
            tau_values.push_back(opt_tau_ptr.value());
        }
        util::sort(tau_values);
        EXPECT_EQ(fvec({10., 20.}), tau_values);
    }
    {
        // Test dynamics:
        // 1. Current at same point on cell 0 at time 10 ms should equal that
        //    on cell 1 at time 20 ms.
        // 2. Current for cell 2 should be sum of currents for cells 0 and 1 at any given time.

        std::vector<double> samples[3];

        sampler_function sampler =
            [&](probe_metadata pm, std::size_t n, const sample_record* records) {
                for (std::size_t i = 0; i<n; ++i) {
                    double v = *util::any_cast<const double*>(records[i].data);
                    samples[pm.id.gid].push_back(v);
                }
            };

        float times[] = {10.f, 20.f};

        auto ctx = make_context(resources);
        auto decomp = partition_load_balance(rec, ctx);
        simulation sim(rec, decomp, ctx);
        sim.add_sampler(all_probes, explicit_schedule(times), sampler);
        sim.run(30.0, 1.f/1024);

        ASSERT_EQ(2u, samples[0].size());
        ASSERT_EQ(2u, samples[1].size());
        ASSERT_EQ(2u, samples[2].size());

        // Integration isn't exact: let's aim for one part in 10'000.
        double relerr = 0.0001;
        EXPECT_TRUE(testing::near_relative(samples[0][0], samples[1][1], relerr));
        EXPECT_TRUE(testing::near_relative(samples[0][0]+samples[1][0], samples[2][0], relerr));
        EXPECT_TRUE(testing::near_relative(samples[0][1]+samples[1][1], samples[2][1], relerr));
    }
}

// Test that ion charge is propagated into mechanism variable.

TEST(fvm_lowered, read_valence) {
    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }

    {
        std::vector<cable_cell> cells(1);

        soma_cell_builder builder(6);
        auto cell = builder.make_cell();
        cell.decorations.paint("\"soma\"", "test_ca_read_valence");
        cable1d_recipe rec(cable_cell{cell});
        rec.catalogue() = make_unit_test_catalogue();

        arb::execution_context context(resources);
        fvm_cell fvcell(context);
        fvcell.initialize({0}, rec);

        // test_ca_read_valence initialization should write ca ion valence
        // to state variable 'record_zca':

        auto mech_ptr = find_mechanism(fvcell, "test_ca_read_valence");
        auto opt_record_z_ptr = util::value_by_key((mech_ptr->*private_field_table_ptr)(), "record_z"s);

        ASSERT_TRUE(opt_record_z_ptr);
        auto& record_z = opt_record_z_ptr.value().first;
        ASSERT_EQ(2.0, record_z[0]);
    }

    {
        // Check ion renaming.
        soma_cell_builder builder(6);
        auto cell = builder.make_cell();
        cell.decorations.paint("\"soma\"", "cr_read_valence");
        cable1d_recipe rec(cable_cell{cell});
        rec.catalogue() = make_unit_test_catalogue();
        rec.catalogue() = make_unit_test_catalogue();

        rec.catalogue().derive("na_read_valence", "test_ca_read_valence", {}, {{"ca", "na"}});
        rec.catalogue().derive("cr_read_valence", "na_read_valence", {}, {{"na", "mn"}});
        rec.add_ion("mn", 7, 0, 0, 0);

        arb::execution_context context(resources);
        fvm_cell fvcell(context);
        fvcell.initialize({0}, rec);

        auto cr_mech_ptr = find_mechanism(fvcell, 0);
        auto cr_opt_record_z_ptr = util::value_by_key((cr_mech_ptr->*private_field_table_ptr)(), "record_z"s);

        ASSERT_TRUE(cr_opt_record_z_ptr);
        auto& cr_record_z = cr_opt_record_z_ptr.value().first;
        ASSERT_EQ(7.0, cr_record_z[0]);
    }
}

// Test correct scaling of ionic currents in reading and writing
TEST(fvm_lowered, ionic_concentrations) {
    auto cat = make_unit_test_catalogue();

    // one cell, one CV:
    fvm_size_type ncell = 1;
    fvm_size_type ncv = 1;
    std::vector<fvm_index_type> cv_to_intdom(ncv, 0);
    std::vector<fvm_value_type> temp(ncv, 23);
    std::vector<fvm_value_type> diam(ncv, 1.);
    std::vector<fvm_value_type> vinit(ncv, -65);
    std::vector<fvm_gap_junction> gj = {};
    std::vector<fvm_index_type> src_to_spike = {};

    fvm_ion_config ion_config;
    mechanism_layout layout;
    mechanism_overrides overrides;

    layout.weight.assign(ncv, 1.);
    for (fvm_size_type i = 0; i<ncv; ++i) {
        layout.cv.push_back(i);
        ion_config.cv.push_back(i);
    }
    ion_config.init_revpot.assign(ncv, 0.);
    ion_config.init_econc.assign(ncv, 0.);
    ion_config.init_iconc.assign(ncv, 0.);
    ion_config.reset_econc.assign(ncv, 0.);
    ion_config.reset_iconc.assign(ncv, 2.3e-4);

    auto read_cai  = cat.instance(backend::kind, "read_cai_init");
    auto write_cai = cat.instance(backend::kind, "write_cai_breakpoint");

    auto& read_cai_mech  = read_cai.mech;
    auto& write_cai_mech = write_cai.mech;

    auto shared_state = std::make_unique<typename backend::shared_state>(
            ncell, ncell, 0, cv_to_intdom, cv_to_intdom, gj, vinit, temp, diam, src_to_spike, read_cai_mech->data_alignment());
    shared_state->add_ion("ca", 2, ion_config);

    shared_state->instantiate(*read_cai_mech, 0, overrides, layout);
    shared_state->instantiate(*write_cai_mech, 1, overrides, layout);

    shared_state->reset();

    // expect 2.3 value in state 's' in read_cai_init after init:
    read_cai_mech->initialize();
    write_cai_mech->initialize();

    std::vector<fvm_value_type> expected_s_values(ncv, 2.3e-4);

    EXPECT_EQ(expected_s_values, mechanism_field(read_cai_mech.get(), "s"));

    // expect 5.2 + 2.3 value in state 's' in read_cai_init after state update:
    read_cai_mech->update_state();
    write_cai_mech->update_state();

    read_cai_mech->update_ions();
    write_cai_mech->update_ions();

    read_cai_mech->update_state();

    expected_s_values.assign(ncv, 7.5e-4);
    EXPECT_EQ(expected_s_values, mechanism_field(read_cai_mech.get(), "s"));
}

TEST(fvm_lowered, ionic_currents) {
    using namespace arb::literals;
    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    soma_cell_builder b(6);

    // Mechanism parameter is in NMODL units, i.e. mA/cm².

    const double jca = 1.5;
    mechanism_desc m1("fixed_ica_current");
    m1["current_density"] = jca;

    // Mechanism models a well-mixed fixed-depth volume without replenishment,
    // giving a linear response to ica over time.
    //
    //     cai' = - coeff · ica
    //
    // with NMODL units: cai' [mM/ms]; ica [mA/cm²], giving coeff in [mol/cm/C].

    const double coeff = 0.5;
    mechanism_desc m2("linear_ca_conc");
    m2["coeff"] = coeff;

    auto c = b.make_cell();
    c.decorations.paint("soma"_lab, m1);
    c.decorations.paint("soma"_lab, m2);

    cable1d_recipe rec({cable_cell{c}});
    rec.catalogue() = make_unit_test_catalogue();

    fvm_cell fvcell(context);
    fvcell.initialize({0}, rec);

    auto& state = *(fvcell.*private_state_ptr).get();
    auto& ion = state.ion_data.at("ca"s);

    // Ionic current should be 15 A/m², and initial concentration zero.
    EXPECT_EQ(15, ion.iX_[0]);
    EXPECT_EQ(0, ion.Xi_[0]);

    // Integration should be (effectively) exact, so check for linear response.
    const double time = 12; // [ms]
    (void)fvcell.integrate(time, 0.1, {}, {});
    double expected_Xi = -time*coeff*jca;
    EXPECT_NEAR(expected_Xi, ion.Xi_[0], 1e-6);
}

// Test correct scaling of an ionic current updated via a point mechanism

TEST(fvm_lowered, point_ionic_current) {
    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    double r = 6.0; // [µm]
    soma_cell_builder b(r);
    auto c = b.make_cell();

    double soma_area_m2 = 4*math::pi<double>*r*r*1e-12; // [m²]

    // Event weight is translated by point_ica_current into a current contribution in nA.
    c.decorations.place(mlocation{0u, 0.5}, "point_ica_current", "syn");

    cable1d_recipe rec({cable_cell{c}});
    rec.catalogue() = make_unit_test_catalogue();

    fvm_cell fvcell(context);
    fvcell.initialize({0}, rec);

    // Only one target, corresponding to our point process on soma.
    double ica_nA = 12.3;
    deliverable_event ev = {0.04, target_handle{0, 0, 0}, (float)ica_nA};

    auto& state = *(fvcell.*private_state_ptr).get();
    auto& ion = state.ion_data.at("ca"s);

    // Ionic current should be 0 A/m² after initialization.
    EXPECT_EQ(0, ion.iX_[0]);

    // Ionic current should be ica_nA/soma_area after integrating past event time.
    const double time = 0.5; // [ms]
    (void)fvcell.integrate(time, 0.01, {ev}, {});

    double expected_iX = ica_nA*1e-9/soma_area_m2;
    EXPECT_FLOAT_EQ(expected_iX, ion.iX_[0]);
}

// Test area-weighted linear combination of ion species concentrations

TEST(fvm_lowered, weighted_write_ion) {
    // Create a cell with 3 branches (same morphopology as in fvm_layout.ion_weights test):
    //   - Soma (part of branch 0) plus three dendrites (d1, d2, d3) meeting at a branch point.
    //   - Dendritic segments are given 1 compartments each.
    //
    //          /
    //         d2
    //        /
    //   s0-d1
    //        \.
    //         d3
    //
    // The CV corresponding to the branch point should comprise the terminal
    // 1/2 of branch 1 and the initial 1/2 of branches 2 and 3.
    //
    // Geometry:
    //   soma 0:  10 µm long, 10 µm diameter cylinder: area = 100π μm²
    //   dend 1: 100 µm long,  1 µm diameter cylinder: area = 100π μm²
    //   dend 2: 200 µm long,  1 µm diameter cylinder: area = 200π μm²
    //   dend 3: 100 µm long,  1 µm diameter cylinder: area = 100π μm²
    //
    // The radius of the soma is chosen such that the surface area of soma is
    // the same as a 100 µm dendrite, which makes it easier to describe the
    // expected weights.

    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    soma_cell_builder b(5);
    b.add_branch(0, 100, 0.5, 0.5, 1, "dend");
    b.add_branch(1, 200, 0.5, 0.5, 1, "dend");
    b.add_branch(1, 100, 0.5, 0.5, 1, "dend");

    auto c = b.make_cell();

    const double con_int = 80;
    const double con_ext = 120;

    // Ca ion reader test_kinlva on CV 2 and 3 via branch 2:
    c.decorations.paint(reg::branch(1), "test_kinlva");

    // Ca ion writer test_ca on CV 2 and 4 via branch 3:
    c.decorations.paint(reg::branch(2), "test_ca");

    cable1d_recipe rec({cable_cell{c}});
    rec.catalogue() = make_unit_test_catalogue();
    rec.add_ion("ca", 2, con_int, con_ext, 0.0);

    fvm_cell fvcell(context);
    fvcell.initialize({0}, rec);

    auto& state = *(fvcell.*private_state_ptr).get();
    auto& ion = state.ion_data.at("ca"s);
    ion.init_concentration();

    std::vector<unsigned> ion_nodes = util::assign_from(ion.node_index_);
    std::vector<unsigned> expected_ion_nodes = {2, 3, 4};
    EXPECT_EQ(expected_ion_nodes, ion_nodes);

    std::vector<double> ion_init_iconc = util::assign_from(ion.init_Xi_);
    std::vector<double> expected_init_iconc = {0.75*con_int, 1.*con_int, 0};
    EXPECT_TRUE(testing::seq_almost_eq<double>(expected_init_iconc, ion_init_iconc));

    auto test_ca = find_mechanism(fvcell, "test_ca");

    auto opt_cai_ptr = util::value_by_key((test_ca->*private_field_table_ptr)(), "cai"s);
    ASSERT_TRUE(opt_cai_ptr);
    auto& test_ca_cai = opt_cai_ptr.value().first;

    auto opt_ca_index_ptr = util::value_by_key((test_ca->*private_ion_table_ptr)(), "ca"s);
    ASSERT_TRUE(opt_ca_index_ptr);
    auto& test_ca_ca_index = opt_ca_index_ptr.value().second;

    double cai_contrib[3] = {200., 0., 300.};
    double test_ca_weight[3] = {0.25, 0., 1.};

    for (int i = 0; i<2; ++i) {
        test_ca_cai[i] = cai_contrib[test_ca_ca_index[i]];
    }

    std::vector<double> expected_iconc(3);
    for (int i = 0; i<3; ++i) {
        expected_iconc[i] = test_ca_weight[i]*cai_contrib[i] + ion_init_iconc[i];
    }

    ion.init_concentration();
    test_ca->update_ions();
    std::vector<double> ion_iconc = util::assign_from(ion.Xi_);
    EXPECT_TRUE(testing::seq_almost_eq<double>(expected_iconc, ion_iconc));
}

TEST(fvm_lowered, gj_coords_simple) {
    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    using pair = std::pair<int, int>;

    class gap_recipe: public recipe {
    public:
        gap_recipe(std::vector<cable_cell> cells) : cells_(cells) {
            gprop_.default_parameters = neuron_parameter_defaults;
        }

        cell_size_type num_cells() const override { return n_; }
        cell_kind get_cell_kind(cell_gid_type) const override { return cell_kind::cable; }
        util::unique_any get_cell_description(cell_gid_type gid) const override {
            return cells_[gid];
        }
        std::vector<arb::gap_junction_connection> gap_junctions_on(cell_gid_type gid) const override{
            std::vector<gap_junction_connection> conns;
            conns.push_back(gap_junction_connection({(gid+1)%2, "gj", lid_selection_policy::assert_univalent}, {"gj", lid_selection_policy::assert_univalent}, 0.5));
            return conns;
        }
        std::any get_global_properties(cell_kind) const override {
            return gprop_;
        }
    protected:
        arb::cable_cell_global_properties gprop_;
        std::vector<cable_cell> cells_;
        cell_size_type n_ = 2;
    };

    fvm_cell fvcell(context);

    std::vector<cable_cell> cells;
    {
        soma_cell_builder b(2.1);
        b.add_branch(0, 10, 0.3, 0.2, 5, "dend");
        auto c = b.make_cell();
        c.decorations.place(b.location({1, 0.8}), gap_junction_site{}, "gj");
        cells.push_back(c);
    }

    {
        soma_cell_builder b(2.4);
        b.add_branch(0, 10, 0.3, 0.2, 2, "dend");
        auto c = b.make_cell();
        c.decorations.place(b.location({1, 1}), gap_junction_site{}, "gj");
        cells.push_back(c);
    }
    gap_recipe rec(cells);

    fvm_cv_discretization D = fvm_cv_discretize(cells, neuron_parameter_defaults, context);

    std::vector<cell_gid_type> gids = {0, 1};
    auto fvm_info = fvcell.initialize(gids, rec);

    auto GJ = fvcell.fvm_gap_junctions(cells, gids, fvm_info.gap_junction_data, rec, D);

    auto weight = [&](fvm_value_type g, fvm_index_type i){
        return g * 1e3 / D.cv_area[i];
    };

    EXPECT_EQ(pair({5,10}), GJ[0].loc);
    EXPECT_EQ(weight(0.5, 5), GJ[0].weight);

    EXPECT_EQ(pair({10,5}), GJ[1].loc);
    EXPECT_EQ(weight(0.5, 10), GJ[1].weight);
}

TEST(fvm_lowered, gj_coords_complex) {
    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    class gap_recipe: public recipe {
    public:
        gap_recipe(std::vector<cable_cell> cells) : cells_(cells) {
            gprop_.default_parameters = neuron_parameter_defaults;
        }

        cell_size_type num_cells() const override { return n_; }
        cell_kind get_cell_kind(cell_gid_type) const override { return cell_kind::cable; }
        util::unique_any get_cell_description(cell_gid_type gid) const override {
            return cells_[gid];
        }
        std::vector<arb::gap_junction_connection> gap_junctions_on(cell_gid_type gid) const override{
            std::vector<gap_junction_connection> conns;
            switch (gid) {
            case 0:
                return {
                    gap_junction_connection({2, "gj0", lid_selection_policy::assert_univalent}, {"gj1", lid_selection_policy::assert_univalent}, 0.01),
                    gap_junction_connection({1, "gj0", lid_selection_policy::assert_univalent}, {"gj0", lid_selection_policy::assert_univalent}, 0.03),
                    gap_junction_connection({1, "gj1", lid_selection_policy::assert_univalent}, {"gj0", lid_selection_policy::assert_univalent}, 0.04)
                };
            case 1:
                return {
                    gap_junction_connection({0, "gj0", lid_selection_policy::assert_univalent}, {"gj0", lid_selection_policy::assert_univalent}, 0.03),
                    gap_junction_connection({0, "gj0", lid_selection_policy::assert_univalent}, {"gj1", lid_selection_policy::assert_univalent}, 0.04),
                    gap_junction_connection({2, "gj1", lid_selection_policy::assert_univalent}, {"gj2", lid_selection_policy::assert_univalent}, 0.02),
                    gap_junction_connection({2, "gj2", lid_selection_policy::assert_univalent}, {"gj3", lid_selection_policy::assert_univalent}, 0.01)
                };
            case 2:
                return {
                    gap_junction_connection({0, "gj1", lid_selection_policy::assert_univalent}, {"gj0", lid_selection_policy::assert_univalent}, 0.01),
                    gap_junction_connection({1, "gj2", lid_selection_policy::assert_univalent}, {"gj1", lid_selection_policy::assert_univalent}, 0.02),
                    gap_junction_connection({1, "gj3", lid_selection_policy::assert_univalent}, {"gj2", lid_selection_policy::assert_univalent}, 0.01)
                };
            default : return {};
            }
            return conns;
        }
        std::any get_global_properties(cell_kind) const override {
            return gprop_;
        }
    protected:
        arb::cable_cell_global_properties gprop_;
        std::vector<cable_cell> cells_;
        cell_size_type n_ = 3;
    };

    // Add 5 gap junctions
    soma_cell_builder b0(2.1);
    b0.add_branch(0, 8, 0.3, 0.2, 4, "dend");

    auto c0 = b0.make_cell();
    mlocation c0_gj[2] = {b0.location({1, 1}), b0.location({1, 0.5})};

    c0.decorations.place(c0_gj[0], gap_junction_site{}, "gj0");
    c0.decorations.place(c0_gj[1], gap_junction_site{}, "gj1");

    soma_cell_builder b1(1.4);
    b1.add_branch(0, 12, 0.3, 0.5, 6, "dend");
    b1.add_branch(1,  9, 0.3, 0.2, 3, "dend");
    b1.add_branch(1,  5, 0.2, 0.2, 5, "dend");

    auto c1 = b1.make_cell();
    mlocation c1_gj[4] = {b1.location({2, 1}), b1.location({1, 1}), b1.location({1, 0.45}), b1.location({1, 0.1})};

    c1.decorations.place(c1_gj[0], gap_junction_site{}, "gj0");
    c1.decorations.place(c1_gj[1], gap_junction_site{}, "gj1");
    c1.decorations.place(c1_gj[2], gap_junction_site{}, "gj2");
    c1.decorations.place(c1_gj[3], gap_junction_site{}, "gj3");


    soma_cell_builder b2(2.9);
    b2.add_branch(0, 4, 0.3, 0.5, 2, "dend");
    b2.add_branch(1, 6, 0.4, 0.2, 2, "dend");
    b2.add_branch(1, 8, 0.1, 0.2, 2, "dend");
    b2.add_branch(2, 4, 0.2, 0.2, 2, "dend");
    b2.add_branch(2, 4, 0.2, 0.2, 2, "dend");

    auto c2 = b2.make_cell();
    mlocation c2_gj[3] = {b2.location({1, 0.5}), b2.location({4, 1}), b2.location({2, 1})};

    c2.decorations.place(c2_gj[0], gap_junction_site{}, "gj0");
    c2.decorations.place(c2_gj[1], gap_junction_site{}, "gj1");
    c2.decorations.place(c2_gj[2], gap_junction_site{}, "gj2");

    std::vector<cable_cell> cells{c0, c1, c2};
    std::vector<cell_gid_type> gids = {0, 1, 2};

    gap_recipe rec(cells);
    fvm_cell fvcell(context);

    auto fvm_info = fvcell.initialize(gids, rec);
    fvcell.fvm_intdom(rec, gids, fvm_info.cell_to_intdom);
    fvm_cv_discretization D = fvm_cv_discretize(cells, neuron_parameter_defaults, context);

    using namespace cv_prefer;
    int c0_gj_cv[2];
    for (int i = 0; i<2; ++i) c0_gj_cv[i] = D.geometry.location_cv(0, c0_gj[i], cv_nonempty);

    int c1_gj_cv[4];
    for (int i = 0; i<4; ++i) c1_gj_cv[i] = D.geometry.location_cv(1, c1_gj[i], cv_nonempty);

    int c2_gj_cv[3];
    for (int i = 0; i<3; ++i) c2_gj_cv[i] = D.geometry.location_cv(2, c2_gj[i], cv_nonempty);

    std::vector<fvm_gap_junction> GJ = fvcell.fvm_gap_junctions(cells, gids, fvm_info.gap_junction_data, rec, D);
    EXPECT_EQ(10u, GJ.size());

    auto weight = [&](fvm_value_type g, fvm_index_type i){
        return g * 1e3 / D.cv_area[i];
    };

    std::vector<fvm_gap_junction> expected = {
        {{c0_gj_cv[0], c1_gj_cv[0]}, weight(0.03, c0_gj_cv[0])},
        {{c0_gj_cv[0], c1_gj_cv[1]}, weight(0.04, c0_gj_cv[0])},
        {{c0_gj_cv[1], c2_gj_cv[0]}, weight(0.01, c0_gj_cv[1])},
        {{c1_gj_cv[0], c0_gj_cv[0]}, weight(0.03, c1_gj_cv[0])},
        {{c1_gj_cv[1], c0_gj_cv[0]}, weight(0.04, c1_gj_cv[1])},
        {{c1_gj_cv[2], c2_gj_cv[1]}, weight(0.02, c1_gj_cv[2])},
        {{c1_gj_cv[3], c2_gj_cv[2]}, weight(0.01, c1_gj_cv[3])},
        {{c2_gj_cv[0], c0_gj_cv[1]}, weight(0.01, c2_gj_cv[0])},
        {{c2_gj_cv[1], c1_gj_cv[2]}, weight(0.02, c2_gj_cv[1])},
        {{c2_gj_cv[2], c1_gj_cv[3]}, weight(0.01, c2_gj_cv[2])}
    };

    using util::sort_by;
    using util::transform_view;

    auto gj_loc = [](const fvm_gap_junction gj) { return gj.loc; };
    auto gj_weight = [](const fvm_gap_junction gj) { return gj.weight; };

    sort_by(GJ, [](fvm_gap_junction gj) { return gj.loc; });
    sort_by(expected, [](fvm_gap_junction gj) { return gj.loc; });

    EXPECT_TRUE(testing::seq_eq(transform_view(expected, gj_loc), transform_view(GJ, gj_loc)));
    EXPECT_TRUE(testing::seq_almost_eq<double>(transform_view(expected, gj_weight), transform_view(GJ, gj_weight)));
}

TEST(fvm_lowered, cell_group_gj) {
    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    }
    else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    using pair = std::pair<int, int>;

    class gap_recipe: public recipe {
    public:
        gap_recipe(const std::vector<cable_cell>& cg0, const std::vector<cable_cell>& cg1) {
            cells_ = cg0;
            cells_.insert(cells_.end(), cg1.begin(), cg1.end());
            gprop_.default_parameters = neuron_parameter_defaults;
        }

        cell_size_type num_cells() const override { return n_; }
        cell_kind get_cell_kind(cell_gid_type) const override { return cell_kind::cable; }
        util::unique_any get_cell_description(cell_gid_type gid) const override {
            return cells_[gid];
        }
        std::vector<arb::gap_junction_connection> gap_junctions_on(cell_gid_type gid) const override{
            std::vector<gap_junction_connection> conns;
            if (gid % 2 == 0) {
                // connect 5 of the first 10 cells in a ring; connect 5 of the second 10 cells in a ring
                auto next_cell = gid == 8 ? 0 : (gid == 18 ? 10 : gid + 2);
                auto prev_cell = gid == 0 ? 8 : (gid == 10 ? 18 : gid - 2);
                conns.push_back(gap_junction_connection({next_cell, "gj", lid_selection_policy::assert_univalent},
                                                             {"gj", lid_selection_policy::assert_univalent}, 0.03));
                conns.push_back(gap_junction_connection({prev_cell, "gj", lid_selection_policy::assert_univalent},
                                                             {"gj", lid_selection_policy::assert_univalent}, 0.03));
            }
            return conns;
        }
        std::any get_global_properties(cell_kind) const override {
            return gprop_;
        }

    protected:
        arb::cable_cell_global_properties gprop_;
        std::vector<cable_cell> cells_;
        cell_size_type n_ = 20;
    };

    std::vector<cable_cell> cell_group0;
    std::vector<cable_cell> cell_group1;

    // Make 20 cells
    for (unsigned i = 0; i < 20; i++) {
        cable_cell_description c = soma_cell_builder(2.1).make_cell();
        if (i % 2 == 0) {
            c.decorations.place(mlocation{0, 1}, gap_junction_site{}, "gj");
        }
        if (i < 10) {
            cell_group0.push_back(c);
        }
        else {
            cell_group1.push_back(c);
        }
    }

    std::vector<cell_gid_type> gids_cg0 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<cell_gid_type> gids_cg1 = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

    gap_recipe rec(cell_group0, cell_group1);

    fvm_cell fvcell0(context);
    fvm_cell fvcell1(context);

    auto fvm_info_0 = fvcell0.initialize(gids_cg0, rec);
    auto fvm_info_1 = fvcell1.initialize(gids_cg1, rec);

    auto num_dom0 = fvcell0.fvm_intdom(rec, gids_cg0, fvm_info_0.cell_to_intdom);
    auto num_dom1 = fvcell1.fvm_intdom(rec, gids_cg1, fvm_info_1.cell_to_intdom);

    fvm_cv_discretization D0 = fvm_cv_discretize(cell_group0, neuron_parameter_defaults, context);
    fvm_cv_discretization D1 = fvm_cv_discretize(cell_group1, neuron_parameter_defaults, context);

    auto GJ0 = fvcell0.fvm_gap_junctions(cell_group0, gids_cg0, fvm_info_0.gap_junction_data, rec, D0);
    auto GJ1 = fvcell1.fvm_gap_junctions(cell_group1, gids_cg1, fvm_info_1.gap_junction_data, rec, D1);

    EXPECT_EQ(10u, GJ0.size());
    EXPECT_EQ(10u, GJ1.size());

    std::vector<pair> expected_loc = {{0, 2}, {0, 8}, {2, 4}, {2, 0}, {4, 6} ,{4, 2}, {6, 8}, {6, 4}, {8, 0}, {8, 6}};

    for (unsigned i = 0; i < GJ0.size(); i++) {
        EXPECT_EQ(expected_loc[i], GJ0[i].loc);
        EXPECT_EQ(expected_loc[i], GJ1[i].loc);
    }

    std::vector<fvm_index_type> expected_doms= {0u, 1u, 0u, 2u, 0u, 3u, 0u, 4u, 0u, 5u};
    EXPECT_EQ(6u, num_dom0);
    EXPECT_EQ(6u, num_dom1);

    EXPECT_EQ(expected_doms, fvm_info_0.cell_to_intdom);
    EXPECT_EQ(expected_doms, fvm_info_1.cell_to_intdom);

}

TEST(fvm_lowered, integration_domains) {
    {
        execution_context context;
        fvm_cell fvcell(context);

        std::vector<cell_gid_type> gids = {11u, 5u, 2u, 3u, 0u, 8u, 7u};
        std::vector<fvm_index_type> cell_to_intdom;

        auto num_dom = fvcell.fvm_intdom(gap_recipe_0(), gids, cell_to_intdom);
        std::vector<fvm_index_type> expected_doms= {0u, 1u, 2u, 2u, 1u, 3u, 2u};

        EXPECT_EQ(4u, num_dom);
        EXPECT_EQ(expected_doms, cell_to_intdom);
    }
    {
        execution_context context;
        fvm_cell fvcell(context);

        std::vector<cell_gid_type> gids = {11u, 5u, 2u, 3u, 0u, 8u, 7u};
        std::vector<fvm_index_type> cell_to_intdom;

        auto num_dom = fvcell.fvm_intdom(gap_recipe_1(), gids, cell_to_intdom);
        std::vector<fvm_index_type> expected_doms= {0u, 1u, 2u, 3u, 4u, 5u, 6u};

        EXPECT_EQ(7u, num_dom);
        EXPECT_EQ(expected_doms, cell_to_intdom);
    }
    {
        execution_context context;
        fvm_cell fvcell(context);

        std::vector<cell_gid_type> gids = {5u, 2u, 3u, 0u};
        std::vector<fvm_index_type> cell_to_intdom;

        auto num_dom = fvcell.fvm_intdom(gap_recipe_2(), gids, cell_to_intdom);
        std::vector<fvm_index_type> expected_doms= {0u, 0u, 0u, 0u};

        EXPECT_EQ(1u, num_dom);
        EXPECT_EQ(expected_doms, cell_to_intdom);
    }
}

TEST(fvm_lowered, post_events_shared_state) {
    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    } else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    class detector_recipe: public arb::recipe {
    public:
        detector_recipe(unsigned ncv, std::vector<unsigned> detectors_per_cell, std::string synapse):
                ncell_(detectors_per_cell.size()),
                ncv_(ncv),
                detectors_per_cell_(detectors_per_cell),
                synapse_(synapse),
                cat_(make_unit_test_catalogue()) {
            const auto default_cat = arb::global_default_catalogue();
            cat_.import(default_cat, "");
        }

        cell_size_type num_cells() const override {
            return ncell_;
        }

        arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
            arb::segment_tree tree;
            tree.append(arb::mnpos, {0, 0, 0.0, 1.0}, {0, 0, 200, 1.0}, 1);

            arb::decor decor;
            decor.set_default(arb::cv_policy_fixed_per_branch(ncv_));

            auto ndetectors = detectors_per_cell_[gid];
            auto offset = 1.0 / ndetectors;
            for (unsigned i = 0; i < ndetectors; ++i) {
                decor.place(arb::mlocation{0, offset * i}, arb::threshold_detector{10}, "detector"+std::to_string(i));
            }
            decor.place(arb::mlocation{0, 0.5}, synapse_, "syanpse");

            return arb::cable_cell(arb::morphology(tree), {}, decor);
        }

        cell_kind get_cell_kind(cell_gid_type gid) const override {
            return cell_kind::cable;
        }

        std::any get_global_properties(arb::cell_kind) const override {
            arb::cable_cell_global_properties gprop;
            gprop.default_parameters = arb::neuron_parameter_defaults;
            gprop.catalogue = &cat_;
            return gprop;
        }

    private:
        unsigned ncell_;
        unsigned ncv_;
        std::vector<unsigned> detectors_per_cell_;
        arb::mechanism_desc synapse_;
        mechanism_catalogue cat_;
    };

    std::vector<unsigned> gids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    const unsigned ncell = gids.size();
    const unsigned cv_per_cell = 10;

    std::vector<std::vector<unsigned>> detectors_per_cell_vec = {
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
            {1, 6, 2, 1, 3, 2, 1, 2, 1, 4},
    };

    for (const auto& detectors_per_cell: detectors_per_cell_vec) {
        detector_recipe rec(cv_per_cell, detectors_per_cell, "post_events_syn");

        fvm_cell fvcell(context);
        fvcell.initialize(gids, rec);

        auto& S = fvcell.*private_state_ptr;

        auto expected_detectors = util::max_value(detectors_per_cell);

        EXPECT_EQ(expected_detectors, S->n_detector);
        EXPECT_EQ(util::sum(detectors_per_cell), S->src_to_spike.size());
        EXPECT_EQ(expected_detectors * ncell, S->time_since_spike.size());

        unsigned detector_id = 0;
        for (unsigned c = 0; c < detectors_per_cell.size(); ++c) {
            for (unsigned d = 0; d < detectors_per_cell[c]; ++d) {
                EXPECT_EQ((int) (c * expected_detectors + d), S->src_to_spike[detector_id++]);
            }
        }
    }
    for (const auto& detectors_per_cell: detectors_per_cell_vec) {
        detector_recipe rec(cv_per_cell, detectors_per_cell, "expsyn");

        fvm_cell fvcell(context);
        fvcell.initialize(gids, rec);

        auto& S = fvcell.*private_state_ptr;

        EXPECT_EQ(0u, S->n_detector);
        EXPECT_EQ(0u, S->src_to_spike.size());
        EXPECT_EQ(0u, S->time_since_spike.size());
    }
}

TEST(fvm_lowered, label_data) {
    arb::proc_allocation resources;
    if (auto nt = arbenv::get_env_num_threads()) {
        resources.num_threads = nt;
    } else {
        resources.num_threads = arbenv::thread_concurrency();
    }
    arb::execution_context context(resources);

    class decorated_recipe: public arb::recipe {
    public:
        decorated_recipe(unsigned ncell): ncell_(ncell) {
            // Cells will have one of 2 decorations:
            //   (1)  1 synapse label with 4 locations; 1 synapse label with 1 location
            //        1 detector label with 1 location
            //        0 gap junctions
            //        if (duplicate) add a duplicate label, check fvm initialize throws
            //   (2)  0 synapse labels
            //        1 detector label with 3 locations; 1 detector label with 2 locations
            //        1 gap-junction label with 2 locations; 1 gap-junction label with 1 location
            using arb::reg::all;
            using arb::ls::uniform;
            arb::segment_tree tree;
            tree.append(arb::mnpos, {0, 0, 0.0, 1.0}, {0, 0, 200, 1.0}, 1);
            {
                arb::decor decor;
                decor.set_default(arb::cv_policy_fixed_per_branch(10));
                decor.place(uniform(all(), 0, 3, 42), "expsyn", "4_synapses");
                decor.place(uniform(all(), 4, 4, 42), "expsyn", "1_synapse");
                decor.place(uniform(all(), 5, 5, 42), arb::threshold_detector{10}, "1_detector");

                cells_.push_back(arb::cable_cell(arb::morphology(tree), {}, decor));
            }
            {
                arb::decor decor;
                decor.set_default(arb::cv_policy_fixed_per_branch(10));
                decor.place(uniform(all(), 0, 2, 24), arb::threshold_detector{10}, "3_detectors");
                decor.place(uniform(all(), 3, 4, 24), arb::threshold_detector{10}, "2_detectors");
                decor.place(uniform(all(), 5, 6, 24), arb::gap_junction_site(), "2_gap_junctions");
                decor.place(uniform(all(), 7, 7, 24), arb::gap_junction_site(), "1_gap_junction");

                cells_.push_back(arb::cable_cell(arb::morphology(tree), {}, decor));
            }
        }

        cell_size_type num_cells() const override {
            return ncell_;
        }

        arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
            if (gid %3 == 0) {
                return cells_[0];
            }
            return cells_[1];
        }

        cell_kind get_cell_kind(cell_gid_type gid) const override {
            return cell_kind::cable;
        }

        std::any get_global_properties(arb::cell_kind) const override {
            arb::cable_cell_global_properties gprop;
            gprop.default_parameters = arb::neuron_parameter_defaults;
            return gprop;
        }

    private:
        unsigned ncell_;
        std::vector<cable_cell> cells_;
    };

    std::vector<unsigned> gids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    const unsigned ncell = gids.size();

    // Check correct synapse, deterctor and gj data
    decorated_recipe rec(ncell);
    std::vector<fvm_index_type> cell_to_intdom;
    std::vector<target_handle> targets;
    probe_association_map probe_map;

    fvm_cell fvcell(context);
    auto fvm_info = fvcell.initialize(gids, rec);

    for (auto gid: gids) {
        if (gid%3 == 0) {
            EXPECT_EQ(5u, fvm_info.num_targets[gid]);
            EXPECT_EQ(1u, fvm_info.num_sources[gid]);
        }
        else {
            EXPECT_EQ(0u, fvm_info.num_targets[gid]);
            EXPECT_EQ(5u, fvm_info.num_sources[gid]);
        }
    }

    // synapses
    {
        auto clg = cell_labels_and_gids(fvm_info.target_data, gids);
        std::vector<cell_size_type> expected_sizes = {2, 0, 0, 2, 0, 0, 2, 0, 0, 2};
        std::vector<std::pair<cell_tag_type, lid_range>> expected_labeled_ranges, actual_labeled_ranges;
        expected_labeled_ranges = {
            {"1_synapse",  {4, 5}}, {"4_synapses", {0, 4}},
            {"1_synapse",  {4, 5}}, {"4_synapses", {0, 4}},
            {"1_synapse",  {4, 5}}, {"4_synapses", {0, 4}},
            {"1_synapse",  {4, 5}}, {"4_synapses", {0, 4}}
        };

        EXPECT_EQ(clg.gids, gids);
        EXPECT_EQ(clg.label_range.sizes(), expected_sizes);
        EXPECT_EQ(clg.label_range.labels().size(), expected_labeled_ranges.size());
        EXPECT_EQ(clg.label_range.ranges().size(), expected_labeled_ranges.size());

        for (unsigned i = 0; i < expected_labeled_ranges.size(); ++i) {
            actual_labeled_ranges.push_back({clg.label_range.labels()[i], clg.label_range.ranges()[i]});
        }

        std::vector<cell_size_type> size_partition;
        auto part = util::make_partition(size_partition, expected_sizes);
        for (const auto& r: part) {
            util::sort(util::subrange_view(actual_labeled_ranges, r));
        }
        EXPECT_EQ(actual_labeled_ranges, expected_labeled_ranges);
    }

    // detectors
    {
        auto clg = cell_labels_and_gids(fvm_info.source_data, gids);
        std::vector<cell_size_type> expected_sizes = {1, 2, 2, 1, 2, 2, 1, 2, 2, 1};
        std::vector<std::pair<cell_tag_type, lid_range>> expected_labeled_ranges, actual_labeled_ranges;
        expected_labeled_ranges = {
            {"1_detector",  {0, 1}},
            {"2_detectors", {3, 5}}, {"3_detectors", {0, 3}},
            {"2_detectors", {3, 5}}, {"3_detectors", {0, 3}},
            {"1_detector",  {0, 1}},
            {"2_detectors", {3, 5}}, {"3_detectors", {0, 3}},
            {"2_detectors", {3, 5}}, {"3_detectors", {0, 3}},
            {"1_detector",  {0, 1}},
            {"2_detectors", {3, 5}}, {"3_detectors", {0, 3}},
            {"2_detectors", {3, 5}}, {"3_detectors", {0, 3}},
            {"1_detector",  {0, 1}}
        };

        EXPECT_EQ(clg.gids, gids);
        EXPECT_EQ(clg.label_range.sizes(), expected_sizes);
        EXPECT_EQ(clg.label_range.labels().size(), expected_labeled_ranges.size());
        EXPECT_EQ(clg.label_range.ranges().size(), expected_labeled_ranges.size());

        for (unsigned i = 0; i < expected_labeled_ranges.size(); ++i) {
            actual_labeled_ranges.push_back({clg.label_range.labels()[i], clg.label_range.ranges()[i]});
        }

        std::vector<cell_size_type> size_partition;
        auto part = util::make_partition(size_partition, expected_sizes);
        for (const auto& r: part) {
            util::sort(util::subrange_view(actual_labeled_ranges, r));
        }
        EXPECT_EQ(actual_labeled_ranges, expected_labeled_ranges);
    }

    // gap_junctions
    {
        auto clg = cell_labels_and_gids(fvm_info.gap_junction_data, gids);
        std::vector<cell_size_type> expected_sizes = {0, 2, 2, 0, 2, 2, 0, 2, 2, 0};
        std::vector<std::pair<cell_tag_type, lid_range>> expected_labeled_ranges, actual_labeled_ranges;
        expected_labeled_ranges = {
            {"1_gap_junction",  {2, 3}}, {"2_gap_junctions", {0, 2}},
            {"1_gap_junction",  {2, 3}}, {"2_gap_junctions", {0, 2}},
            {"1_gap_junction",  {2, 3}}, {"2_gap_junctions", {0, 2}},
            {"1_gap_junction",  {2, 3}}, {"2_gap_junctions", {0, 2}},
            {"1_gap_junction",  {2, 3}}, {"2_gap_junctions", {0, 2}},
            {"1_gap_junction",  {2, 3}}, {"2_gap_junctions", {0, 2}},
        };

        EXPECT_EQ(clg.gids, gids);
        EXPECT_EQ(clg.label_range.sizes(), expected_sizes);
        EXPECT_EQ(clg.label_range.labels().size(), expected_labeled_ranges.size());
        EXPECT_EQ(clg.label_range.ranges().size(), expected_labeled_ranges.size());

        for (unsigned i = 0; i < expected_labeled_ranges.size(); ++i) {
            actual_labeled_ranges.push_back({clg.label_range.labels()[i], clg.label_range.ranges()[i]});
        }

        std::vector<cell_size_type> size_partition;
        auto part = util::make_partition(size_partition, expected_sizes);
        for (const auto& r: part) {
            util::sort(util::subrange_view(actual_labeled_ranges, r));
        }
        EXPECT_EQ(actual_labeled_ranges, expected_labeled_ranges);
    }
}