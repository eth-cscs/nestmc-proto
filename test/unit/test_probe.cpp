#include "../gtest.h"

#include <cmath>
#include <map>
#include <vector>

#include <arbor/cable_cell.hpp>
#include <arbor/common_types.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/mechanism.hpp>
#include <arbor/mechcat.hpp>
#include <arbor/mechinfo.hpp>
#include <arbor/sampling.hpp>
#include <arbor/schedule.hpp>
#include <arbor/simple_sampler.hpp>
#include <arbor/simulation.hpp>
#include <arbor/util/any_cast.hpp>
#include <arbor/util/any_ptr.hpp>
#include <arbor/util/pp_util.hpp>
#include <arbor/version.hpp>
#include <arborenv/gpu_env.hpp>
#include <arbor/mechanism.hpp>

#include "backends/event.hpp"
#include "backends/multicore/fvm.hpp"
#ifdef ARB_GPU_ENABLED
#include "backends/gpu/fvm.hpp"
#endif
#include "fvm_lowered_cell_impl.hpp"
#include "memory/gpu_wrappers.hpp"
#include "util/filter.hpp"
#include "util/rangeutil.hpp"

#include "common.hpp"
#include "common_morphologies.hpp"
#include "unit_test_catalogue.hpp"
#include "../common_cells.hpp"
#include "../simple_recipes.hpp"

using namespace arb;
using util::any_cast;

using multicore_fvm_cell = fvm_lowered_cell_impl<multicore::backend>;
using multicore_shared_state = multicore::backend::shared_state;
ACCESS_BIND(std::unique_ptr<multicore_shared_state> multicore_fvm_cell::*, multicore_fvm_state_ptr, &multicore_fvm_cell::state_);

template <typename Backend>
struct backend_access {
    using fvm_cell = multicore_fvm_cell;

    static multicore_shared_state& state(fvm_cell& cell) {
        return *(cell.*multicore_fvm_state_ptr).get();
    }

    template <typename fvm_type>
    static fvm_type deref(const fvm_type* p) { return *p; }
};

#ifdef ARB_GPU_ENABLED

using gpu_fvm_cell = fvm_lowered_cell_impl<gpu::backend>;
using gpu_shared_state = gpu::backend::shared_state;
ACCESS_BIND(std::unique_ptr<gpu_shared_state> gpu_fvm_cell::*, gpu_fvm_state_ptr, &gpu_fvm_cell::state_);

template <>
struct backend_access<gpu::backend> {
    using fvm_cell = gpu_fvm_cell;

    static gpu_shared_state& state(fvm_cell& cell) {
        return *(cell.*gpu_fvm_state_ptr).get();
    }

    template <typename fvm_type>
    static fvm_type deref(const fvm_type* p) {
        fvm_type r;
        memory::gpu_memcpy_d2h(&r, p, sizeof(r));
        return r;
    }
};

#endif

// Used in a number of tests below, produce a simple Y-shaped 3-branch cell morphology with
// linearly tapered branches.

static morphology make_y_morphology() {
    segment_tree tree;
    tree.append(mnpos, {0., 0., 0., 1.}, {100., 0., 0., 0.8}, 1);
    tree.append(0, {100., 100., 0., 0.5}, 0);
    tree.append(0, {100., 0,  100., 0.4}, 0);
    return {tree};
}

static morphology make_stick_morphology() {
    segment_tree tree;
    tree.append(mnpos, {0., 0., 0., 1.}, {100., 0., 0., 1.0}, 1);
    return {tree};
}

template <typename Backend>
void run_v_i_probe_test(const context& ctx) {
    using fvm_cell = typename backend_access<Backend>::fvm_cell;
    auto deref = [](const fvm_value_type* p) { return backend_access<Backend>::deref(p); };

    soma_cell_builder builder(12.6157/2.0);
    builder.add_branch(0, 200, 1.0/2, 1.0/2, 1, "dend");
    builder.add_branch(0, 200, 1.0/2, 1.0/2, 1, "dend");
    auto bs = builder.make_cell();

    bs.decorations.set_default(cv_policy_fixed_per_branch(1));

    auto stim = i_clamp::box(0, 100, 0.3);
    bs.decorations.place(mlocation{1, 1}, stim, "clamp");

    cable1d_recipe rec((cable_cell(bs)));

    mlocation loc0{0, 0};
    mlocation loc1{1, 1};
    mlocation loc2{1, 0.3};

    rec.add_probe(0, 10, cable_probe_membrane_voltage{loc0});
    rec.add_probe(0, 20, cable_probe_membrane_voltage{loc1});
    rec.add_probe(0, 30, cable_probe_total_ion_current_density{loc2});

    fvm_cell lcell(*ctx);
    auto fvm_info = lcell.initialize({0}, rec);

    const auto& probe_map = fvm_info.probe_map;
    EXPECT_EQ(3u, rec.get_probes(0).size());
    EXPECT_EQ(3u, probe_map.size());

    EXPECT_EQ(10, probe_map.tag.at({0, 0}));
    EXPECT_EQ(20, probe_map.tag.at({0, 1}));
    EXPECT_EQ(30, probe_map.tag.at({0, 2}));

    auto get_probe_raw_handle = [&](cell_member_type x, unsigned i = 0) {
        return probe_map.data_on(x).front().raw_handle_range()[i];
    };

    // Voltage probes are interpolated, so expect fvm_probe_info
    // to wrap an fvm_probe_interpolated; ion current density is
    // implemented as an interpolated probe too, in order to account
    // for stimulus contributions.

    ASSERT_TRUE(std::get_if<fvm_probe_interpolated>(&probe_map.data_on({0, 0}).front().info));
    ASSERT_TRUE(std::get_if<fvm_probe_interpolated>(&probe_map.data_on({0, 1}).front().info));
    ASSERT_TRUE(std::get_if<fvm_probe_interpolated>(&probe_map.data_on({0, 2}).front().info));

    probe_handle p0a = get_probe_raw_handle({0, 0}, 0);
    probe_handle p0b = get_probe_raw_handle({0, 0}, 1);
    probe_handle p1a = get_probe_raw_handle({0, 1}, 0);
    probe_handle p1b = get_probe_raw_handle({0, 1}, 1);
    probe_handle p2a = get_probe_raw_handle({0, 2}, 0);
    probe_handle p2b = get_probe_raw_handle({0, 2}, 1);

    // Ball-and-stick cell with default discretization policy should
    // have three CVs, one for branch 0, one trivial one covering the
    // branch point, and one for branch 1.
    //
    // Consequently, expect the interpolated voltage probe handles
    // to be on CVs 0 and 1 for probe 0,0 on branch 0,
    // and on CVs 1 and 2 for probe 0,1 on branch 1.

    auto& state = backend_access<Backend>::state(lcell);
    auto& voltage = state.voltage;

    EXPECT_EQ(voltage.data(),   p0a);
    EXPECT_EQ(voltage.data()+1, p0b);
    EXPECT_EQ(voltage.data()+1, p1a);
    EXPECT_EQ(voltage.data()+2, p1b);

    // Expect initial raw probe handle values to be the resting potential for
    // the voltage probes (cell membrane potential should be constant), and
    // zero for the current probe (including stimulus component).

    fvm_value_type resting = voltage[0];
    EXPECT_NE(0.0, resting);

    EXPECT_EQ(resting, deref(p0a));
    EXPECT_EQ(resting, deref(p0b));
    EXPECT_EQ(resting, deref(p1a));
    EXPECT_EQ(resting, deref(p1a));
    EXPECT_EQ(0.0, deref(p2a));
    EXPECT_EQ(0.0, deref(p2b));

    // After an integration step, expect voltage probe values
    // to differ from resting, and for there to be a non-zero current.

    lcell.integrate(0.01, 0.0025, {}, {});

    EXPECT_NE(resting, deref(p0a));
    EXPECT_NE(resting, deref(p0b));
    EXPECT_NE(resting, deref(p1a));
    EXPECT_NE(resting, deref(p1b));
    EXPECT_NE(0.0, deref(p2a));
    EXPECT_NE(0.0, deref(p2b));
}

template <typename Backend>
void run_v_cell_probe_test(const context& ctx) {
    using fvm_cell = typename backend_access<Backend>::fvm_cell;

    // Take the per-cable voltage over a Y-shaped cell with and without
    // interior forks in the discretization. The metadata will be used
    // to determine the corresponding CVs for each cable, and the raw
    // pointer to backend data checked against the expected CV offset.

    auto m = make_y_morphology();

    std::pair<const char*, cv_policy> test_policies[] = {
        {"trivial fork", cv_policy_fixed_per_branch(3, cv_policy_flag::none)},
        {"interior fork", cv_policy_fixed_per_branch(3, cv_policy_flag::interior_forks)},
    };

    for (auto& testcase: test_policies) {
        SCOPED_TRACE(testcase.first);
        decor d;
        d.set_default(testcase.second);

        cable_cell cell(m, {}, d);

        cable1d_recipe rec(cell, false);
        rec.add_probe(0, 0, cable_probe_membrane_voltage_cell{});

        fvm_cell lcell(*ctx);
        auto fvm_info = lcell.initialize({0}, rec);

        ASSERT_EQ(1u, fvm_info.probe_map.size());

        const fvm_probe_multi* h_ptr = std::get_if<fvm_probe_multi>(&fvm_info.probe_map.data_on({0, 0}).front().info);
        ASSERT_TRUE(h_ptr);
        auto& h = *h_ptr;

        const mcable_list* cl_ptr = std::get_if<mcable_list>(&h_ptr->metadata);
        ASSERT_TRUE(cl_ptr);
        auto& cl = *cl_ptr;

        ASSERT_EQ(h.raw_handles.size(), cl.size());

        // Independetly discretize the cell so we can follow cable–CV relationship.

        cv_geometry geom = cv_geometry_from_ends(cell, testcase.second.cv_boundary_points(cell));

        // For each cable in metadata, get CV from geom and confirm raw handle is
        // state voltage + CV.

        auto& state = backend_access<Backend>::state(lcell);
        auto& voltage = state.voltage;

        for (auto i: util::count_along(*cl_ptr)) {
            mlocation cable_mid{cl[i].branch, 0.5*(cl[i].prox_pos+cl[i].dist_pos)};
            auto cv = geom.location_cv(0, cable_mid, cv_prefer::cv_empty);

            EXPECT_EQ(voltage.data()+cv, h.raw_handles[i]);
        }
    }
}

template <typename Backend>
void run_expsyn_g_probe_test(const context& ctx) {
    using fvm_cell = typename backend_access<Backend>::fvm_cell;
    auto deref = [](const fvm_value_type* p) { return backend_access<Backend>::deref(p); };

    const double tau = 2.0;
    EXPECT_EQ(tau, global_default_catalogue()["expsyn"].parameters.at("tau").default_value);

    // Ball-and-stick cell, two synapses, both in same CV.
    mlocation loc0{1, 0.8};
    mlocation loc1{1, 1.0};

    soma_cell_builder builder(12.6157/2.0);
    builder.add_branch(0, 200, 1.0/2, 1.0/2, 1, "dend");
    builder.add_branch(0, 200, 1.0/2, 1.0/2, 1, "dend");
    auto bs = builder.make_cell();
    bs.decorations.place(loc0, "expsyn", "syn0");
    bs.decorations.place(loc1, "expsyn", "syn1");
    bs.decorations.set_default(cv_policy_fixed_per_branch(2));

    auto run_test = [&](bool coalesce_synapses) {
        cable1d_recipe rec(cable_cell(bs), coalesce_synapses);
        rec.add_probe(0, 10, cable_probe_point_state{0u, "expsyn", "g"});
        rec.add_probe(0, 20, cable_probe_point_state{1u, "expsyn", "g"});

        fvm_cell lcell(*ctx);
        auto fvm_info = lcell.initialize({0}, rec);
        const auto& probe_map = fvm_info.probe_map;
        const auto& targets = fvm_info.target_handles;

        EXPECT_EQ(2u, rec.get_probes(0).size());
        EXPECT_EQ(2u, probe_map.size());
        ASSERT_EQ(1u, probe_map.data.count({0, 0}));
        ASSERT_EQ(1u, probe_map.data.count({0, 1}));

        EXPECT_EQ(10, probe_map.tag.at({0, 0}));
        EXPECT_EQ(20, probe_map.tag.at({0, 1}));

        auto get_probe_raw_handle = [&](cell_member_type x, unsigned i = 0) {
            return probe_map.data_on(x).front().raw_handle_range()[i];
        };

        probe_handle p0 = get_probe_raw_handle({0, 0});
        probe_handle p1 = get_probe_raw_handle({0, 1});

        // Expect initial probe values to be intial synapse g == 0.

        EXPECT_EQ(0.0, deref(p0));
        EXPECT_EQ(0.0, deref(p1));

        if (coalesce_synapses) {
            // Should be the same raw pointer!
            EXPECT_EQ(p0, p1);
        }

        // Integrate to 3 ms, with one event at 1ms to first expsyn weight 0.5,
        // and another at 2ms to second, weight 1.

        std::vector<deliverable_event> evs = {
            {1.0, targets[0], 0.5},
            {2.0, targets[1], 1.0}
        };
        const double tfinal = 3.;
        const double dt = 0.001;
        lcell.integrate(tfinal, dt, evs, {});

        fvm_value_type g0 = deref(p0);
        fvm_value_type g1 = deref(p1);

        // Expected value: weight*exp(-(t_final-t_event)/tau).
        double expected_g0 = 0.5*std::exp(-(tfinal-1.0)/tau);
        double expected_g1 = 1.0*std::exp(-(tfinal-2.0)/tau);

        const double rtol = 1e-6;
        if (coalesce_synapses) {
            EXPECT_TRUE(testing::near_relative(expected_g0+expected_g1, g0, rtol));
            EXPECT_TRUE(testing::near_relative(expected_g0+expected_g1, g1, rtol));
        }
        else {
            EXPECT_TRUE(testing::near_relative(expected_g0, g0, rtol));
            EXPECT_TRUE(testing::near_relative(expected_g1, g1, rtol));
        }
    };

    {
        SCOPED_TRACE("uncoalesced synapses");
        run_test(false);
    }

    {
        SCOPED_TRACE("coalesced synapses");
        run_test(true);
    }
}

template <typename Backend>
void run_expsyn_g_cell_probe_test(const context& ctx) {
    using fvm_cell = typename backend_access<Backend>::fvm_cell;
    auto deref = [](const auto* p) { return backend_access<Backend>::deref(p); };

    // Place a mixture of expsyn and exp2syn synapses over two cells.
    // Confirm that a whole-cell expsyn g state probe gives:
    //  * A metadata element for each placed expsyn probe.
    //  * A multiplicity that matches the number of probes sharing a CV if synapses
    //    were coalesced.
    //  * A raw handle that sees an event sent to the corresponding target,

    cv_policy policy = cv_policy_fixed_per_branch(3);

    auto m  = make_y_morphology();
    arb::decor d;
    d.set_default(policy);

    std::unordered_map<cell_lid_type, mlocation> expsyn_target_loc_map;

    unsigned n_expsyn = 0;
    for (unsigned bid = 0; bid<3u; ++bid) {
        for (unsigned j = 0; j<10; ++j) {
            auto idx = (bid*10+j)*2;
            mlocation expsyn_loc{bid, 0.1*j};
            d.place(expsyn_loc, "expsyn", "syn"+std::to_string(idx));
            expsyn_target_loc_map[2*n_expsyn] = expsyn_loc;
            d.place(mlocation{bid, 0.1*j+0.05}, "exp2syn", "syn"+std::to_string(idx+1));
            ++n_expsyn;
        }
    }

    std::vector<cable_cell> cells(2, arb::cable_cell(m, {}, d));

    auto run_test = [&](bool coalesce_synapses) {
        cable1d_recipe rec(cells, coalesce_synapses);

        rec.add_probe(0, 0, cable_probe_point_state_cell{"expsyn", "g"});
        rec.add_probe(1, 0, cable_probe_point_state_cell{"expsyn", "g"});

        fvm_cell lcell(*ctx);
        auto fvm_info = lcell.initialize({0, 1}, rec);
        const auto& probe_map = fvm_info.probe_map;
        const auto& targets = fvm_info.target_handles;

        // Send an event to each expsyn synapse with a weight = target+100*cell_gid, and
        // integrate for a tiny time step.

        std::vector<deliverable_event> events;
        for (unsigned i: {0u, 1u}) {
            // Cells have the same number of targets, so the offset for cell 1 is exactly...
            cell_local_size_type cell_offset = i==0? 0: targets.size()/2;

            for (auto target_id: util::keys(expsyn_target_loc_map)) {
                deliverable_event ev{0., targets.at(target_id+cell_offset), float(target_id+100*i)};
                events.push_back(ev);
            }
        }
        (void)lcell.integrate(1e-5, 1e-5, events, {});

        // Independently get cv geometry to compute CV indices.

        cv_geometry geom = cv_geometry_from_ends(cells[0], policy.cv_boundary_points(cells[0]));
        append(geom, cv_geometry_from_ends(cells[1], policy.cv_boundary_points(cells[1])));

        ASSERT_EQ(2u, probe_map.size());
        for (unsigned i: {0u, 1u}) {
            const auto* h_ptr = std::get_if<fvm_probe_multi>(&probe_map.data_on({i, 0}).front().info);
            ASSERT_TRUE(h_ptr);

            const auto* m_ptr = std::get_if<std::vector<cable_probe_point_info>>(&h_ptr->metadata);
            ASSERT_TRUE(m_ptr);

            const fvm_probe_multi& h = *h_ptr;
            const std::vector<cable_probe_point_info> m = *m_ptr;

            ASSERT_EQ(h.raw_handles.size(), m.size());
            ASSERT_EQ(n_expsyn, m.size());

            std::vector<double> expected_coalesced_cv_value(geom.size());
            std::vector<double> expected_uncoalesced_value(targets.size());

            std::vector<double> target_cv(targets.size(), (unsigned)-1);
            std::unordered_map<fvm_size_type, unsigned> cv_expsyn_count;

            for (unsigned j = 0; j<n_expsyn; ++j) {
                ASSERT_EQ(1u, expsyn_target_loc_map.count(m[j].target));
                EXPECT_EQ(expsyn_target_loc_map.at(m[j].target), m[j].loc);

                auto cv = geom.location_cv(i, m[j].loc, cv_prefer::cv_nonempty);
                target_cv[j] = cv;
                ++cv_expsyn_count[cv];

                double event_weight = m[j].target+100*i;
                expected_uncoalesced_value[j] = event_weight;
                expected_coalesced_cv_value[cv] += event_weight;
            }

            for (unsigned j = 0; j<n_expsyn; ++j) {
                if (coalesce_synapses) {
                    EXPECT_EQ(cv_expsyn_count.at(target_cv[j]), m[j].multiplicity);
                }
                else {
                    EXPECT_EQ(1u, m[j].multiplicity);
                }
            }

            for (unsigned j = 0; j<n_expsyn; ++j) {
                double expected_value = coalesce_synapses?
                    expected_coalesced_cv_value[target_cv[j]]:
                    expected_uncoalesced_value[j];

                double value = deref(h.raw_handles[j]);

                EXPECT_NEAR(expected_value, value, 0.01); // g values will have decayed a little.
            }
        }
    };

    {
        SCOPED_TRACE("uncoalesced synapses");
        run_test(false);
    }

    {
        SCOPED_TRACE("coalesced synapses");
        run_test(true);
    }
}

template <typename Backend>
void run_ion_density_probe_test(const context& ctx) {
    using fvm_cell = typename backend_access<Backend>::fvm_cell;
    auto deref = [](const fvm_value_type* p) { return backend_access<Backend>::deref(p); };

    // Use test mechanism write_Xi_Xo to check ion concentration probes and
    // density mechanism state probes.

    auto cat = make_unit_test_catalogue();
    cat.derive("write_ca1", "write_Xi_Xo", {{"xi0", 1.25}, {"xo0", 1.5}, {"s0", 1.75}}, {{"x", "ca"}});
    cat.derive("write_ca2", "write_Xi_Xo", {{"xi0", 2.25}, {"xo0", 2.5}, {"s0", 2.75}}, {{"x", "ca"}});
    cat.derive("write_na3", "write_Xi_Xo", {{"xi0", 3.25}, {"xo0", 3.5}, {"s0", 3.75}}, {{"x", "na"}});

    // Simple constant diameter cable, 3 CVs.

    auto m = make_stick_morphology();
    decor d;
    d.set_default(cv_policy_fixed_per_branch(3));

    // Calcium ions everywhere, half written by write_ca1, half by write_ca2.
    // Sodium ions only on distal half.

    d.paint(mcable{0, 0., 0.5}, "write_ca1");
    d.paint(mcable{0, 0.5, 1.}, "write_ca2");
    d.paint(mcable{0, 0.5, 1.}, "write_na3");

    // Place probes in each CV.

    mlocation loc0{0, 0.1};
    mlocation loc1{0, 0.5};
    mlocation loc2{0, 0.9};

    cable1d_recipe rec(cable_cell(m, {}, d));
    rec.catalogue() = cat;

    // Probe (0, 0): ca internal on CV 0.
    rec.add_probe(0, 0, cable_probe_ion_int_concentration{loc0, "ca"});
    // Probe (0, 1): ca internal on CV 1.
    rec.add_probe(0, 0, cable_probe_ion_int_concentration{loc1, "ca"});
    // Probe (0, 2): ca internal on CV 2.
    rec.add_probe(0, 0, cable_probe_ion_int_concentration{loc2, "ca"});

    // Probe (0, 3): ca external on CV 0.
    rec.add_probe(0, 0, cable_probe_ion_ext_concentration{loc0, "ca"});
    // Probe (0, 4): ca external on CV 1.
    rec.add_probe(0, 0, cable_probe_ion_ext_concentration{loc1, "ca"});
    // Probe (0, 5): ca external on CV 2.
    rec.add_probe(0, 0, cable_probe_ion_ext_concentration{loc2, "ca"});

    // Probe (0, 6): na internal on CV 0.
    rec.add_probe(0, 0, cable_probe_ion_int_concentration{loc0, "na"});
    // Probe (0, 7): na internal on CV 2.
    rec.add_probe(0, 0, cable_probe_ion_int_concentration{loc2, "na"});

    // Probe (0, 8): write_ca2 state 's' in CV 0.
    rec.add_probe(0, 0, cable_probe_density_state{loc0, "write_ca2", "s"});
    // Probe (0, 9): write_ca2 state 's' in CV 1.
    rec.add_probe(0, 0, cable_probe_density_state{loc1, "write_ca2", "s"});
    // Probe (0, 10): write_ca2 state 's' in CV 2.
    rec.add_probe(0, 0, cable_probe_density_state{loc2, "write_ca2", "s"});

    // Probe (0, 11): na internal on whole cell.
    rec.add_probe(0, 0, cable_probe_ion_int_concentration_cell{"na"});
    // Probe (0, 12): ca external on whole cell.
    rec.add_probe(0, 0, cable_probe_ion_ext_concentration_cell{"ca"});

    fvm_cell lcell(*ctx);
    auto fvm_info = lcell.initialize({0}, rec);
    const auto& probe_map = fvm_info.probe_map;

    // Should be no sodium ion instantiated on CV 0, so probe (0, 6) should
    // have been silently discared. Similarly, write_ca2 is not instantiated on
    // CV 0, and so probe (0, 8) should have been discarded. All other probes
    // should be in the map.

    EXPECT_EQ(13u, rec.get_probes(0).size());
    EXPECT_EQ(11u, probe_map.size());

    auto get_probe_raw_handle = [&](cell_member_type x, unsigned i = 0) {
        return probe_map.data_on(x).front().raw_handle_range()[i];
    };

    probe_handle ca_int_cv0 = get_probe_raw_handle({0, 0});
    probe_handle ca_int_cv1 = get_probe_raw_handle({0, 1});
    probe_handle ca_int_cv2 = get_probe_raw_handle({0, 2});
    probe_handle ca_ext_cv0 = get_probe_raw_handle({0, 3});
    probe_handle ca_ext_cv1 = get_probe_raw_handle({0, 4});
    probe_handle ca_ext_cv2 = get_probe_raw_handle({0, 5});
    EXPECT_EQ(0u, probe_map.data.count({0, 6}));
    probe_handle na_int_cv2 = get_probe_raw_handle({0, 7});
    EXPECT_EQ(0u, probe_map.data.count({0, 8}));
    probe_handle write_ca2_s_cv1 = get_probe_raw_handle({0, 9});
    probe_handle write_ca2_s_cv2 = get_probe_raw_handle({0, 10});

    // Ion concentrations should have been written in initialization.
    // For CV 1, calcium concentration should be mean of the two values
    // from write_ca1 and write_ca2.

    EXPECT_EQ(1.25, deref(ca_int_cv0));
    EXPECT_DOUBLE_EQ((1.25+2.25)/2., deref(ca_int_cv1));
    EXPECT_EQ(2.25, deref(ca_int_cv2));

    EXPECT_EQ(1.5, deref(ca_ext_cv0));
    EXPECT_DOUBLE_EQ((1.5+2.5)/2., deref(ca_ext_cv1));
    EXPECT_EQ(2.5, deref(ca_ext_cv2));

    EXPECT_EQ(3.25, deref(na_int_cv2));

    // State variable in write_ca2 should be the same in both CV 1 and 2.
    // The raw handles should be different addresses, however.

    EXPECT_EQ(2.75, deref(write_ca2_s_cv1));
    EXPECT_EQ(2.75, deref(write_ca2_s_cv2));
    EXPECT_NE(write_ca2_s_cv1, write_ca2_s_cv2);

    // For the all cell sodium concentration probe, check that metadata
    // cables comprise a cable for all of CV 1 and a cable for all of CV 2:
    // half coverage of CV 1 by sodium ion-requiring mechanism implies
    // ion values are valid across whole CV.
    //
    // Implementation detail has that the cables (and corresponding raw handles) are
    // sorted by CV in the fvm_probe_weighted_multi object; this is assumed
    // below.

    auto* p_ptr = std::get_if<fvm_probe_multi>(&probe_map.data_on({0, 11}).front().info);
    ASSERT_TRUE(p_ptr);
    const fvm_probe_multi& na_int_all_info = *p_ptr;

    auto* m_ptr = std::get_if<mcable_list>(&na_int_all_info.metadata);
    ASSERT_TRUE(m_ptr);
    mcable_list na_int_all_metadata = *m_ptr;

    ASSERT_EQ(2u, na_int_all_metadata.size());
    ASSERT_EQ(2u, na_int_all_info.raw_handles.size());

    EXPECT_DOUBLE_EQ(1./3., na_int_all_metadata[0].prox_pos);
    EXPECT_DOUBLE_EQ(2./3., na_int_all_metadata[0].dist_pos);
    EXPECT_DOUBLE_EQ(2./3., na_int_all_metadata[1].prox_pos);
    EXPECT_DOUBLE_EQ(1.,    na_int_all_metadata[1].dist_pos);
    EXPECT_EQ(na_int_cv2,   na_int_all_info.raw_handles[1]);
    EXPECT_EQ(na_int_cv2-1, na_int_all_info.raw_handles[0]);

    p_ptr = std::get_if<fvm_probe_multi>(&probe_map.data_on({0, 12}).front().info);
    ASSERT_TRUE(p_ptr);
    const fvm_probe_multi& ca_ext_all_info = *p_ptr;

    m_ptr = std::get_if<mcable_list>(&ca_ext_all_info.metadata);
    ASSERT_TRUE(m_ptr);
    mcable_list ca_ext_all_metadata = *m_ptr;

    ASSERT_EQ(3u, ca_ext_all_metadata.size());
    ASSERT_EQ(3u, ca_ext_all_info.raw_handles.size());

    EXPECT_DOUBLE_EQ(0.,    ca_ext_all_metadata[0].prox_pos);
    EXPECT_DOUBLE_EQ(1./3., ca_ext_all_metadata[1].prox_pos);
    EXPECT_DOUBLE_EQ(2./3., ca_ext_all_metadata[2].prox_pos);
    EXPECT_EQ(ca_ext_cv0,   ca_ext_all_info.raw_handles[0]);
    EXPECT_EQ(ca_ext_cv1,   ca_ext_all_info.raw_handles[1]);
    EXPECT_EQ(ca_ext_cv2,   ca_ext_all_info.raw_handles[2]);
}

template <typename Backend>
void run_partial_density_probe_test(const context& ctx) {
    using fvm_cell = typename backend_access<Backend>::fvm_cell;
    auto deref = [](const fvm_value_type* p) { return backend_access<Backend>::deref(p); };

    // Use test mechanism param_as_state to query averaged state values in CVs with
    // partial coverage by the mechanism.

    auto cat = make_unit_test_catalogue();

    cable_cell cells[2];

    // Each cell is a simple constant diameter cable, with 3 CVs each.

    auto m = make_stick_morphology();

    decor d0, d1;
    d0.set_default(cv_policy_fixed_per_branch(3));
    d1.set_default(cv_policy_fixed_per_branch(3));

    // Paint the mechanism on every second 10% interval of each cell.
    // Expected values on a CV are the weighted mean of the parameter values
    // over the intersections of the support and the CV.
    //
    // Cell 0:       [0.0, 0.1], [0.2, 0.3] [0.4, 0.5] [0.6, 0.7] [0.8, 0.9].
    // Param values:     2           3          4          5          6
    // Expected values:
    //    CV 0:   2.5
    //    CV 1:   4.4
    //    CV 2:   5.75
    //
    // Cell 1:       [0.1, 0.2], [0.3, 0.4] [0.5, 0.6] [0.7, 0.8] [0.9, 1.0].
    // Param values:     7           8          9         10         11
    // Expected values:
    //    CV 3:   7.25
    //    CV 4:   8.6
    //    CV 5:  10.5

    auto mk_mech = [](double param) { return mechanism_desc("param_as_state").set("p", param); };

    d0.paint(mcable{0, 0.0, 0.1}, mk_mech(2));
    d0.paint(mcable{0, 0.2, 0.3}, mk_mech(3));
    d0.paint(mcable{0, 0.4, 0.5}, mk_mech(4));
    d0.paint(mcable{0, 0.6, 0.7}, mk_mech(5));
    d0.paint(mcable{0, 0.8, 0.9}, mk_mech(6));

    d1.paint(mcable{0, 0.1, 0.2}, mk_mech(7));
    d1.paint(mcable{0, 0.3, 0.4}, mk_mech(8));
    d1.paint(mcable{0, 0.5, 0.6}, mk_mech(9));
    d1.paint(mcable{0, 0.7, 0.8}, mk_mech(10));
    d1.paint(mcable{0, 0.9, 1.0}, mk_mech(11));

    cells[0] = cable_cell(m, {}, d0);
    cells[1] = cable_cell(m, {}, d1);

    // Place probes in the middle of each 10% interval, i.e. at 0.05, 0.15, etc.
    struct test_probe {
        double pos;
        double expected[2]; // Expected value in each cell, NAN => n/a
    } test_probes[] = {
        { 0.05, { 2.5,   NAN  }},  // CV 0, 3
        { 0.15, { NAN,   7.25 }},  // CV 0, 3
        { 0.25, { 2.5,   NAN  }},  // CV 0, 3
        { 0.35, { NAN,   8.6  }},  // CV 1, 4
        { 0.45, { 4.4,   NAN  }},  // CV 1, 4
        { 0.55, { NAN,   8.6  }},  // CV 1, 4
        { 0.65, { 4.4,   NAN  }},  // CV 1, 4
        { 0.75, { NAN,  10.5  }},  // CV 2, 5
        { 0.85, { 5.75,  NAN  }},  // CV 2, 5
        { 0.95, { NAN,  10.5  }}   // CV 2, 5
    };

    cable1d_recipe rec(cells);
    rec.catalogue() = cat;

    for (auto tp: test_probes) {
        rec.add_probe(0, 0, cable_probe_density_state{mlocation{0, tp.pos}, "param_as_state", "s"});
        rec.add_probe(1, 0, cable_probe_density_state{mlocation{0, tp.pos}, "param_as_state", "s"});
    }

    fvm_cell lcell(*ctx);
    auto fvm_info = lcell.initialize({0, 1}, rec);
    const auto& probe_map = fvm_info.probe_map;

    // There should be 10 probes on each cell, but only 10 in total in the probe map,
    // as only those probes that are in the mechanism support should have an entry.

    EXPECT_EQ(10u, rec.get_probes(0).size());
    EXPECT_EQ(10u, rec.get_probes(1).size());
    EXPECT_EQ(10u, probe_map.size());

    auto get_probe_raw_handle = [&](cell_member_type x, unsigned i = 0) {
        return probe_map.data_on(x).front().raw_handle_range()[i];
    };

    // Check probe values against expected values.

    cell_lid_type probe_lid = 0;
    for (auto tp: test_probes) {
        for (cell_gid_type gid: {0, 1}) {
            cell_member_type probe_id{gid, probe_lid};
            if (std::isnan(tp.expected[gid])) {
                EXPECT_EQ(0u, probe_map.data.count(probe_id));
            }
            else {
                probe_handle h = get_probe_raw_handle(probe_id);
                EXPECT_DOUBLE_EQ(tp.expected[gid], deref(h));
            }
        }
        ++probe_lid;
    }
}

template <typename Backend>
void run_axial_and_ion_current_sampled_probe_test(const context& ctx) {
    // On a passive cable in steady-state, the capacitive membrane current will be zero,
    // and the axial currents should balance the stimulus and ionic membrane currents in any CV.
    //
    // Membrane currents not associated with an ion are not visible, so we use a custom
    // mechanism 'ca_linear' that presents a passive, constant conductance membrane current
    // as a calcium ion current.

    auto cat = make_unit_test_catalogue();

    // Cell is a tapered cable with 3 CVs.

    auto m = make_stick_morphology();
    arb::decor d;

    const unsigned n_cv = 3;
    cv_policy policy = cv_policy_fixed_per_branch(n_cv);
    d.set_default(policy);

    d.place(mlocation{0, 0}, i_clamp(0.3), "clamp");

    // The time constant will be membrane capacitance / membrane conductance.
    // For τ = 0.1 ms, set conductance to 0.01 S/cm² and membrance capacitance
    // to 0.01 F/m².

    d.paint(reg::all(), mechanism_desc("ca_linear").set("g", 0.01)); // [S/cm²]
    d.set_default(membrane_capacitance{0.01}); // [F/m²]
    const double tau = 0.1; // [ms]

    cable1d_recipe rec(cable_cell(m, {}, d));
    rec.catalogue() = cat;

    cable_cell cell(m, {}, d);

    // Place axial current probes at CV boundaries and make cell-wide probes for
    // total ionic membrane current and stimulus currents.

    mlocation_list cv_boundaries;
    util::assign(cv_boundaries,
        util::filter(thingify(policy.cv_boundary_points(cell), cell.provider()),
            [](mlocation loc) { return loc.pos!=0 && loc.pos!=1; }));

    ASSERT_EQ(n_cv-1, cv_boundaries.size());
    const unsigned n_axial_probe = n_cv-1;

    for (mlocation loc: cv_boundaries) {
        // Probe ids will be be (0, 0), (0, 1), etc.
        rec.add_probe(0, 0, cable_probe_axial_current{loc});
    }

    // Use tags 1 and 2 for the whole-cell probes.
    rec.add_probe(0, 1, cable_probe_total_ion_current_cell{});
    rec.add_probe(0, 2, cable_probe_stimulus_current_cell{});

    partition_hint_map phints = {
       {cell_kind::cable, {partition_hint::max_size, partition_hint::max_size, true}}
    };
    simulation sim(rec, partition_load_balance(rec, ctx, phints), ctx);

    // Take a sample at 20 tau, and run sim for just a bit longer.

    std::vector<double> i_axial(n_axial_probe);
    std::vector<double> i_memb, i_stim;

    sim.add_sampler(all_probes, explicit_schedule({20*tau}),
        [&](probe_metadata pm, std::size_t n_sample, const sample_record* samples)
        {
            // Expect exactly one sample.
            ASSERT_EQ(1u, n_sample);

            if (pm.tag!=0) { // (whole cell probe)
                const mcable_list* m = any_cast<const mcable_list*>(pm.meta);
                ASSERT_NE(nullptr, m);
                // Metadata should comprise one cable per CV.
                ASSERT_EQ(n_cv, m->size());

                const cable_sample_range* s = any_cast<const cable_sample_range*>(samples[0].data);
                ASSERT_NE(nullptr, s);
                ASSERT_EQ(s->first+n_cv, s->second);

                auto* i_var = pm.tag==1? &i_memb: &i_stim;
                for (const double* p = s->first; p!=s->second; ++p) {
                    i_var->push_back(*p);
                }
            }
            else { // axial current probe
                // Probe id tells us which axial current this is.
                ASSERT_LT(pm.id.index, n_axial_probe);

                const mlocation* m = any_cast<const mlocation*>(pm.meta);
                ASSERT_NE(nullptr, m);

                const double* s = any_cast<const double*>(samples[0].data);
                ASSERT_NE(nullptr, s);

                i_axial.at(pm.id.index) = *s;
            }
        });

    const double dt = 0.025; // [ms]
    sim.run(20*tau+dt, dt);

    ASSERT_EQ(n_cv, i_memb.size());

    for (unsigned i = 0; i<n_cv; ++i) {
        // Axial currents are in the distal (increasing CV index) direction,
        // while membrane currents are from intra- to extra-cellular medium.
        //
        // Net outward flux from CV is axial current on distal side, minus
        // axial current on proximal side, plus membrane current, and should
        // sum to zero.

        double net_axial_flux = i<n_axial_probe? i_axial[i]: 0;
        net_axial_flux -= i>0? i_axial[i-1]: 0;

        EXPECT_TRUE(testing::near_relative(net_axial_flux, -i_memb[i]-i_stim[i], 1e-6));
    }
}

// Run given cells taking samples from the provied probes on one of the cells.
//
// Use the default mechanism catalogue augmented by unit test specific mechanisms.
// (Timestep fixed at 0.025 ms).

template <typename SampleData, typename SampleMeta = void>
auto run_simple_samplers(
    const arb::context& ctx,
    double t_end,
    const std::vector<cable_cell>& cells,
    cell_gid_type probe_cell,
    const std::vector<std::any>& probe_addrs,
    const std::vector<double>& when)
{
    cable1d_recipe rec(cells, false);
    rec.catalogue() = make_unit_test_catalogue(global_default_catalogue());
    unsigned n_probe = probe_addrs.size();

    for (auto& addr: probe_addrs) {
        rec.add_probe(probe_cell, 0, addr);
    }

    partition_hint_map phints = {
       {cell_kind::cable, {partition_hint::max_size, partition_hint::max_size, true}}
    };
    simulation sim(rec, partition_load_balance(rec, ctx, phints), ctx);

    std::vector<trace_vector<SampleData, SampleMeta>> traces(n_probe);
    for (unsigned i = 0; i<n_probe; ++i) {
        sim.add_sampler(one_probe({probe_cell, i}), explicit_schedule(when), make_simple_sampler(traces[i]));
    }

    sim.run(t_end, 0.025);
    return traces;
}

template <typename SampleData, typename SampleMeta = void>
auto run_simple_sampler(
    const arb::context& ctx,
    double t_end,
    const std::vector<cable_cell>& cells,
    cell_gid_type probe_cell,
    const std::any& probe_addr,
    const std::vector<double>& when)
{
    return run_simple_samplers<SampleData, SampleMeta>(ctx, t_end, cells, probe_cell, {probe_addr}, when).at(0);
}

template <typename Backend>
void run_multi_probe_test(const context& ctx) {
    // Construct and run thorugh simple sampler a probe defined over
    // cell terminal points; check metadata and values.

    // m_mlt_b6 has terminal branches 1, 2, 4, and 5.
    auto m = common_morphology::m_mlt_b6;
    decor d;

    // Paint mechanism on branches 1, 2, and 5, omitting branch 4.
    d.paint(reg::branch(1), mechanism_desc("param_as_state").set("p", 10.));
    d.paint(reg::branch(2), mechanism_desc("param_as_state").set("p", 20.));
    d.paint(reg::branch(5), mechanism_desc("param_as_state").set("p", 50.));

    auto tracev = run_simple_sampler<double, mlocation>(ctx, 0.1, {cable_cell{m, {}, d}}, 0, cable_probe_density_state{ls::terminal(), "param_as_state", "s"}, {0.});

    // Expect to have received a sample on each of the terminals of branches 1, 2, and 5.
    ASSERT_EQ(3u, tracev.size());

    std::vector<std::pair<mlocation, double>> vals;
    for (auto& trace: tracev) {
        ASSERT_EQ(1u, trace.size());
        vals.push_back({trace.meta, trace[0].v});
    }

    util::sort(vals);
    EXPECT_EQ((mlocation{1, 1.}), vals[0].first);
    EXPECT_EQ((mlocation{2, 1.}), vals[1].first);
    EXPECT_EQ((mlocation{5, 1.}), vals[2].first);
    EXPECT_EQ(10., vals[0].second);
    EXPECT_EQ(20., vals[1].second);
    EXPECT_EQ(50., vals[2].second);
}

template <typename Backend>
void run_v_sampled_probe_test(const context& ctx) {
    soma_cell_builder builder(12.6157/2.0);
    builder.add_branch(0, 200, 1.0/2, 1.0/2, 1, "dend");
    builder.add_branch(0, 200, 1.0/2, 1.0/2, 1, "dend");

    auto bs = builder.make_cell();
    bs.decorations.set_default(cv_policy_fixed_per_branch(1));
    auto d0 = bs.decorations;
    auto d1 = bs.decorations;

    // Add stims, up to 0.5 ms on cell 0, up to 1.0 ms on cell 1, so that
    // samples at the same point on each cell will give the same value at
    // 0.3 ms, but different at 0.6 ms.

    d0.place(mlocation{1, 1}, i_clamp::box(0, 0.5, 1.), "clamp0");
    d1.place(mlocation{1, 1}, i_clamp::box(0, 1.0, 1.), "clamp1");
    mlocation probe_loc{1, 0.2};

    std::vector<cable_cell> cells = {{bs.morph, bs.labels, d0}, {bs.morph, bs.labels, d1}};

    const double t_end = 1.; // [ms]
    std::vector<double> when = {0.3, 0.6}; // Sample at 0.3 and 0.6 ms.

    auto trace0 = run_simple_sampler<double, mlocation>(ctx, t_end, cells, 0, cable_probe_membrane_voltage{probe_loc}, when).at(0);
    ASSERT_TRUE(trace0);

    EXPECT_EQ(probe_loc, trace0.meta);
    EXPECT_EQ(2u, trace0.size());

    auto trace1 = run_simple_sampler<double, mlocation>(ctx, t_end, cells, 1, cable_probe_membrane_voltage{probe_loc}, when).at(0);
    ASSERT_TRUE(trace1);

    EXPECT_EQ(probe_loc, trace1.meta);
    EXPECT_EQ(2u, trace1.size());

    EXPECT_EQ(trace0[0].t, trace1[0].t);
    EXPECT_EQ(trace0[0].v, trace1[0].v);

    EXPECT_EQ(trace0[1].t, trace1[1].t);
    EXPECT_NE(trace0[1].v, trace1[1].v);
}


template <typename Backend>
void run_total_current_probe_test(const context& ctx) {
    // Model two passive Y-shaped cells with a similar but not identical
    // time constant τ.
    //
    // Sample each cell's total membrane currents at the same time,
    // approximately equal to the time constants τ.
    //
    // Net current flux in each cell should be zero, but currents should
    // differ between the cells.

    auto m = make_y_morphology();
    decor d0;

    const unsigned n_cv_per_branch = 3;
    const unsigned n_branch = 3;

    // The time constant will be membrane capacitance / membrane conductance.
    // For τ = 0.1 ms, set conductance to 0.01 S/cm² and membrance capacitance
    // to 0.01 F/m².

    const double tau = 0.1;     // [ms]
    d0.place(mlocation{0, 0}, i_clamp(0.3), "clamp0");

    d0.paint(reg::all(), mechanism_desc("ca_linear").set("g", 0.01)); // [S/cm²]
    d0.set_default(membrane_capacitance{0.01}); // [F/m²]
    // Tweak membrane capacitance on cells[1] so as to change dynamics a bit.
    auto d1 = d0;
    d1.set_default(membrane_capacitance{0.009}); // [F/m²]

    // We'll run each set of tests twice: once with a trivial (zero-volume) CV
    // at the fork points, and once with a non-trivial CV centred on the fork
    // point.

    trace_data<std::vector<double>, mcable_list> traces[2];
    trace_data<std::vector<double>, mcable_list> ion_traces[2];
    trace_data<std::vector<double>, mcable_list> stim_traces[2];

    // Run the cells sampling at τ and 20τ for both total membrane
    // current and total membrane ionic current.

    auto run_cells = [&](bool interior_forks) {
        auto flags = interior_forks? cv_policy_flag::interior_forks: cv_policy_flag::none;
        cv_policy policy = cv_policy_fixed_per_branch(n_cv_per_branch, flags);
        d0.set_default(policy);
        d1.set_default(policy);
        std::vector<cable_cell> cells = {{m, {}, d0}, {m, {}, d1}};


        for (unsigned i = 0; i<2; ++i) {
            SCOPED_TRACE(i);

            const double t_end = 21*tau; // [ms]

            traces[i] = run_simple_sampler<std::vector<double>, mcable_list>(ctx, t_end, cells, i,
                    cable_probe_total_current_cell{}, {tau, 20*tau}).at(0);

            ion_traces[i] = run_simple_sampler<std::vector<double>, mcable_list>(ctx, t_end, cells, i,
                    cable_probe_total_ion_current_cell{}, {tau, 20*tau}).at(0);

            stim_traces[i] = run_simple_sampler<std::vector<double>, mcable_list>(ctx, t_end, cells, i,
                    cable_probe_stimulus_current_cell{}, {tau, 20*tau}).at(0);

            ASSERT_EQ(2u, traces[i].size());
            ASSERT_EQ(2u, ion_traces[i].size());

            // Check metadata size:
            //  * With trivial forks, should have n_cv_per_branch*n_branch cables; zero-length cables
            //    associated with the trivial CVs at forks should not be included.
            //  * With nontrivial forks, we'll have an extra cable at the head of each branch, which
            //    for all but the root branch will be a component cable of the CV on the fork.
            //
            // Total membrane current and total ionic mebrane current should have the
            // same support and same metadata.

            ASSERT_EQ((n_cv_per_branch+(int)interior_forks)*n_branch, traces[i].meta.size());
            ASSERT_EQ(ion_traces[i].meta, traces[i].meta);
            EXPECT_EQ(ion_traces[i][0].v.size(), traces[i][0].v.size());
            EXPECT_EQ(ion_traces[i][1].v.size(), traces[i][1].v.size());

            // Check total membrane currents + stimulus currents are individually non-zero, but sum is,
            // both at t=τ (j=0) and t=20τ (j=1).

            ASSERT_EQ(ion_traces[i].meta, stim_traces[i].meta);

            for (unsigned j: {0u, 1u}) {
                double max_abs_current = 0;
                double sum_current = 0;
                for (auto k: util::count_along(traces[i].meta)) {
                    double current = traces[i][j].v[k] + stim_traces[i][j].v[k];

                    EXPECT_NE(0.0, current);
                    max_abs_current = std::max(max_abs_current, std::abs(current));
                    sum_current += current;
                }

                ASSERT_NEAR(0.0, sum_current, 1e-6*max_abs_current);
            }

            // Confirm that total (non-stim) and ion currents differ at τ but are close at 20τ.

            for (unsigned k = 0; k<traces[i].size(); ++k) {
                const double rtol_large = 1e-3;
                EXPECT_FALSE(testing::near_relative(traces[i][0].v.at(k), ion_traces[i][0].v.at(k), rtol_large));
            }

            for (unsigned k = 0; k<traces[i].size(); ++k) {
                const double rtol_small = 1e-6;
                EXPECT_TRUE(testing::near_relative(traces[i][1].v.at(k), ion_traces[i][1].v.at(k), rtol_small));
            }

        }

        // Total membrane currents should differ between the two cells at t=τ.

        for (unsigned k = 0; k<traces[0][0].v.size(); ++k) {
            EXPECT_NE(traces[0][0].v.at(k), traces[1][0].v.at(k));
        }
    };

    {
        SCOPED_TRACE("trivial fork CV");
        run_cells(false);
    }

    {
        SCOPED_TRACE("non-trival fork CV");
        run_cells(true);
    }
}


template <typename Backend>
void run_stimulus_probe_test(const context& ctx) {
    // Model two simple stick cable cells, 3 CVs each, and stimuli on cell 0, cv 1
    // and cell 1, cv 2. Run both cells in the same cell group.

    const double stim_until = 1.; // [ms]
    auto m = make_stick_morphology();
    cv_policy policy = cv_policy_fixed_per_branch(3);

    decor d0, d1;
    d0.set_default(policy);
    d0.place(mlocation{0, 0.5}, i_clamp::box(0., stim_until, 10.), "clamp0");
    d0.place(mlocation{0, 0.5}, i_clamp::box(0., stim_until, 20.), "clamp1");
    double expected_stim0 = 30;

    d1.set_default(policy);
    d1.place(mlocation{0, 1}, i_clamp::box(0., stim_until, 30.), "clamp0");
    d1.place(mlocation{0, 1}, i_clamp::box(0., stim_until, -10.), "clamp1");
    double expected_stim1 = 20;

    std::vector<cable_cell> cells = {{m, {}, d0}, {m, {}, d1}};

    // Sample the cells during the stimulus, and after.

    trace_data<std::vector<double>, mcable_list> traces[2];

    for (unsigned i: {0u, 1u}) {
        traces[i] = run_simple_sampler<std::vector<double>, mcable_list>(ctx, 2.5*stim_until, cells, i,
                cable_probe_stimulus_current_cell{}, {stim_until/2, 2*stim_until}).at(0);

        ASSERT_EQ(3u, traces[i].meta.size());
        for (unsigned cv: {0u, 1u, 2u}) {
            ASSERT_EQ(2u, traces[i].size());
        }
    }

    // Every sample in each trace should be zero _except_ the first sample for cell 0, cv 1
    // and the first sample for cell 1, cv 2.

    EXPECT_EQ((std::vector<double>{0, expected_stim0, 0}), traces[0][0]);
    EXPECT_EQ((std::vector<double>{0, 0, expected_stim1}), traces[1][0]);
    EXPECT_EQ((std::vector<double>(3)), traces[0][1]);
    EXPECT_EQ((std::vector<double>(3)), traces[1][1]);
}

template <typename Backend>
void run_exact_sampling_probe_test(const context& ctx) {
    // As the exact sampling implementation interacts with the event delivery
    // implementation within in cable cell groups, construct a somewhat
    // elaborate model with 4 cells and a gap junction between cell 1 and 3.

    struct adhoc_recipe: recipe {
        std::vector<cable_cell> cells_;
        cable_cell_global_properties gprop_;

        adhoc_recipe() {
            gprop_.default_parameters = neuron_parameter_defaults;

            soma_cell_builder builder(12.6157/2.0);
            builder.add_branch(0, 200, 1.0/2, 1.0/2, 1, "dend");
            builder.add_branch(0, 200, 1.0/2, 1.0/2, 1, "dend");

            std::vector<cable_cell_description> cd;
            cd.assign(4, builder.make_cell());

            cd[0].decorations.place(mlocation{1, 0.1}, "expsyn", "syn");
            cd[1].decorations.place(mlocation{1, 0.1}, "exp2syn", "syn");
            cd[2].decorations.place(mlocation{1, 0.9}, "expsyn", "syn");
            cd[3].decorations.place(mlocation{1, 0.9}, "exp2syn", "syn");

            cd[1].decorations.place(mlocation{1, 0.2}, gap_junction_site{}, "gj");
            cd[3].decorations.place(mlocation{1, 0.2}, gap_junction_site{}, "gj");

            for (auto& d: cd) cells_.push_back(d);
        }

        cell_size_type num_cells() const override { return cells_.size(); }

        util::unique_any get_cell_description(cell_gid_type gid) const override {
            return cells_.at(gid);
        }

        cell_kind get_cell_kind(cell_gid_type) const override {
            return cell_kind::cable;
        }

        std::vector<probe_info> get_probes(cell_gid_type) const override {
            return {cable_probe_membrane_voltage{mlocation{1, 0.5}}};
        }

        std::vector<gap_junction_connection> gap_junctions_on(cell_gid_type gid) const override {
            switch (gid) {
            case 1:
                return {gap_junction_connection({3, "gj", lid_selection_policy::assert_univalent}, {"gj", lid_selection_policy::assert_univalent}, 1.)};
            case 3:
                return {gap_junction_connection({1, "gj", lid_selection_policy::assert_univalent}, {"gj", lid_selection_policy::assert_univalent}, 1.)};
            default:
                return {};
            }
        }

        std::vector<event_generator> event_generators(cell_gid_type gid) const override {
            // Send a single event to cell i at 0.1*i milliseconds.
            explicit_generator::lse_vector spikes = {{{"syn"}, 0.1*gid, 1.f}};
            return {explicit_generator(spikes)};
        }

        std::any get_global_properties(cell_kind k) const override {
            return k==cell_kind::cable? gprop_: std::any{};
        }
    };

    // Check two things:
    // 1. Membrane voltage is similar with and without exact sampling.
    // 2. Sample times are in fact exact with exact sampling.

    std::vector<trace_vector<double>> lax_traces(4), exact_traces(4);

    const double max_dt = 0.001;
    const double t_end = 1.;
    std::vector<time_type> sched_times{1./7., 3./7., 4./7., 6./7.};
    schedule sample_sched = explicit_schedule(sched_times);

    adhoc_recipe rec;
    unsigned n_cell = rec.num_cells();
    unsigned n_sample_time = sched_times.size();

    partition_hint_map phints = {
       {cell_kind::cable, {partition_hint::max_size, partition_hint::max_size, true}}
    };
    domain_decomposition one_cell_group = partition_load_balance(rec, ctx, phints);

    simulation lax_sim(rec, one_cell_group, ctx);
    for (unsigned i = 0; i<n_cell; ++i) {
        lax_sim.add_sampler(one_probe({i, 0}), sample_sched, make_simple_sampler(lax_traces.at(i)), sampling_policy::lax);
    }
    lax_sim.run(t_end, max_dt);

    simulation exact_sim(rec, one_cell_group, ctx);
    for (unsigned i = 0; i<n_cell; ++i) {
        exact_sim.add_sampler(one_probe({i, 0}), sample_sched, make_simple_sampler(exact_traces.at(i)), sampling_policy::exact);
    }
    exact_sim.run(t_end, max_dt);

    for (unsigned i = 0; i<n_cell; ++i) {
        ASSERT_EQ(1u, lax_traces.at(i).size());
        ASSERT_EQ(n_sample_time, lax_traces.at(i).at(0).size());
        ASSERT_EQ(1u, exact_traces.at(i).size());
        ASSERT_EQ(n_sample_time, exact_traces.at(i).at(0).size());
    }

    for (unsigned i = 0; i<n_cell; ++i) {
        auto& lax_trace = lax_traces[i][0];
        auto& exact_trace = exact_traces[i][0];

        for (unsigned j = 0; j<n_sample_time; ++j) {
            EXPECT_NE(sched_times.at(j), lax_trace.at(j).t);
            EXPECT_EQ(sched_times.at(j), exact_trace.at(j).t);

            EXPECT_TRUE(testing::near_relative(lax_trace.at(j).v, exact_trace.at(j).v, 0.01));
        }
    }
}

// Generate unit tests multicore_X and gpu_X for each entry X in PROBE_TESTS,
// which establish the appropriate arbor context and then call run_X_probe_test.

#undef PROBE_TESTS
#define PROBE_TESTS \
    v_i, v_cell, v_sampled, expsyn_g, expsyn_g_cell, ion_density, \
    axial_and_ion_current_sampled, partial_density, exact_sampling, \
    multi, total_current

#undef RUN_MULTICORE
#define RUN_MULTICORE(x) \
TEST(probe, multicore_##x) { \
    context ctx = make_context(); \
    run_##x##_probe_test<multicore::backend>(ctx); \
}

ARB_PP_FOREACH(RUN_MULTICORE, PROBE_TESTS)

#ifdef ARB_GPU_ENABLED

#undef RUN_GPU
#define RUN_GPU(x) \
TEST(probe, gpu_##x) { \
    context ctx = make_context(proc_allocation{1, arbenv::default_gpu()}); \
    if (has_gpu(ctx)) { \
        run_##x##_probe_test<gpu::backend>(ctx); \
    } \
}

ARB_PP_FOREACH(RUN_GPU, PROBE_TESTS)

#endif // def ARB_GPU_ENABLED

// Test simulator `get_probe_metadata` interface.
// (No need to run this on GPU back-end as well.)

TEST(probe, get_probe_metadata) {
    // Reuse multiprobe test set-up to confirm simulator::get_probe_metadata returns
    // correct vector of metadata.

    auto m = common_morphology::m_mlt_b6;
    decor d;

    // Paint mechanism on branches 1, 2, and 5, omitting branch 4.
    d.paint(reg::branch(1), mechanism_desc("param_as_state").set("p", 10.));
    d.paint(reg::branch(2), mechanism_desc("param_as_state").set("p", 20.));
    d.paint(reg::branch(5), mechanism_desc("param_as_state").set("p", 50.));

    cable1d_recipe rec(cable_cell{m, {}, d}, false);
    rec.catalogue() = make_unit_test_catalogue(global_default_catalogue());
    rec.add_probe(0, 7, cable_probe_density_state{ls::terminal(), "param_as_state", "s"});

    context ctx = make_context();
    partition_hint_map phints = {
       {cell_kind::cable, {partition_hint::max_size, partition_hint::max_size, true}}
    };
    simulation sim(rec, partition_load_balance(rec, ctx, phints), ctx);

    std::vector<probe_metadata> mm = sim.get_probe_metadata({0, 0});
    ASSERT_EQ(3u, mm.size());

    EXPECT_EQ((cell_member_type{0, 0}), mm[0].id);
    EXPECT_EQ((cell_member_type{0, 0}), mm[1].id);
    EXPECT_EQ((cell_member_type{0, 0}), mm[2].id);

    EXPECT_EQ(0u, mm[0].index);
    EXPECT_EQ(1u, mm[1].index);
    EXPECT_EQ(2u, mm[2].index);

    EXPECT_EQ(7, mm[0].tag);
    EXPECT_EQ(7, mm[1].tag);
    EXPECT_EQ(7, mm[2].tag);

    const mlocation* l0 = any_cast<const mlocation*>(mm[0].meta);
    const mlocation* l1 = any_cast<const mlocation*>(mm[1].meta);
    const mlocation* l2 = any_cast<const mlocation*>(mm[2].meta);

    ASSERT_TRUE(l0);
    ASSERT_TRUE(l1);
    ASSERT_TRUE(l2);

    std::vector<mlocation> locs = {*l0, *l1, *l2};
    util::sort(locs);
    EXPECT_EQ((mlocation{1, 1.}), locs[0]);
    EXPECT_EQ((mlocation{2, 1.}), locs[1]);
    EXPECT_EQ((mlocation{5, 1.}), locs[2]);
}
