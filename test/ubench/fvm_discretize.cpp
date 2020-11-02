// Discretization can be too slow: use benchmark infrastrucure to hone in
// on the issue.

#include <fstream>

#include <arbor/cable_cell.hpp>
#include <arbor/morph/morphology.hpp>

#include <arborio/swcio.hpp>

#include <benchmark/benchmark.h>

#include "event_queue.hpp"
#include "fvm_layout.hpp"

#ifndef DATADIR
#define DATADIR "."
#endif

#define STRING(s) #s

#undef SWCFILE
#define SWCFILE STRING(DATADIR) "/motoneuron.swc"

using namespace arb;

arb::morphology from_swc(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("could not open "+path);

    return morphology(arborio::load_swc_arbor(arborio::parse_swc(in, arborio::swc_mode::strict)));
}

void run_cv_geom(benchmark::State& state) {
    auto gdflt = neuron_parameter_defaults;
    const std::size_t ncv_per_branch = state.range(0);

    cable_cell c(from_swc(SWCFILE));
    auto ends = cv_policy_fixed_per_branch(ncv_per_branch).cv_boundary_points(c);

    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(cv_geometry_from_ends(c, ends));
    }
}

void run_cv_geom_every_segment(benchmark::State& state) {
    auto gdflt = neuron_parameter_defaults;

    cable_cell c(from_swc(SWCFILE));
    auto ends = cv_policy_every_segment().cv_boundary_points(c);

    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(cv_geometry_from_ends(c, ends));
    }
}

void run_cv_geom_explicit(benchmark::State& state) {
    auto gdflt = neuron_parameter_defaults;

    cable_cell c(from_swc(SWCFILE));

    while (state.KeepRunning()) {
        auto ends = cv_policy_every_segment().cv_boundary_points(c);
        auto ends2 = cv_policy_explicit(std::move(ends)).cv_boundary_points(c);

        benchmark::DoNotOptimize(cv_geometry_from_ends(c, ends2));
    }
}

void run_discretize(benchmark::State& state) {
    auto gdflt = neuron_parameter_defaults;
    const std::size_t ncv_per_branch = state.range(0);

    cable_cell c(from_swc(SWCFILE));
    c.default_parameters.discretization = cv_policy_fixed_per_branch(ncv_per_branch);

    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(fvm_cv_discretize(c, gdflt));
    }
}

void run_discretize_every_segment(benchmark::State& state) {
    auto gdflt = neuron_parameter_defaults;

    cable_cell c(from_swc(SWCFILE));
    c.default_parameters.discretization = cv_policy_every_segment();

    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(fvm_cv_discretize(c, gdflt));
    }
}

void run_discretize_explicit(benchmark::State& state) {
    auto gdflt = neuron_parameter_defaults;

    cable_cell c(from_swc(SWCFILE));

    while (state.KeepRunning()) {
        auto ends = cv_policy_every_segment().cv_boundary_points(c);
        c.default_parameters.discretization = cv_policy_explicit(std::move(ends));

        benchmark::DoNotOptimize(fvm_cv_discretize(c, gdflt));
    }
}

BENCHMARK(run_cv_geom)->RangeMultiplier(2)->Range(1,32)->Unit(benchmark::kMicrosecond);
BENCHMARK(run_discretize)->RangeMultiplier(2)->Range(1,32)->Unit(benchmark::kMicrosecond);

BENCHMARK(run_cv_geom_every_segment)->Unit(benchmark::kMillisecond);
BENCHMARK(run_discretize_every_segment)->Unit(benchmark::kMillisecond);

BENCHMARK(run_cv_geom_explicit)->Unit(benchmark::kMillisecond);
BENCHMARK(run_discretize_explicit)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
