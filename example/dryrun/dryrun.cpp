/*
 * A miniapp that demonstrates how to use dry_run mode
 *
 */

#include <any>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <nlohmann/json.hpp>

#include <arborio/label_parse.hpp>

#include <arbor/common_types.hpp>
#include <arbor/context.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/profile/meter_manager.hpp>
#include <arbor/profile/profiler.hpp>
#include <arbor/simple_sampler.hpp>
#include <arbor/simulation.hpp>
#include <arbor/symmetric_recipe.hpp>
#include <arbor/recipe.hpp>
#include <arbor/version.hpp>

#include <sup/ioutil.hpp>
#include <sup/json_meter.hpp>
#include <sup/json_params.hpp>

#include "branch_cell.hpp"

#ifdef ARB_MPI_ENABLED
#include <mpi.h>
#include <arborenv/with_mpi.hpp>
#endif

using namespace arborio::literals;

struct run_params {
    std::string name = "default";
    bool dry_run = false;
    unsigned num_cells_per_rank = 10;
    unsigned num_ranks = 1;
    double min_delay = 10;
    double duration = 100;
    cell_parameters cell;
    bool defaulted = true;
};

void write_trace_json(const arb::trace_data<double>& trace);
run_params read_options(int argc, char** argv);

using arb::cell_gid_type;
using arb::cell_lid_type;
using arb::cell_size_type;
using arb::cell_member_type;
using arb::cell_kind;
using arb::time_type;

// Generate a cell.
arb::cable_cell branch_cell(arb::cell_gid_type gid, const cell_parameters& params);

class tile_desc: public arb::tile {
public:
    tile_desc(unsigned num_cells, unsigned num_tiles, cell_parameters params, unsigned min_delay):
            num_cells_(num_cells),
            num_tiles_(num_tiles),
            cell_params_(params),
            min_delay_(min_delay)
    {}

    cell_size_type num_cells() const override {
        return num_cells_;
    }

    cell_size_type num_tiles() const override {
        return num_tiles_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        return branch_cell(gid, cell_params_);
    }

    cell_kind get_cell_kind(cell_gid_type gid) const override {
        return cell_kind::cable;
    }

    std::any get_global_properties(arb::cell_kind) const override {
        arb::cable_cell_global_properties gprop;
        gprop.default_parameters = arb::neuron_parameter_defaults;
        return gprop;
    }

    // Each cell has one incoming connection, from any cell in the network spanning all ranks:
    // src gid in {0, ..., num_cells_*num_tiles_ - 1}.
    std::vector<arb::cell_connection> connections_on(cell_gid_type gid) const override {
        std::uniform_int_distribution<cell_gid_type> source_distribution(0, num_cells_*num_tiles_ - 2);

        auto src_gen = std::mt19937(gid);
        auto src = source_distribution(src_gen);
        if (src>=gid) ++src;

        return {arb::cell_connection({src, "detector"}, {"synapse"}, event_weight_, min_delay_)};
    }

    // Return an event generator on every 20th gid. This function needs to generate events
    // for ALL cells on ALL ranks. This is because the symmetric recipe can not easily
    // translate the src gid of an event generator
    std::vector<arb::event_generator> event_generators(cell_gid_type gid) const override {
        std::vector<arb::event_generator> gens;
        if (gid%20 == 0) {
            gens.push_back(arb::explicit_generator({{{"synapse"}, 1.0, event_weight_}}));
        }
        return gens;
    }

    std::vector<arb::probe_info> get_probes(cell_gid_type gid) const override {
        // One probe per cell, sampling membrane voltage at end of soma.
        return {arb::cable_probe_membrane_voltage{arb::mlocation{0, 0.0}}};
    }

private:
    cell_size_type num_cells_;
    cell_size_type num_tiles_;
    cell_parameters cell_params_;
    double min_delay_;
    float event_weight_ = 0.05;
};

int main(int argc, char** argv) {
    try {
#ifdef ARB_MPI_ENABLED
        arbenv::with_mpi guard(argc, argv, false);
#endif
        bool root = true;
        auto params = read_options(argc, argv);

        auto resources = arb::proc_allocation();
        auto ctx = arb::make_context(resources);

        if (params.dry_run) {
            ctx = arb::make_context(resources, arb::dry_run_info(params.num_ranks, params.num_cells_per_rank));
        }
#ifdef ARB_MPI_ENABLED
        else {
            ctx = arb::make_context(resources, MPI_COMM_WORLD);
            if (params.defaulted) params.num_ranks = arb::num_ranks(ctx);
            root = arb::rank(ctx)==0;
        }
#endif

#ifdef ARB_PROFILE_ENABLED
        arb::profile::profiler_initialize(ctx);
#endif
        std::cout << sup::mask_stream(root);

        // Print a banner with information about hardware configuration
        std::cout << "gpu:      " << (has_gpu(ctx)? "yes": "no") << "\n";
        std::cout << "threads:  " << num_threads(ctx) << "\n";
        std::cout << "mpi:      " << (has_mpi(ctx)? "yes": "no") << "\n";
        std::cout << "ranks:    " << num_ranks(ctx) << "(" << params.num_ranks << ")\n" << std::endl;
        std::cout << "run mode: " << distribution_type(ctx) << "\n";

        assert(arb::num_ranks(ctx)==params.num_ranks);

        arb::profile::meter_manager meters;
        meters.start(ctx);

        // Create an instance of our tile and use it to make a symmetric_recipe.
        auto tile = std::make_unique<tile_desc>(params.num_cells_per_rank,
                params.num_ranks, params.cell, params.min_delay);
        arb::symmetric_recipe recipe(std::move(tile));

        auto decomp = arb::partition_load_balance(recipe, ctx);

        // Construct the model.
        arb::simulation sim(recipe, decomp, ctx);

        // The id of the only probe on the cell: the cell_member type points to (cell 0, probe 0)
        auto probe_id = cell_member_type{0, 0};
        // The schedule for sampling is 10 samples every 1 ms.
        auto sched = arb::regular_schedule(1);
        // This is where the voltage samples will be stored as (time, value) pairs
        arb::trace_vector<double> voltage;
        // Now attach the sampler at probe_id, with sampling schedule sched, writing to voltage
        sim.add_sampler(arb::one_probe(probe_id), sched, arb::make_simple_sampler(voltage));

        // Set up recording of spikes to a vector on the root process.
        std::vector<arb::spike> recorded_spikes;
        if (root) {
            sim.set_global_spike_callback(
                    [&recorded_spikes](const std::vector<arb::spike>& spikes) {
                        recorded_spikes.insert(recorded_spikes.end(), spikes.begin(), spikes.end());
                    });
        }

        meters.checkpoint("model-init", ctx);

        // Run the simulation for 100 ms, with time steps of 0.025 ms.
        sim.run(params.duration, 0.025);

        meters.checkpoint("model-run", ctx);

        auto ns = sim.num_spikes();
        std::cout << "\n" << ns << " spikes generated at rate of "
                  << params.duration/ns << " ms between spikes\n\n";

        // Write spikes to file
        if (root) {
            std::ofstream fid("spikes.gdf");
            if (!fid.good()) {
                std::cerr << "Warning: unable to open file spikes.gdf for spike output\n";
            }
            else {
                char linebuf[45];
                for (auto spike: recorded_spikes) {
                    auto n = std::snprintf(
                            linebuf, sizeof(linebuf), "%u %.4f\n",
                            unsigned{spike.source.gid}, float(spike.time));
                    fid.write(linebuf, n);
                }
            }
            // Write the samples to a json file.
            write_trace_json(voltage.at(0));
        }

        auto profile = arb::profile::profiler_summary();
        std::cout << profile << "\n";

        auto report = arb::profile::make_meter_report(meters, ctx);
        std::cout << report;
    }
    catch (std::exception& e) {
        std::cerr << "exception caught in ring miniapp:\n" << e.what() << "\n";
        return 1;
    }

    return 0;
}

void write_trace_json(const arb::trace_data<double>& trace) {
    std::string path = "./voltages.json";

    nlohmann::json json;
    json["name"] = "ring demo";
    json["units"] = "mV";
    json["cell"] = "0.0";
    json["probe"] = "0";

    auto& jt = json["data"]["time"];
    auto& jy = json["data"]["voltage"];

    for (const auto& sample: trace) {
        jt.push_back(sample.t);
        jy.push_back(sample.v);
    }

    std::ofstream file(path);
    file << std::setw(1) << json << "\n";
}

run_params read_options(int argc, char** argv) {
    using sup::param_from_json;

    run_params params;
    if (argc<2) {
        std::cout << "Using default parameters.\n";
        return params;
    }
    else {
        params.defaulted = false;
    }
    if (argc>2) {
        throw std::runtime_error("More than one command line option is not permitted.");
    }

    std::string fname = argv[1];
    std::cout << "Loading parameters from file: " << fname << "\n";
    std::ifstream f(fname);

    if (!f.good()) {
        throw std::runtime_error("Unable to open input parameter file: "+fname);
    }

    nlohmann::json json;
    f >> json;

    param_from_json(params.name, "name", json);
    param_from_json(params.dry_run, "dry-run", json);
    param_from_json(params.num_cells_per_rank, "num-cells-per-rank", json);
    param_from_json(params.num_ranks, "num-ranks", json);
    param_from_json(params.duration, "duration", json);
    param_from_json(params.min_delay, "min-delay", json);
    params.cell = parse_cell_parameters(json);

    if (!json.empty()) {
        for (auto it=json.begin(); it!=json.end(); ++it) {
            std::cout << "  Warning: unused input parameter: \"" << it.key() << "\"\n";
        }
        std::cout << "\n";
    }

    return params;
}
