#pragma once

#include <cmath>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <arbor/assert.hpp>
#include <arbor/common_types.hpp>
#include <arbor/fvm_types.hpp>
#include <arbor/simd/simd.hpp>

#include "backends/event.hpp"
#include "util/padded_alloc.hpp"
#include "util/rangeutil.hpp"

#include "matrix_state.hpp"
#include "multi_event_stream.hpp"
#include "threshold_watcher.hpp"
#include "fvm_layout.hpp"
#include "multicore_common.hpp"
#include "partition_by_constraint.hpp"

namespace arb {
namespace multicore {

/*
 * Ion state fields correspond to NMODL ion variables, where X
 * is replaced with the name of the ion. E.g. for calcium 'ca':
 *
 *     Field   NMODL variable   Meaning
 *     -------------------------------------------------------
 *     iX_     ica              calcium ion current density
 *     eX_     eca              calcium ion channel reversal potential
 *     Xi_     cai              internal calcium concentration
 *     Xo_     cao              external calcium concentration
 */

struct ion_state {
    unsigned alignment = 1; // Alignment and padding multiple.

    iarray node_index_;     // Instance to CV map.
    array iX_;              // (A/m²) current density
    array eX_;              // (mV) reversal potential
    array Xi_;              // (mM) internal concentration
    array Xo_;              // (mM) external concentration

    array init_Xi_;         // (mM) area-weighted initial internal concentration
    array init_Xo_;         // (mM) area-weighted initial external concentration
    array reset_Xi_;        // (mM) area-weighted user-set internal concentration
    array reset_Xo_;        // (mM) area-weighted user-set internal concentration
    array init_eX_;         // (mV) initial reversal potential

    array charge;           // charge of ionic species (global value, length 1)

    ion_state() = default;

    ion_state(
        int charge,
        const fvm_ion_config& ion_data,
        unsigned align
    );

    // Set ion concentrations to weighted proportion of default concentrations.
    void init_concentration();

    // Set ionic current density to zero.
    void zero_current();

    // Zero currents, reset concentrations, and reset reversal potential from initial values.
    void reset();
};

struct istim_state {
    unsigned alignment = 1; // Alignment and padding multiple.

    // Immutable data (post initialization):
    iarray accu_index_;     // Instance to accumulator index (accu_stim_ index) map.
    iarray accu_to_cv_;     // Accumulator index to CV map.

    array frequency_;       // (kHz) stimulus frequency per instance.
    array phase_;           // (rad) stimulus waveform phase at t=0.
    array envl_amplitudes_; // (A/m²) stimulus envelope amplitudes, partitioned by instance.
    array envl_times_;      // (A/m²) stimulus envelope timepoints, partitioned by instance.
    iarray envl_divs_;      // Partition divisions for envl_ arrays,

    // Mutable data:
    array accu_stim_;       // (A/m²) accumulated stim current / CV area, one per CV with a stimulus.
    iarray envl_index_;     // Per instance index into envl_ arrays, corresponding to last sample time.

    // Zero stim current.
    void zero_current();

    // Zero stim current, reset indices.
    void reset();

    // Contribute to current density:
    void add_current(const array& time, const iarray& cv_to_intdom, array& current_density);

    // Construct state from i_clamp data:
    istim_state(const fvm_stimulus_config& stim_data, unsigned align);

    istim_state() = default;
};

struct shared_state {
    struct mech_storage {
        array data_;
        iarray indices_;
        constraint_partition constraints_;
    };

    unsigned alignment = 1;   // Alignment and padding multiple.
    util::padded_allocator<> alloc;  // Allocator with corresponging alignment/padding.

    fvm_size_type n_intdom = 0; // Number of integration domains.
    fvm_size_type n_detector = 0; // Max number of detectors on all cells.
    fvm_size_type n_cv = 0;   // Total number of CVs.
    fvm_size_type n_gj = 0;   // Total number of GJs.

    iarray cv_to_intdom;      // Maps CV index to integration domain index.
    iarray cv_to_cell;        // Maps CV index to the first spike
    gjarray gap_junctions;   // Stores gap_junction info.
    array time;               // Maps intdom index to integration start time [ms].
    array time_to;            // Maps intdom index to integration stop time [ms].
    array dt_intdom;          // Maps  index to (stop time) - (start time) [ms].
    array dt_cv;              // Maps CV index to dt [ms].
    array voltage;            // Maps CV index to membrane voltage [mV].
    array current_density;    // Maps CV index to membrane current density contributions [A/m²].
    array conductivity;       // Maps CV index to membrane conductivity [kS/m²].

    array init_voltage;       // Maps CV index to initial membrane voltage [mV].
    array temperature_degC;   // Maps CV to local temperature (read only) [°C].
    array diam_um;            // Maps CV to local diameter (read only) [µm].

    array time_since_spike;   // Stores time since last spike on any detector, organized by cell.
    iarray src_to_spike;      // Maps spike source index to spike index

    arb_value_type* time_ptr;

    istim_state stim_data;
    std::unordered_map<std::string, ion_state> ion_data;
    deliverable_event_stream deliverable_events;
    std::unordered_map<unsigned, mech_storage> storage;

    shared_state() = default;

    shared_state(
        fvm_size_type n_intdom,
        fvm_size_type n_cell,
        fvm_size_type n_detector,
        const std::vector<fvm_index_type>& cv_to_intdom_vec,
        const std::vector<fvm_index_type>& cv_to_cell_vec,
        const std::vector<fvm_gap_junction>& gj_vec,
        const std::vector<fvm_value_type>& init_membrane_potential,
        const std::vector<fvm_value_type>& temperature_K,
        const std::vector<fvm_value_type>& diam,
        const std::vector<fvm_index_type>& src_to_spike,
        unsigned align
    );

    void instantiate(mechanism&, unsigned, const mechanism_overrides&, const mechanism_layout&);

    void set_parameter(mechanism&, const std::string&, const std::vector<arb_value_type>&);

    void add_ion(
        const std::string& ion_name,
        int charge,
        const fvm_ion_config& ion_data);

    void configure_stimulus(const fvm_stimulus_config&);

    void zero_currents();

    void ions_init_concentration();

    void ions_nernst_reversal_potential(fvm_value_type temperature_K);

    // Set time_to to earliest of time+dt_step and tmax.
    void update_time_to(fvm_value_type dt_step, fvm_value_type tmax);

    // Set the per-integration domain and per-compartment dt from time_to - time.
    void set_dt();

    // Update gap_junction state
    void add_gj_current();

    // Update stimulus state and add current contributions.
    void add_stimulus_current();

    // Return minimum and maximum time value [ms] across cells.
    std::pair<fvm_value_type, fvm_value_type> time_bounds() const;

    // Return minimum and maximum voltage value [mV] across cells.
    // (Used for solution bounds checking.)
    std::pair<fvm_value_type, fvm_value_type> voltage_bounds() const;

    // Take samples according to marked events in a sample_event_stream.
    void take_samples(
        const sample_event_stream::state& s,
        array& sample_time,
        array& sample_value);

    void reset();
};

// For debugging only:
std::ostream& operator<<(std::ostream& o, const shared_state& s);


} // namespace multicore
} // namespace arb
