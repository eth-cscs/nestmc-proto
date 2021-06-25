#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include <arbor/mechinfo.hpp>
#include <arbor/fvm_types.hpp>
#include <arbor/mechanism_abi.h>

namespace arb {

class mechanism;
using mechanism_ptr = std::unique_ptr<mechanism>;

struct ion_state_view {
    fvm_value_type* current_density;
    fvm_value_type* reversal_potential;
    fvm_value_type* internal_concentration;
    fvm_value_type* external_concentration;
    fvm_value_type* ionic_charge;
};

// Generated mechanism field, global and ion table lookup types.
// First component is name, second is pointer to corresponing member in
// the mechanism's parameter pack, or for field_default_table,
// the scalar value used to initialize the field.
using global_table_entry = std::pair<const char*, fvm_value_type>;
using mechanism_global_table = std::vector<global_table_entry>;

using state_table_entry = std::pair<const char*, std::pair<fvm_value_type*, fvm_value_type>>;
using mechanism_state_table = std::vector<state_table_entry>;

using field_table_entry = std::pair<const char*, std::pair<fvm_value_type*, fvm_value_type>>;
using mechanism_field_table = std::vector<field_table_entry>;

using ion_state_entry = std::pair<const char*, std::pair<ion_state_view, fvm_index_type*>>;
using mechanism_ion_table = std::vector<ion_state_entry>;

class mechanism {
public:
    using value_type = fvm_value_type;
    using index_type = fvm_index_type;
    using size_type  = fvm_size_type;

    mechanism(const arb_mechanism_type m, const arb_mechanism_interface& i): mech_{m}, iface_{i} {}
    mechanism() = default;
    mechanism(const mechanism&) = delete;

    // Return fingerprint of mechanism dynamics source description for validation/replication.
    const mechanism_fingerprint fingerprint() const { return mech_.fingerprint; };

    // Name as given in mechanism source.
    std::string internal_name() const { return mech_.name; }

    // Density or point mechanism?
    arb_mechanism_kind kind() const { return mech_.kind; };

    // Does the implementation require padding and alignment of shared data structures?
    unsigned data_alignment() const { return iface_.alignment; }

    // Make a new object of the derived mechanism type, but does not copy any state.
    mechanism_ptr clone() const { return std::make_unique<mechanism>(mech_, iface_); }

    // Non-global parameters can be set post-instantiation:
    void set_parameter(const std::string&, const std::vector<arb_value_type>&);

    // Forward to interface methods
    void initialize()     { ppack_.vec_t = *time_ptr_ptr; iface_.init_mechanism(&ppack_); }
    void update_current() { ppack_.vec_t = *time_ptr_ptr; iface_.compute_currents(&ppack_); }
    void update_state()   { ppack_.vec_t = *time_ptr_ptr; iface_.advance_state(&ppack_); }
    void update_ions()    { ppack_.vec_t = *time_ptr_ptr; iface_.write_ions(&ppack_); }
    void post_event()     { ppack_.vec_t = *time_ptr_ptr; iface_.post_event(&ppack_); }
    void deliver_events(arb_deliverable_event_stream& stream) {
        ppack_.vec_t  = *time_ptr_ptr;
        ppack_.events = stream;
        iface_.apply_events(&ppack_);
    }

    // Peek into state variable
    fvm_value_type* field_data(const std::string& var);

    mechanism_field_table field_table();
    mechanism_global_table global_table();
    mechanism_state_table state_table();
    mechanism_ion_table ion_table();

    virtual ~mechanism() = default;

    // Per-cell group identifier for an instantiated mechanism.
    unsigned mechanism_id() const { return ppack_.mechanism_id; }

    arb_mechanism_type  mech_;
    arb_mechanism_interface iface_;
    arb_mechanism_ppack ppack_;
    bool mult_in_place_;                             // perform multiplication in place?
    size_type num_ions_ = 0;                         // Ion count

    size_type width_padded_;

    std::vector<arb_value_type>  globals_;
    std::vector<arb_value_type*> parameters_;
    std::vector<arb_value_type*> state_vars_;
    std::vector<arb_ion_state>   ion_states_;

   arb_value_type** time_ptr_ptr;
};

struct mechanism_layout {
    // Maps in-instance index to CV index.
    std::vector<fvm_index_type> cv;

    // Maps in-instance index to compartment contribution.
    std::vector<fvm_value_type> weight;

    // Number of logical point processes at in-instance index;
    // if empty, point processes are not coalesced and all multipliers are 1.
    std::vector<fvm_index_type> multiplicity;
};

struct mechanism_overrides {
    // Global scalar parameters (any value down-conversion to fvm_value_type is the
    // responsibility of the concrete mechanism).
    std::unordered_map<std::string, double> globals;

    // Ion renaming: keys are ion dependency names as
    // reported by the mechanism info.
    std::unordered_map<std::string, std::string> ion_rebind;
};

} // namespace arb
