#pragma once

#include <arbor/common_types.hpp>
#include <arbor/fvm_types.hpp>

// Structures for the representation of event delivery targets and
// staged events.

namespace arb {

// Post-synaptic spike events

struct target_handle {
    cell_local_size_type mech_id;    // mechanism type identifier (per cell group).
    cell_local_size_type mech_index; // instance of the mechanism
    cell_size_type intdom_index;     // which integration domain (acts as index into e.g. vec_t)

    target_handle() {}
    target_handle(cell_local_size_type mech_id, cell_local_size_type mech_index, cell_size_type intdom_index):
        mech_id(mech_id), mech_index(mech_index), intdom_index(intdom_index) {}
};

struct deliverable_event {
    time_type time;
    float weight;
    target_handle handle;

    deliverable_event() {}
    deliverable_event(time_type time, target_handle handle, float weight):
        time(time), weight(weight), handle(handle) {}
};

// Stream index accessor function for multi_event_stream:
inline cell_size_type event_index(const deliverable_event& ev) {
    return ev.handle.intdom_index;
}

// Subset of event information required for mechanism delivery.
struct deliverable_event_data {
    cell_local_size_type mech_id;    // same as target_handle::mech_id
    cell_local_size_type mech_index; // same as target_handle::mech_index
    float weight;
};

// Delivery data accessor function for multi_event_stream:
inline deliverable_event_data event_data(const deliverable_event& ev) {
    return {ev.handle.mech_id, ev.handle.mech_index, ev.weight};
}


// Sample events (raw values from back-end state).

using probe_handle = const fvm_value_type*;

struct raw_probe_info {
    probe_handle handle;      // where the to-be-probed value sits
    sample_size_type offset;  // offset into array to store raw probed value
};

struct sample_event {
    time_type time;
    cell_size_type intdom_index;  // which integration domain probe is on
    raw_probe_info raw;           // event payload: what gets put where on sample
};

inline raw_probe_info event_data(const sample_event& ev) {
    return ev.raw;
}

inline cell_size_type event_index(const sample_event& ev) {
    return ev.intdom_index;
}


} // namespace arb
