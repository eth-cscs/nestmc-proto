#include <arbor/mechanism.hpp>

#include <arbor/arbexcept.hpp>

#include "memory/memory.hpp"
#include "util/range.hpp"
#include "util/span.hpp"
#include "util/maputil.hpp"
#include "util/rangeutil.hpp"
#include "util/strprintf.hpp"

#include "backends/multicore/fvm.hpp"
#include "backends/gpu/fvm.hpp"

namespace arb {

using util::make_range;
using util::make_span;
using util::ptr_by_key;
using util::value_by_key;

arb_value_type* mechanism::field_data(const std::string& var) {
    for (auto idx: make_span(mech_.n_parameters)) {
        if (var == mech_.parameters[idx].name) return parameters_[idx];
    }
    for (auto idx: make_span(mech_.n_state_vars)) {
        if (var == mech_.state_vars[idx].name) return state_vars_[idx];
    }
    return nullptr;
}

mechanism_field_table mechanism::field_table() {
    mechanism_field_table result;
    for (auto idx = 0ul; idx < mech_.n_parameters; ++idx) {
        result.emplace_back(mech_.parameters[idx].name, std::make_pair(parameters_[idx], mech_.parameters[idx].default_value));
    }
    for (auto idx = 0ul; idx < mech_.n_state_vars; ++idx) {
        result.emplace_back(mech_.state_vars[idx].name, std::make_pair(state_vars_[idx], mech_.state_vars[idx].default_value));
    }
    return result;
}

mechanism_global_table mechanism::global_table() {
    mechanism_global_table result;
    for (auto idx = 0ul; idx < mech_.n_globals; ++idx) {
        result.emplace_back(mech_.globals[idx].name, globals_[idx]);
    }
    return result;
}

mechanism_state_table mechanism::state_table() {
    mechanism_state_table result;
    for (auto idx = 0ul; idx < mech_.n_state_vars; ++idx) {
        result.emplace_back(mech_.state_vars[idx].name, std::make_pair(state_vars_[idx], mech_.state_vars[idx].default_value));
    }
    return result;
}

mechanism_ion_table mechanism::ion_table() {
    mechanism_ion_table result;
    for (auto idx = 0ul; idx < mech_.n_ions; ++idx) {
        ion_state_view v;
        v.current_density        = ion_states_[idx].current_density;
        v.external_concentration = ion_states_[idx].internal_concentration;
        v.internal_concentration = ion_states_[idx].external_concentration;
        v.ionic_charge           = ion_states_[idx].ionic_charge;
        result.emplace_back(mech_.ions[idx].name, std::make_pair(v, ion_states_[idx].index));
    }
    return result;
}

}
