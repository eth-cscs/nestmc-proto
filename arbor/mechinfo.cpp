#include <arbor/mechinfo.hpp>

#include "util/span.hpp"

namespace arb {
mechanism_info::mechanism_info(const arb_mechanism_type& m) {
    post_events = m.has_post_events;
    linear      = m.is_linear;
    fingerprint = m.fingerprint;
    for (auto idx: util::make_span(m.n_globals)) {
        const auto& v = m.globals[idx];
        globals[v.name] = { mechanism_field_spec::field_kind::global, v.unit, v.default_value, v.range_low, v.range_high };
    }
    for (auto idx: util::make_span(m.n_parameters)) {
        const auto& v = m.parameters[idx];
        parameters[v.name] = { mechanism_field_spec::field_kind::parameter, v.unit, v.default_value, v.range_low, v.range_high };
    }
    for (auto idx: util::make_span(m.n_state_vars)) {
        const auto& v = m.state_vars[idx];
        state[v.name] = { mechanism_field_spec::field_kind::state, v.unit, v.default_value, v.range_low, v.range_high };
    }
    for (auto idx: util::make_span(m.n_ions)) {
        const auto& v = m.ions[idx];
        ions[v.name] = { v.write_int_concentration,
        v.write_ext_concentration,
        v.read_rev_potential,
        v.write_rev_potential,
        v.read_valence,
        v.verify_valence,
        v.expected_valence };
    }
}

mechanism_info::mechanism_info(const std::unordered_map<std::string, mechanism_field_spec>& g,
                               const std::unordered_map<std::string, mechanism_field_spec>& p,
                               const std::unordered_map<std::string, mechanism_field_spec>& s,
                               const std::unordered_map<std::string, ion_dependency>& i,
                               const mechanism_fingerprint& f,
                               bool l,
                               bool q): globals{g}, parameters{p}, state{s}, ions{i}, fingerprint{f}, linear{l}, post_events{q} {}

}
