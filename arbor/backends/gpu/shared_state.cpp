#include <cstddef>
#include <vector>

#include <arbor/constants.hpp>
#include <arbor/fvm_types.hpp>
#include <arbor/math.hpp>

#include "backends/event.hpp"
#include "backends/gpu/gpu_store_types.hpp"
#include "backends/gpu/shared_state.hpp"
#include "backends/multi_event_stream_state.hpp"
#include "memory/copy.hpp"
#include "memory/wrappers.hpp"
#include "util/index_into.hpp"
#include "util/rangeutil.hpp"
#include "util/maputil.hpp"
#include "util/range.hpp"
#include "util/strprintf.hpp"

using arb::memory::make_const_view;

namespace arb {
namespace gpu {

// CUDA implementation entry points:

void update_time_to_impl(
    std::size_t n, fvm_value_type* time_to, const fvm_value_type* time,
    fvm_value_type dt, fvm_value_type tmax);

void update_time_to_impl(
    std::size_t n, fvm_value_type* time_to, const fvm_value_type* time,
    fvm_value_type dt, fvm_value_type tmax);

void set_dt_impl(
    fvm_size_type nintdom, fvm_size_type ncomp, fvm_value_type* dt_intdom, fvm_value_type* dt_comp,
    const fvm_value_type* time_to, const fvm_value_type* time, const fvm_index_type* cv_to_intdom);

void add_gj_current_impl(
    fvm_size_type n_gj, const fvm_gap_junction* gj, const fvm_value_type* v, fvm_value_type* i);

void take_samples_impl(
    const multi_event_stream_state<raw_probe_info>& s,
    const fvm_value_type* time, fvm_value_type* sample_time, fvm_value_type* sample_value);

void add_scalar(std::size_t n, fvm_value_type* data, fvm_value_type v);

// GPU-side minmax: consider CUDA kernel replacement.
std::pair<fvm_value_type, fvm_value_type> minmax_value_impl(fvm_size_type n, const fvm_value_type* v) {
    auto v_copy = memory::on_host(memory::const_device_view<fvm_value_type>(v, n));
    return util::minmax_value(v_copy);
}

// Ion state methods:

ion_state::ion_state(
    int charge,
    const fvm_ion_config& ion_data,
    unsigned // alignment/padding ignored.
):
    node_index_(make_const_view(ion_data.cv)),
    iX_(ion_data.cv.size(), NAN),
    eX_(ion_data.cv.size(), NAN),
    Xi_(ion_data.cv.size(), NAN),
    Xo_(ion_data.cv.size(), NAN),
    init_Xi_(make_const_view(ion_data.init_iconc)),
    init_Xo_(make_const_view(ion_data.init_econc)),
    reset_Xi_(make_const_view(ion_data.reset_iconc)),
    reset_Xo_(make_const_view(ion_data.reset_econc)),
    init_eX_(make_const_view(ion_data.init_revpot)),
    charge(1u, charge)
{
    arb_assert(node_index_.size()==init_Xi_.size());
    arb_assert(node_index_.size()==init_Xo_.size());
    arb_assert(node_index_.size()==init_eX_.size());
}

void ion_state::init_concentration() {
    memory::copy(init_Xi_, Xi_);
    memory::copy(init_Xo_, Xo_);
}

void ion_state::zero_current() {
    memory::fill(iX_, 0);
}

void ion_state::reset() {
    zero_current();
    memory::copy(reset_Xi_, Xi_);
    memory::copy(reset_Xo_, Xo_);
    memory::copy(init_eX_, eX_);
}

// istim_state methods:

istim_state::istim_state(const fvm_stimulus_config& stim) {
    using util::assign;

    // Translate instance-to-CV index from stim to istim_state index vectors.
    std::vector<fvm_index_type> accu_index_stage;
    assign(accu_index_stage, util::index_into(stim.cv, stim.cv_unique));

    std::size_t n = accu_index_stage.size();
    std::vector<fvm_value_type> envl_a, envl_t;
    std::vector<fvm_index_type> edivs;

    frequency_ = make_const_view(stim.frequency);
    phase_ = make_const_view(stim.phase);

    arb_assert(n==frequency_.size());
    arb_assert(n==stim.envelope_time.size());
    arb_assert(n==stim.envelope_amplitude.size());

    edivs.reserve(n+1);
    edivs.push_back(0);

    for (auto i: util::make_span(n)) {
        arb_assert(stim.envelope_time[i].size()==stim.envelope_amplitude[i].size());
        arb_assert(util::is_sorted(stim.envelope_time[i]));

        util::append(envl_a, stim.envelope_amplitude[i]);
        util::append(envl_t, stim.envelope_time[i]);
        edivs.push_back(fvm_index_type(envl_t.size()));
    }

    accu_index_ = make_const_view(accu_index_stage);
    accu_to_cv_ = make_const_view(stim.cv_unique);
    accu_stim_ = array(accu_index_.size());
    envl_amplitudes_ = make_const_view(envl_a);
    envl_times_ = make_const_view(envl_t);
    envl_divs_ = make_const_view(edivs);

    // Initial indices into envelope match partition divisions; ignore last (index n) element.
    envl_index_ = envl_divs_;

    // Initialize ppack pointers.
    ppack_.accu_index = accu_index_.data();
    ppack_.accu_to_cv = accu_to_cv_.data();
    ppack_.frequency = frequency_.data();
    ppack_.phase = phase_.data();
    ppack_.envl_amplitudes = envl_amplitudes_.data();
    ppack_.envl_times = envl_times_.data();
    ppack_.envl_divs = envl_divs_.data();
    ppack_.accu_stim = accu_stim_.data();
    ppack_.envl_index = envl_index_.data();
    // The following ppack fields must be set in add_current() before queuing kernel.
    ppack_.time = nullptr;
    ppack_.cv_to_intdom = nullptr;
    ppack_.current_density = nullptr;
}

std::size_t istim_state::size() const {
    return frequency_.size();
}

void istim_state::zero_current() {
    memory::fill(accu_stim_, 0.);
}

void istim_state::reset() {
    zero_current();
    memory::copy(envl_divs_, envl_index_);
}

void istim_state::add_current(const array& time, const iarray& cv_to_intdom, array& current_density) {
    ppack_.time = time.data();
    ppack_.cv_to_intdom = cv_to_intdom.data();
    ppack_.current_density = current_density.data();
    istim_add_current_impl((int)size(), ppack_);
}

// Shared state methods:

shared_state::shared_state(
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
    unsigned // alignment parameter ignored.
    ):
    n_intdom(n_intdom),
    n_detector(n_detector),
    n_cv(cv_to_intdom_vec.size()),
    n_gj(gj_vec.size()),
    cv_to_intdom(make_const_view(cv_to_intdom_vec)),
    cv_to_cell(make_const_view(cv_to_cell_vec)),
    gap_junctions(make_const_view(gj_vec)),
    time(n_intdom),
    time_to(n_intdom),
    dt_intdom(n_intdom),
    dt_cv(n_cv),
    voltage(n_cv),
    current_density(n_cv),
    conductivity(n_cv),
    init_voltage(make_const_view(init_membrane_potential)),
    temperature_degC(make_const_view(temperature_K)),
    diam_um(make_const_view(diam)),
    time_since_spike(n_cell*n_detector),
    src_to_spike(make_const_view(src_to_spike)),
    deliverable_events(n_intdom)
{
    memory::fill(time_since_spike, -1.0);
    add_scalar(temperature_degC.size(), temperature_degC.data(), -273.15);
}

namespace {
template<typename T>
memory::const_host_view<T> std_view(const std::vector<T>& v) { return memory::const_host_view<T>{v.data(), v.size()}; }

template<typename T>
void append_chunk(const std::vector<T>& in, T*& out, T*& ptr, size_t n) {
    memory::copy(std_view(in), memory::device_view<T>(ptr, n));
    out = ptr;
    ptr += n;
};
  
template<typename T>
void append_const(T in, T*& out, T*& ptr, size_t n) {
    memory::fill(memory::device_view<T>(ptr, n), in);
    out = ptr;
    ptr += n;
};
  
}

void shared_state::set_parameter(mechanism& m, const std::string& key, const std::vector<arb_value_type>& values) {
    if (values.size()!=m.ppack_.width) throw arbor_internal_error("mechanism parameter size mismatch");
    auto field_ptr = m.field_data(key);
    if (!field_ptr) throw arbor_internal_error(util::pprintf("no such mechanism parameter '{}'", key));
    if (!m.ppack_.width) return;
    assert (m.iface_.backend);
    memory::copy(memory::make_const_view(values), memory::device_view<arb_value_type>(field_ptr, m.ppack_.width));
}

void shared_state::instantiate(mechanism& m, unsigned id, const mechanism_overrides& overrides, const mechanism_layout& pos_data) {
    assert(m.iface_.backend == arb_backend_kind_gpu);
    using util::make_range;
    using util::make_span;
    using util::ptr_by_key;
    using util::value_by_key;

    // Set internal variables
    m.mult_in_place_ = !pos_data.multiplicity.empty();
    m.num_ions_      = m.mech_.n_ions;
    m.time_ptr_ptr   = &time_ptr;

    auto alignment    = std::max(array::alignment(), iarray::alignment());
    auto width        = pos_data.cv.size();
    auto width_padded = math::round_up(pos_data.cv.size(), alignment);

    // Assign non-owning views onto shared state:
    m.ppack_ = {0};
    m.ppack_.width            = width;
    m.ppack_.mechanism_id     = id;
    m.ppack_.vec_ci           = cv_to_cell.data();
    m.ppack_.vec_di           = cv_to_intdom.data();
    m.ppack_.vec_dt           = dt_cv.data();
    m.ppack_.vec_v            = voltage.data();
    m.ppack_.vec_i            = current_density.data();
    m.ppack_.vec_g            = conductivity.data();
    m.ppack_.temperature_degC = temperature_degC.data();
    m.ppack_.diam_um          = diam_um.data();
    m.ppack_.time_since_spike = time_since_spike.data();
    m.ppack_.n_detectors      = n_detector;

    // Allocate view pointers
    m.state_vars_ = std::vector<arb_value_type*>(m.mech_.n_state_vars);
    m.parameters_ = std::vector<arb_value_type*>(m.mech_.n_parameters);
    m.ion_states_ = std::vector<arb_ion_state>(m.mech_.n_ions);
    m.globals_    = std::vector<arb_value_type>(m.mech_.n_globals);

    if (storage.find(id) != storage.end()) throw arb::arbor_internal_error("Duplicate mech id in shared state");
    auto& store = storage[id];

    // Set ion views
    for (auto idx: make_span(m.mech_.n_ions)) {
        auto ion = m.mech_.ions[idx].name;
        auto ion_binding = value_by_key(overrides.ion_rebind, ion).value_or(ion);
        ion_state* oion = ptr_by_key(ion_data, ion_binding);
        if (!oion) throw arbor_internal_error("gpu/mechanism: mechanism holds ion with no corresponding shared state");
        m.ion_states_[idx] = { oion->iX_.data(), oion->eX_.data(), oion->Xi_.data(), oion->Xo_.data(), oion->charge.data(), nullptr };
    }

    // If there are no sites (is this ever meaningful?) there is nothing more to do.
    if (width==0) return;

    // Allocate and initialize state and parameter vectors with default values.
    {
        // Allocate bulk storage
        auto count    = (m.mech_.n_state_vars + m.mech_.n_parameters + 1)*width_padded + m.mech_.n_globals;
        store.data_   = array(count, NAN);
        auto base_ptr = store.data_.data();
        // First sub-array of data_ is used for weight_
        append_chunk(pos_data.weight, m.ppack_.weight, base_ptr, width);
        // Set fields
        for (auto idx: make_span(m.mech_.n_parameters)) {
            append_const(m.mech_.parameters[idx].default_value, m.parameters_[idx], base_ptr, width);
        }
        for (auto idx: make_span(m.mech_.n_state_vars)) {
            append_const(m.mech_.state_vars[idx].default_value, m.state_vars_[idx], base_ptr, width);
        }
        // Assign global scalar parameters. NB: Last chunk, since it breaks the width_padded alignment
        for (auto idx: make_span(m.mech_.n_globals)) m.globals_[idx] = m.mech_.globals[idx].default_value;
        for (auto& [k, v]: overrides.globals) {
            auto found = false;
            for (auto idx: make_span(m.mech_.n_globals)) {
                if (m.mech_.globals[idx].name == k) {
                    m.globals_[idx] = v;
                    found = true;
                    break;
                }
            }
            if (!found) throw arbor_internal_error(util::pprintf("gpu/mechanism: no such mechanism global '{}'", k));
        }
        memory::copy(std_view(m.globals_), memory::device_view<arb_value_type>(base_ptr, m.mech_.n_globals));
        m.ppack_.globals = base_ptr;
        base_ptr += m.mech_.n_globals;
    }

    // Allocate and initialize index vectors, viz. node_index_ and any ion indices.
    {
        // Allocate bulk storage
        auto count     = (m.mult_in_place_ ? 1 : 0) + m.mech_.n_ions + 1;
        store.indices_ = iarray(count*width_padded);
        auto base_ptr  = store.indices_.data();
        // Setup node indices
        append_chunk(pos_data.cv, m.ppack_.node_index, base_ptr, width);
        // Create ion indices
        for (auto idx: make_span(m.mech_.n_ions)) {
            auto  ion = m.mech_.ions[idx].name;
            auto& index_ptr = m.ion_states_[idx].index;
            // Index into shared_state respecting ion rebindings
            auto ion_binding = value_by_key(overrides.ion_rebind, ion).value_or(ion);
            ion_state* oion = ptr_by_key(ion_data, ion_binding);
            if (!oion) throw arbor_internal_error("gpu/mechanism: mechanism holds ion with no corresponding shared state");
            // Obtain index and move data
            auto ni = memory::on_host(oion->node_index_);
            auto indices = util::index_into(pos_data.cv, ni);
            std::vector<arb_index_type> mech_ion_index(indices.begin(), indices.end());
            append_chunk(mech_ion_index, index_ptr, base_ptr, width);
        }

        if (m.mult_in_place_) append_chunk(pos_data.multiplicity, m.ppack_.multiplicity, base_ptr, width);
    }

    // Shift data to GPU, set up pointers
    store.parameters_d_ = memory::device_vector<arb_value_type*>(m.parameters_.size());
    memory::copy(std_view(m.parameters_), store.parameters_d_);
    m.ppack_.parameters = store.parameters_d_.data();
    store.state_vars_d_ = memory::device_vector<arb_value_type*>(m.state_vars_.size());
    memory::copy(std_view(m.state_vars_), store.state_vars_d_);
    m.ppack_.state_vars = store.state_vars_d_.data();
    store.ion_states_d_ = memory::device_vector<arb_ion_state>(m.ion_states_.size());
    memory::copy(std_view(m.ion_states_), store.ion_states_d_);
    m.ppack_.ion_states = store.ion_states_d_.data();
}

void shared_state::add_ion(
    const std::string& ion_name,
    int charge,
    const fvm_ion_config& ion_info)
{
    ion_data.emplace(std::piecewise_construct,
        std::forward_as_tuple(ion_name),
        std::forward_as_tuple(charge, ion_info, 1u));
}

void shared_state::configure_stimulus(const fvm_stimulus_config& stims) {
    stim_data = istim_state(stims);
}

void shared_state::reset() {
    memory::copy(init_voltage, voltage);
    memory::fill(current_density, 0);
    memory::fill(conductivity, 0);
    memory::fill(time, 0);
    memory::fill(time_to, 0);
    memory::fill(time_since_spike, -1.0);

    for (auto& i: ion_data) {
        i.second.reset();
    }
    stim_data.reset();
}

void shared_state::zero_currents() {
    memory::fill(current_density, 0);
    memory::fill(conductivity, 0);
    for (auto& i: ion_data) {
        i.second.zero_current();
    }
    stim_data.zero_current();
}

void shared_state::ions_init_concentration() {
    for (auto& i: ion_data) {
        i.second.init_concentration();
    }
}

void shared_state::update_time_to(fvm_value_type dt_step, fvm_value_type tmax) {
    update_time_to_impl(n_intdom, time_to.data(), time.data(), dt_step, tmax);
}

void shared_state::set_dt() {
    set_dt_impl(n_intdom, n_cv, dt_intdom.data(), dt_cv.data(), time_to.data(), time.data(), cv_to_intdom.data());
}

void shared_state::add_gj_current() {
    add_gj_current_impl(n_gj, gap_junctions.data(), voltage.data(), current_density.data());
}

void shared_state::add_stimulus_current() {
    stim_data.add_current(time, cv_to_intdom, current_density);
}

std::pair<fvm_value_type, fvm_value_type> shared_state::time_bounds() const {
    return minmax_value_impl(n_intdom, time.data());
}

std::pair<fvm_value_type, fvm_value_type> shared_state::voltage_bounds() const {
    return minmax_value_impl(n_cv, voltage.data());
}

void shared_state::take_samples(const sample_event_stream::state& s, array& sample_time, array& sample_value) {
    take_samples_impl(s, time.data(), sample_time.data(), sample_value.data());
}

// Debug interface
std::ostream& operator<<(std::ostream& o, shared_state& s) {
    o << " cv_to_intdom " << s.cv_to_intdom << "\n";
    o << " time         " << s.time << "\n";
    o << " time_to      " << s.time_to << "\n";
    o << " dt_intdom    " << s.dt_intdom << "\n";
    o << " dt_cv        " << s.dt_cv << "\n";
    o << " voltage      " << s.voltage << "\n";
    o << " init_voltage " << s.init_voltage << "\n";
    o << " temperature  " << s.temperature_degC << "\n";
    o << " diameter     " << s.diam_um << "\n";
    o << " current      " << s.current_density << "\n";
    o << " conductivity " << s.conductivity << "\n";
    for (auto& ki: s.ion_data) {
        auto& kn = ki.first;
        auto& i = ki.second;
        o << " " << kn << "/current_density        " << i.iX_ << "\n";
        o << " " << kn << "/reversal_potential     " << i.eX_ << "\n";
        o << " " << kn << "/internal_concentration " << i.Xi_ << "\n";
        o << " " << kn << "/external_concentration " << i.Xo_ << "\n";
        o << " " << kn << "/intconc_initial        " << i.init_Xi_ << "\n";
        o << " " << kn << "/extconc_initial        " << i.init_Xo_ << "\n";
        o << " " << kn << "/revpot_initial         " << i.init_eX_ << "\n";
        o << " " << kn << "/node_index             " << i.node_index_ << "\n";
    }
    return o;
}

} // namespace gpu
} // namespace arb
