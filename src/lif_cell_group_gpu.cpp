#include <lif_cell_group_gpu.hpp>

using namespace nest::mc;

// Constructor containing gid of first cell in a group and a container of all cells.
lif_cell_group_gpu::lif_cell_group_gpu(cell_gid_type first_gid, const std::vector<util::unique_any>& cells):
gid_base_(first_gid)
{
    cells_.reserve(cells.size());
    lambda_.resize(cells.size());
    generator_.resize(cells.size());
    next_poiss_time_.resize(cells.size());

    last_time_updated_.resize(cells.size());

    tau_m.reserve(cells.size());
    V_th.reserve(cells.size());
    C_m.reserve(cells.size());
    E_L.reserve(cells.size());
    V_m.reserve(cells.size());
    V_reset.reserve(cells.size());
    t_ref.reserve(cells.size());

    n_poiss.reserve(cells.size());
    w_poiss.reserve(cells.size());
    d_poiss.reserve(cells.size());

    for (const auto& cell : cells) {
        cells_.push_back(cell);

        tau_m.push_back(cell.tau_m);
        V_th.push_back(cell.V_th);
        C_m.push_back(cell.C_m);
        E_L.push_back(cell.E_L);
        V_m.push_back(cell.V_m);
        V_reset.push_back(cell.V_reset);
        t_ref.push_back(cell.t_ref);

        n_poiss.push_back(cell.n_poiss);
        w_poiss.push_back(cell.w_poiss);
        d_poiss.push_back(cell.d_poiss);
    }

    // Initialize variables for the external Poisson input.
    for (auto lid : util::make_span(0, cells_.size())) {
        // If a cell receives some external Poisson input then initialize the corresponding variables.
        if (cells_[lid].n_poiss > 0) {
            lambda_[lid] = (1.0/(cells_[lid].rate * cells_[lid].n_poiss));
            generator_[lid].seed(1000 + first_gid + lid);
            sample_next_poisson(lid);
        }
    }
}

cell_kind lif_cell_group_gpu::get_cell_kind() const {
    return cell_kind::lif_neuron;
}

void lif_cell_group_gpu::advance(time_type tfinal, time_type dt) {
    PE("lif");
    for (size_t lid = 0; lid < cells_.size(); ++lid) {
        // Advance each cell independently.
        advance_cell(tfinal, dt, lid);
    }
    PL();
}

void lif_cell_group_gpu::enqueue_events(const std::vector<postsynaptic_spike_event>& events) {
    // Distribute incoming events to individual cells.
    for (auto& e: events) {
        cell_events_[e.target.gid - gid_base_].push(e);
    }
}

const std::vector<spike>& lif_cell_group_gpu::spikes() const {
    return spikes_;
}

void lif_cell_group_gpu::clear_spikes() {
    spikes_.clear();
}

// TODO: implement sampler
void lif_cell_group_gpu::add_sampler(cell_member_type probe_id, sampler_function s, time_type start_time) {}

// TODO: implement binner_
void lif_cell_group_gpu::set_binning_policy(binning_kind policy, time_type bin_interval) {
}

// no probes in single-compartment cells
std::vector<probe_record> lif_cell_group_gpu::probes() const {
    return {};
}

void lif_cell_group_gpu::reset() {
    spikes_.clear();
    // Clear all the event queues.
    for (auto& queue : cell_events_) {
        queue.clear();
    }
}

// Returns the next most recent event that is yet to be processed.
// It can be either Poisson event or the queue event.
// Only events that happened before tfinal are considered.
__device__
bool lif_cell_group_gpu::next_event(
       cell_gid_type gid_base_,
       time_type* next_poiss_time_,
       float* w_poiss,
       float * d_poiss,
       cell_gid_type lid,
       time_type tfinal,
       double* lambda_,
       unsigned* cell_begin,
       unsigned* cell_end,
       postsynaptic_spike_event* event_buffer,
       postsynaptic_spike_event& event)
{
    auto t_poiss = next_poiss_time_[lid] + d_poiss[lid];
    auto poiss_ev = postsynaptic_spike_event{{cell_lid_type(gid_base_ + lid), 0}, t_poiss, w_poiss[lid]};

    // There are still queue events.
    if (cell_begin[lid] < cell_end[lid]) {
        postsynaptic_spike_event q_ev = event_buffer[cell_begin[lid]];
        // Queue event is the most recent one.
        if (q_ev.time < std::min(tfinal, t_poiss)) {
            cell_begin[lid]++;
            event = q_ev;
            return true;
        }

        // Poisson event is the most recent one.
        if (t_poiss < tfinal) {
            // Sample next Poisson event.
            next_poiss_time_[lid] += exp_dist_(generator_[lid]) * lambda_[lid];
            event = poiss_ev;
            return true;
        }

        // There are events, but not before tfinal
        return false;
    }

    // There are no more queue events but possibly Poisson events.
    if (t_poiss < tfinal) {
        // Sample next Poisson event.
        next_poiss_time_[lid] += exp_dist_(generator_[lid]) * lambda_[lid];
        event = poiss_ev;
        return true;
    }

    // There are neither Poisson nor queue events before tfinal.
    return false;
}

__global__
void advance_kernel (cell_gid_type gid_base_,
            time_type tfinal,
            unsigned num_cells,
            double* tau_m,
            double* V_th,
            double* C_m,
            double* E_L,
            double* V_m,
            double* V_reset,
            double* t_ref,
            unsigned* n_poiss,
            float* rate,
            float* w_poiss,
            float* d_poiss,
            double* lambda_,
            time_type* last_time_updated_,
            time_type* next_poiss_time_,
            unsigned* cell_begin,
            unsigned* cell_end,
            postsynaptic_spike_event* event_buffer)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx < num_cells)
        // Current time of last update.
        auto t = last_time_updated_[lid];
        postsynaptic_spike_event ev;
        // If a neuron was in the refractory period,
        // ignore any new events that happened before t,
        // including poisson events as well.
        while (next_event(gid_base_,
                        next_poiss_time_,
                        w_poiss,
                        d_poiss,
                        lid,
                        t,
                        lambda_,
                        cell_begin,
                        cell_end,
                        postsynaptic_spike_event* event_buffer,
                        ev)) {};

        // Integrate until tfinal using the exact solution of membrane voltage differential equation.
        while (next_event(gid_base_,
                        next_poiss_time_,
                        w_poiss,
                        d_poiss,
                        lid,
                        tfinal,
                        lambda_,
                        cell_begin,
                        cell_end,
                        postsynaptic_spike_event* event_buffer,
                        ev))  {
            auto weight = ev->weight;
            auto event_time = ev->time;

            // If a neuron is in refractory period, ignore this event.
            if (event_time < t) {
                continue;
            }

            // Let the membrane potential decay.
            V_m[lid] *= exp(-(event_time - t) / tau_m[lid]);
            // Add jump due to spike.
            V_m[lid] += weight/C_m[lid];

            t = event_time;

            // If crossing threshold occurred
            if (V_m[lid] >= V_th[lid]) {
                cell_member_type spike_neuron_gid = {gid_base_ + lid, 0};
                spike s = {spike_neuron_gid, t};

                spikes_.push_back(s);

                // Advance last_time_updated.
                t += cell.t_ref;

                // Reset the voltage to the resting potential.
                V_m[lid] = E_L[lid];
            }
            // This is the last time a cell was updated.
            last_time_updated_[lid] = t;
        }
    }
}

// Advances a single cell (lid) with the exact solution (jumps can be arbitrary).
// Parameter dt is ignored, since we make jumps between two consecutive spikes.
void lif_cell_group_gpu::advance_cell(time_type tfinal, time_type dt, cell_gid_type lid) {
    event_buffer.clear();
    cell_begin.clear();
    cell_end.clear();

    for (unsigned i : util:make_span(0, cells_.size())) {
        cell_begin.push_back(event_buffer.size());
        auto& q = cell_events_[i];
        while(auto e = q.pop_if_before(tfinal)){
            event_buffer.push_back(*e);
        }
        cell_end.push_back(event_buffer.size());
    }

    unsigned block_dim = 128;

    unsigned grid_dim = (cells_.size() - 1) / block_dim + 1;

    managed_vector<unsigned> spike_buffer(event_buffer.size());

    advance_kernel<<<grid_dim, block_dim>>>(gid_base_,
                                            tfinal,
                                            cells_.size(),
                                            tau_m.data(),
                                            V_th.data(),
                                            C_m.data(),
                                            E_L.data(),
                                            V_m.data(),
                                            V_reset.data(),
                                            t_ref.data(),
                                            n_poiss.data(),
                                            rate.data(),
                                            w_poiss.data(),
                                            d_poiss.data(),
                                            lambda_.data(),
                                            last_time_updated_.data(),
                                            next_poiss_time_.data(),
                                            cell_begin.data(),
                                            cell_end.data(),
                                            event_buffer.data());
    /*
    // TODO: Wait for GPU to finish, then process the spikes
    unsigned cell_lid = 0;
    for (unsigned i = 0; i < spike_buffer.size(); ++i) {
        if ( i >= cell_end[i])
    }
    */
}
