#include <cmath>

#include <arbor/fvm_types.hpp>
#include <arbor/gpu/gpu_api.hpp>
#include <arbor/gpu/gpu_common.hpp>
#include <arbor/gpu/math_cu.hpp>

#include "backends/gpu/stimulus.hpp"

namespace arb {
namespace gpu {

// TODO: Implement version with reproducibility-friendly accumulations.
// See Arbor issue #1059.

namespace kernel {

__global__
void istim_add_current_impl(int n, istim_pp pp) {
    constexpr double two_pi = 2*pi;

    auto i = threadIdx.x + blockDim.x*blockIdx.x;
    if (i>=n) return;

    fvm_index_type ei_left = pp.envl_divs[i];
    fvm_index_type ei_right = pp.envl_divs[i+1];

    fvm_index_type ai = pp.accu_index[i];
    fvm_index_type cv = pp.accu_to_cv[ai];
    double t = pp.time[pp.cv_to_intdom[cv]];

    if (ei_left==ei_right || t<pp.envl_times[ei_left]) return;

    fvm_index_type& ei = pp.envl_index[i];
    while (ei+1<ei_right && pp.envl_times[ei+1]<=t) ++ei;

    double J = pp.envl_amplitudes[ei]; // current density (A/m²)
    if (ei+1<ei_right) {
        // linearly interpolate:
        double J1 = pp.envl_amplitudes[ei+1];
        double u = (t-pp.envl_times[ei])/(pp.envl_times[ei+1]-pp.envl_times[ei]);
        J = lerp(J, J1, u);
    }

    if (double f = pp.frequency[i]) {
         J *= std::sin(two_pi*f*t + pp.phase[i]);
    }

    gpu_atomic_add(&pp.accu_stim[ai], J);
    gpu_atomic_sub(&pp.current_density[cv], J);
}

} // namespace kernel

void istim_add_current_impl(int n, const istim_pp& pp) {
    constexpr unsigned block_dim = 128;
    const unsigned grid_dim = impl::block_count(n, block_dim);
    if (!grid_dim) return;
    kernel::istim_add_current_impl<<<grid_dim, block_dim>>>(n, pp);
}

} // namespace gpu
} // namespace arb
