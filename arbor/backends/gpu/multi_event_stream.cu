#include <arbor/common_types.hpp>
#include <arbor/gpu/gpu_common.hpp>

#include "backends/event.hpp"

namespace arb {
namespace gpu {

namespace kernels {
    template <typename T, typename I>
    __global__ void mark_until_after(
        unsigned n,
        I* __restrict__ const mark,
        const I* __restrict__ const span_end,
        const T* __restrict__ const ev_time,
        const T* __restrict__ const t_until)
    {
        unsigned i = threadIdx.x+blockIdx.x*blockDim.x;
        if (i<n) {
            auto t = t_until[i];
            auto end = span_end[i];
            auto &m = mark[i];

            while (m!=end && !(ev_time[m]>t)) {
                ++m;
            }
        }
    }

    template <typename T, typename I>
    __global__ void mark_until(
        unsigned n,
        I* __restrict__ const mark,
        const I* __restrict__ const span_end,
        const T* __restrict__ const ev_time,
        const T* __restrict__ const t_until)
    {
        unsigned i = threadIdx.x+blockIdx.x*blockDim.x;
        if (i<n) {
            auto t = t_until[i];
            auto end = span_end[i];
            auto &m = mark[i];

            while (m!=end && t>ev_time[m]) {
                ++m;
            }
        }
    }

    template <typename I>
    __global__ void drop_marked_events(
        unsigned n,
        I* __restrict__ const n_nonempty,
        I* __restrict__ const span_begin,
        const I* __restrict__ const span_end,
        const I* __restrict__ const mark)
    {
        unsigned i = threadIdx.x+blockIdx.x*blockDim.x;
        if (i<n) {
            bool emptied = (span_begin[i]<span_end[i] && mark[i]==span_end[i]);
            span_begin[i] = mark[i];
            if (emptied) {
                atomicAdd(n_nonempty, (I)-1);
            }
        }
    }

    template <typename T, typename I>
    __global__ void event_time_if_before(
        unsigned n,
        const I* __restrict__ const span_begin,
        const I* __restrict__ const span_end,
        const T* __restrict__ const ev_time,
        T* __restrict__ const t_until)
    {
        unsigned i = threadIdx.x+blockIdx.x*blockDim.x;
        if (i<n) {
            if (span_begin[i]<span_end[i]) {
                auto ev_t = ev_time[span_begin[i]];
                if (t_until[i]>ev_t) {
                    t_until[i] = ev_t;
                }
            }
        }
    }
} // namespace kernels

void mark_until_after_w(unsigned n,
        fvm_index_type* mark,
        fvm_index_type* span_end,
        fvm_value_type* ev_time,
        const fvm_value_type* t_until)
{
    const int nblock = impl::block_count(n, 128);
    kernels::mark_until_after
        <<<nblock, 128>>>
        (n, mark, span_end, ev_time, t_until);
}

void mark_until_w(unsigned n,
        fvm_index_type* mark,
        fvm_index_type* span_end,
        fvm_value_type* ev_time,
        const fvm_value_type* t_until)
{
    const int nblock = impl::block_count(n, 128);
    kernels::mark_until
        <<<nblock, 128>>>
        (n, mark, span_end, ev_time, t_until);
}

void drop_marked_events_w(unsigned n,
        fvm_index_type* n_nonempty_stream,
        fvm_index_type* span_begin,
        fvm_index_type* span_end,
        fvm_index_type* mark)
{
    const int nblock = impl::block_count(n, 128);
    kernels::drop_marked_events
        <<<nblock, 128>>>
        (n, n_nonempty_stream, span_begin, span_end, mark);

}

void event_time_if_before_w(unsigned n,
        fvm_index_type* span_begin,
        fvm_index_type* span_end,
        fvm_value_type* ev_time,
        fvm_value_type* t_until)
{
    const int nblock = impl::block_count(n, 128);
    kernels::event_time_if_before
        <<<nblock, 128>>>
        (n, span_begin, span_end, ev_time, t_until);
}

} // namespace gpu
} // namespace arb
