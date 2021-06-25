#pragma once

#include <arbor/gpu/gpu_common.hpp>
#include "stack_storage.hpp"

namespace arb {
namespace gpu {

template <typename T>
__device__
void push_back(stack_storage<T>& s, const T& value) {
    // Atomically increment the stores counter. The atomicAdd returns
    // the value of stores before the increment, which is the location
    // at which this thread can store value.
    unsigned position = atomicAdd(&(s.stores), 1u);

    // It is possible that stores>capacity. In this case, only capacity
    // entries are stored, and additional values are lost. The stores
    // contains the total number of attempts to push.
    if (position<s.capacity) {
        s.data[position] = value;
    }

    // Note: there are no guards against s.stores overflowing: in which
    // case the stores counter would start again from 0, and values would
    // be overwritten from the front of the stack.
}

} // namespace gpu
} // namespace arb
