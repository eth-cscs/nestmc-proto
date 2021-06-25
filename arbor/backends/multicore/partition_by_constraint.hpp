#pragma once

#include <vector>

#include <arbor/simd/simd.hpp>

namespace arb {
namespace multicore {

namespace S = ::arb::simd;
using S::index_constraint;

struct constraint_partition {
    using iarray = arb::multicore::iarray;

    iarray contiguous;
    iarray constant;
    iarray independent;
    iarray none;
};

template <typename It>
bool is_contiguous_n(It first, unsigned width) {
    while (--width) {
        It next = first;
        if ((*first) +1 != *++next) {
            return false;
        }
        first = next;
    }
    return true;
}

template <typename It>
bool is_constant_n(It first, unsigned width) {
    while (--width) {
        It next = first;
        if (*first != *++next) {
            return false;
        }
        first = next;
    }
    return true;
}

template <typename It>
bool is_independent_n(It first, unsigned width) {
    while (--width) {
        It next = first;
        if (*first == *++next) {
            return false;
        }
        first = next;
    }
    return true;
}

template <typename It>
index_constraint idx_constraint(It it, unsigned simd_width) {
    if (is_contiguous_n(it, simd_width)) {
        return index_constraint::contiguous;
    }
    else if (is_constant_n(it, simd_width)) {
        return index_constraint::constant;
    }
    else if (is_independent_n(it, simd_width)) {
        return index_constraint::independent;
    }
    else {
        return index_constraint::none;
    }
}

template <typename T>
constraint_partition make_constraint_partition(const T& node_index, unsigned width, unsigned simd_width) {
    constraint_partition part;
    if (simd_width) {
        for (unsigned i = 0; i < width; i+= simd_width) {
            auto ptr = &node_index[i];
            if (is_contiguous_n(ptr, simd_width)) {
                part.contiguous.push_back(i);
            }
            else if (is_constant_n(ptr, simd_width)) {
                part.constant.push_back(i);
            }
            else if (is_independent_n(ptr, simd_width)) {
                part.independent.push_back(i);
            }
            else {
                part.none.push_back(i);
            }
        }
    }
    return part;
}

bool constexpr is_constraint_stronger(index_constraint a, index_constraint b) {
    return a==b ||
           a==index_constraint::none ||
           (a==index_constraint::independent && b==index_constraint::contiguous);
}

template <typename T, typename U>
bool compatible_index_constraints(const T& node_index, const U& ion_index, unsigned simd_width){
    for (unsigned i = 0; i < node_index.size(); i+= simd_width) {
        auto nc = idx_constraint(&node_index[i], simd_width);
        auto ic = idx_constraint(&ion_index[i], simd_width);

        if (!is_constraint_stronger(nc, ic)) {
            return false;
        }
    }
    return true;
}

} // namespace util
} // namespace arb



