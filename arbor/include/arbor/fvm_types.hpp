#pragma once

#include <arbor/common_types.hpp>

#include <arbor/arb_types.h>

namespace arb {

using fvm_value_type = arb_value_type;
using fvm_size_type  = arb_size_type;
using fvm_index_type = arb_index_type;

struct fvm_gap_junction {
    using value_type = fvm_value_type;
    using index_type = fvm_index_type;

    std::pair<index_type, index_type> loc;
    value_type weight;

    fvm_gap_junction() {}
    fvm_gap_junction(std::pair<index_type, index_type> l, value_type w): loc(l), weight(w) {}

};

} // namespace arb
