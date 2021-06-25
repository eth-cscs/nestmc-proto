#include <string>

#include <arbor/mechanism.hpp>

#include "fvm.hpp"
#include "util/span.hpp"

// Provides implementation of backend::mechanism_field_data.

namespace arb {
namespace multicore {

fvm_value_type* backend::mechanism_field_data(arb::mechanism* mptr, const std::string& field) {
    return mptr ? mptr->field_data(field): nullptr;
}

} // namespace multicore
} // namespace arb
