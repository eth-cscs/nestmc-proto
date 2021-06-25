#include <string>

#include <arbor/mechanism.hpp>
#include "fvm.hpp"

// Provides implementation of backend::mechanism_field_data.

namespace arb {
namespace gpu {

fvm_value_type* backend::mechanism_field_data(arb::mechanism* m, const std::string& field) {
    return m? m->field_data(field): nullptr;
}

void multiply_in_place_(fvm_value_type* s, const fvm_index_type* p, int n);

void backend::multiply_in_place(fvm_value_type* s, const fvm_index_type* p, int n) {
    multiply_in_place_(s, p, n);
}


} // namespace gpu
} // namespace arb
