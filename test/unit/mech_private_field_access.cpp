#include <arbor/version.hpp>

#include <arbor/mechanism.hpp>
#include "backends/multicore/fvm.hpp"
#include "util/maputil.hpp"

#ifdef ARB_GPU_ENABLED
#include "backends/gpu/fvm.hpp"
#include "memory/gpu_wrappers.hpp"
#endif

#include "common.hpp"
#include "mech_private_field_access.hpp"

using namespace arb;

// Multicore mechanisms:

ACCESS_BIND(mechanism_field_table (mechanism::*)(), field_table_ptr, &mechanism::field_table);
ACCESS_BIND(arb_mechanism_type  mechanism::*, mech_type_ptr,  &mechanism::mech_);
ACCESS_BIND(arb_mechanism_ppack mechanism::*, mech_ppack_ptr, &mechanism::ppack_);

std::vector<fvm_value_type> mc_mechanism_field(mechanism* m, const std::string& key) {
    auto opt_ptr = util::value_by_key((m->*field_table_ptr)(), key);
    if (!opt_ptr) throw std::logic_error("internal error: no such field in mechanism");

    const fvm_value_type* field_data = opt_ptr.value().first;
    auto ppack = m->*mech_ppack_ptr;
    return std::vector<fvm_value_type>(field_data, field_data+ppack.width);
}

// GPU mechanisms:

#ifdef ARB_GPU_ENABLED

std::vector<fvm_value_type> gpu_mechanism_field(mechanism* m, const std::string& key) {
    auto opt_ptr = util::value_by_key((m->*field_table_ptr)(), key);
    if (!opt_ptr) throw std::logic_error("internal error: no such field in mechanism");

    const fvm_value_type* field_data = opt_ptr.value().first;
    std::vector<fvm_value_type> values(m->ppack_.width);

    memory::gpu_memcpy_d2h(values.data(), field_data, sizeof(fvm_value_type)*m->ppack_.width);
    return values;
}
#endif

// Generic access:

std::vector<fvm_value_type> mechanism_field(mechanism* m, const std::string& key) {
    if (m->iface_.backend == arb_backend_kind_cpu) {
        return mc_mechanism_field(m, key);
    }

#ifdef ARB_GPU_ENABLED
    if (m->iface_.backend == arb_backend_kind_gpu) {
        return gpu_mechanism_field(m, key);
    }
#endif

    throw std::logic_error("internal error: mechanism instantiated on unknown backend");
}

