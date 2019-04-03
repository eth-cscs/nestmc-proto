#include <string>
#include <sstream>

#include <arbor/common_types.hpp>
#include <arbor/context.hpp>

#include "strings.hpp"

namespace pyarb {

std::string cell_member_string(const arb::cell_member_type& m) {
    std::stringstream s;
    s << "<cell_member: gid " << m.gid
      << ", index " << m.index << ">";
    return s.str();
}

std::string context_string(const arb::context& c) {
    std::stringstream s;
    const bool gpu = arb::has_gpu(c);
    const bool mpi = arb::has_mpi(c);
    s << "<context: threads " << arb::num_threads(c)
      << ", gpu " << (gpu? "yes": "None")
      << ", distributed " << (mpi? "MPI": "Local")
      << " ranks " << arb::num_ranks(c)
      << ">";
    return s.str();
}

std::string proc_allocation_string(const arb::proc_allocation& a) {
    std::stringstream s;
    s << "<hardware resource allocation: threads " << a.num_threads << ", gpu ";
    if (a.has_gpu()) {
        s << a.gpu_id;
    }
    else {
        s << "None";
    }
    s << ">";
    return s.str();
}

} // namespace pyarb
