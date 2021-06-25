#include <ostream>
#include <set>
#include <string>
#include <regex>

#define FMT_HEADER_ONLY YES
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/compile.h>

#include "blocks.hpp"
#include "infoprinter.hpp"
#include "module.hpp"
#include "printerutil.hpp"

#include "io/ostream_wrappers.hpp"
#include "io/prefixbuf.hpp"

using io::quote;

std::string build_info_header(const Module& m, const printer_options& opt, bool cpu, bool gpu) {
    using io::indent;
    using io::popindent;

    std::string name = m.module_name();

    io::pfxstringstream out;

    std::string fingerprint = "<placeholder>";

    out << fmt::format("#pragma once\n\n"
                       "#include <cmath>\n"
                       "#include <{}mechanism_abi.h>\n\n",
                       arb_header_prefix());

    auto vars = local_module_variables(m);
    auto ion_deps = m.ion_deps();


    std::unordered_map<std::string, Id> name2id;
    for (const auto& id: m.parameter_block().parameters) name2id[id.name()] = id;
    for (const auto& id: m.state_block().state_variables) name2id[id.name()] = id;

    auto fmt_var = [&](const auto& v) {
        auto kv = name2id.find(v->name());
        auto lo = std::numeric_limits<double>::lowest();
        auto hi = std::numeric_limits<double>::max();
        std::string unit = "";
        if (kv != name2id.end()) {
            auto id = kv->second;
            unit = id.unit_string();
            if (id.has_range()) {
                auto lo = id.range.first;
                auto hi = id.range.second;
            }
        }
        return fmt::format("{{ \"{}\", \"{}\", {}, {}, {} }}",
                           v->name(),
                           unit,
                           std::isnan(v->value()) ? "NAN" : std::to_string(v->value()),
                           lo, hi);
    };

    auto fmt_ion = [](const auto& i) { return fmt::format(FMT_COMPILE("{{ \"{}\", {}, {}, {}, {}, {}, {}, {}}}"),
                                                          i.name,
                                                          i.writes_concentration_int() ? "true" : "false",
                                                          i.writes_concentration_ext() ? "true" : "false",
                                                          i.writes_rev_potential() ? "true" : "false",
                                                          i.uses_rev_potential() ? "true" : "false",
                                                          i.uses_valence() ? "true" : "false",
                                                          i.verifies_valence() ? "true" : "false",
                                                          i.expected_valence); };


    out << fmt::format("extern \"C\" {{\n"
                       "  arb_mechanism_type make_{0}_{1}() {{\n",
                       std::regex_replace(opt.cpp_namespace, std::regex{"::"}, "_"),
                       name);

    out << "    // Tables\n";
    {
        auto n = 0ul;
        io::separator sep("", ",\n                                        ");
        out << "    static arb_field_info globals[] = { ";
        for (const auto& var: vars.scalars) {
            out << sep << fmt_var(var);
            ++n;
        }
        out << " };\n"
            << "    static arb_size_type n_globals = " << n << ";\n";
    }

    {
        auto n = 0ul;
        io::separator sep("", ",\n                                           ");
        out << "    static arb_field_info state_vars[] = { ";
        for (const auto& var: vars.arrays) {
            if(var->is_state()) {
                out << sep << fmt_var(var);
                ++n;
            }
        }
        out << " };\n"
            << "    static arb_size_type n_state_vars = " << n << ";\n";
    }
    {
        auto n = 0ul;
        io::separator sep("", ",\n                                           ");
        out << "    static arb_field_info parameters[] = { ";
        for (const auto& var: vars.arrays) {
            if(!var->is_state()) {
                out << sep << fmt_var(var);
                ++n;
            }
        }
        out << " };\n"
            << "    static arb_size_type n_parameters = " << n << ";\n";
    }

    {
        io::separator sep("", ",\n");
        out << "    static arb_ion_info ions[] = { ";
        auto n = 0ul;
        for (const auto& var: ion_deps) {
            out << sep << fmt_ion(var);
            ++n;
        }
        out << " };\n"
            << "    static arb_size_type n_ions = " << n << ";\n";
    }

    out << fmt::format(FMT_COMPILE("\n"
                                   "    arb_mechanism_type result;\n"
                                   "    result.abi_version=ARB_MECH_ABI_VERSION;\n"
                                   "    result.fingerprint=\"{1}\";\n"
                                   "    result.name=\"{0}\";\n"
                                   "    result.kind={2};\n"
                                   "    result.is_linear={3};\n"
                                   "    result.has_post_events={4};\n"
                                   "    result.globals=globals;\n"
                                   "    result.n_globals=n_globals;\n"
                                   "    result.ions=ions;\n"
                                   "    result.n_ions=n_ions;\n"
                                   "    result.state_vars=state_vars;\n"
                                   "    result.n_state_vars=n_state_vars;\n"
                                   "    result.parameters=parameters;\n"
                                   "    result.n_parameters=n_parameters;\n"
                                   "    return result;\n"
                                   "  }}\n"
                                   "\n"),
                       name,
                       fingerprint,
                       module_kind_str(m),
                       m.is_linear(),
                       m.has_post_events())
        << fmt::format("  arb_mechanism_interface* make_{0}_{1}_interface_multicore(){2}\n"
                       "  arb_mechanism_interface* make_{0}_{1}_interface_gpu(){3}\n"
                       "}}\n",
                       std::regex_replace(opt.cpp_namespace, std::regex{"::"}, "_"),
                       name,
                       cpu ? ";" : " { return nullptr; }",
                       gpu ? ";" : " { return nullptr; }");
    return out.str();
}
