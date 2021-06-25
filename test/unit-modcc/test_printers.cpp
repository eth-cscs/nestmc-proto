#include <memory>
#include <regex>
#include <string>
#include <sstream>

#include "common.hpp"
#include "io/bulkio.hpp"

#include "printer/cexpr_emit.hpp"
#include "printer/cprinter.hpp"
#include "printer/gpuprinter.hpp"
#include "expression.hpp"
#include "symdiff.hpp"

struct testcase {
    const char* source;
    const char* expected;
};

static std::string strip(std::string text) {
    // Strip all spaces, except when between two minus symbols, where instead
    // they should be replaced by a single space:
    //
    // 1. Replace whitespace with two spaces.
    // 2. Replace '-  -' with '- -'.
    // 3. Replace '  ' with ''.

    static std::regex rx1("\\s+");
    static std::regex rx2("-  -");
    static std::regex rx3("  ");

    text = std::regex_replace(text, rx1, "  ");
    text = std::regex_replace(text, rx2, "- -");
    text = std::regex_replace(text, rx3, "");

    return text;
}

TEST(scalar_printer, constants) {
    testcase testcases[] = {
        {"1./0.",      "INFINITY"},
        {"-1./0.",     "-INFINITY"},
        {"(-1)^0.5",   "NAN"},
        {"1/(-1./0.)", "-0."},
        {"1-1",        "0."},
    };

    for (const auto& tc: testcases) {
        auto expr = constant_simplify(parse_expression(tc.source));
        ASSERT_TRUE(expr && expr->is_number());

        std::stringstream s;
        s << as_c_double(expr->is_number()->value());

        EXPECT_EQ(std::string(tc.expected), s.str());
    }
}

TEST(scalar_printer, statement) {
    std::vector<testcase> testcases = {
        {"y=x+3",            "y=x+3.0"},
        {"y=y^z",            "y=pow(y,z)"},
        {"y=exp((x/2) + 3)", "y=exp(x/2.0+3.0)"},
        {"z=a/b/c",          "z=a/b/c"},
        {"z=a/(b/c)",        "z=a/(b/c)"},
        {"z=(a*b)/c",        "z=a*b/c"},
        {"z=a-(b+c)",        "z=a-(b+c)"},
        {"z=(a>0)<(b>0)",    "z=a>0.<(b>0.)"},
        {"z=a- -2",          "z=a- -2.0"},
        {"z=fabs(x-z)",      "z=abs(x-z)"},
        {"z=min(x,y)",       "z=min(x,y)"},
        {"z=min(max(a,b),y)","z=min(max(a,b),y)"},
    };

    // create a scope that contains the symbols used in the tests
    Scope<Symbol>::symbol_map globals;
    auto scope = std::make_shared<Scope<Symbol>>(globals);

    for (auto var: {"x", "y", "z", "a", "b", "c"}) {
        scope->add_local_symbol(var, make_symbol<LocalVariable>(Location(), var, localVariableKind::local));
    }

    for (const auto& tc: testcases) {
        auto e = parse_line_expression(tc.source);
        ASSERT_TRUE(e);

        e->semantic(scope);
        if(e->has_error()) {
            std::cerr << e->error_message() << std::endl;
            FAIL();
        }
        {
            SCOPED_TRACE("CPrinter");
            std::stringstream out;
            auto printer = std::make_unique<CPrinter>(out);
            e->accept(printer.get());
            std::string text = out.str();

            verbose_print(e->to_string(), " :--: ", text);
            EXPECT_EQ(strip(tc.expected), strip(text));
        }

        {
            SCOPED_TRACE("GpuPrinter");
            std::stringstream out;
            auto printer = std::make_unique<GpuPrinter>(out);
            e->accept(printer.get());
            std::string text = out.str();

            verbose_print(e->to_string(), " :--: ", text);
            EXPECT_EQ(strip(tc.expected), strip(text));
        }
    }
}

TEST(CPrinter, proc_body) {
    std::vector<testcase> testcases = {
        {
            "PROCEDURE trates(v) {\n"
            "    LOCAL k\n"
            "    minf = 1-1/(1+exp((v-k)/k))\n"
            "    hinf = 1/(1+exp((v-k)/k))\n"
            "    mtau = 0.5\n"
            "    htau = 1500\n"
            "}"
            ,
            "arb_value_type k;\n"
            "_pp_var_minf[i_] = 1.0-1.0/(1.0+exp((v-k)/k));\n"
            "_pp_var_hinf[i_] = 1.0/(1.0+exp((v-k)/k));\n"
            "_pp_var_mtau[i_] = 0.5;\n"
            "_pp_var_htau[i_] = 1500.0;\n"
        }
    };

    // create a scope that contains the symbols used in the tests
    Scope<Symbol>::symbol_map globals;
    globals["minf"] = make_symbol<VariableExpression>(Location(), "minf");
    globals["hinf"] = make_symbol<VariableExpression>(Location(), "hinf");
    globals["mtau"] = make_symbol<VariableExpression>(Location(), "mtau");
    globals["htau"] = make_symbol<VariableExpression>(Location(), "htau");
    globals["v"]    = make_symbol<VariableExpression>(Location(), "v");

    for (const auto& tc: testcases) {
        expression_ptr e = parse_procedure(tc.source);
        ASSERT_TRUE(e->is_symbol());

        auto procname = e->is_symbol()->name();
        auto& proc = (globals[procname] = symbol_ptr(e.release()->is_symbol()));

        proc->semantic(globals);
        std::stringstream out;
        auto v = std::make_unique<CPrinter>(out);
        proc->is_procedure()->body()->accept(v.get());
        std::string text = out.str();

        verbose_print(proc->is_procedure()->body()->to_string());
        verbose_print(" :--: ", text);

        EXPECT_EQ(strip(tc.expected), strip(text));
    }
}

TEST(CPrinter, proc_body_const) {
    std::vector<testcase> testcases = {
            {
                    "PROCEDURE trates(v) {\n"
                    "    mtau = 0.5 - t0 + t1\n"
                    "}"
                    ,
                    "_pp_var_mtau[i_] = 0.5 - -0.5 + 1.2;\n"
            }
    };

    // create a scope that contains the symbols used in the tests
    Scope<Symbol>::symbol_map globals;
    globals["mtau"] = make_symbol<VariableExpression>(Location(), "mtau");

    for (const auto& tc: testcases) {
        Parser p(tc.source);
        p.constants_map_.insert({"t0","-0.5"});
        p.constants_map_.insert({"t1","1.2"});
        expression_ptr e = p.parse_procedure();
        ASSERT_TRUE(e->is_symbol());

        auto procname = e->is_symbol()->name();
        auto& proc = (globals[procname] = symbol_ptr(e.release()->is_symbol()));

        proc->semantic(globals);
        std::stringstream out;
        auto v = std::make_unique<CPrinter>(out);
        proc->is_procedure()->body()->accept(v.get());
        std::string text = out.str();

        verbose_print(proc->is_procedure()->body()->to_string());
        verbose_print(" :--: ", text);

        EXPECT_EQ(strip(tc.expected), strip(text));
    }
}

TEST(CPrinter, proc_body_inlined) {
    const char* expected =
        "ll0_ = 0.;\n"
        "r_6_ = 0.;\n"
        "r_7_ = 0.;\n"
        "r_8_ = 0.;\n"
        "r_9_=_pp_var_s2[i_]*0.33333333333333331;\n"
        "r_8_=_pp_var_s1[i_]+2.0;\n"
        "if(_pp_var_s1[i_]==3.0){\n"
        "   r_7_=2.0*r_8_;\n"
        "}\n"
        "else{\n"
        "   if(_pp_var_s1[i_]==4.0){\n"
        "       r_11_ = 0.;\n"
        "       r_12_ = 0.;\n"
        "       r_12_=6.0+_pp_var_s1[i_];\n"
        "       r_11_=r_12_;\n"
        "       r_7_=r_8_*r_11_;\n"
        "   }\n"
        "   else{\n"
        "       r_10_=exp(r_8_);\n"
        "       r_7_=r_10_*_pp_var_s1[i_];\n"
        "   }\n"
        "}\n"
        "r_13_=0.;\n"
        "r_14_=0.;\n"
        "r_14_=r_9_/_pp_var_s2[i_];\n"
        "r_15_=log(r_14_);\n"
        "r_13_=42.0*r_15_;\n"
        "r_6_=r_9_*r_13_;\n"
        "t0=r_7_*r_6_;\n"
        "t1=exprelr(t0);\n"
        "ll0_=t1+2.0;\n"
        "if(ll0_==3.0){\n"
        "   t2=10.0;\n"
        "}\n"
        "else{\n"
        "   if(ll0_==4.0){\n"
        "       r_17_=0.;\n"
        "       r_18_=0.;\n"
        "       r_18_=6.0+ll0_;\n"
        "       r_17_=r_18_;\n"
        "       t2=5.0*r_17_;\n"
        "   }\n"
        "   else{\n"
        "       r_16_=148.4131591025766;\n"
        "       t2=r_16_*ll0_;\n"
        "   }\n"
        "}\n"
        "_pp_var_s2[i_]=t2+4.0;\n";

    Module m(io::read_all(DATADIR "/mod_files/test6.mod"), "test6.mod");
    Parser p(m, false);
    p.parse();
    m.semantic();

    auto& proc_rates = m.symbols().at("rates");

    ASSERT_TRUE(proc_rates->is_symbol());

    std::stringstream out;
    auto v = std::make_unique<CPrinter>(out);
    proc_rates->is_procedure()->body()->accept(v.get());
    std::string text = out.str();

    verbose_print(proc_rates->is_procedure()->body()->to_string());
    verbose_print(" :--: ", text);

    // Remove the first statement that declares the locals
    // Their print order is not fixed
    auto proc_with_locals = strip(text);
    proc_with_locals.erase(0, proc_with_locals.find(";") + 1);

    EXPECT_EQ(strip(expected), proc_with_locals);
}

TEST(SimdPrinter, simd_if_else) {
    std::vector<const char*> expected_procs = {
            "simd_value u;\n"
            "simd_mask mask_0_ = S::cmp_gt(i, (double)2.0);\n"
            "S::where(mask_0_,u) = (double)7.0;\n"
            "S::where(S::logical_not(mask_0_),u) = (double)5.0;\n"
            "indirect(_pp_var_s+i_, simd_width_) = S::where(S::logical_not(mask_0_),simd_cast<simd_value>((double)42.0));\n"
            "indirect(_pp_var_s+i_, simd_width_) = u;"
            ,
            "simd_value u;\n"
            "simd_mask mask_1_ = S::cmp_gt(i, (double)2.0);\n"
            "S::where(mask_1_,u) = (double)7.0;\n"
            "S::where(S::logical_not(mask_1_),u) = (double)5.0;\n"
            "indirect(_pp_var_s+i_, simd_width_) = S::where(S::logical_and(S::logical_not(mask_1_), mask_input_),simd_cast<simd_value>((double)42.0));\n"
            "indirect(_pp_var_s+i_, simd_width_) = S::where(mask_input_, u);"
            ,
            "simd_mask mask_2_ = S::cmp_gt(simd_cast<simd_value>(indirect(_pp_var_g+i_, simd_width_)), (double)2.0);\n"
            "simd_mask mask_3_ = S::cmp_gt(simd_cast<simd_value>(indirect(_pp_var_g+i_, simd_width_)), (double)3.0);\n"
            "S::where(S::logical_and(mask_2_,mask_3_),i) = (double)0.;\n"
            "S::where(S::logical_and(mask_2_,S::logical_not(mask_3_)),i) = (double)1.0;\n"
            "simd_mask mask_4_ = S::cmp_lt(simd_cast<simd_value>(indirect(_pp_var_g+i_, simd_width_)), (double)1.0);\n"
            "indirect(_pp_var_s+i_, simd_width_) = S::where(S::logical_and(S::logical_not(mask_2_),mask_4_),simd_cast<simd_value>((double)2.0));\n"
            "rates(i_, S::logical_and(S::logical_not(mask_2_),S::logical_not(mask_4_)), i);"
    };

    Module m(io::read_all(DATADIR "/mod_files/test7.mod"), "test7.mod");
    Parser p(m, false);
    p.parse();
    m.semantic();

    struct proc {
        std::string name;
        bool masked;
    };

    std::vector<proc> procs = {{"rates", false}, {"rates", true}, {"foo", false}};
    for (unsigned i = 0; i < procs.size(); i++) {
        auto p = procs[i];
        std::stringstream out;
        auto &proc = m.symbols().at(p.name);
        ASSERT_TRUE(proc->is_symbol());

        auto v = std::make_unique<SimdPrinter>(out);
        if (p.masked) {
            v->set_input_mask("mask_input_");
        }
        proc->is_procedure()->body()->accept(v.get());
        std::string text = out.str();

        verbose_print(proc->is_procedure()->body()->to_string());
        verbose_print(" :--: ", text);

        auto proc_with_locals = strip(text);
        EXPECT_EQ(strip(expected_procs[i]), proc_with_locals);

    }
}
