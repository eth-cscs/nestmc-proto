#include <cmath>
#include <iomanip>
#include <ostream>
#include <unordered_map>

#include "printerutil.hpp"
#include "cexpr_emit.hpp"
#include "error.hpp"
#include "lexer.hpp"
#include "io/ostream_wrappers.hpp"
#include "astmanip.hpp"
#include "io/prefixbuf.hpp"

std::ostream& operator<<(std::ostream& out, as_c_double wrap) {
    bool neg = std::signbit(wrap.value);

    switch (std::fpclassify(wrap.value)) {
    case FP_INFINITE:
        return out << (neg? "-": "") << "INFINITY";
    case FP_NAN:
        return out << "NAN";
    case FP_ZERO:
        return out << (neg? "-0.": "0.");
    default:
        double val;
        std::stringstream s;
        // If wrap.value is an int print it as X.0, this is needed for std::max and std::min
        if (std::modf(wrap.value, &val) == 0) {
            s << io::classic << std::fixed << std::setprecision(1) << wrap.value;
        } else {
            s << io::classic << std::setprecision(17) << wrap.value;
        }
        return out << s.rdbuf();
    }
}

void CExprEmitter::emit_as_call(const char* sub, Expression* e) {
    out_ << sub << '(';
    e->accept(this);
    out_ << ')';
}

void CExprEmitter::emit_as_call(const char* sub, Expression* e1, Expression* e2) {
    out_ << sub << '(';
    e1->accept(this);
    out_ << ", ";
    e2->accept(this);
    out_ << ')';
}

void CExprEmitter::visit(NumberExpression* e) {
    out_ << " " << as_c_double(e->value());
}

void CExprEmitter::visit(UnaryExpression* e) {
    // Place a space in front of minus sign to avoid invalid
    // expressions of the form: (v[i]--67)
    static std::unordered_map<tok, const char*> unaryop_tbl = {
        {tok::minus,   " -"},
        {tok::exp,     "exp"},
        {tok::cos,     "cos"},
        {tok::sin,     "sin"},
        {tok::log,     "log"},
        {tok::abs,     "abs"},
        {tok::exprelr, "exprelr"},
        {tok::safeinv, "safeinv"}
    };

    if (!unaryop_tbl.count(e->op())) {
        throw compiler_exception(
            "CExprEmitter: unsupported unary operator "+token_string(e->op()), e->location());
    }

    const char* op_spelling = unaryop_tbl.at(e->op());
    Expression* inner = e->expression();

    // No need to use parenthesis for unary minus if inner expression is
    // not binary.
    if (e->op()==tok::minus && !inner->is_binary()) {
        out_ << op_spelling;
        inner->accept(this);
    }
    else {
        emit_as_call(op_spelling, inner);
    }
}

void CExprEmitter::visit(AssignmentExpression* e) {
    e->lhs()->accept(this);
    out_ << " = ";
    e->rhs()->accept(this);
}

void CExprEmitter::visit(PowBinaryExpression* e) {
    emit_as_call("pow", e->lhs(), e->rhs());
}

void CExprEmitter::visit(BinaryExpression* e) {
    static std::unordered_map<tok, const char*> binop_tbl = {
        {tok::minus,    "-"},
        {tok::plus,     "+"},
        {tok::times,    "*"},
        {tok::divide,   "/"},
        {tok::lt,       "<"},
        {tok::lte,      "<="},
        {tok::gt,       ">"},
        {tok::gte,      ">="},
        {tok::equality, "=="},
        {tok::land,     "&&"},
        {tok::lor,      "||"},
        {tok::ne,       "!="},
        {tok::min,      "min"},
        {tok::max,      "max"},
    };

    if (!binop_tbl.count(e->op())) {
        throw compiler_exception(
            "CExprEmitter: unsupported binary operator "+token_string(e->op()), e->location());
    }

    auto rhs = e->rhs();
    auto lhs = e->lhs();
    const char* op_spelling = binop_tbl.at(e->op());

    if (e->is_infix()) {
        associativityKind assoc = Lexer::operator_associativity(e->op());
        int op_prec = Lexer::binop_precedence(e->op());

        auto need_paren = [op_prec](Expression* subexpr, bool assoc_side) -> bool {
            if (auto b = subexpr->is_binary()) {
                int sub_prec = Lexer::binop_precedence(b->op());
                return sub_prec<op_prec || (!assoc_side && sub_prec==op_prec);
            }
            return false;
        };

        if (need_paren(lhs, assoc==associativityKind::left)) {
            emit_as_call("", lhs);
        }
        else {
            lhs->accept(this);
        }

        out_ << op_spelling;

        if (need_paren(rhs, assoc==associativityKind::right)) {
            emit_as_call("", rhs);
        }
        else {
            rhs->accept(this);
        }
    }
    else {
        emit_as_call(op_spelling, lhs, rhs);
    }
}

void CExprEmitter::visit(IfExpression* e) {
    out_ << "if (";
    e->condition()->accept(this);
    out_ << ") {\n" << io::indent;
    e->true_branch()->accept(this);
    out_ << io::popindent << "}\n";

    if (auto fb = e->false_branch()) {
        out_ << "else ";
        if (fb->is_if()) {
            fb->accept(this);
        }
        else {
            out_ << "{\n" << io::indent;
            fb->accept(this);
            out_ << io::popindent << "}\n";
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::unordered_set<std::string> SimdExprEmitter::mask_names_;

void SimdExprEmitter::visit(PowBinaryExpression* e) {
    out_ << "S::pow(";
    e->lhs()->accept(this);
    out_ << ", ";
    e->rhs()->accept(this);
    out_ << ')';
}

void SimdExprEmitter::visit(NumberExpression* e) {
    out_ << " (double)" << as_c_double(e->value());
} 

void SimdExprEmitter::visit(UnaryExpression* e) {
    static std::unordered_map<tok, const char*> unaryop_tbl = {
        {tok::minus,   "S::neg"},
        {tok::exp,     "S::exp"},
        {tok::cos,     "S::cos"},
        {tok::sin,     "S::sin"},
        {tok::log,     "S::log"},
        {tok::abs,     "S::abs"},
        {tok::exprelr, "S::exprelr"},
        {tok::safeinv, "safeinv"}
    };

    if (!unaryop_tbl.count(e->op())) {
        throw compiler_exception(
            "CExprEmitter: unsupported unary operator "+token_string(e->op()), e->location());
    }

    const char* op_spelling = unaryop_tbl.at(e->op());
    Expression* inner = e->expression();

    auto iden = inner->is_identifier();
    bool is_scalar = iden && scalars_.count(iden->name()); 
    if (e->op()==tok::minus && is_scalar) {
        out_ << "simd_cast<simd_value>(-";
        inner->accept(this);
        out_ << ")";
    }
    else {
        emit_as_call(op_spelling, inner);
    }
}

std::string id_prefix(IdentifierExpression* id) {
    if (id) {
        if (auto symbol = id->symbol()->is_symbol()) {
            if (auto var = symbol->is_variable()) {
                if (!var->is_local_variable()) {
                    return pp_var_pfx + id->name();
                }
            }
        }
    }
    return id->name();
}


void SimdExprEmitter::visit(BinaryExpression* e) {
    static std::unordered_map<tok, const char *> func_tbl = {
            {tok::minus,    "S::sub"},
            {tok::plus,     "S::add"},
            {tok::times,    "S::mul"},
            {tok::divide,   "S::div"},
            {tok::lt,       "S::cmp_lt"},
            {tok::lte,      "S::cmp_leq"},
            {tok::gt,       "S::cmp_gt"},
            {tok::gte,      "S::cmp_geq"},
            {tok::equality, "S::cmp_eq"},
            {tok::land,     "S::logical_and"},
            {tok::lor,      "S::logical_or"},
            {tok::ne,       "S::cmp_neq"},
            {tok::min,      "S::min"},
            {tok::max,      "S::max"},
    };

    static std::unordered_map<tok, const char *> binop_tbl = {
            {tok::minus,    "-"},
            {tok::plus,     "+"},
            {tok::times,    "*"},
            {tok::divide,   "/"},
            {tok::lt,       "<"},
            {tok::lte,      "<="},
            {tok::gt,       ">"},
            {tok::gte,      ">="},
            {tok::equality, "=="},
            {tok::land,     "&&"},
            {tok::lor,      "||"},
            {tok::ne,       "!="},
            {tok::min,      "min"},
            {tok::max,      "max"},
    };


    if (!binop_tbl.count(e->op())) {
        throw compiler_exception(
                "CExprEmitter: unsupported binary operator " + token_string(e->op()), e->location());
    }

    std::string rhs_name, lhs_name, rhs_pfxd, lhs_pfxd;

    auto rhs = e->rhs();
    auto lhs = e->lhs();

    const char *op_spelling = binop_tbl.at(e->op());
    const char *func_spelling = func_tbl.at(e->op());

    if (auto id = rhs->is_identifier()) {
        rhs_name = id->name();
        rhs_pfxd = id_prefix(id);
    }
    if (auto id = lhs->is_identifier()) {
        lhs_name = id->name();
        lhs_pfxd = id_prefix(id);
    }

    if (scalars_.count(rhs_name) && scalars_.count(lhs_name)) {
        if (e->is_infix()) {
            associativityKind assoc = Lexer::operator_associativity(e->op());
            int op_prec = Lexer::binop_precedence(e->op());

            auto need_paren = [op_prec](Expression *subexpr, bool assoc_side) -> bool {
                if (auto b = subexpr->is_binary()) {
                    int sub_prec = Lexer::binop_precedence(b->op());
                    return sub_prec < op_prec || (!assoc_side && sub_prec == op_prec);
                }
                return false;
            };

            out_ << "simd_cast<simd_value>(";
            if (need_paren(lhs, assoc == associativityKind::left)) {
                emit_as_call("", lhs);
            } else {
                lhs->accept(this);
            }

            out_ << op_spelling;

            if (need_paren(rhs, assoc == associativityKind::right)) {
                emit_as_call("", rhs);
            } else {
                rhs->accept(this);
            }
            out_ << ")";
        } else {
            out_ << "simd_cast<simd_value>(";
            emit_as_call(op_spelling, lhs, rhs);
            out_ << ")";
        }
    } else if (scalars_.count(rhs_name) && !scalars_.count(lhs_name)) {
        out_ << func_spelling << '(';
        lhs->accept(this);
        out_ << ", simd_cast<simd_value>(" << rhs_pfxd;
        out_ << "))";
    } else if (!scalars_.count(rhs_name) && scalars_.count(lhs_name)) {
        out_ << func_spelling << "(simd_cast<simd_value>(" << lhs_pfxd << "), ";
        rhs->accept(this);
        out_ << ")";
    } else {
        out_ << func_spelling << '(';
        lhs->accept(this);
        out_ << ", ";
        rhs->accept(this);
        out_ << ')';
    }
}

void SimdExprEmitter::visit(BlockExpression* block) {
    for (auto& stmt: block->statements()) {
        if (!stmt->is_local_declaration()) {
            stmt->accept(this);
            if (!stmt->is_if() && !stmt->is_block()) {
                out_ << ";\n";
            }
        }
    }
}

void SimdExprEmitter::visit(CallExpression* e) {
    if(is_indirect_)
        out_ << e->name() << "(index_";
    else
        out_ << e->name() << "(i_";

    if (processing_true_ && !current_mask_.empty()) {
        out_ << ", " << current_mask_;
    } else if (!processing_true_ && !current_mask_bar_.empty()) {
        out_ << ", " << current_mask_bar_;
    }
    for (auto& arg: e->args()) {
        out_ << ", ";
        arg->accept(this);
    }
    out_ << ")";
}

void SimdExprEmitter::visit(AssignmentExpression* e) {
    ENTER(out_);
    if (!e->lhs() || !e->lhs()->is_identifier() || !e->lhs()->is_identifier()->symbol()) {
        throw compiler_exception("Expect symbol on lhs of assignment: "+e->to_string());
    }

    auto mask = processing_true_ ? current_mask_ : current_mask_bar_;
    Symbol* lhs = e->lhs()->is_identifier()->symbol();

    auto lhs_pfxd = id_prefix(e->lhs()->is_identifier());

    if (lhs->is_variable() && lhs->is_variable()->is_range()) {
        if (!input_mask_.empty()) {
            mask = "S::logical_and(" + mask + ", " + input_mask_ + ")";
        }
        if(is_indirect_)
            out_ << "indirect(" << lhs_pfxd << "+index_, simd_width_) = ";
        else
            out_ << "indirect(" << lhs_pfxd << "+i_, simd_width_) = ";

        out_ << "S::where(" << mask << ", ";

        bool cast = e->rhs()->is_number();
        if (cast) out_ << "simd_cast<simd_value>(";
        e->rhs()->accept(this);

        out_ << ")";

        if (cast) out_ << ")";
    } else {
        out_ << "S::where(" << mask << ", ";
        e->lhs()->accept(this);
        out_ << ") = ";
        e->rhs()->accept(this);
    }
    EXIT(out_);
}

void SimdExprEmitter::visit(IfExpression* e) {

    // Save old masks
    auto old_mask     = current_mask_;
    auto old_mask_bar = current_mask_bar_;
    auto old_branch   = processing_true_;

    // Create new mask name
    auto new_mask = make_unique_var(e->scope(), "mask_");

    // Set new masks
    out_ << "simd_mask " << new_mask << " = ";
    e->condition()->accept(this);
    out_ << ";\n";

    if (!current_mask_.empty()) {
        auto base_mask = processing_true_ ? current_mask_ : current_mask_bar_;
        current_mask_bar_ = "S::logical_and(" + base_mask + ", S::logical_not(" + new_mask + "))";
        current_mask_     = "S::logical_and(" + base_mask + ", " + new_mask + ")";

    } else {
        current_mask_bar_ = "S::logical_not(" + new_mask + ")";
        current_mask_ = new_mask;
    }

    processing_true_ = true;
    e->true_branch()->accept(this);

    processing_true_ = false;
    if (auto fb = e->false_branch()) {
        fb->accept(this);
    }

    // Reset old masks
    current_mask_     = old_mask;
    current_mask_bar_ = old_mask_bar;
    processing_true_  = old_branch;

}
