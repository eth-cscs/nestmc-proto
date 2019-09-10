#include <iostream>

#include "astmanip.hpp"
#include "error.hpp"
#include "functioninliner.hpp"
#include "errorvisitor.hpp"

expression_ptr inline_function_call(const expression_ptr& e)
{
    auto assign_to_func = e->is_assignment();
    auto ret_name = assign_to_func->lhs()->is_identifier()->clone();

    if(auto f = assign_to_func->rhs()->is_function_call()) {
        auto body = f->function()->body()->clone();

        FunctionInliner func_inliner(f->name(), ret_name, f->function()->args(), f->args(), e->scope());

        body->accept(&func_inliner);
        return body;
    }
    return {};
}
///////////////////////////////////////////////////////////////////////////////
//  function inliner
///////////////////////////////////////////////////////////////////////////////
void FunctionInliner::replace_with_args(Expression* e) {
    for(auto i=0u; i<fargs_.size(); ++i) {
        if(auto id = cargs_[i]->is_identifier()) {
            VariableReplacer v(fargs_[i], id->spelling());
            e->accept(&v);
        }
        else if(auto value = cargs_[i]->is_number()) {
            ValueInliner v(fargs_[i], value->value());
            e->accept(&v);
        }
        else {
            throw compiler_exception("can't inline functions with expressions as arguments", e->location());
        }
    }
    e->semantic(scope_);

    ErrorVisitor v("");
    e->accept(&v);
    if(v.num_errors()) {
        throw compiler_exception("something went wrong with inlined function call ", e->location());
    }
}

void FunctionInliner::visit(Expression* e) {
    throw compiler_exception(
            "I don't know how to function inlining for this statement : "
            + e->to_string(), e->location());
}

void FunctionInliner::visit(LocalDeclaration* e) {
    auto loc = e->location();

    std::map<std::string, Token> new_vars;
    for (auto& var: e->variables()) {
        auto unique_decl = make_unique_local_decl(scope_, loc, "r_");
        auto unique_name = unique_decl.id->is_identifier()->spelling();

        fargs_.push_back(var.first);
        cargs_.push_back(unique_decl.id->clone());

        auto e_tok = var.second;
        e_tok.spelling = unique_name;
        new_vars[unique_name] =  e_tok;
    }
    e->variables().swap(new_vars);
}

void FunctionInliner::visit(BlockExpression* e) {
    for (auto& expr: e->statements()) {
        expr->accept(this);
    }
}

void FunctionInliner::visit(UnaryExpression* e) {
    replace_with_args(e);
}

void FunctionInliner::visit(BinaryExpression* e) {
    replace_with_args(e);
}

void FunctionInliner::visit(AssignmentExpression* e) {
    std::cout << e->to_string() << std::endl;
    e->rhs()->accept(this);
    if (auto elhs = e->lhs()->is_identifier()) {
        if (elhs->spelling() == func_name_) {
            e->replace_lhs(lhs_->clone());
        } else {
            for (unsigned i = 0;  i < fargs_.size(); i++) {
                if (fargs_[i] == elhs->spelling()) {
                    e->replace_lhs(cargs_[i]->clone());
                    break;
                }
            }
        }
    }
}

void FunctionInliner::visit(IfExpression* e) {
    e->condition()->accept(this);
    e->true_branch()->accept(this);
    e->false_branch()->accept(this);
}

///////////////////////////////////////////////////////////////////////////////
//  variable replacer
///////////////////////////////////////////////////////////////////////////////

void VariableReplacer::visit(Expression *e) {
    throw compiler_exception(
            "I don't know how to variable inlining for this statement : "
            + e->to_string(), e->location());
}

void VariableReplacer::visit(UnaryExpression *e) {
    auto exp = e->expression()->is_identifier();
    if(exp && exp->spelling()==source_) {
        e->replace_expression(
            make_expression<IdentifierExpression>(exp->location(), target_)
        );
    }
    else if(!exp) {
        e->expression()->accept(this);
    }
}

void VariableReplacer::visit(BinaryExpression *e) {
    auto lhs = e->lhs()->is_identifier();
    if(lhs && lhs->spelling()==source_) {
        e->replace_lhs(
            make_expression<IdentifierExpression>(lhs->location(), target_)
        );
    }
    else if(!lhs){ // only inspect subexpressions that are not themselves identifiers
        e->lhs()->accept(this);
    }

    auto rhs = e->rhs()->is_identifier();
    if(rhs && rhs->spelling()==source_) {
        e->replace_rhs(
            make_expression<IdentifierExpression>(rhs->location(), target_)
        );
    }
    else if(!rhs){ // only inspect subexpressions that are not themselves identifiers
        e->rhs()->accept(this);
    }
}

///////////////////////////////////////////////////////////////////////////////
//  value inliner
///////////////////////////////////////////////////////////////////////////////

void ValueInliner::visit(Expression *e) {
    throw compiler_exception(
            "I don't know how to value inlining for this statement : "
            + e->to_string(), e->location());
}

void ValueInliner::visit(UnaryExpression *e) {
    auto exp = e->expression()->is_identifier();
    if(exp && exp->spelling()==source_) {
        e->replace_expression(
            make_expression<NumberExpression>(exp->location(), value_)
        );
    }
    else if(!exp){
        e->expression()->accept(this);
    }
}

void ValueInliner::visit(BinaryExpression *e) {
    auto lhs = e->lhs()->is_identifier();
    if(lhs && lhs->spelling()==source_) {
        e->replace_lhs(
            make_expression<NumberExpression>(lhs->location(), value_)
        );
    }
    else if(!lhs) {
        e->lhs()->accept(this);
    }

    auto rhs = e->rhs()->is_identifier();
    if(rhs && rhs->spelling()==source_) {
        e->replace_rhs(
            make_expression<NumberExpression>(rhs->location(), value_)
        );
    }
    else if(!rhs){
        e->rhs()->accept(this);
    }
}
