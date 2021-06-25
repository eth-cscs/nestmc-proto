#pragma once

#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include "error.hpp"
#include "identifier.hpp"
#include "scope.hpp"
#include "token.hpp"

#include "io/pprintf.hpp"

class Visitor;

class Expression;
class CallExpression;
class BlockExpression;
class IfExpression;
class LocalDeclaration;
class ArgumentExpression;
class FunctionExpression;
class DerivativeExpression;
class PrototypeExpression;
class ProcedureExpression;
class IdentifierExpression;
class NumberExpression;
class IntegerExpression;
class BinaryExpression;
class UnaryExpression;
class AssignmentExpression;
class ConserveExpression;
class LinearExpression;
class ReactionExpression;
class StoichExpression;
class StoichTermExpression;
class CompartmentExpression;
class ConditionalExpression;
class InitialBlock;
class SolveExpression;
class Symbol;
class ConductanceExpression;
class PDiffExpression;
class VariableExpression;
class NetReceiveExpression;
class PostEventExpression;
class APIMethod;
class IndexedVariable;
class LocalVariable;

using expression_ptr = std::unique_ptr<Expression>;
using symbol_ptr = std::unique_ptr<Symbol>;
using scope_type = Scope<Symbol>;
using scope_ptr = std::shared_ptr<scope_type>;
using expr_list_type = std::list<expression_ptr>;

template <typename T, typename... Args>
expression_ptr make_expression(Args&&... args) {
    return expression_ptr(new T(std::forward<Args>(args)...));
}

template <typename T, typename... Args>
symbol_ptr make_symbol(Args&&... args) {
    return symbol_ptr(new T(std::forward<Args>(args)...));
}

// helper functions for generating unary and binary expressions
expression_ptr unary_expression(Location, tok, expression_ptr&&);
expression_ptr unary_expression(tok, expression_ptr&&);
expression_ptr binary_expression(Location, tok, expression_ptr&&, expression_ptr&&);
expression_ptr binary_expression(tok, expression_ptr&&, expression_ptr&&);

/// specifies special properties of a ProcedureExpression
enum class procedureKind {
    normal,      ///< PROCEDURE
    api,         ///< API PROCEDURE
    initial,     ///< INITIAL
    net_receive, ///< NET_RECEIVE
    post_event,  ///< POST_EVENT
    breakpoint,  ///< BREAKPOINT
    kinetic,     ///< KINETIC
    derivative,  ///< DERIVATIVE
    linear,      ///< LINEAR
};
std::string to_string(procedureKind k);

/// classification of different symbol kinds
enum class symbolKind {
    function,         ///< function call
    procedure,        ///< procedure call
    variable,         ///< variable at module scope
    indexed_variable, ///< a variable that is indexed
    local_variable,   ///< variable at local scope
};
std::string to_string(symbolKind k);

/// methods for time stepping state
enum class solverMethod {
    cnexp, // for diagonal linear ODE systems.
    sparse, // for non-diagonal linear ODE systems.
    none
};

enum class solverVariant {
    regular,
    steadystate
};

static std::string to_string(solverMethod m) {
    switch(m) {
        case solverMethod::cnexp:  return std::string("cnexp");
        case solverMethod::sparse: return std::string("sparse");
        case solverMethod::none:   return std::string("none");
    }
    return std::string("<error : undefined solverMethod>");
}

class Expression {
public:
    explicit Expression(Location location)
    :   location_(location)
    {}

    virtual ~Expression() {};

    // This printer should be implemented with a visitor pattern
    // expressions must provide a method for stringification
    virtual std::string to_string() const = 0;

    Location const& location() const { return location_; }

    scope_ptr scope() { return scope_; }

    // set scope explicitly
    void scope(scope_ptr s) { scope_ = s; }

    void error(std::string const& str) {
        error_        = true;
        error_string_ += str;
    }
    void warning(std::string const& str) {
        warning_        = true;
        warning_string_ += str;
    }
    bool has_error()   { return error_; }
    bool has_warning() { return warning_; }
    std::string const& error_message()   const { return error_string_;   }
    std::string const& warning_message() const { return warning_string_; }

    // perform semantic analysis
    virtual void semantic(scope_ptr);
    virtual void semantic(scope_type::symbol_map&) {
        throw compiler_exception("unable to perform semantic analysis for " + this->to_string(), location_);
    };

    // clone an expression
    virtual expression_ptr clone() const;

    // easy lookup of properties
    virtual CallExpression*        is_call()              {return nullptr;}
    virtual CallExpression*        is_function_call()     {return nullptr;}
    virtual CallExpression*        is_procedure_call()    {return nullptr;}
    virtual BlockExpression*       is_block()             {return nullptr;}
    virtual IfExpression*          is_if()                {return nullptr;}
    virtual LocalDeclaration*      is_local_declaration() {return nullptr;}
    virtual ArgumentExpression*    is_argument()          {return nullptr;}
    virtual FunctionExpression*    is_function()          {return nullptr;}
    virtual DerivativeExpression*  is_derivative()        {return nullptr;}
    virtual PrototypeExpression*   is_prototype()         {return nullptr;}
    virtual IdentifierExpression*  is_identifier()        {return nullptr;}
    virtual NumberExpression*      is_number()            {return nullptr;}
    virtual IntegerExpression*     is_integer()           {return nullptr;}
    virtual BinaryExpression*      is_binary()            {return nullptr;}
    virtual UnaryExpression*       is_unary()             {return nullptr;}
    virtual AssignmentExpression*  is_assignment()        {return nullptr;}
    virtual ConserveExpression*    is_conserve()          {return nullptr;}
    virtual LinearExpression*      is_linear()            {return nullptr;}
    virtual ReactionExpression*    is_reaction()          {return nullptr;}
    virtual StoichExpression*      is_stoich()            {return nullptr;}
    virtual StoichTermExpression*  is_stoich_term()       {return nullptr;}
    virtual ConditionalExpression* is_conditional()       {return nullptr;}
    virtual CompartmentExpression* is_compartment()       {return nullptr;}
    virtual InitialBlock*          is_initial_block()     {return nullptr;}
    virtual SolveExpression*       is_solve_statement()   {return nullptr;}
    virtual Symbol*                is_symbol()            {return nullptr;}
    virtual ConductanceExpression* is_conductance_statement() {return nullptr;}
    virtual PDiffExpression*       is_pdiff()             {return nullptr;}

    virtual bool is_lvalue() const {return false;}
    virtual bool is_global_lvalue() const {return false;}

    // force all derived classes to implement visitor
    // this might be a bad idea
    virtual void accept(Visitor *v) = 0;

protected:
    // these are used to flag errors when performing semantic analysis
    // we might want to extend these to an additional "contaminated" flag
    // which marks whether an error was found in a subnode of a node.
    bool error_=false;
    bool warning_=false;
    std::string error_string_;
    std::string warning_string_;

    Location location_;
    scope_ptr scope_;
};

class Symbol : public Expression {
public :
    Symbol(Location loc, std::string name, symbolKind kind)
    :   Expression(std::move(loc)),
        name_(std::move(name)),
        kind_(kind)
    {}

    std::string const& name() const {
        return name_;
    }

    symbolKind kind() const {
        return kind_;
    }

    Symbol* is_symbol() override {
        return this;
    }

    std::string to_string() const override;
    void accept(Visitor *v) override;

    virtual VariableExpression*   is_variable()          {return nullptr;}
    virtual ProcedureExpression*  is_procedure()         {return nullptr;}
    virtual NetReceiveExpression* is_net_receive()       {return nullptr;}
    virtual PostEventExpression*  is_post_event()        {return nullptr;}
    virtual APIMethod*            is_api_method()        {return nullptr;}
    virtual IndexedVariable*      is_indexed_variable()  {return nullptr;}
    virtual LocalVariable*        is_local_variable()    {return nullptr;}

private :
    std::string name_;

    symbolKind kind_;
};

enum class localVariableKind {
    local, argument
};

// an identifier
class IdentifierExpression : public Expression {
public:
    IdentifierExpression(Location loc, std::string const& spelling)
    :   Expression(loc), spelling_(spelling)
    {}

    IdentifierExpression(IdentifierExpression const& other)
    :   Expression(other.location()), spelling_(other.spelling())
    {}

    IdentifierExpression(IdentifierExpression const* other)
    :   Expression(other->location()), spelling_(other->spelling())
    {}

    std::string const& spelling() const {
        return spelling_;
    }

    std::string to_string() const override {
        return yellow(pprintf("%", spelling_));
    }

    expression_ptr clone() const override;

    void semantic(scope_ptr scp) override;

    void symbol(Symbol* sym) { symbol_ = sym; }
    Symbol* symbol() { return symbol_; };

    void accept(Visitor *v) override;

    IdentifierExpression* is_identifier() override {return this;}

    bool is_lvalue() const override;
    bool is_global_lvalue() const override;

    ~IdentifierExpression() {}

    std::string const& name() const {
        if(symbol_) return symbol_->name();
        throw compiler_exception(
            " attempt to look up name of identifier for which no symbol_ yet defined",
            location_);
    }

protected:
    Symbol* symbol_ = nullptr;

    // there has to be some pointer to a table of identifiers
    std::string spelling_;
};

// an identifier for a derivative
class DerivativeExpression : public IdentifierExpression {
public:
    DerivativeExpression(Location loc, std::string const& name)
    :   IdentifierExpression(loc, name)
    {}

    std::string to_string() const override {
        return blue("diff") + "(" + yellow(spelling()) + ")";
    }

    expression_ptr clone() const override;

    void semantic(scope_ptr scp) override;

    DerivativeExpression* is_derivative() override { return this; }

    ~DerivativeExpression() {}

    void accept(Visitor *v) override;
};

// a number
class NumberExpression : public Expression {
public:
    NumberExpression(Location loc, std::string const& value)
        : Expression(loc), value_(std::stold(value))
    {}

    NumberExpression(Location loc, double value)
        : Expression(loc), value_(value)
    {}

    virtual double value() const {return value_;};

    std::string to_string() const override {
        return purple(pprintf("%", value_));
    }

    // do nothing for number semantic analysis
    void semantic(scope_ptr scp) override {};
    expression_ptr clone() const override;

    NumberExpression* is_number() override {return this;}

    ~NumberExpression() {}

    void accept(Visitor *v) override;
private:
    double value_;
};

// an integral number
class IntegerExpression : public NumberExpression {
public:
    IntegerExpression(Location loc, std::string const& value)
        : NumberExpression(loc, value), integer_(std::stoll(value))
    {}

    IntegerExpression(Location loc, long long integer)
        : NumberExpression(loc, static_cast<double>(integer)), integer_(integer)
    {}

    long long integer_value() const {return integer_;}

    std::string to_string() const override {
        return purple(pprintf("%", integer_));
    }

    // do nothing for number semantic analysis
    void semantic(scope_ptr scp) override {};
    expression_ptr clone() const override;

    IntegerExpression* is_integer() override {return this;}

    ~IntegerExpression() {}

    void accept(Visitor *v) override;
private:
    long long integer_;
};



// declaration of a LOCAL variable
class LocalDeclaration : public Expression {
public:
    LocalDeclaration(Location loc)
    :   Expression(loc)
    {}
    LocalDeclaration(Location loc, std::string const& name)
    :   Expression(loc)
    {
        Token tok(tok::identifier, name, loc);
        add_variable(tok);
    }

    std::string to_string() const override;

    bool add_variable(Token name);
    LocalDeclaration* is_local_declaration() override {return this;}
    void semantic(scope_ptr scp) override;
    std::vector<Symbol*>& symbols() {return symbols_;}
    std::map<std::string, Token>& variables() {return vars_;}
    expression_ptr clone() const override;
    ~LocalDeclaration() {}
    void accept(Visitor *v) override;
private:
    std::vector<Symbol*> symbols_;
    // there has to be some pointer to a table of identifiers
    std::map<std::string, Token> vars_;
};

// declaration of an argument
class ArgumentExpression : public Expression {
public:
    ArgumentExpression(Location loc, Token const& tok)
    :   Expression(loc),
        token_(tok),
        name_(tok.spelling)
    {}

    std::string to_string() const override;

    bool add_variable(Token name);
    ArgumentExpression* is_argument() override {return this;}
    void semantic(scope_ptr scp) override;
    Token   token()  {return token_;}
    std::string const& name()  {return name_;}
    void set_name(std::string const& n) {
        name_ = n;
    }
    const std::string& spelling() const {
        return token_.spelling;
    }
    expression_ptr clone() const override;
    ~ArgumentExpression() {}
    void accept(Visitor *v) override;
private:
    Token token_;
    std::string name_;
};

// variable definition
class VariableExpression : public Symbol {
public:
    VariableExpression(Location loc, std::string name)
    :   Symbol(loc, std::move(name), symbolKind::variable)
    {}

    std::string to_string() const override;

    void access(accessKind a) {
        access_ = a;
    }
    void visibility(visibilityKind v) {
        visibility_ = v;
    }
    void linkage(linkageKind l) {
        linkage_ = l;
    }
    void range(rangeKind r) {
        range_kind_ = r;
    }
    void ion_channel(std::string i) {
        ion_channel_ = std::move(i);
    }
    void state(bool s) {
        is_state_ = s;
    }
    void shadows(Symbol* s) {
        shadows_ = s;
    }

    accessKind access() const {
        return access_;
    }
    visibilityKind visibility() const {
        return visibility_;
    }
    linkageKind linkage() const {
        return linkage_;
    }
    const std::string& ion_channel() const {
        return ion_channel_;
    }

    Symbol* shadows() const {
        return shadows_;
    }

    bool is_ion()       const {return !ion_channel_.empty();}
    bool is_state()     const {return is_state_;}
    bool is_range()     const {return range_kind_  == rangeKind::range;}
    bool is_scalar()    const {return !is_range();}

    bool is_readable()  const {return    access_==accessKind::read
                                      || access_==accessKind::readwrite;}

    bool is_writeable() const {return    access_==accessKind::write
                                      || access_==accessKind::readwrite;}

    double value()       const {return value_;}
    void value(double v) {value_ = v;}

    void accept(Visitor *v) override;
    VariableExpression* is_variable() override {return this;}

    ~VariableExpression() {}
protected:

    bool           is_state_    = false;
    accessKind     access_      = accessKind::readwrite;
    visibilityKind visibility_  = visibilityKind::local;
    linkageKind    linkage_     = linkageKind::external;
    rangeKind      range_kind_  = rangeKind::range;
    std::string    ion_channel_ = "";
    double         value_       = std::numeric_limits<double>::quiet_NaN();
    Symbol*        shadows_     = nullptr;
};

// Indexed variables refer to data held in the shared simulation state.
// Printers will rewrite reads from or assignments from indexed variables
// according to its data source and ion channel.

class IndexedVariable : public Symbol {
public:
    IndexedVariable(Location loc,
                    std::string lookup_name,
                    sourceKind data_source,
                    accessKind acc,
                    std::string channel="")
    :   Symbol(std::move(loc), std::move(lookup_name), symbolKind::indexed_variable),
        access_(acc),
        ion_channel_(std::move(channel)),
        data_source_(data_source)
    {
        // external symbols are either read or write only
        if(access()==accessKind::readwrite) {
            throw compiler_exception(
                pprintf("attempt to generate an index % with readwrite access", yellow(lookup_name)),
                location_);
        }
    }

    std::string to_string() const override;

    accessKind access() const { return access_; }
    std::string ion_channel() const { return ion_channel_; }
    sourceKind data_source() const { return data_source_; }
    void data_source(sourceKind k) { data_source_ = k; }

    bool is_ion()   const { return !ion_channel_.empty(); }
    bool is_read()  const { return access_ == accessKind::read;   }
    bool is_write() const { return access_ == accessKind::write;  }

    void accept(Visitor *v) override;
    IndexedVariable* is_indexed_variable() override {return this;}

    ~IndexedVariable() {}
protected:
    accessKind  access_;
    std::string ion_channel_;
    sourceKind  data_source_;
};

class LocalVariable : public Symbol {
public :
    LocalVariable(Location loc,
                  std::string name,
                  localVariableKind kind=localVariableKind::local)
    :   Symbol(std::move(loc), std::move(name), symbolKind::local_variable),
        kind_(kind)
    {}

    LocalVariable* is_local_variable() override {
        return this;
    }

    localVariableKind kind() const {
        return kind_;
    }

    bool is_indexed() const {
        return external_!=nullptr && external_->data_source()!=sourceKind::no_source;
    }

    std::string ion_channel() const {
        return external_? external_->ion_channel(): "";
    }

    bool is_read() const {
        if(is_indexed()) return external_->is_read();
        return true;
    }

    bool is_write() const {
        if(is_indexed()) return external_->is_write();
        return true;
    }

    bool is_local() const {
        return kind_==localVariableKind::local;
    }

    bool is_arg() const {
        return kind_==localVariableKind::argument;
    }

    IndexedVariable* external_variable() {
        return external_;
    }

    void external_variable(IndexedVariable *i) {
        external_ = i;
    }

    std::string to_string() const override;
    void accept(Visitor *v) override;

private :
    IndexedVariable *external_=nullptr;
    localVariableKind kind_;
};


// a SOLVE statement
class SolveExpression : public Expression {
public:
    SolveExpression(
            Location loc,
            std::string name,
            solverMethod method,
            solverVariant variant)
    :   Expression(loc), name_(std::move(name)), method_(method), variant_(variant), procedure_(nullptr)
    {}

    std::string to_string() const override {
        return blue("solve") + "(" + yellow(name_) + ", "
            + green(::to_string(method_)) + ")";
    }

    std::string const& name() const {
        return name_;
    }

    solverMethod method() const {
        return method_;
    }

    solverVariant variant() const {
        return variant_;
    }

    ProcedureExpression* procedure() const {
        return procedure_;
    }

    void procedure(ProcedureExpression *e) {
        procedure_ = e;
    }

    SolveExpression* is_solve_statement() override {
        return this;
    }

    expression_ptr clone() const override;

    void semantic(scope_ptr scp) override;
    void accept(Visitor *v) override;

    ~SolveExpression() {}
private:
    /// pointer to the variable symbol for the state variable to be solved for
    std::string name_;
    solverMethod method_;
    solverVariant variant_;

    ProcedureExpression* procedure_;
};

// a CONDUCTANCE statement
class ConductanceExpression : public Expression {
public:
    ConductanceExpression(
            Location loc,
            std::string name,
            std::string channel)
    :   Expression(loc), name_(std::move(name)), ion_channel_(std::move(channel))
    {}

    std::string to_string() const override {
        return blue("conductance") + "(" + yellow(name_) + ", "
            + green(ion_channel_.empty()? "none": ion_channel_) + ")";
    }

    std::string const& name() const {
        return name_;
    }

    std::string const& ion_channel() const {
        return ion_channel_;
    }

    ConductanceExpression* is_conductance_statement() override {
        return this;
    }

    expression_ptr clone() const override;

    void semantic(scope_ptr scp) override;
    void accept(Visitor *v) override;

    ~ConductanceExpression() {}
private:
    /// pointer to the variable symbol for the state variable to be solved for
    std::string name_;
    std::string ion_channel_;
};

////////////////////////////////////////////////////////////////////////////////
// recursive if statement
// requires a BlockExpression that is a simple wrapper around a std::list
// of Expressions surrounded by {}
////////////////////////////////////////////////////////////////////////////////

class BlockExpression : public Expression {
protected:
    expr_list_type statements_;
    bool is_nested_ = false;

public:
    BlockExpression(
        Location loc,
        expr_list_type&& statements,
        bool is_nested)
    :   Expression(loc),
        statements_(std::move(statements)),
        is_nested_(is_nested)
    {}

    BlockExpression* is_block() override {
        return this;
    }

    expr_list_type& statements() {
        return statements_;
    }

    expression_ptr clone() const override;

    // provide iterators for easy iteration over statements
    auto begin() -> decltype(statements_.begin()) {
        return statements_.begin();
    }
    auto end() -> decltype(statements_.end()) {
        return statements_.end();
    }
    auto back() -> decltype(statements_.back()) {
        return statements_.back();
    }
    auto front() -> decltype(statements_.front()) {
        return statements_.front();
    }
    bool is_nested() const {
        return is_nested_;
    }

    void semantic(scope_ptr scp) override;
    void accept(Visitor* v) override;

    std::string to_string() const override;
};

class IfExpression : public Expression {
public:
    IfExpression(Location loc, expression_ptr&& con, expression_ptr&& tb, expression_ptr&& fb)
    :   Expression(loc), condition_(std::move(con)), true_branch_(std::move(tb)), false_branch_(std::move(fb))
    {}

    IfExpression* is_if() override {
        return this;
    }
    Expression* condition() {
        return condition_.get();
    }
    Expression* true_branch() {
        return true_branch_.get();
    }
    Expression* false_branch() {
        return false_branch_.get();
    }

    expression_ptr clone() const override;

    std::string to_string() const override;

    void replace_condition(expression_ptr&& other);
    void semantic(scope_ptr scp) override;

    void accept(Visitor* v) override;
private:
    expression_ptr condition_;
    expression_ptr true_branch_;
    expression_ptr false_branch_;
};

// a proceduce prototype
class PrototypeExpression : public Expression {
public:
    PrototypeExpression(
            Location loc,
            std::string const& name,
            std::vector<expression_ptr>&& args)
    :   Expression(loc), name_(name), args_(std::move(args))
    {}

    std::string const& name() const {return name_;}

    std::vector<expression_ptr>&      args()       {return args_;}
    std::vector<expression_ptr>const& args() const {return args_;}
    PrototypeExpression* is_prototype() override {return this;}

    // TODO: printing out the vector of unique pointers is an unsolved problem...
    std::string to_string() const override {
        return name_; //+ pprintf("(% args : %)", args_.size(), args_);
    }

    ~PrototypeExpression() {}

    void accept(Visitor *v) override;
private:
    std::string name_;
    std::vector<expression_ptr> args_;
};

class ReactionExpression : public Expression {
public:
    ReactionExpression(Location loc,
                       expression_ptr&& lhs_terms,
                       expression_ptr&& rhs_terms,
                       expression_ptr&& fwd_rate_expr,
                       expression_ptr&& rev_rate_expr)
    : Expression(loc),
      lhs_(std::move(lhs_terms)), rhs_(std::move(rhs_terms)),
      fwd_rate_(std::move(fwd_rate_expr)), rev_rate_(std::move(rev_rate_expr))
    {}

    ReactionExpression* is_reaction() override {return this;}

    std::string to_string() const override;
    void semantic(scope_ptr scp) override;
    expression_ptr clone() const override;
    void accept(Visitor *v) override;

    expression_ptr& lhs() { return lhs_; }
    const expression_ptr& lhs() const { return lhs_; }

    expression_ptr& rhs() { return rhs_; }
    const expression_ptr& rhs() const { return rhs_; }

    expression_ptr& fwd_rate() { return fwd_rate_; }
    const expression_ptr& fwd_rate() const { return fwd_rate_; }

    expression_ptr& rev_rate() { return rev_rate_; }
    const expression_ptr& rev_rate() const { return rev_rate_; }

private:
    expression_ptr lhs_;
    expression_ptr rhs_;
    expression_ptr fwd_rate_;
    expression_ptr rev_rate_;
};

class CompartmentExpression : public Expression {
public:
    CompartmentExpression(Location loc,
                          expression_ptr&& scale_factor,
                          std::vector<expression_ptr>&& state_vars)
    : Expression(loc), scale_factor_(std::move(scale_factor)), state_vars_(std::move(state_vars)) {}

    CompartmentExpression* is_compartment() override {return this;}

    std::string to_string() const override;
    void semantic(scope_ptr scp) override;
    expression_ptr clone() const override;
    void accept(Visitor *v) override;

    expression_ptr& scale_factor() { return scale_factor_; }
    const expression_ptr& scale_factor() const { return scale_factor_; }

    std::vector<expression_ptr>& state_vars() { return state_vars_; }
    const std::vector<expression_ptr>& state_vars() const { return state_vars_; }

    ~CompartmentExpression() {}

private:
    expression_ptr scale_factor_;
    std::vector<expression_ptr> state_vars_;
};

class StoichTermExpression : public Expression {
public:
    StoichTermExpression(Location loc,
                         expression_ptr&& coeff,
                         expression_ptr&& ident)
    : Expression(loc),
      coeff_(std::move(coeff)), ident_(std::move(ident))
    {}

    StoichTermExpression* is_stoich_term() override {return this;}

    std::string to_string() const override {
        return pprintf("% %", coeff()->to_string(), ident()->to_string());
    }
    void semantic(scope_ptr scp) override;
    expression_ptr clone() const override;
    void accept(Visitor *v) override;

    expression_ptr& coeff() { return coeff_; }
    const expression_ptr& coeff() const { return coeff_; }

    expression_ptr& ident() { return ident_; }
    const expression_ptr& ident() const { return ident_; }

    bool negative() const {
        auto iexpr = coeff_->is_integer();
        return iexpr && iexpr->integer_value()<0;
    }

private:
    expression_ptr coeff_;
    expression_ptr ident_;
};

class StoichExpression : public Expression {
public:
    StoichExpression(Location loc, std::vector<expression_ptr>&& terms)
    : Expression(loc), terms_(std::move(terms))
    {}

    StoichExpression(Location loc)
    : Expression(loc)
    {}

    StoichExpression* is_stoich() override {return this;}

    std::string to_string() const override;
    void semantic(scope_ptr scp) override;
    expression_ptr clone() const override;
    void accept(Visitor *v) override;

    std::vector<expression_ptr>& terms() { return terms_; }
    const std::vector<expression_ptr>& terms() const { return terms_; }

private:
    std::vector<expression_ptr> terms_;
};

// marks a call site in the AST
// is used to mark both function and procedure calls
class CallExpression : public Expression {
public:
    CallExpression(Location loc, std::string spelling, std::vector<expression_ptr>&& args)
    :   Expression(loc), spelling_(std::move(spelling)), args_(std::move(args))
    {}

    std::vector<expression_ptr>& args() { return args_; }
    std::vector<expression_ptr> const& args() const { return args_; }

    std::string& name() { return spelling_; }
    std::string const& name() const { return spelling_; }

    void semantic(scope_ptr scp) override;
    expression_ptr clone() const override;

    std::string to_string() const override;

    void accept(Visitor *v) override;

    CallExpression* is_call() override {
        return this;
    }
    CallExpression* is_function_call()  override {
        return symbol_->kind() == symbolKind::function ? this : nullptr;
    }
    CallExpression* is_procedure_call() override {
        return symbol_->kind() == symbolKind::procedure ? this : nullptr;
    }

    FunctionExpression* function() {
        return symbol_->kind() == symbolKind::function
            ? symbol_->is_function() : nullptr;
    }

    ProcedureExpression* procedure() {
        return symbol_->kind() == symbolKind::procedure
            ? symbol_->is_procedure() : nullptr;
    }

private:
    Symbol* symbol_;

    std::string spelling_;
    std::vector<expression_ptr> args_;
};

class ProcedureExpression : public Symbol {
public:
    ProcedureExpression( Location loc,
                         std::string name,
                         std::vector<expression_ptr>&& args,
                         expression_ptr&& body,
                         procedureKind k=procedureKind::normal)
    :   Symbol(loc, std::move(name), symbolKind::procedure), args_(std::move(args)), kind_(k)
    {
        if(!body->is_block()) {
            throw compiler_exception(
                " attempt to initialize ProcedureExpression with non-block expression, i.e.\n"
                + body->to_string(),
                location_);
        }
        body_ = std::move(body);
    }

    std::vector<expression_ptr>& args() {
        return args_;
    }
    BlockExpression* body() {
        return body_.get()->is_block();
    }
    void body(expression_ptr&& new_body) {
        if(!new_body->is_block()) {
            Location loc = new_body? new_body->location(): Location{};
            throw compiler_exception(
                " attempt to set ProcedureExpression body with non-block expression, i.e.\n"
                + new_body->to_string(),
                loc);
        }
        body_ = std::move(new_body);
    }

    void semantic(scope_ptr scp) override;
    void semantic(scope_type::symbol_map &scp) override;
    ProcedureExpression* is_procedure() override {return this;}
    std::string to_string() const override;
    void accept(Visitor *v) override;

    /// can be used to determine whether the procedure has been lowered
    /// from a special block, e.g. BREAKPOINT, INITIAL, NET_RECEIVE, etc
    procedureKind kind() const {return kind_;}

protected:
    Symbol* symbol_;

    std::vector<expression_ptr> args_;
    expression_ptr body_;
    procedureKind kind_ = procedureKind::normal;
};

class APIMethod : public ProcedureExpression {
public:
    APIMethod( Location loc,
               std::string name,
               std::vector<expression_ptr>&& args,
               expression_ptr&& body)
        :   ProcedureExpression(loc, std::move(name), std::move(args), std::move(body), procedureKind::api)
    {}

    using ProcedureExpression::semantic;
    void semantic(scope_type::symbol_map &scp) override;
    APIMethod* is_api_method() override {return this;}
    void accept(Visitor *v) override;

    std::string to_string() const override;
};

/// stores the INITIAL block in a NET_RECEIVE block, if there is one
/// should not be used anywhere but NET_RECEIVE
class InitialBlock : public BlockExpression {
public:
    InitialBlock(
        Location loc,
        expr_list_type&& statements)
    :   BlockExpression(loc, std::move(statements), true)
    {}

    std::string to_string() const override;

    // currently we use the semantic for a BlockExpression
    // this could be overriden to check for statements
    // specific to initialization of a net_receive block
    //void semantic() override;

    void accept(Visitor *v) override;
    InitialBlock* is_initial_block() override {return this;}
};

/// handle NetReceiveExpressions as a special case of ProcedureExpression
class NetReceiveExpression : public ProcedureExpression {
public:
    NetReceiveExpression( Location loc,
                          std::string name,
                          std::vector<expression_ptr>&& args,
                          expression_ptr&& body)
    :   ProcedureExpression(loc, std::move(name), std::move(args), std::move(body), procedureKind::net_receive)
    {}

    void semantic(scope_type::symbol_map &scp) override;
    NetReceiveExpression* is_net_receive() override {return this;}
    /// hard code the kind
    procedureKind kind() {return procedureKind::net_receive;}

    void accept(Visitor *v) override;
    InitialBlock* initial_block() {return initial_block_;}
protected:
    InitialBlock* initial_block_ = nullptr;
};

/// handle PostEventExpression as a special case of ProcedureExpression
class PostEventExpression : public ProcedureExpression {
public:
    PostEventExpression( Location loc,
                          std::string name,
                          std::vector<expression_ptr>&& args,
                          expression_ptr&& body)
            :   ProcedureExpression(loc, std::move(name), std::move(args), std::move(body), procedureKind::post_event)
    {}

    void semantic(scope_type::symbol_map &scp) override;
    PostEventExpression* is_post_event() override {return this;}
    /// hard code the kind
    procedureKind kind() {return procedureKind::post_event;}

    void accept(Visitor *v) override;
};

class FunctionExpression : public Symbol {
public:
    FunctionExpression( Location loc,
                        std::string name,
                        std::vector<expression_ptr>&& args,
                        expression_ptr&& body)
    :   Symbol(loc, std::move(name), symbolKind::function),
        args_(std::move(args))
    {
        if(!body->is_block()) {
            throw compiler_exception(
                " attempt to initialize FunctionExpression with non-block expression, i.e.\n"
                + body->to_string(),
                location_);
        }
        body_ = std::move(body);
    }

    std::vector<expression_ptr>& args() {
        return args_;
    }
    BlockExpression* body() {
        return body_->is_block();
    }
    void body(expression_ptr&& new_body) {
        if(!new_body->is_block()) {
            Location loc = new_body? new_body->location(): Location{};
            throw compiler_exception(
                    " attempt to set FunctionExpression body with non-block expression, i.e.\n"
                    + new_body->to_string(),
                    loc);
        }
        body_ = std::move(new_body);
    }

    FunctionExpression* is_function() override {return this;}
    void semantic(scope_type::symbol_map&) override;
    std::string to_string() const override;
    void accept(Visitor *v) override;

private:
    Symbol* symbol_;

    std::vector<expression_ptr> args_;
    expression_ptr body_;
};

////////////////////////////////////////////////////////////
// unary expressions
////////////////////////////////////////////////////////////

/// Unary expression
class UnaryExpression : public Expression {
protected:
    expression_ptr expression_;
    tok op_;
public:
    UnaryExpression(Location loc, tok op, expression_ptr&& e)
    :   Expression(loc),
        expression_(std::move(e)),
        op_(op)
    {}

    std::string to_string() const override {
        return pprintf("(% %)", green(token_string(op_)), expression_->to_string());
    }

    expression_ptr clone() const override;

    tok op() const {return op_;}
    UnaryExpression* is_unary() override {return this;};
    Expression* expression() {return expression_.get();}
    const Expression* expression() const {return expression_.get();}
    void semantic(scope_ptr scp) override;
    void accept(Visitor *v) override;
    void replace_expression(expression_ptr&& other);
};

/// negation unary expression, i.e. -x
class NegUnaryExpression : public UnaryExpression {
public:
    NegUnaryExpression(Location loc, expression_ptr e)
    :   UnaryExpression(loc, tok::minus, std::move(e))
    {}

    void accept(Visitor *v) override;
};

/// exponential unary expression, i.e. e^x or exp(x)
class ExpUnaryExpression : public UnaryExpression {
public:
    ExpUnaryExpression(Location loc, expression_ptr e)
    :   UnaryExpression(loc, tok::exp, std::move(e))
    {}

    void accept(Visitor *v) override;
};

// logarithm unary expression, i.e. log_10(x)
class LogUnaryExpression : public UnaryExpression {
public:
    LogUnaryExpression(Location loc, expression_ptr e)
    :   UnaryExpression(loc, tok::log, std::move(e))
    {}

    void accept(Visitor *v) override;
};

// absolute value unary expression, i.e. abs(x)
class AbsUnaryExpression : public UnaryExpression {
public:
    AbsUnaryExpression(Location loc, expression_ptr e)
    :   UnaryExpression(loc, tok::abs, std::move(e))
    {}

    void accept(Visitor *v) override;
};

class SafeInvUnaryExpression : public UnaryExpression {
public:
    SafeInvUnaryExpression(Location loc, expression_ptr e)
    :   UnaryExpression(loc, tok::safeinv, std::move(e))
    {}

    void accept(Visitor *v) override;
};

// exprel reciprocal unary expression,
// i.e. x/(exp(x)-1)=x/expm1(x) with exprelr(0)=1
class ExprelrUnaryExpression : public UnaryExpression {
public:
    ExprelrUnaryExpression(Location loc, expression_ptr e)
    :   UnaryExpression(loc, tok::exprelr, std::move(e))
    {}

    void accept(Visitor *v) override;
};

// cosine unary expression, i.e. cos(x)
class CosUnaryExpression : public UnaryExpression {
public:
    CosUnaryExpression(Location loc, expression_ptr e)
    :   UnaryExpression(loc, tok::cos, std::move(e))
    {}

    void accept(Visitor *v) override;
};

// sin unary expression, i.e. sin(x)
class SinUnaryExpression : public UnaryExpression {
public:
    SinUnaryExpression(Location loc, expression_ptr e)
    :   UnaryExpression(loc, tok::sin, std::move(e))
    {}

    void accept(Visitor *v) override;
};

////////////////////////////////////////////////////////////
// binary expressions

////////////////////////////////////////////////////////////

/// binary expression base class
/// never used directly in the AST, instead the specializations that derive from
/// it are inserted into the AST.
class BinaryExpression : public Expression {
protected:
    expression_ptr lhs_;
    expression_ptr rhs_;
    tok op_;
public:
    BinaryExpression(Location loc, tok op, expression_ptr&& lhs, expression_ptr&& rhs)
    :   Expression(loc),
        lhs_(std::move(lhs)),
        rhs_(std::move(rhs)),
        op_(op)
    {}

    tok op() const {return op_;}
    virtual bool is_infix() const {return true;}
    Expression* lhs() {return lhs_.get();}
    Expression* rhs() {return rhs_.get();}
    const Expression* lhs() const {return lhs_.get();}
    const Expression* rhs() const {return rhs_.get();}
    BinaryExpression* is_binary() override {return this;}
    void semantic(scope_ptr scp) override;
    expression_ptr clone() const override;
    void replace_rhs(expression_ptr&& other);
    void replace_lhs(expression_ptr&& other);
    std::string to_string() const override;
    void accept(Visitor *v) override;
};

class AssignmentExpression : public BinaryExpression {
public:
    AssignmentExpression(Location loc, expression_ptr&& lhs, expression_ptr&& rhs)
    :   BinaryExpression(loc, tok::eq, std::move(lhs), std::move(rhs))
    {}

    AssignmentExpression* is_assignment() override {return this;}

    void semantic(scope_ptr scp) override;

    void accept(Visitor *v) override;
};

class ConserveExpression : public BinaryExpression {
public:
    ConserveExpression(Location loc, expression_ptr&& lhs, expression_ptr&& rhs)
    :   BinaryExpression(loc, tok::eq, std::move(lhs), std::move(rhs))
    {}

    ConserveExpression* is_conserve() override {return this;}
    expression_ptr clone() const override;

    void semantic(scope_ptr scp) override;

    void accept(Visitor *v) override;
};

class LinearExpression : public BinaryExpression {
public:
    LinearExpression(Location loc, expression_ptr&& lhs, expression_ptr&& rhs)
            :   BinaryExpression(loc, tok::eq, std::move(lhs), std::move(rhs))
    {}

    LinearExpression* is_linear() override {return this;}
    expression_ptr clone() const override;

    void semantic(scope_ptr scp) override;

    void accept(Visitor *v) override;
};

class AddBinaryExpression : public BinaryExpression {
public:
    AddBinaryExpression(Location loc, expression_ptr&& lhs, expression_ptr&& rhs)
    :   BinaryExpression(loc, tok::plus, std::move(lhs), std::move(rhs))
    {}

    void accept(Visitor *v) override;
};

class SubBinaryExpression : public BinaryExpression {
public:
    SubBinaryExpression(Location loc, expression_ptr&& lhs, expression_ptr&& rhs)
    :   BinaryExpression(loc, tok::minus, std::move(lhs), std::move(rhs))
    {}

    void accept(Visitor *v) override;
};

class MulBinaryExpression : public BinaryExpression {
public:
    MulBinaryExpression(Location loc, expression_ptr&& lhs, expression_ptr&& rhs)
    :   BinaryExpression(loc, tok::times, std::move(lhs), std::move(rhs))
    {}

    void accept(Visitor *v) override;
};

class DivBinaryExpression : public BinaryExpression {
public:
    DivBinaryExpression(Location loc, expression_ptr&& lhs, expression_ptr&& rhs)
    :   BinaryExpression(loc, tok::divide, std::move(lhs), std::move(rhs))
    {}

    void accept(Visitor *v) override;
};

class MinBinaryExpression : public BinaryExpression {
public:
    MinBinaryExpression(Location loc, expression_ptr&& lhs, expression_ptr&& rhs)
    :   BinaryExpression(loc, tok::min, std::move(lhs), std::move(rhs))
    {}

    // min is a prefix binop
    bool is_infix() const override {return false;}

    void accept(Visitor *v) override;
};

class MaxBinaryExpression : public BinaryExpression {
public:
    MaxBinaryExpression(Location loc, expression_ptr&& lhs, expression_ptr&& rhs)
    :   BinaryExpression(loc, tok::max, std::move(lhs), std::move(rhs))
    {}

    // max is a prefix binop
    bool is_infix() const override {return false;}

    void accept(Visitor *v) override;
};

class PowBinaryExpression : public BinaryExpression {
public:
    PowBinaryExpression(Location loc, expression_ptr&& lhs, expression_ptr&& rhs)
    :   BinaryExpression(loc, tok::pow, std::move(lhs), std::move(rhs))
    {}

    void accept(Visitor *v) override;
};

class ConditionalExpression : public BinaryExpression {
public:
    ConditionalExpression(Location loc, tok op, expression_ptr&& lhs, expression_ptr&& rhs)
    :   BinaryExpression(loc, op, std::move(lhs), std::move(rhs))
    {}

    ConditionalExpression* is_conditional() override {return this;}

    void accept(Visitor *v) override;
};

class PDiffExpression : public Expression {
public:
    PDiffExpression(Location loc, expression_ptr&& var, expression_ptr&& arg)
    :  Expression(loc), var_(std::move(var)), arg_(std::move(arg))
    {}

    std::string to_string() const override;
    void accept(Visitor *v) override;
    void semantic(scope_ptr scp) override;
    expression_ptr clone() const override;

    PDiffExpression* is_pdiff()  override { return this; }

    Expression* var() { return var_.get(); }
    Expression* arg() { return arg_.get(); }

private:
    expression_ptr var_;
    expression_ptr arg_;
};

