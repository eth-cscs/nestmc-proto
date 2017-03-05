#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>

#include "errorvisitor.hpp"
#include "expressionclassifier.hpp"
#include "functionexpander.hpp"
#include "functioninliner.hpp"
#include "kineticrewriter.hpp"
#include "module.hpp"
#include "parser.hpp"
#include "solvers.hpp"
#include "symdiff.hpp"
#include "visitor.hpp"

class NrnCurrentRewriter: public BlockRewriterBase {
    expression_ptr id(const std::string& name, Location loc) {
        return make_expression<IdentifierExpression>(loc, name);
    }

    expression_ptr id(const std::string& name) {
        return id(name, loc_);
    }

    static ionKind is_ion_update(Expression* e) {
        if(auto a = e->is_assignment()) {
            if(auto sym = a->lhs()->is_identifier()->symbol()) {
                if(auto var = sym->is_local_variable()) {
                    return var->ion_channel();
                }
            }
        }
        return ionKind::none;
    }

    moduleKind kind_;
    bool has_current_update_ = false;

public:
    using BlockRewriterBase::visit;

    explicit NrnCurrentRewriter(moduleKind kind): kind_(kind) {}

    virtual void finalize() override {
        if (has_current_update_) {
            // Initialize current_ as first statement.
            statements_.push_front(make_expression<AssignmentExpression>(loc_,
                    id("current_"),
                    make_expression<NumberExpression>(loc_, 0.0)));

            if (kind_==moduleKind::density) {
                statements_.push_back(make_expression<AssignmentExpression>(loc_,
                    id("current_"),
                    make_expression<MulBinaryExpression>(loc_,
                        id("weights_"),
                        id("current_"))));
            }
        }
    }

    virtual void visit(SolveExpression *e) override {}
    virtual void visit(ConductanceExpression *e) override {}
    virtual void visit(AssignmentExpression *e) override {
        statements_.push_back(e->clone());
        auto loc = e->location();

        if (is_ion_update(e)!=ionKind::none) {
            has_current_update_ = true;

            if (!linear_test(e->rhs(), {"v"}).is_linear) {
                error({"current update expressions must be linear in v: "+e->rhs()->to_string(),
                       e->location()});
                return;
            }
            else {
                statements_.push_back(make_expression<AssignmentExpression>(loc,
                    id("current_", loc),
                    make_expression<AddBinaryExpression>(loc,
                        id("current_", loc),
                        e->lhs()->clone())));
            }
        }
    }
};

Module::Module(std::string const& fname)
: fname_(fname)
{
    // open the file at the end
    std::ifstream fid;
    fid.open(fname.c_str(), std::ios::binary | std::ios::ate);
    if(!fid.is_open()) { // return if no file opened
        return;
    }

    // determine size of file
    std::size_t size = fid.tellg();
    fid.seekg(0, std::ios::beg);

    // allocate space for storage and read
    buffer_.resize(size+1);
    fid.read(buffer_.data(), size);
    buffer_[size] = 0; // append \0 to terminate string
}

Module::Module(std::vector<char> const& buffer) {
    buffer_ = buffer;

    // add \0 to end of buffer if not already present
    if (buffer_[buffer_.size()-1] != 0)
        buffer_.push_back(0);
}

Module::Module(const char* buffer, size_t count) {
    auto size = std::distance(buffer, std::find(buffer, buffer+count, '\0'));
    buffer_.reserve(size+1);
    buffer_.insert(buffer_.end(), buffer, buffer+size);
    buffer_.push_back(0);
}

std::vector<Module::symbol_ptr>&
Module::procedures() {
    return procedures_;
}

std::vector<Module::symbol_ptr>const&
Module::procedures() const {
    return procedures_;
}

std::vector<Module::symbol_ptr>&
Module::functions() {
    return functions_;
}

std::vector<Module::symbol_ptr>const&
Module::functions() const {
    return functions_;
}

Module::symbol_map&
Module::symbols() {
    return symbols_;
}

Module::symbol_map const&
Module::symbols() const {
    return symbols_;
}

std::string Module::error_string() const {
    std::string str;
    for (const error_entry& entry: errors()) {
        if (!str.empty()) str += '\n';
        str += red("error   ");
        str += white(pprintf("%:% ", file_name(), entry.location));
        str += entry.message;
    }
    return str;
}

std::string Module::warning_string() const {
    std::string str;
    for (const error_entry& entry: errors()) {
        if (!str.empty()) str += '\n';
        str += purple("error   ");
        str += white(pprintf("%:% ", file_name(), entry.location));
        str += entry.message;
    }
    return str;
}

bool Module::semantic() {
    ////////////////////////////////////////////////////////////////////////////
    // create the symbol table
    // there are three types of symbol to look up
    //  1. variables
    //  2. function calls
    //  3. procedure calls
    // the symbol table is generated, then we can traverse the AST and verify
    // that all symbols are correctly used
    ////////////////////////////////////////////////////////////////////////////

    // first add variables defined in the NEURON, ASSIGNED and PARAMETER
    // blocks these symbols have "global" scope, i.e. they are visible to all
    // functions and procedurs in the mechanism
    add_variables_to_symbols();

    // Helper which iterates over a vector of Symbols, moving them into the
    // symbol table.
    // Returns false if a symbol name clashes with the name of a symbol that
    // is already in the symbol table.
    auto move_symbols = [this] (std::vector<symbol_ptr>& symbol_list) {
        for(auto& symbol: symbol_list) {
            bool is_found = (symbols_.find(symbol->name()) != symbols_.end());
            if(is_found) {
                error(
                    pprintf("'%' clashes with previously defined symbol",
                            symbol->name()),
                    symbol->location()
                );
                return false;
            }
            // move symbol to table
            symbols_[symbol->name()] = std::move(symbol);
        }
        return true;
    };

    // move functions and procedures to the symbol table
    if(!move_symbols(functions_))  return false;
    if(!move_symbols(procedures_)) return false;

    ////////////////////////////////////////////////////////////////////////////
    // now iterate over the functions and procedures and perform semantic
    // analysis on each. This includes
    //  -   variable, function and procedure lookup
    //  -   generate local variable table for each function/procedure
    //  -   inlining function calls
    ////////////////////////////////////////////////////////////////////////////
#ifdef LOGGING
    std::cout << white("===================================\n");
    std::cout << cyan("        Function Inlining\n");
    std::cout << white("===================================\n");
#endif
    int errors = 0;
    for(auto& e : symbols_) {
        auto& s = e.second;

        if(    s->kind() == symbolKind::function
            || s->kind() == symbolKind::procedure)
        {
#ifdef LOGGING
            std::cout << "\nfunction inlining for " << s->location() << "\n"
                      << s->to_string() << "\n"
                      << green("\n-call site lowering-\n\n");
#endif
            // first perform semantic analysis
            s->semantic(symbols_);

            // then use an error visitor to print out all the semantic errors
            ErrorVisitor v(file_name());
            s->accept(&v);
            errors += v.num_errors();

            // inline function calls
            // this requires that the symbol table has already been built
            if(v.num_errors()==0) {
                auto &b = s->kind()==symbolKind::function ?
                    s->is_function()->body()->statements() :
                    s->is_procedure()->body()->statements();

                // lower function call sites so that all function calls are of
                // the form : variable = call(<args>)
                // e.g.
                //      a = 2 + foo(2+x, y, 1)
                // becomes
                //      ll0_ = foo(2+x, y, 1)
                //      a = 2 + ll0_
                for(auto e=b.begin(); e!=b.end(); ++e) {
                    b.splice(e, lower_function_calls((*e).get()));
                }
#ifdef LOGGING
                std::cout << "body after call site lowering\n";
                for(auto& l : b) std::cout << "  " << l->to_string() << " @ " << l->location() << "\n";
                std::cout << green("\n-argument lowering-\n\n");
#endif

                // lower function arguments that are not identifiers or literals
                // e.g.
                //      ll0_ = foo(2+x, y, 1)
                //      a = 2 + ll0_
                // becomes
                //      ll1_ = 2+x
                //      ll0_ = foo(ll1_, y, 1)
                //      a = 2 + ll0_
                for(auto e=b.begin(); e!=b.end(); ++e) {
                    if(auto be = (*e)->is_binary()) {
                        // only apply to assignment expressions where rhs is a
                        // function call because the function call lowering step
                        // above ensures that all function calls are of this form
                        if(auto rhs = be->rhs()->is_function_call()) {
                            b.splice(e, lower_function_arguments(rhs->args()));
                        }
                    }
                }

#ifdef LOGGING
                std::cout << "body after argument lowering\n";
                for(auto& l : b) std::cout << "  " << l->to_string() << " @ " << l->location() << "\n";
                std::cout << green("\n-inlining-\n\n");
#endif

                // Do the inlining, which currently only works for functions
                // that have a single statement in their body
                // e.g. if the function foo in the examples above is defined as follows
                //
                //  function foo(a, b, c) {
                //      foo = a*(b + c)
                //  }
                //
                // the full inlined example is
                //      ll1_ = 2+x
                //      ll0_ = ll1_*(y + 1)
                //      a = 2 + ll0_
                for(auto e=b.begin(); e!=b.end(); ++e) {
                    if(auto ass = (*e)->is_assignment()) {
                        if(ass->rhs()->is_function_call()) {
                            ass->replace_rhs(inline_function_call(ass->rhs()));
                        }
                    }
                }

#ifdef LOGGING
                std::cout << "body after inlining\n";
                for(auto& l : b) std::cout << "  " << l->to_string() << " @ " << l->location() << "\n";
#endif
            }
        }
    }

    if(errors) {
        error("There were "+std::to_string(errors)+" errors in the semantic analysis");
        return false;
    }

    // All API methods are generated from statements in one of the special procedures
    // defined in NMODL, e.g. the nrn_init() API call is based on the INITIAL block.
    // When creating an API method, the first task is to look up the source procedure,
    // i.e. the INITIAL block for nrn_init(). This lambda takes care of this repetative
    // lookup work, with error checking.
    auto make_empty_api_method = [this]
            (std::string const& name, std::string const& source_name)
            -> std::pair<APIMethod*, ProcedureExpression*>
    {
        if( !has_symbol(source_name, symbolKind::procedure) ) {
            error(pprintf("unable to find symbol '%'", yellow(source_name)),
                   Location());
            return std::make_pair(nullptr, nullptr);
        }

        auto source = symbols_[source_name]->is_procedure();
        auto loc = source->location();

        if( symbols_.find(name)!=symbols_.end() ) {
            error(pprintf("'%' clashes with reserved name, please rename it",
                          yellow(name)),
                  symbols_.find(name)->second->location());
            return std::make_pair(nullptr, source);
        }

        symbols_[name] = make_symbol<APIMethod>(
                          loc, name,
                          std::vector<expression_ptr>(), // no arguments
                          make_expression<BlockExpression>
                            (loc, expr_list_type(), false)
                         );

        auto proc = symbols_[name]->is_api_method();
        return std::make_pair(proc, source);
    };

    //.........................................................................
    // nrn_init : based on the INITIAL block (i.e. the 'initial' procedure
    //.........................................................................
    auto initial_api = make_empty_api_method("nrn_init", "initial");
    auto api_init  = initial_api.first;
    auto proc_init = initial_api.second;

    if(api_init)
    {
        auto& body = api_init->body()->statements();

        for(auto& e : *proc_init->body()) {
            body.emplace_back(e->clone());
        }

        api_init->semantic(symbols_);
    }
    else {
        if(!proc_init) {
            error("an INITIAL block is required", Location());
        }
        return false;
    }

    // Look in the symbol table for a procedure with the name "breakpoint".
    // This symbol corresponds to the BREAKPOINT block in the .mod file
    // There are two APIMethods generated from BREAKPOINT.
    // The first is nrn_state, which is the first case handled below.
    // The second is nrn_current, which is handled after this block
    auto state_api  = make_empty_api_method("nrn_state", "breakpoint");
    auto api_state  = state_api.first;
    auto breakpoint = state_api.second;

    api_state->semantic(symbols_);
    scope_ptr nrn_state_scope = api_state->scope();

    if(breakpoint) {
        //..........................................................
        // nrn_state : The temporal integration of state variables
        //..........................................................

        // grab SOLVE statements, put them in `nrn_state` after translation.
        bool found_solve = false;
        bool found_non_solve = false;
        std::set<std::string> solved_ids;

        for(auto& e: (breakpoint->body()->statements())) {
            SolveExpression* solve_expression = e->is_solve_statement();
            if(!solve_expression) {
                found_non_solve = true;
                continue;
            }
            if(found_non_solve) {
                error("SOLVE statements must come first in BREAKPOINT block",
                    e->location());
                return false;
            }

            found_solve = true;
            std::unique_ptr<SolverVisitorBase> solver;

            switch(solve_expression->method()) {
            case solverMethod::cnexp:
                solver = make_unique<CnexpSolverVisitor>();
                break;
            case solverMethod::sparse:
                solver = make_unique<SparseSolverVisitor>();
                break;
            case solverMethod::none:
                solver = make_unique<DirectSolverVisitor>();
                break;
            }

            // If the derivative block is a kinetic block, perform kinetic
            // rewrite first.

            auto deriv = solve_expression->procedure();

            if (deriv->kind()==procedureKind::kinetic) {
                kinetic_rewrite(deriv->body())->accept(solver.get());
            }
            else {
                deriv->body()->accept(solver.get());
            }

            if (auto solve_block = solver->as_block(false)) {
                // Check that we didn't solve an already solved variable.
                for (const auto& id: solver->solved_identifiers()) {
                    if (solved_ids.count(id)>0) {
                        error("Variable "+id+" solved twice!", e->location());
                        return false;
                    }
                    solved_ids.insert(id);
                }

                // May have now redundant local variables; remove these first.
                solve_block = remove_unused_locals(solve_block->is_block());

                // Copy body into nrn_state.
                for (auto& stmt: solve_block->is_block()->statements()) {
                    api_state->body()->statements().push_back(std::move(stmt));
                }
            }
            else {
                // Something went wrong: copy errors across.
                append_errors(solver->errors());
                return false;
            }
        }

        // handle the case where there is no SOLVE in BREAKPOINT
        if(!found_solve) {
            warning(" there is no SOLVE statement, required to update the"
                    " state variables, in the BREAKPOINT block",
                    breakpoint->location());
        }
        else {
            // redo semantic pass in order to elimate any removed local symbols.
            api_state->semantic(symbols_);
        }

        //..........................................................
        // nrn_current : update contributions to currents
        //..........................................................
        NrnCurrentRewriter nrn_current_rewriter(kind());
        breakpoint->accept(&nrn_current_rewriter);
        auto nrn_current_block = nrn_current_rewriter.as_block();
        if (!nrn_current_block) {
            append_errors(nrn_current_rewriter.errors());
            return false;
        }

        symbols_["nrn_current"] =
            make_symbol<APIMethod>(
                    breakpoint->location(), "nrn_current",
                    std::vector<expression_ptr>(),
                    constant_simplify(nrn_current_block));
        symbols_["nrn_current"]->semantic(symbols_);
    }
    else {
        error("a BREAKPOINT block is required");
        return false;
    }

    return !has_error();
}

/// populate the symbol table with class scope variables
void Module::add_variables_to_symbols() {
    // add reserved symbols (not v, because for some reason it has to be added
    // by the user)
    auto create_variable = [this] (const char* name, rangeKind rng, accessKind acc) {
        auto t = new VariableExpression(Location(), name);
        t->state(false);
        t->linkage(linkageKind::local);
        t->ion_channel(ionKind::none);
        t->range(rng);
        t->access(acc);
        t->visibility(visibilityKind::global);
        symbols_[name] = symbol_ptr{t};
    };

    create_variable("t",  rangeKind::scalar, accessKind::read);
    create_variable("dt", rangeKind::scalar, accessKind::read);
    // density mechanisms use a vector of weights from current densities to
    // units of nA
    if (kind()==moduleKind::density) {
        create_variable("weights_", rangeKind::range, accessKind::read);
    }

    // add indexed variables to the table
    auto create_indexed_variable = [this]
        (std::string const& name, std::string const& indexed_name,
         tok op, accessKind acc, ionKind ch, Location loc)
    {
        if(symbols_.count(name)) {
            throw compiler_exception(
                "trying to insert a symbol that already exists",
                loc);
        }
        symbols_[name] =
            make_symbol<IndexedVariable>(loc, name, indexed_name, acc, op, ch);
    };

    create_indexed_variable("current_", "vec_i", tok::plus,
                            accessKind::write, ionKind::none, Location());
    create_indexed_variable("v", "vec_v", tok::eq,
                            accessKind::read,  ionKind::none, Location());

    // add state variables
    for(auto const &var : state_block()) {
        VariableExpression *id = new VariableExpression(Location(), var.name());

        id->state(true);    // set state to true
        // state variables are private
        //      what about if the state variables is an ion concentration?
        id->linkage(linkageKind::local);
        id->visibility(visibilityKind::local);
        id->ion_channel(ionKind::none);    // no ion channel
        id->range(rangeKind::range);       // always a range
        id->access(accessKind::readwrite);

        symbols_[var.name()] = symbol_ptr{id};
    }

    // add the parameters
    for(auto const& var : parameter_block()) {
        auto name = var.name();
        if(name == "v") { // global voltage values
            // ignore voltage, which is added as an indexed variable by default
            continue;
        }
        VariableExpression *id = new VariableExpression(Location(), name);

        id->state(false);           // never a state variable
        id->linkage(linkageKind::local);
        // parameters are visible to Neuron
        id->visibility(visibilityKind::global);
        id->ion_channel(ionKind::none);
        // scalar by default, may later be upgraded to range
        id->range(rangeKind::scalar);
        id->access(accessKind::read);

        // check for 'special' variables
        if(name == "celcius") { // global celcius parameter
            id->linkage(linkageKind::external);
        }

        // set default value if one was specified
        if(var.value.size()) {
            id->value(std::stod(var.value));
        }

        symbols_[name] = symbol_ptr{id};
    }

    // add the assigned variables
    for(auto const& var : assigned_block()) {
        auto name = var.name();
        if(name == "v") { // global voltage values
            // ignore voltage, which is added as an indexed variable by default
            continue;
        }
        VariableExpression *id = new VariableExpression(var.token.location, name);

        id->state(false);           // never a state variable
        id->linkage(linkageKind::local);
        // local visibility by default
        id->visibility(visibilityKind::local);
        id->ion_channel(ionKind::none); // can change later
        // ranges because these are assigned to in loop
        id->range(rangeKind::range);
        id->access(accessKind::readwrite);

        symbols_[name] = symbol_ptr{id};
    }

    ////////////////////////////////////////////////////
    // parse the NEURON block data, and use it to update
    // the variables in symbols_
    ////////////////////////////////////////////////////
    // first the ION channels
    // add ion channel variables
    auto update_ion_symbols = [this, create_indexed_variable]
            (Token const& tkn, accessKind acc, ionKind channel)
    {
        auto const& var = tkn.spelling;

        // add the ion variable's indexed shadow
        if(has_symbol(var)) {
            auto sym = symbols_[var].get();

            // has the user declared a range/parameter with the same name?
            if(sym->kind()!=symbolKind::indexed_variable) {
                warning(
                    pprintf("the symbol % clashes with the ion channel variable,"
                            " and will be ignored", yellow(var)),
                    sym->location()
                );
                // erase symbol
                symbols_.erase(var);
            }
        }

        create_indexed_variable(var, "ion_"+var,
                                acc==accessKind::read ? tok::eq : tok::plus,
                                acc, channel, tkn.location);
    };

    // check for nonspecific current
    if( neuron_block().has_nonspecific_current() ) {
        auto const& i = neuron_block().nonspecific_current;
        update_ion_symbols(i, accessKind::write, ionKind::nonspecific);
    }


    for(auto const& ion : neuron_block().ions) {
        for(auto const& var : ion.read) {
            update_ion_symbols(var, accessKind::read, ion.kind());
        }
        for(auto const& var : ion.write) {
            update_ion_symbols(var, accessKind::write, ion.kind());
        }
    }

    // then GLOBAL variables
    for(auto const& var : neuron_block().globals) {
        if(!symbols_[var.spelling]) {
            error( yellow(var.spelling) +
                   " is declared as GLOBAL, but has not been declared in the" +
                   " ASSIGNED block",
                   var.location);
            return;
        }
        auto& sym = symbols_[var.spelling];
        if(auto id = sym->is_variable()) {
            id->visibility(visibilityKind::global);
        }
        else if (!sym->is_indexed_variable()){
            throw compiler_exception(
                "unable to find symbol " + yellow(var.spelling) + " in symbols",
                Location());
        }
    }

    // then RANGE variables
    for(auto const& var : neuron_block().ranges) {
        if(!symbols_[var.spelling]) {
            error( yellow(var.spelling) +
                   " is declared as RANGE, but has not been declared in the" +
                   " ASSIGNED or PARAMETER block",
                   var.location);
            return;
        }
        auto& sym = symbols_[var.spelling];
        if(auto id = sym->is_variable()) {
            id->range(rangeKind::range);
        }
        else if (!sym->is_indexed_variable()){
            throw compiler_exception(
                "unable to find symbol " + yellow(var.spelling) + " in symbols",
                var.location);
        }
    }
}

bool Module::optimize() {
    // how to structure the optimizer
    // loop over APIMethods
    //      - apply optimization to each in turn
    ConstantFolderVisitor folder;
    for(auto &symbol : symbols_) {
        auto kind = symbol.second->kind();
        BlockExpression* body;
        if(kind == symbolKind::procedure) {
            // we are only interested in true procedures and APIMethods
            auto proc = symbol.second->is_procedure();
            auto pkind = proc->kind();
            if(pkind == procedureKind::normal || pkind == procedureKind::api )
                body = symbol.second->is_procedure()->body();
            else
                continue;
        }
        // for now don't look at functions
        //else if(kind == symbolKind::function) {
        //    body = symbol.second.expression->is_function()->body();
        //}
        else {
            continue;
        }

        /////////////////////////////////////////////////////////////////////
        // loop over folding and propogation steps until there are no changes
        /////////////////////////////////////////////////////////////////////

        // perform constant folding
        for(auto& line : *body) {
            line->accept(&folder);
        }

        // preform expression simplification
        // i.e. removing zeros/refactoring reciprocals/etc

        // perform constant propogation

        /////////////////////////////////////////////////////////////////////
        // remove dead local variables
        /////////////////////////////////////////////////////////////////////
    }

    return true;
}
