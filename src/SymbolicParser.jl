#= 
# A set of tools for parsing the variational forms into symEngine expressions.
=#
module SymbolicParser
export sp_parse, add_custom_op, add_custom_op_file, sym_var

# See finch_import_symbols.jl for a list of all imported symbols.
import ..Finch: @import_finch_symbols
@import_finch_symbols()
import ..Finch: build_symexpressions

using SymEngine, LinearAlgebra

# An operator struct that ties a symbol to a function.
struct SymOperator
    symbol::Symbol      # The name used in the input expression
    op                  # Function handle for the operator
end

struct CallbackPlaceholder
    name::String
    args::Array
end
Base.show(io::IO, x::CallbackPlaceholder) = print(io, to_string(x));
function to_string(x::CallbackPlaceholder)
    argstr = "";
    for i=1:length(x.args)
        if i>1
            argstr *= ", ";
        end
        if typeof(x.args[i]) <: Array
            if length(x.args[i]) == 1
                argstr *= string(x.args[i][1]);
            else
                argstr *= "[";
                for j=1:length(x.args[i])
                    if j>1
                        argstr *= ", ";
                    end
                    argstr *= string(x.args[i][j]);
                end
                argstr *= "]";
            end
        else
            argstr *= string(x.args[i]);
        end
    end
    return x.name*"("*argstr*")"
end

#### globals ########################
# Basic symbolic operators that are included automatically and have precedence
op_names = [];
ops = [];

# Custom operators that can be added by a user
custom_op_names = [];
custom_ops = [];
#####################################

include("basic_ops.jl");
include("symbolic_stepper_reformat.jl");

# These are special function symbols that need to be defined.
# They can be used in operators.
# They are defined in symexpression.jl for now
@funs(conditional)
@funs(symbolmin)
@funs(symbolmax)
@funs(isgreaterthan)
@funs(islessthan)
@funs(indexing_operator)

# a special operator for dealing with scalar multiplication when the scalar is in an array
import Base.*
function *(a::Array{Basic,1}, b::Array{Basic,1})
    if size(a) == (1,)
        return b.*a[1];
    elseif size(b) == (1,)
        return a.*b[1];
    elseif size(a) == size(b)
        # This should be an error, but I will treat it as a dot product for lazy people
        return [transpose(a) * b];
    else
        return SymEngine.*(a,b);
    end
end

# Adds a single custom operator
function add_custom_op(s, handle)
    global custom_op_names;
    global custom_ops;
    push!(custom_op_names, s); # add the Symbol s to the list
    push!(custom_ops, SymOperator(s, handle));
    log_entry("Added custom operator: "*string(s), 2);
end

# Includes a set of custom operators
# The file will include an array of symbols and an array of function handles
function add_custom_op_file(file)
    include(file);
    for i=1:length(_names)
        add_custom_op(_names[i], _handles[i]);
    end
end

# Builds an array of symbolic layer symbols
function sym_var(name, type, components=1)
    if type == VAR_ARRAY
        # These must be indexed
        symvar = [symbols("_"*name*"_INDEXED")];
        
    else
        symvar = [symbols("_"*name*"_1")];
        for i=2:components
            push!(symvar, symbols("_"*name*"_"*string(i)));
        end
    end
    
    return symvar;
end

# Parses a variational form expression into a SymEngine Basic expression 
# Takes an Expr, variable symbol
# The symbols are only for determining lhs vs. rhs
# Returns a tuple of SymEngine expressions (lhs, rhs) for the equation lhs = rhs
# lhs contains terms including the unknown variable
# rhs contains terms without it
function sp_parse(ex, var; is_FV=false, is_flux=false)
    lhs = nothing;
    rhs = nothing;
    varcount = 1;
    timederiv = false;
    placeholders = [];
    
    log_entry("SP expr = "*string(ex), 3);
    
    # Check that there are as many vars as exs
    if typeof(var) <: Array
        if !(typeof(ex) <:Array && length(var) == length(ex))
            printerr("Error: Need same # of unknowns and equations");
            return (nothing,nothing);
        end
        varcount = length(var);
    end
    
    # Work with a copy of ex
    if typeof(ex) == Symbol
        symex = ex;
    else
        symex = copy(ex);
    end
    
    # Insert parameters
    symex = insert_parameters(symex);
    log_entry("SP insert parameters -> "*string(symex), 3);
    
    # Replace symbols for variables, coefficients, test functions, and special operators
    symex = replace_symbols(symex);
    log_entry("SP replace symbols -> "*string(symex), 3);
    
    # Look for indexes like [i] where i is an indexer
    if length(indexers) > 0
        symex = handle_indexers(symex);
        log_entry("SP symbolic indicies -> "*string(symex), 3);
    end
    
    # Look for callbacks like CALLBACK_f(a,b,c), process their arguments, and replace them with placeholder vars
    symex = handle_callbacks(symex, placeholders);
    log_entry("SP handle callbacks -> "*string(symex), 3);
    for i=1:length(placeholders)
        log_entry("    PLACEHOLDER_"*string(i)*": "*to_string(placeholders[i]), 3);
    end
    
    # change some operators like ^ and / to broadcast versions .^ ./
    symex = broadcast_ops(symex);
    
    # Evaluate the expression to apply symbolic operators
    symex = apply_ops(symex);
    log_entry("SP apply ops -> "*string(symex), 3);
    
    # If the result of this is not in the proper array structure, fix it
    if typeof(var) <: Array
        if !(typeof(symex) <: Array)
            symex = [[symex]];
        else
            for i=1:length(symex)
                if !(typeof(symex[i]) <: Array)
                    symex[i] = [symex[i]];
                end
            end
        end
    else
        if !(typeof(symex) <: Array)
            symex = [symex];
        end
    end
    
    
    # Expand the expression and separate terms
    sterms = get_sym_terms(symex);
    log_entry("SP sterms = "*string(sterms), 3);
    
    # Check for time derivatives
    timederiv = check_for_dt(sterms);
    
    # Check for surface integrals
    has_surface = check_for_surface(sterms);
    
    # Each element has an lhs and rhs
    lhs = copy(sterms); # set up the container right
    rhs = copy(sterms);
    if varcount > 1
        for i=1:length(symex)
            sz = size(symex[i]);
            (lhs[i],rhs[i]) = split_left_right(sterms[i],sz,var);
            if length(rhs[i]) == 0
                rhs[i] = [Basic(0)];
            end
        end
    else
        sz = size(symex);
        (lhs,rhs) = split_left_right(sterms,sz,var);
        if length(rhs[1]) == 0
            rhs[1] = [Basic(0)];
        end
    end
    
    # If there was a time derivative, separate those terms as well
    if timederiv
        dtlhs = copy(lhs);
        dtrhs = copy(rhs);
        if varcount > 1
            for i=1:length(symex)
                sz = size(symex[i]);
                (dtlhs[i],lhs[i]) = split_dt(lhs[i],sz);
                (dtrhs[i],rhs[i]) = split_dt(rhs[i],sz);
            end
        else
            sz = size(symex);
            (dtlhs,lhs) = split_dt(lhs,sz);
            (dtrhs,rhs) = split_dt(rhs,sz);
        end
    end
    
    # If there was a surface integral, separate those terms as well
    if has_surface
        surflhs = copy(lhs);
        surfrhs = copy(rhs);
        if varcount > 1
            for i=1:length(symex)
                sz = size(symex[i]);
                (surflhs[i],lhs[i]) = split_surf(lhs[i],sz);
                (surfrhs[i],rhs[i]) = split_surf(rhs[i],sz);
            end
        else
            sz = size(symex);
            (surflhs,lhs) = split_surf(lhs,sz);
            (surfrhs,rhs) = split_surf(rhs,sz);
        end
    end
    
    log_entry("SP LHS = "*string(lhs), 3);
    log_entry("SP RHS = "*string(rhs), 3);
    if timederiv
        log_entry("SP dtLHS = "*string(dtlhs), 3);
        log_entry("SP dtRHS = "*string(dtrhs), 3);
        if has_surface
            log_entry("SP surfLHS = "*string(surflhs), 3);
            log_entry("SP surfRHS = "*string(surfrhs), 3);
        end
        
    else
        if has_surface
            log_entry("SP surfLHS = "*string(surflhs), 3);
            log_entry("SP surfRHS = "*string(surfrhs), 3);
        end
    end
    
    # If needed, reformat for time stepper
    if timederiv || is_FV
        stepper_type = config.stepper;
        if is_FV
            if stepper_type == MIXED_STEPPER && typeof(time_stepper) <: Array
                stepper_type = time_stepper[2].type;
            end
            # There is an assumed Dt(u) added which is the only time derivative.
            if is_flux
                log_entry("flux, before modifying for time: "*string(lhs)*" - "*string(rhs));
                (newlhs, newrhs) = reformat_for_stepper_fv_flux(lhs, rhs, stepper_type);
                log_entry("flux, modified for time stepping: "*string(newlhs)*" + "*string(newrhs));
            else # source
                log_entry("source, before modifying for time: "*string(lhs)*" - "*string(rhs));
                (newlhs, newrhs) = reformat_for_stepper_fv_source(lhs, rhs, stepper_type);
                log_entry("source, modified for time stepping: "*string(newlhs)*" + "*string(newrhs));
            end
            # Parse Basic->Expr and insert placeholders
            newlhs = basic_to_expr_and_place(newlhs, placeholders)
            newrhs = basic_to_expr_and_place(newrhs, placeholders)
            # Build a SymExpression for each of the pieces. 
            (lhs_symexpr, rhs_symexpr) = build_symexpressions(var, newlhs, newrhs, remove_zeros=true);
            lhs = lhs_symexpr;
            rhs = rhs_symexpr;
            
        else # FE
            if stepper_type == MIXED_STEPPER && typeof(time_stepper) <: Array
                stepper_type = time_stepper[1].type;
            end
            if has_surface
                log_entry("Weak form, before modifying for time: Dt("*string(dtlhs)*") + "*string(lhs)*" + surface("*string(surflhs)*") = "*string(rhs)*" + surface("*string(surfrhs)*")");
                (newlhs, newrhs, newsurflhs, newsurfrhs) = reformat_for_stepper((dtlhs, lhs), rhs, surflhs, surfrhs, stepper_type);
                log_entry("Weak form, modified for time stepping: "*string(newlhs)*" + surface("*string(newsurflhs)*") = "*string(newrhs)*" + surface("*string(newsurfrhs)*")");
                # Parse Basic->Expr and insert placeholders
                newlhs = basic_to_expr_and_place(newlhs, placeholders)
                newrhs = basic_to_expr_and_place(newrhs, placeholders)
                newsurflhs = basic_to_expr_and_place(newsurflhs, placeholders)
                newsurfrhs = basic_to_expr_and_place(newsurfrhs, placeholders)
                # Build a SymExpression for each of the pieces. 
                (lhs_symexpr, rhs_symexpr, lhs_surf_symexpr, rhs_surf_symexpr) = build_symexpressions(var, newlhs, newrhs, newsurflhs, newsurfrhs, remove_zeros=true);
                lhs = lhs_symexpr;
                rhs = rhs_symexpr;
                surflhs = lhs_surf_symexpr;
                surfrhs = rhs_surf_symexpr;
            else # no surface
                log_entry("Weak form, before modifying for time: Dt("*string(dtlhs)*") + "*string(lhs)*" = "*string(rhs));
                (newlhs, newrhs) = reformat_for_stepper((dtlhs, lhs), rhs, stepper_type);
                log_entry("Weak form, modified for time stepping: "*string(newlhs)*" = "*string(newrhs));
                # Parse Basic->Expr and insert placeholders
                newlhs = basic_to_expr_and_place(newlhs, placeholders)
                newrhs = basic_to_expr_and_place(newrhs, placeholders)
                # Build a SymExpression for each of the pieces. 
                (lhs_symexpr, rhs_symexpr) = build_symexpressions(var, newlhs, newrhs, remove_zeros=true);
                lhs = lhs_symexpr;
                rhs = rhs_symexpr;
            end
        end
        
    else #FE with no time derivative or surface
        # Parse Basic->Expr and insert placeholders
        lhs = basic_to_expr_and_place(lhs, placeholders)
        rhs = basic_to_expr_and_place(rhs, placeholders)
        # Build a SymExpression for each of the pieces. 
        (lhs_symexpr, rhs_symexpr) = build_symexpressions(var, lhs, rhs, remove_zeros=true);
        lhs = lhs_symexpr;
        rhs = rhs_symexpr;
    end
    
    # Returns
    if has_surface
        return (lhs, rhs, surflhs, surfrhs);
    else
        return (lhs, rhs);
    end
end

# Replaces parameter symbols with their values
function insert_parameters(ex)
    if typeof(ex) == Symbol
        # parameter?
        for p in parameters
            if ex === p.symbol
                if p.type == SCALAR
                    return insert_parameters(p.value[1]);
                else
                    return insert_parameters(p.value);
                end
            end
        end
        # nope
        return ex;
        
    elseif typeof(ex) == Expr && length(ex.args) > 0
        for i=1:length(ex.args)
            ex.args[i] = insert_parameters(ex.args[i]); # recursively replace on args if ex
        end
        return ex;
    elseif typeof(ex) <:Array
        result = copy(ex);
        for i=1:length(ex)
            result[i] = insert_parameters(ex[i]);
        end
        return result;
    else
        return ex;
    end
end

# Replaces variable, coefficient and operator symbols in the expression
function replace_symbols(ex)
    if typeof(ex) <: Number
        return Basic(ex);
    elseif typeof(ex) == Symbol
        # variable?
        for v in variables
            if ex === v.symbol
                return v.symvar;
            end
        end
        # coefficient?
        for c in coefficients
            if ex === c.symbol
                # constant coefficients are entered as numbers
                if c.type == SCALAR && typeof(c.value[1]) <: Number
                    return [Basic(c.value[1])];
                end
                return c.symvar;
            end
        end
        # operator?
        for i=1:length(ops)
            if ex === ops[i].symbol
                #return Symbol(string(ops[i].op));
                return :(ops[$i].op);
            end
        end
        for i=1:length(custom_ops)
            if ex === custom_ops[i].symbol
                #return Symbol(string(custom_ops[i].op));
                return :(custom_ops[$i].op);
            end
        end
        # test function?
        for c in test_functions
            if ex === c.symbol
                return c.symvar;
            end
        end
        # callback functions
        for c in callback_functions
            if string(ex) == c.name
                # Make a special SymEngine function for this and process it later
                cname = Symbol("CALLBACK_" * c.name);
                cstr = string(cname);
                glop = Expr(:global, Expr(:(=), cname, Expr(:call, :(SymEngine.SymFunction), cstr)));
                eval(glop);
                return cname;
            end
        end
        # none of them?
        return ex;
    elseif typeof(ex) == Expr && length(ex.args) > 0
        for i=1:length(ex.args)
            ex.args[i] = replace_symbols(ex.args[i]); # recursively replace on args if ex
        end
        
        return ex;
    elseif typeof(ex) <:Array
        result = copy(ex);
        for i=1:length(ex)
            if typeof(ex[i]) == SymEngine.Basic
                result[i] = replace_symbols(ex[i])[1];
            else
                result[i] = replace_symbols(ex[i]);
            end
            
        end
        return result;
    else
        return ex;
    end
end

# Turns u[i] into u_INDEXEDBYi
function handle_indexers(ex)
    # Traverse the Expr looking for ref heads
    if typeof(ex) <:Array
        result = copy(ex);
        for i=1:length(ex)
            result[i] = handle_indexers(ex[i]);
        end
        return result;
        
    elseif typeof(ex) == Expr
        if ex.head === :ref
            needs_indexing = false;
            for i=2:length(ex.args)
                if !(typeof(ex.args[i]) <: Number)
                    needs_indexing = true;
                end
            end
            if needs_indexing
                
                var = ex.args[1];
                if typeof(var) <: Array && length(var) == 1 && occursin("INDEXED", string(var[1]))
                    newname = string(var[1]);
                    for i=2:length(ex.args)
                        newname *= "BY"*string(ex.args[i]);
                    end
                    return [symbols(newname)];
                else
                    printerr("Expected an array type variable when indexing: "*string(ex));
                end
                
                
                # newex = Expr(:call, :indexing_operator, ex.args[1]);
                # # the index symbols need to be turned into symengine things
                # for i=2:length(ex.args)
                #     if typeof(ex.args[i]) == Symbol
                #         ex.args[i] = symbols(string(ex.args[i]))
                #     end
                # end
                # append!(newex.args, ex.args[2:end]);
                # return newex;
            end
            
        else
            for i=1:length(ex.args)
                ex.args[i] = handle_indexers(ex.args[i]);
            end
        end
    end
    
    return ex;
end

# Turns callbackfun([a],[b],2) into PLACEHOLDER_i, processes the arguments ([a],[b],2), and stores in placeholders[i]
function handle_callbacks(ex, placeholders)
    # Traverse the Expr looking for :call with "CALLBACK_"
    if typeof(ex) <:Array
        result = copy(ex);
        for i=1:length(ex)
            result[i] = handle_callbacks(ex[i], placeholders);
        end
        return result;
        
    elseif typeof(ex) == Expr
        if ex.head === :call && occursin("CALLBACK_", string(ex.args[1]))
            
            # build a placeholder var
            pind = length(placeholders)+1;
            phvar = symbols("PLACEHOLDER_"*string(pind));
            # process the arguments
            for i=2:length(ex.args)
                ex.args[i] = apply_ops(ex.args[i]);
            end
            # build a callbackplaceholder struct
            push!(placeholders, CallbackPlaceholder(string(ex.args[1]), ex.args[2:end]))
            
            #swap in the phvar
            ex = phvar;
        else
            for i=1:length(ex.args)
                ex.args[i] = handle_callbacks(ex.args[i], placeholders);
            end
        end
    end
    
    return ex;
    
end

function broadcast_ops(ex)
    if typeof(ex) <:Array
        result = copy(ex);
        for i=1:length(ex)
            result[i] = broadcast_ops(ex[i]);
        end
        return result;
    elseif typeof(ex) == Expr
        for i=1:length(ex.args)
            ex.args[i] = broadcast_ops(ex.args[i]);
        end
        return ex;
    elseif typeof(ex) == Symbol
        # replace ^,/,+,- with .^,./,.+,.-
        if ex === :^ ex = :.^
        elseif ex === :/ ex = :./
        elseif ex === :* ex = :.*
        #elseif ex === :+ ex = :.+
        #elseif ex === :- ex = :.-
        end
    end
    return ex;
end

# Eval to apply the sym_*_op ops to create a SymEngine expression
function apply_ops(ex)
    # try
    #     if typeof(ex) <:Array
    #         result = [];
    #         for i=1:length(ex)
    #             #println("evaluating: "*string(ex[i]));
    #             push!(result, eval(ex[i]));
    #         end
    #         return result;
    #     else
    #         return eval(ex);
    #     end
    # catch e
    #     printerr("Problem evaluating the symbolic expression: "*string(ex));
    #     println(string(e));
    #     return 0;
    # end
    if typeof(ex) <:Array
        result = [];
        for i=1:length(ex)
            #println("evaluating: "*string(ex[i]));
            push!(result, eval(ex[i]));
        end
        return result;
    else
        return eval(ex);
    end
end

function has_unknown(ex, var)
    str = string(ex);
    result = false;
    if typeof(var) <: Array
        for i=1:length(var)
            vs = "_"*string(var[i])*"_";
            if occursin(vs, str)
                result = true;
            end
        end
    else
        vs = "_"*string(var)*"_";
        result = occursin(vs, str);
    end
    
    return result;
end

# Separate terms for each variable ex.
function get_sym_terms(ex)
    # Recursively work on each term of the array
    if typeof(ex) <: Array
        result = [get_sym_terms(ex[1])];
        for i=2:length(ex)
            push!(result,  get_sym_terms(ex[i]));
        end
        
        return result;
    end
    
    # ex is a symbolic expression(not array of them)
    # if ex is just a number, turn it into a Basic for expand
    if typeof(ex) <: Number
       ex = Basic(ex); 
    end
    # First expand it
    newex = expand(ex);
    log_entry("expanded: "*string(newex), 3);
    
    # Then separate the terms into an array
    return get_all_terms(newex);
end

function get_all_terms(ex)
    if typeof(ex) == Basic
        # convert to Expr, separate, convert to Basic
        expr = Meta.parse(string(ex));
        log_entry("Basic->Expr: "*string(expr), 3);
        
        terms = get_all_terms(expr);
        #println("Exprterms = "*string(terms));
        bterms = Array{Basic,1}(undef,0);
        for i=1:length(terms)
            if !(terms[i] == 0)
                push!(bterms, Basic(terms[i]));
            end
        end
        return bterms;
    end
    
    # At this point ex must be an Expr or symbol or number
    terms = [];
    if !(typeof(ex) == Expr)
        push!(terms, ex);
        return terms;
    end
    if ex.head === :call
        if ex.args[1] === :+ || ex.args[1] === :.+
            for i=2:length(ex.args)
                terms = append!(terms, get_all_terms(ex.args[i]));
            end
        elseif ex.args[1] === :-
            # apply -() to minused terms
            # Remember that the result of this will only be added terms, no minuses
            terms = get_all_terms(ex.args[2]);
            if length(ex.args) < 3
                for j=1:length(terms)
                    terms[j] = apply_negative(terms[j]);
                end
            else
                for i=3:length(ex.args)
                    tmp = get_all_terms(ex.args[i]);
                    for j=1:length(tmp)
                        tmp[j] = apply_negative(tmp[j]);
                    end
                    append!(terms, tmp);
                end
            end
        else
            push!(terms, ex);
        end
    else
        push!(terms, ex);
    end
    return terms;
end

function check_for_dt(terms)
    result = false;
    if typeof(terms) <: Array
        for i=1:length(terms)
            result = result || check_for_dt(terms[i]);
        end
    else
        result = occursin("TIMEDERIV", string(terms));
    end
    return result;
end

function check_for_surface(terms)
    result = false;
    if typeof(terms) <: Array
        for i=1:length(terms)
            result = result || check_for_surface(terms[i]);
        end
    else
        result = occursin("SURFACEINTEGRAL", string(terms));
    end
    return result;
end

function split_left_right(sterms,sz,var)
    lhs = copy(sterms); # set up the container right
    rhs = copy(sterms);
    if length(sz) == 0 # probably just a number
        rhs = sterms;
        lhs = [];
    elseif length(sz) == 1 # vector or scalar
        for i=1:sz[1]
            lhs[i] = Array{Basic,1}(undef,0);
            rhs[i] = Array{Basic,1}(undef,0);
            for ti=1:length(sterms[i])
                if has_unknown(sterms[i][ti], var)
                    #println("lhs: "*string(sterms[i][ti]));
                    push!(lhs[i], sterms[i][ti]);
                else
                    #println("rhs: "*string(sterms[i][ti]));
                    # switch sign to put on RHS
                    push!(rhs[i], -sterms[i][ti]);
                end
            end
        end
    elseif length(sz) == 2 # matrix
        for j=1:sz[2]
            for i=1:sz[1]
                lhs[i,j] = Basic(0);
                rhs[i,j] = Basic(0);
                for ti=1:length(sterms[i,j])
                    if has_unknown(sterms[i,j][ti], var)
                        #println("lhs: "*string(sterms[i,j][ti]));
                        push!(lhs[i,j], sterms[i,j][ti]);
                    else
                        #println("rhs: "*string(sterms[i,j][ti]));
                        push!(rhs[i,j], -sterms[i,j][ti]);
                    end
                end
            end
        end
    elseif length(sz) == 3 # rank 3
        for k=1:sz[3]
            for j=1:sz[2]
                for i=1:sz[1]
                    lhs[i,j,k] = Basic(0);
                    rhs[i,j,k] = Basic(0);
                    for ti=1:length(sterms[i,j,k])
                        if has_unknown(sterms[i,j,k][ti], var)
                            #println("lhs: "*string(sterms[i,j,k][ti]));
                            push!(lhs[i,j,k], sterms[i,j,k][ti]);
                        else
                            #println("rhs: "*string(sterms[i,j,k][ti]));
                            push!(rhs[i,j,k], -sterms[i,j,k][ti]);
                        end
                    end
                end
            end
        end
    end
    return (lhs,rhs);
end

function split_dt(terms, sz)
    hasdt = copy(terms); # set up the container right
    nodt = copy(terms);
    TIMEDERIV = symbols("TIMEDERIV"); # will be removed from terms
    if length(sz) == 1 # vector or scalar
        for i=1:sz[1]
            hasdt[i] = Array{Basic,1}(undef,0);
            nodt[i] = Array{Basic,1}(undef,0);
            for ti=1:length(terms[i])
                if check_for_dt(terms[i][ti])
                    #println("hasdt: "*string(terms[i][ti]));
                    terms[i][ti] = subs(terms[i][ti], TIMEDERIV=>1);
                    push!(hasdt[i], terms[i][ti]);
                else
                    #println("nodt: "*string(terms[i][ti]));
                    push!(nodt[i], terms[i][ti]);
                end
            end
        end
    elseif length(sz) == 2 # matrix
        for j=1:sz[2]
            for i=1:sz[1]
                hasdt[i,j] = Basic(0);
                nodt[i,j] = Basic(0);
                for ti=1:length(terms[i,j])
                    if check_for_dt(terms[i,j][ti])
                        #println("hasdt: "*string(terms[i,j][ti]));
                        terms[i,j][ti] = subs(terms[i,j][ti], TIMEDERIV=>1);
                        push!(hasdt[i,j], terms[i,j][ti]);
                    else
                        #println("nodt: "*string(terms[i,j][ti]));
                        push!(nodt[i,j], terms[i,j][ti]);
                    end
                end
            end
        end
    elseif length(sz) == 3 # rank 3
        for k=1:sz[3]
            for j=1:sz[2]
                for i=1:sz[1]
                    hasdt[i,j,k] = Basic(0);
                    nodt[i,j,k] = Basic(0);
                    for ti=1:length(terms[i,j,k])
                        if check_for_dt(terms[i,j,k][ti])
                            #println("hasdt: "*string(terms[i,j,k][ti]));
                            terms[i,j,k][ti] = subs(terms[i,j,k][ti], TIMEDERIV=>1);
                            push!(hasdt[i,j,k], terms[i,j,k][ti]);
                        else
                            #println("nodt: "*string(terms[i,j,k][ti]));
                            push!(nodt[i,j,k], terms[i,j,k][ti]);
                        end
                    end
                end
            end
        end
    end
    
    return (hasdt, nodt);
end

function split_surf(terms, sz)
    hassurf = copy(terms); # set up the container right
    nosurf = copy(terms);
    SURFACEINTEGRAL = symbols("SURFACEINTEGRAL"); # will be removed from terms
    if length(sz) == 1 # vector or scalar
        for i=1:sz[1]
            hassurf[i] = Array{Basic,1}(undef,0);
            nosurf[i] = Array{Basic,1}(undef,0);
            for ti=1:length(terms[i])
                if check_for_surface(terms[i][ti])
                    terms[i][ti] = subs(terms[i][ti], SURFACEINTEGRAL=>1);
                    push!(hassurf[i], terms[i][ti]);
                else
                    push!(nosurf[i], terms[i][ti]);
                end
            end
        end
    elseif length(sz) == 2 # matrix
        for j=1:sz[2]
            for i=1:sz[1]
                hassurf[i,j] = Basic(0);
                nosurf[i,j] = Basic(0);
                for ti=1:length(terms[i,j])
                    if check_for_surface(terms[i,j][ti])
                        terms[i,j][ti] = subs(terms[i,j][ti], SURFACEINTEGRAL=>1);
                        push!(hassurf[i,j], terms[i,j][ti]);
                    else
                        push!(nosurf[i,j], terms[i,j][ti]);
                    end
                end
            end
        end
    elseif length(sz) == 3 # rank 3
        for k=1:sz[3]
            for j=1:sz[2]
                for i=1:sz[1]
                    hassurf[i,j,k] = Basic(0);
                    nosurf[i,j,k] = Basic(0);
                    for ti=1:length(terms[i,j,k])
                        if check_for_surface(terms[i,j,k][ti])
                            terms[i,j,k][ti] = subs(terms[i,j,k][ti], SURFACEINTEGRAL=>1);
                            push!(hassurf[i,j,k], terms[i,j,k][ti]);
                        else
                            push!(nosurf[i,j,k], terms[i,j,k][ti]);
                        end
                    end
                end
            end
        end
    end
    
    return (hassurf, nosurf);
end

# Parse Basic into Expr and insert placeholders
# This works recursively so each piece in the array is parsed, then the Expr is traversed to swap placeholders.
function basic_to_expr_and_place(ex, placeholders)
    if typeof(ex) <: Array
        exarray = [];
        for i=1:length(ex)
            push!(exarray, basic_to_expr_and_place(ex[i], placeholders));
        end
        return exarray;
        
    elseif typeof(ex) <: SymEngine.Basic
        # ex is a SymEngine.Basic object
        # First parse it into an Expr
        jex = Meta.parse(string(ex));
        # Then process it
        newex = basic_to_expr_and_place(jex, placeholders);
        
        return newex;
        
    elseif typeof(ex) == Expr
        for i=1:length(ex.args)
            ex.args[i] = basic_to_expr_and_place(ex.args[i], placeholders)
        end
        return ex;
        
    elseif typeof(ex) == Symbol
        # Here's where the replacement happens
        str = string(ex);
        if occursin("PLACEHOLDER_", str)
            phind = parse(Int, str[13:end]);
            ex = Meta.parse(to_string(placeholders[phind]));
            
            if typeof(ex) == Expr
                for i=1:length(ex.args)
                    ex.args[i] = basic_to_expr_and_place(ex.args[i], placeholders)
                end
            end
        end
        
        return ex;
        
    else # could be a number or string or something
        return ex;
    end
end

function apply_negative(ex)
    negex = :(-a);
    if typeof(ex) == Symbol
        negex.args[2] = ex;
        return negex;
    elseif typeof(ex) <: Number
        return -ex;
    elseif typeof(ex) == Expr
        if ex.args[1] == :- && length(ex.args) == 2
            return ex.args[2];
        else
            negex.args[2] = ex;
            return negex;
        end
    end
end

# Applies a string flag to all free_symbols in ex
function apply_flag_to_all_symbols(flag, ex)
    if typeof(ex) <: Array
        for i=1:length(ex)
            ex[i] = apply_flag_to_all_symbols(flag, ex[i]);
        end
        
    elseif typeof(ex) <: Basic
        symbs = free_symbols(ex);
        for i=1:length(symbs)
            ex = subs(ex, symbs[i], symbols(flag*"_"*string(symbs[i])));
        end
    end
    
    return ex;
end

end # module