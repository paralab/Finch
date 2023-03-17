# Takes a SymExpression like 2*a+b and makes an IR treating all symbols
# as named coefficient arrays.
# This works recursively
function arithmetic_expr_to_IR(ex)
    IRtypes = IR_entry_types();
    if typeof(ex) === nothing
        return nothing;
    elseif typeof(ex) <: Number
        return ex;
    elseif typeof(ex) == Symbol
        return ex;
    elseif typeof(ex) == SymEntity
        if ex.index == -1 # This could be a constant number or special symbol
            val = tryparse(Float64, ex.name);
            if val === nothing
                return Symbol(ex.name);
            else
                return val;
            end
        end
        # These are special symbols that don't need any modifiers
        specials = ["ELEMENTDIAMETER", "ELEMENTVOLUME", "BOUNDARYVALUE", "FACENORMAL1", "FACENORMAL2", "TRUENORMAL", "DIST2BDRY"];
        
        # It is a named coefficient
        tag = "";
        type_label = "value_";
        for i=1:length(ex.flags)
            tag = ex.flags[i] * tag;
            if ex.flags[i] == "NEIGHBORHOOD"
                tmp = type_label[1];
                type_label = "noBroadcast_";
            end
        end
        for i=1:length(ex.derivs)
            tag = "D"*string(ex.derivs[i]) * tag;
        end
        
        if string(ex.name) in specials
            str = string(ex.name);
        else
            str = type_label*tag*"_"*string(ex.name);
        end
        
        if typeof(ex.index) == Int
            str *= "_"*string(ex.index);
        else
            str *= "_";
            for i=1:length(ex.index)
                str *= string(ex.index[i]);
            end
        end
        
        return IR_data_node(IRtypes.float_data, Symbol(str), [1], []);
        
    elseif typeof(ex) == Expr && ex.head === :call
        args = [];
        for i=1:length(ex.args)
            push!(args, arithmetic_expr_to_IR(ex.args[i]));
        end
        return IR_operation_node(IRtypes.math_op, args);
        
    elseif typeof(ex) == SymExpression
        # turn the tree expr into IR
        return arithmetic_expr_to_IR(ex.tree);
    else # What else could it be?
        return ex;
    end
end

# For an array type IR_data_node this adds an index
# (A, i) -> A[i]
function apply_indexed_access(IR, index, IRtypes::Union{IR_entry_types, Nothing} = nothing)
    if typeof(IR) == IR_data_node && length(IR.size) > 0
        return IR_data_node(IR.type, IR.label, IR.size, index);
        
    elseif typeof(IR) == IR_operation_node
        # make a copy
        newargs = copy(IR.args)
        # apply to all args
        for i=1:length(newargs)
            newargs[i] = apply_indexed_access(newargs[i], index, IRtypes)
        end
        return IR_operation_node(IR.type, newargs);
        
    else
        return IR;
    end
end

# Traverses the IR looking for data nodes with this label.
# Remove size and index info to make it a scalar.
# This is useful for GPU things.
function scalarize_data_nodes(IR, symbol)
    if typeof(IR) == IR_data_node && IR.label === symbol
        IR = IR_data_node(IR.type, IR.label, [], []);
        
    elseif typeof(IR) == IR_operation_node
        # apply to all args
        for i=1:length(IR.args)
            IR.args[i] = scalarize_data_nodes(IR.args[i], symbol)
        end
        
    elseif typeof(IR) == IR_block_node
        # apply to all IR_parts
        for i=1:length(IR.parts)
            IR.parts[i] = scalarize_data_nodes(IR.parts[i], symbol)
        end
        
    elseif typeof(IR) == IR_loop_node
        # apply to body
        IR.body = scalarize_data_nodes(IR.body, symbol)
        
    elseif typeof(IR) == IR_conditional_node
        # apply to body
        IR.body = scalarize_data_nodes(IR.body, symbol)
        IR.elsepart = scalarize_data_nodes(IR.elsepart, symbol)
    end
    
    return IR;
end

# Wraps some IR in a timeroutputs timer with a given label
# It is up to the caller to keep labels unique
function wrap_in_timer(label::Symbol, content::IR_part)
    IRtypes = IR_entry_types();
    return IR_operation_node(IRtypes.named_op, [:TIMER, label, content]);
end

function extract_entities(symex::Array)
    # symex is an array of arrays of SymExpressions which are Expr trees with SymEntities as leaves. (array for variable components, terms)
    # In the case of muliple variables, it's an array of arrays of arrays. (variables, components, terms)
    # The Symexpression contains a list of all leaves that need to be evaluated before combining.
    # First collect all of them and eliminate identical ones.
    entities = []
    for vi=1:length(symex)
        for ci=1:length(symex[vi])
            for ti=1:length(symex[vi][ci])
                for i=1:length(symex[vi][ci][ti].entities) # loop over entities for this variable/component
                    entity_present = false;
                    for j=1:length(entities) # check against existing entities
                        if is_same_entity(symex[vi][ci][ti].entities[i], entities[j])
                            entity_present = true;
                            break;
                        end
                    end
                    if !entity_present
                        push!(entities, symex[vi][ci][ti].entities[i]); # add it if unique
                    end
                end
            end
        end
    end
    
    return entities;
end

# Check if this is just an empty set of arrays
function is_empty_expression(ex)
    if typeof(ex) <: Array
        for i=1:length(ex)
            if !(is_empty_expression(ex[i]))
                return false;
            end
        end
    else # this must be some object in the array
        return false;
    end
    # If we get here, ex was an empty array
    return true;
end

# They are the same if all parts of them are the same.
function is_same_entity(a::SymEntity, b::SymEntity)
    return a.name == b.name && a.index == b.index && a.derivs == b.derivs && a.flags == b.flags;
end

function is_test_function(ent)
    for t in finch_state.test_functions
        if string(t.symbol) == ent.name
            return true;
        end
    end
    return false;
end

function is_unknown_var(ent, vars)
    if typeof(vars) <: Variable
        return ent.name == string(vars.symbol);
    elseif typeof(vars) <:Array
        for v in vars
            if string(v.symbol) == ent.name
                return true;
            end
        end
    end
    return false;
end

# Makes a name for this value based on coefficient or variable name, derivative and other modifiers, and component index.
# c is a SymEntity
function make_entity_name(c::SymEntity)
    # Special case for dt or other special symbols
    if c.index == -1
        return c.name;
    end
    # These are special symbols that don't need any modifiers
    specials = ["ELEMENTDIAMETER", "ELEMENTVOLUME", "BOUNDARYVALUE", "FACENORMAL1", "FACENORMAL2", "TRUENORMAL", "DIST2BDRY"];
    
    tag = "";
    type_label = "value_";
    for i=1:length(c.flags)
        tag = c.flags[i] * tag;
        if c.flags[i] == "NEIGHBORHOOD"
            tmp = type_label[1];
            type_label = "noBroadcast_";
        end
    end
    for i=1:length(c.derivs)
        tag = "D"*string(c.derivs[i]) * tag;
    end
    
    if string(c.name) in specials
        str = string(c.name);
    else
        str = type_label*tag*"_"*string(c.name);
    end
    
    if typeof(c.index) == Int
        str *= "_"*string(c.index);
    else
        str *= "_";
        for i=1:length(c.index)
            str *= string(c.index[i]);
        end
    end
    
    return str;
end

# Checks if the coefficient has constant value.
# If so, also returns the value.
function is_constant_coef(c)
    isit = false;
    val = 0;
    for i=1:length(finch_state.coefficients)
        if c === finch_state.coefficients[i].symbol
            isit = (typeof(finch_state.coefficients[i].value[1]) <: Number);
            if isit
                val = finch_state.coefficients[i].value[1];
            else
                val = finch_state.coefficients[i].value[1].name;
            end
        end
    end
    
    return (isit, val);
end

# Checks the type of coefficient: constant, genfunction, or variable
# Returns: (type, val)
# special: type=-1, val=0
# number: type=0, val=number
# constant coef: type=1, val=number
# genfunction: type=2, val= index in genfunctions array
# variable: type=3, val=index in variables array
# indexed coefficient: type=4, val=index in coefficients array
function get_coef_val(c)
    if c.index == -1
        # It is a special symbol like dt
        return (-1, 0);
    elseif c.index == 0
        # It is a number
        return(0, c.name);
    end
    
    type = 0;
    val = 0;
    for i=1:length(finch_state.coefficients)
        if c.name == string(finch_state.coefficients[i].symbol)
            if typeof(c.index) == Int
                isit = (typeof(finch_state.coefficients[i].value[c.index]) <: Number);
                if isit
                    type = 1; # a constant wrapped in a coefficient ... not ideal
                    val = finch_state.coefficients[i].value[c.index];
                else
                    type = 2; # a function
                    name = finch_state.coefficients[i].value[c.index].name;
                    for j=1:length(finch_state.genfunctions)
                        if name == finch_state.genfunctions[j].name
                            val = j;
                        end
                    end
                end
            else
                # indexed coefficient
                type = 4;
                val = i;
            end
        end
    end
    if type == 0 # a variable
        for i=1:length(finch_state.variables)
            if c.name == string(finch_state.variables[i].symbol)
                type = 3;
                val = finch_state.variables[i].index;
            end
        end
    end
    if type == 0 # something special with a nonzero index
        return (-1,0);
    end
    
    return (type, val);
end

# Returns the index in finch's coefficient array
# or -1 if it is not there.
function get_coef_index(c)
    ind = -1;
    for i=1:length(finch_state.coefficients)
        if c.name == string(finch_state.coefficients[i].symbol)
            ind = finch_state.coefficients[i].index;
        end
    end
    
    return ind;
end

# Which side of a face is this entity on?
# Return a number: 0 = no side specified, 1/2 = side 1/2 in face2element order, 3 = average both, 4 = neighborhood
function get_face_side_info(ent)
    side = 0; # 0 means no side flag
    if typeof(ent) == SymEntity
        for flagi=1:length(ent.flags)
            if occursin("DGSIDE1", ent.flags[flagi]) || occursin("CELL1", ent.flags[flagi])
                return 1;
            elseif occursin("DGSIDE2", ent.flags[flagi]) || occursin("CELL2", ent.flags[flagi])
                return 2;
            elseif occursin("CENTRAL", ent.flags[flagi])
                return 3;
            elseif occursin("NEIGHBORHOOD", ent.flags[flagi])
                return 4;
            end
        end
        
    elseif typeof(ent) == Expr && length(ent.args) == 2 # This should only be allowed for negative signs :(-ent)
        side = get_face_side_info(ent.args[2]);
        
    # else number or ?? -> 0
    end
    
    return side;
end

# Determine if this expression has an unknown variable entity
# If so, return it's index.
# If not, return -1.
function get_first_var_index(ex, var)
    result = -1;
    if typeof(ex) == Expr
        # recursively check
        for i=1:length(ex.args)
            tmp = get_first_var_index(ex.args[i], var);
            if tmp > 0
                return tmp;
            end
        end
        
    elseif typeof(ex) == SymEntity
        numind = ex.index;
        if typeof(numind) <: Array
            numind = 1;
        end
        if typeof(var) <: Array
            tmpind = 0;
            for vi=1:length(var)
                if is_unknown_var(ex,var[vi])
                    result = tmpind + numind;
                else
                    tmpind += length(var[vi].symvar);
                end
            end
        else
            if is_unknown_var(ex,var)
                result = numind;
            end
        end
    end
    
    return result;
end

# Assume the term looks like a*b*c or a*(b*c) or something similar.
function separate_factors(ex, var=nothing)
    test_part = nothing;
    trial_part = nothing;
    coef_part = nothing;
    test_ind = 0;
    trial_ind = 0;
    
    if typeof(ex) == Expr
        if (ex.args[1] === :.- || ex.args[1] === :-) && length(ex.args)==2
            # a negative sign. Stick it on one of the factors
            (test_part, trial_part, coef_part, test_ind, trial_ind) = separate_factors(ex.args[2], var);
            negex = :(-a);
            if !(coef_part === nothing)
                negex.args[2] = coef_part;
                coef_part = negex;
            elseif !(test_part === nothing)
                negex.args[2] = test_part;
                test_part = negex;
            elseif !(trial_part === nothing)
                negex.args[2] = trial_part;
                trial_part = negex;
            end
            
        elseif ex.args[1] === :.* || ex.args[1] === :*
            tmpcoef = [];
            # Recursively separate each arg to handle a*(b*(c*...))
            for i=2:length(ex.args)
                (testi, triali, coefi, testindi, trialindi) = separate_factors(ex.args[i], var);
                if !(testi === nothing)
                    test_part = testi;
                    test_ind = testindi;
                end
                if !(triali === nothing)
                    trial_part = triali;
                    trial_ind = trialindi;
                end
                if !(coefi === nothing)
                    push!(tmpcoef, coefi);
                end
            end
            if length(tmpcoef) > 1
                coef_part = :(a*b);
                coef_part.args = [:*];
                append!(coef_part.args, tmpcoef);
            elseif length(tmpcoef) == 1
                coef_part = tmpcoef[1];
            end
            
        else# It is an Expr, but not -() or * (note: at this point a/b is a*(1/b) )
            # It could be a function like conditional_function.
            # Search for unknown variables within and use the first one found for the index.
            # If non are found, treat as coefficient part
            tmp_index = get_first_var_index(ex, var)
            if tmp_index > 0
                trial_ind = tmp_index;
                trial_part = ex;
            else
                coef_part = ex;
            end
        end
        
    elseif typeof(ex) == SymEntity
        numind = ex.index;
        if typeof(numind) <: Array
            numind = 1;
        end
        if is_test_function(ex)
            test_part = ex;
            test_ind = numind;
        elseif !(var === nothing) && is_unknown_var(ex, var)
            trial_part = ex;
            # find the trial index
            if typeof(var) <: Array
                tmpind = 0;
                for vi=1:length(var)
                    if is_unknown_var(ex,var[vi])
                        trial_ind = tmpind + numind;
                    else
                        tmpind += length(var[vi].symvar);
                    end
                end
            else
                trial_ind = numind;
            end
        else
            coef_part = ex;
        end
    else # a number?
        coef_part = ex;
    end
    
    return (test_part, trial_part, coef_part, test_ind, trial_ind);
end

# symex could be an array. If so, make a similar array containing an array of terms.
# This does any final processing of the terms. At this point symex should be a single term 
# without a top level + or -. for example, a*b/c or a^2*b but not a*b+c
# What processing is done:
#   a/b -> a*(1/b)
#   Broadcast all ops: +-*/^ -> .+.-.*./.^
#   If it's a SymExpression, return the proccessed Expr in symex.tree
function process_terms(symex)
    if typeof(symex) <: Array
        terms = [];
        for i=1:length(symex)
            push!(terms, process_terms(symex[i]));
        end
        
    else
        if typeof(symex) == SymExpression
            newex = copy(symex.tree);
        elseif typeof(symex) == Expr
            newex = copy(symex);
        else
            newex = symex;
        end
        
        #println(string(newex) * " : " * string(typeof(newex)))
        
        if typeof(newex) == Expr && newex.head === :call
            if (newex.args[1] === :+ || newex.args[1] === :.+ || 
                (newex.args[1] === :- || newex.args[1] === :.-) && length(newex.args) > 2 )
                println("unexpected term in process terms: "*string(newex));
            else
                if newex.args[1] === :/ || newex.args[1] === :./
                    # change a/b to a*1/b
                    mulex = :(a.*b);
                    mulex.args[2] = newex.args[2];
                    newex.args[2] = 1;
                    mulex.args[3] = newex;
                    newex = mulex;
                end
            end
            
            # insert indices indexing_operator(u, i, j) -> u[index_val_i, index_val_j]
            newex = insert_indices(newex);
            
            # broadcast ops
            # newex = broadcast_all_ops(newex);
            
            terms = newex;
        else
            # There is just one factor
            terms = newex;
        end
    end
    
    return terms;
end

# Turn entities with index=[i,j] into u[index_val_i, index_val_j]
function insert_indices(ex)
    if typeof(ex) <: Array
        result = [];
        for i=1:length(ex)
            push!(result, insert_indices(ex[i]));
        end
        return result;
        
    elseif typeof(ex) == SymEntity
        if typeof(ex.index) <: Array
            indices = copy(ex.index);
            
            for i=1:length(indices)
                indices[i] = Symbol("index_val_"*string(indices[i]));
            end
            ex = Expr(:ref, indexed_thing);
            append!(ex.args, indices);
            
        else
            for i=1:length(ex.args)
                ex.args[i] = insert_indices(ex.args[i]);
            end
        end
        
    end
    
    return ex;
end

# Combine terms that have the same factor
# a * b * c + a * d * c -> a * (c * (b + d))
# Assume that each term will look like *(a,b,c) or *(a,b)
function combine_shared_factors(terms)
    nt = length(terms);
    newterms = Vector{IR_part}(undef,0);
    test_factors = [];
    coef_factors = [];
    trial_factors = [];
    test_indices = [zeros(Int,0)]; # indices in terms for matching factors in test_factors
    coef_indices = [zeros(Int,0)];
    trial_indices = [zeros(Int,0)];
    term2test = zeros(Int,nt);
    term2coef = zeros(Int,nt);
    term2trial = zeros(Int,nt);
    test_symbols = [];
    coef_symbols = [];
    trial_symbols = [];
    
    symbolic_expr = Basic(0);
    
    # collect all unique factors of each type and make lists of terms that include them
    # also set up a symengine expression with variables like t1, c1, r1
    for i=1:nt
        if typeof(terms[i]) == IR_operation_node && terms[i].args[1] === :*
            nf = length(terms[i].args) - 1;
            test = terms[i].args[2];
            coef = terms[i].args[3];
            test_sym = Basic(1);
            coef_sym = Basic(1);
            trial_sym = Basic(1);
            
            # Check if they have been found already
            found = false;
            for j=1:length(test_factors)
                if string(test) == string(test_factors[j])
                    found = true;
                    push!(test_indices[j],i);
                    term2test[i] = j;
                    test_sym = test_symbols[j];
                end
            end
            if !found
                # add it to test_factors
                push!(test_factors, test);
                push!(test_indices, [i]);
                term2test[i] = length(test_factors);
                test_sym = symbols("t"*string(length(test_factors)))
                push!(test_symbols, test_sym);
            end
            
            found = false;
            for j=1:length(coef_factors)
                if string(coef) == string(coef_factors[j])
                    found = true;
                    push!(coef_indices[j],i);
                    term2coef[i] = j;
                    coef_sym = coef_symbols[j];
                end
            end
            if !found
                # add it to coef_factors
                push!(coef_factors, coef);
                push!(coef_indices, [i]);
                term2coef[i] = length(coef_factors);
                coef_sym = symbols("c"*string(length(coef_factors)))
                push!(coef_symbols, coef_sym);
            end
            
            if nf == 4 # there is a trial factor
                trial = terms[i].args[3];
                found = false;
                for j=1:length(trial_factors)
                    if string(trial) == string(trial_factors[j])
                        found = true;
                        push!(trial_indices[j],i);
                        term2trial[i] = j;
                        trial_sym = trial_symbols[j];
                    end
                end
                if !found
                    # add it to trial_factors
                    push!(trial_factors, trial);
                    push!(trial_indices, [i]);
                    term2trial[i] = length(trial_factors);
                    trial_sym = symbols("c"*string(length(trial_factors)));
                    push!(trial_symbols, trial_sym);
                end
            end
            
            # At this point we have all the symbols.
            # Add the term to the expression
            symbolic_expr = symbolic_expr + test_sym * coef_sym * trial_sym;
        end
    end
    
    # The symbolic expression has been built. Use symengine to simplify.
    # simple_expr = 
    
    return newterms;
end

# Given the factor info, group into simplified terms
# output is [[t,i,[[t,j,[[t,k], ...]], ...]], ...]
# for fi*(fj*(fk+...)+...)+...
# t = 1,2,3 for test,coef,trial
# i,j,k are indices in those arrays
function get_factor_groups(test, coef, trial, term2test, term2coef, term2trial)
    remaining = fill(false,nt);
    num_remaining = nt;
    
    # Find the most common from remaining factors
    largest = 0;
    largest_ind = 0;
    largest_type = 0;
    for i=1:nt
        if remaining[i] # This term has not been included yet
            num_common = 0;
            f_ind = term2trial[i]
            for j=1:length(test_indices[f_ind])
                if remaining[test_indices[f_ind][j]]
                    num_common+=1;
                end
            end
            if num_common > largest
                largest = num_common;
                largest_ind = f_ind;
                largest_type = 1;
            end
            
            num_common = 0;
            f_ind = term2coef[i];
            for j=1:length(coef_indices[f_ind])
                if remaining[coef_indices[f_ind][j]]
                    num_common+=1;
                end
            end
            if num_common > largest
                largest = num_common;
                largest_ind = f_ind;
                largest_type = 2;
            end
            
            num_common = 0;
            f_ind = term2trial[i];
            if f_ind > 0
                for j=1:length(trial_indices[f_ind])
                    if remaining[trial_indices[f_ind][j]]
                        num_common+=1;
                    end
                end
                if num_common > largest
                    largest = num_common;
                    largest_ind = f_ind;
                    largest_type = 3;
                end
            end
        end
    end
    
    # That most common factor will be put here -> a*(...) + ...
    
    this_term = [largest_type, largest_ind, []]
    
end