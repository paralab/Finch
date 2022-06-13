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
        
        if typeof(ex.index) == Int
            str = type_label*tag*"_"*string(ex.name)*"_"*string(ex.index);
        else
            str = type_label*tag*"_"*string(ex.name)*"_";
            for i=1:length(ex.index)
                str *= string(ex.index[i]);
            end
        end
        
        return IR_data_node(IRtypes.array_data, Symbol(str));
        
    elseif typeof(ex) == Expr && ex.head === :call
        args = [];
        for i=2:length(ex.args)
            push!(args, arithmetic_expr_to_IR(ex.args[i]));
        end
        return IR_eval_node(IRtypes.math_eval, ex.head, args);
        
    elseif typeof(ex) == SymExpression
        # turn the tree expr into IR
        return arithmetic_expr_to_IR(ex.tree);
    else # What else could it be?
        return ex;
    end
end

function extract_entities(symex::Array, multivar::Bool)
    # symex is an array of arrays of SymExpressions which are Expr trees with SymEntities as leaves. (array for variable components, terms)
    # In the case of muliple variables, it's an array of arrays of arrays. (variables, components, terms)
    # The Symexpression contains a list of all leaves that need to be evaluated before combining.
    # First collect all of them and eliminate identical ones.
    entities = []
    if multivar
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
    else # same thing as above, but for symex rather than symex[vi]
        for ci=1:length(symex)
            for ti=1:length(symex[ci])
                for i=1:length(symex[ci][ti].entities) # loop over entities for this variable/component
                    entity_present = false;
                    for j=1:length(entities) # check against existing entities
                        if is_same_entity(symex[ci][ti].entities[i], entities[j])
                            entity_present = true;
                            break;
                        end
                    end
                    if !entity_present
                        push!(entities, symex[ci][ti].entities[i]); # add it if unique
                    end
                end
            end
        end
    end
    
    return entities;
end

# They are the same if all parts of them are the same.
function is_same_entity(a::SymEntity, b::SymEntity)
    return a.name == b.name && a.index == b.index && a.derivs == b.derivs && a.flags == b.flags;
end

function is_test_function(ent)
    for t in test_functions
        if string(t.symbol) == ent.name
            return true;
        end
    end
    return false;
end

function is_unknown_var(ent, vars)
    if typeof(vars) == Variable
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
    
    if typeof(c.index) == Int
        str = type_label*tag*"_"*string(c.name)*"_"*string(c.index);
    else
        str = type_label*tag*"_"*string(c.name)*"_";
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
    for i=1:length(coefficients)
        if c === coefficients[i].symbol
            isit = (typeof(coefficients[i].value[1]) <: Number);
            if isit
                val = coefficients[i].value[1];
            else
                val = coefficients[i].value[1].name;
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
    for i=1:length(coefficients)
        if c.name == string(coefficients[i].symbol)
            if typeof(c.index) == Int
                isit = (typeof(coefficients[i].value[c.index]) <: Number);
                if isit
                    type = 1; # a constant wrapped in a coefficient ... not ideal
                    val = coefficients[i].value[c.index];
                else
                    type = 2; # a function
                    name = coefficients[i].value[c.index].name;
                    for j=1:length(genfunctions)
                        if name == genfunctions[j].name
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
        for i=1:length(variables)
            if c.name == string(variables[i].symbol)
                type = 3;
                val = variables[i].index;
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
    for i=1:length(coefficients)
        if c.name == string(coefficients[i].symbol)
            ind = coefficients[i].index;
        end
    end
    
    return ind;
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
                coef_part = :(a.*b);
                coef_part.args = [:.*];
                append!(coef_part.args, tmpcoef);
            elseif length(tmpcoef) == 1
                coef_part = tmpcoef[1];
            end
            
        else # It is an Expr, but not -() or * (note: at this point a/b is a*(1/b) )
            # This simplification may cause trouble somewhere. Revisit if needed.
            coef_part = ex;
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
    
    # println(ex)
    # println(test_part)
    # println(trial_part)
    # println(coef_part)
    # println("")
    
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