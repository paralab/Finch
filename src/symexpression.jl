#=
# The final symbolic form will be an Expr expression tree.
# The leaf nodes are all SymEntitys as described below.
# The parents will be math operators. 
# All used entities will also be kept in a list.
#
# A SymEntity has
# - name, the name of this variable, coefficient, etc. (or number for constants)
# - index, the component of this, such as a vector component (-1 for constants)
# - derivs
# - flags
=#
export build_symexpressions

# struct SymExpression
#     tree                # The expression tree
#     entities::Array     # The SymEntitys corresponding to the leaf nodes of the tree
# end

# struct SymEntity
#     name::Union{Float64, String}    # The string or number 
#     index::Union{Int64, Array}     # vector component, 1 for scalars, -1 for numbers, "INDEXEDBY..." for indexed vars
#     derivs::Array                   # any derivatives applied
#     flags::Array                    # any non-derivative modifiers or flags
# end

Base.copy(x::SymEntity) = SymEntity(x.name, x.index, copy(x.derivs), copy(x.flags));

# For printing, write the entity in a more readable format
Base.show(io::IO, x::SymEntity) = print(io, symentity_string(x));
function symentity_string(a::SymEntity)
    result = "";
    for i=1:length(a.derivs)
        result *= "D"*string(a.derivs[i]) * "(";
    end
    
    for i=1:length(a.flags)
        result *= a.flags[i];
    end
    
    result *= a.name;
    if typeof(a.index) == Int && a.index > 0
        result *= "_"*string(a.index);
    elseif typeof(a.index) <: Array && length(a.index) > 0
        if length(a.index) > 1
            result *= "_{"*a.index[1];
            for i=2:length(a.index)
                result *= ","*a.index[i];
            end
            result *= "}";
        elseif length(a.index) == 1
            result *= "_"*a.index[1];
        end
    end
    
    for i=1:length(a.derivs)
        result *= ")";
    end
    
    return result;
end

# Make a Latex string for the expression
# Recursively traverses the tree turning each piece into latex.
function symexpression_to_latex(ex)
    result = "";
    
    if typeof(ex) <: Array
        if length(ex) < 1
            return "";
        end
        
        if typeof(ex[1]) <: Array
            result = "\\left[";
            for i=1:length(ex)
                if i > 1
                    result *= ", ";
                end
                result *= symexpression_to_latex(ex[i]);
            end
            result *= "\\right]";
            
        else # The innermost array is an array of added terms
            result = symexpression_to_latex(ex[1]);
            for i=2:length(ex)
                result *= " + " * symexpression_to_latex(ex[i]);
            end
        end
        
    elseif typeof(ex) == SymExpression
        result = symexpression_to_latex(ex.tree);
        
    elseif typeof(ex) == Expr
        if ex.head === :call
            if ex.args[1] in [:+, :.+, :-, :.-, :*, :.*] && length(ex.args) > 2
                if typeof(ex.args[2]) == Expr
                    result *= "\\left(";
                end
                result *= symexpression_to_latex(ex.args[2]);
                if typeof(ex.args[2]) == Expr
                    result *= "\\right)";
                end
                for i=3:length(ex.args)
                    result *= " "*string(ex.args[1])[end] * " ";
                    if typeof(ex.args[i]) == Expr
                        result *= "\\left(";
                    end
                    result *= symexpression_to_latex(ex.args[i]);
                    if typeof(ex.args[i]) == Expr
                        result *= "\\right)";
                    end
                end
                
            elseif ex.args[1] in [:-, :.-] && length(ex.args) == 2 # negative
                result *= "-";
                if typeof(ex.args[2]) == Expr
                    result *= "\\left(";
                end
                result *= symexpression_to_latex(ex.args[2]);
                if typeof(ex.args[2]) == Expr
                    result *= "\\right)";
                end
                
            elseif ex.args[1] in [:/, :./]
                top = symexpression_to_latex(ex.args[2]);
                if length(ex.args) > 3
                    # multiply them first
                    newex = :(a*b);
                    newex.args = [:*];
                    append!(newex.args, ex.args[3:end]);
                    bottom = symexpression_to_latex(newex);
                else
                    bottom = symexpression_to_latex(ex.args[3]);
                end
                
                result = "\\frac{" * top * "}{" * bottom * "}";
                
            elseif ex.args[1] in [:^, :.^]
                bottom = symexpression_to_latex(ex.args[2]);
                power = symexpression_to_latex(ex.args[3]);
                result = bottom * "^{" * power * "}";
                
            elseif ex.args[1] in [:sqrt]
                result = "\\sqrt{" * symexpression_to_latex(ex.args[2]) * "}";
                
            elseif ex.args[1] in [:sin, :cos, :tan]
                result = "\\"*string(ex.args[1]) * "\\left(" * symexpression_to_latex(ex.args[2]) * "\\right)";
                
            elseif ex.args[1] === :abs
                result = "\\left|" * symexpression_to_latex(ex.args[2]) * "\\right|";
                
            elseif ex.args[1] === :conditional_function
                result =  "conditional\\left(" * symexpression_to_latex(ex.args[2]);
                for i=3:length(ex.args)
                    result *= ", " * symexpression_to_latex(ex.args[i]);
                end
                result *= "\\right)";
                
            else
                result = string(ex.args[1]) * "\\left(" * symexpression_to_latex(ex.args[2]);
                for i=3:length(ex.args)
                    result *= ", " * symexpression_to_latex(ex.args[i]);
                end
                result *= "\\right)";
            end
            
        elseif ex.head === :.
            # This is a broadcast symbol. args[1] is the operator, args[2] is a tuple for the operator
            # eg. change :(sin.(a)) to :(sin(a))
            newex = Expr(:call, ex.args[1]);
            append!(newex.args, ex.args[2].args);
            result = symexpression_to_latex(newex);
            
        else
            # There are some other possible symbols, but they could be ignored
            result = symexpression_to_latex(ex.args[1]);
        end
        
    elseif typeof(ex) == SymEntity
        # Apply derivatives
        dimchars = ["x","y","z"];
        name = ex.name;
        if name == "FACENORMAL1"
            name = "normal";
        elseif name == "FACENORMAL2"
            name = "-normal";
        end
        if typeof(ex.index) == Int && ex.index > 0
            result = name * "_{"*string(ex.index);
            if "DGSIDE1" in ex.flags || "CELL1" in ex.flags
                result *= "+";
            elseif "DGSIDE2" in ex.flags || "CELL2" in ex.flags
                result *= "-";
            end
            if length(ex.derivs) > 0
                result *= ","
                for i=1:length(ex.derivs)
                    result *= dimchars[ex.derivs[i]];
                end
            end
            result *= "}";
        else
            result = name;
        end
        
        
    elseif typeof(ex) <: Number
        result = string(ex);
        
    else
        result = string(ex);
    end
    
    return result;
end

# Builds all of them one at a time.
function build_symexpressions(var, expressions; remove_zeros=false)
    symexpressions = [];
    for e in expressions
        push!(symexpressions, build_symexpression(var, e, remove_zero_in_array=remove_zeros));
    end
    
    return symexpressions;
end

# Builds one symexpression from a SymEngine object or an Expr
function build_symexpression(var, ex; remove_zero_in_array=false)
    if ex === nothing
        return nothing
    end
    
    if typeof(ex) <: Array
        exarray = [];
        for i=1:length(ex)
            if !(remove_zero_in_array && ex[i] == 0)
                push!(exarray, build_symexpression(var, ex[i], remove_zero_in_array=remove_zero_in_array));
            end
        end
        return exarray;
        
    elseif typeof(ex) <: SymEngine.Basic
        # ex is a SymEngine.Basic object
        # First parse it into an Expr
        jex = Meta.parse(string(ex));
        
        # Then traverse it to find SymTerm leaves and build a SymExpression
        (newex, sents) = extract_symentities(var, jex);
        
        # Replace some special operators with their real expressions
        newex = replace_special_ops(newex);
        
        # If using FEM, reorder factors for easier generation
        if finch_state.config.solver_type == CG || finch_state.config.solver_type == DG
            newex = order_expression_for_fem(var, finch_state.test_functions, newex);
        end
        
        symex = SymExpression(newex, sents);
        
        return symex;
        
    elseif typeof(ex) == Expr || typeof(ex) == Symbol || typeof(ex) <: Number
        # Traverse it to find SymTerm leaves and build a SymExpression
        (newex, sents) = extract_symentities(var, ex);
        
        # Replace some special operators with their real expressions
        newex = replace_special_ops(newex);
        
        # If using FEM, reorder factors for easier generation
        if finch_state.config.solver_type == CG || finch_state.config.solver_type == DG
            newex = order_expression_for_fem(var, finch_state.test_functions, newex);
        end
        
        symex = SymExpression(newex, sents);
        
        return symex;
        
    end
end

# Traverse the tree and replace appropriate symbols with SymEntitys.
# Also, keep a list of symentities used.
function extract_symentities(var, ex)
    sents = [];
    if typeof(ex) == Expr
        # Recursively traverse the tree
        newex = copy(ex);
        if ex.head === :call # an operation like +(a,b)
            for i=2:length(ex.args)
                (subex, subsents) = extract_symentities(var, ex.args[i]);
                newex.args[i] = subex;
                append!(sents, subsents);
            end
        elseif ex.head === :if # a conditional like if a<b c else d end
            # This has three expr: one for a<b and a block for c and a block for d
            for i=1:length(ex.args)
                (subex, subsents) = extract_symentities(var, ex.args[i]);
                newex.args[i] = subex;
                append!(sents, subsents);
            end
            
        else # There are a few other things it could be. Just process the args.
            for i=1:length(ex.args)
                (subex, subsents) = extract_symentities(var, ex.args[i]);
                newex.args[i] = subex;
                append!(sents, subsents);
            end
        end
        
    elseif typeof(ex) == Symbol
        sen = build_symentity(var, ex);
        push!(sents, sen);
        newex = sen;
        
    else # probably a number
        newex = ex;
    end
    
    return (newex, sents);
end

# Builds a SymEntity struct from a symbolic name such as:
# DGSIDE1_D2__u_1 -> name=u, index=1, derivs=D2, flags=DGSIDE1
function build_symentity(var, ex)
    if typeof(ex) <: Number
        return SymEntity(ex, 0, [], []);
    elseif typeof(ex) == Symbol
        # Extract the index, derivs, flags
        (ind, symb, mods) = extract_entity_parts(ex);
        # Separate the derivative mods from the other flags
        (derivs, others) = get_derivative_mods(mods);
        
        return SymEntity(symb, ind, derivs, others);
        
    else
        printerr("unexpected type in build_symentity: "*string(typeof(ex))*" : "*string(ex));
    end
end

# Extract the indexing symbols from the index string
function extract_index_symbols_from_index_string(ind_str)
    indices = [];
    b = 10;
    e = 10;
    for i=10:length(ind_str)
        if ind_str[i-1:i] == "BY"
            e = i-2;
            push!(indices, ind_str[b:e]);
            b = i+1;
        elseif i == length(ind_str)
            push!(indices, ind_str[b:end]);
        end
    end
    
    return indices;
end

# Extract meaning from the symbolic object name
# The format of the input symbol should look like this
#   MOD1_MOD2_..._var_n
# There could be any number of mods, _var_n is the symvar symbol (n is the vector index, or 1 for scalar, or INDEXEDBY... for array type)
# Returns (n, var, [MOD1,MOD2,...]) all as strings
function extract_entity_parts(ex)
    str = string(ex);
    index = 0;
    var = "";
    mods = [];
    l = lastindex(str);
    e = l; # end of variable name
    b = l; # beginning of variable name
    
    # dt is a special symbol that will be passed as a number value in the generated function.
    if str == "dt"
        return(-1, str, []);
    end
    
    # Extract the index
    for i=l:-1:1
        if str[i] == '_' && i < l
            index_str = str[i+1:l];
            # This could be an integer or INDEXEDBY...
            if occursin("INDEXED", index_str)
                index = extract_index_symbols_from_index_string(index_str);
            else
                try
                    index = parse(Int, index_str);
                catch
                    return (-1,str,[]); # an unexpected index, perhaps it is a special symbol
                end
            end
            
            e = i-1; # end of variable name
            break;
        end
    end
    
    # Extract the name
    for i=e:-1:1
        if str[i] == '_'
            var = str[i+1:e];
            b = i-1;
            break;
        end
    end
    
    # extract the modifiers like D1_ etc. separated by underscores
    if b>1
        e = b-1;
        for i=e:-1:1
            if str[i] == '_'
                b = i+1;
                push!(mods, str[b:e]);
                e = b-2;
                
            elseif i == 1
                push!(mods, str[1:e]);
            end
        end
    end
    
    return (index, var, mods);
end

# Separate the derivative mods(Dn_) from the other mods and return an array 
# of derivative indices and the remaining mods.
function get_derivative_mods(mods)
    derivs = [];
    others = [];
    
    for i=1:length(mods)
        if length(mods[i]) == 2 && mods[i][1] == 'D'
            try
                index = parse(Int, mods[i][2])
                push!(derivs, index);
            catch
                printerr("Unexpected modifier: "*mods[i]*" see get_derivative_mods() in symexpression.jl")
                push!(others, mods[i]);
            end
        else
            push!(others, mods[i]);
        end
    end
    
    return (derivs, others);
end

# changes special operators into the intended expressions
function replace_special_ops(ex)
    if typeof(ex) == Expr
        newex = copy(ex);
        for i=1:length(ex.args)
            newex.args[i] = replace_special_ops(newex.args[i]);
        end
        
    elseif typeof(ex) == Symbol
        # check against these special ones that need to be replaced
        # This will eventually contain more complex things.
        specials = [:symbolmin, :symbolmax, :isgreaterthan, :islessthan, :conditional, 
                    :exp, :sin, :cos, :tan, :sinh, :cosh, :tanh, :abs];
        replacements = Dict([(:symbolmin, :min), (:symbolmax, :max), (:isgreaterthan, :>), 
                            (:islessthan, :<), (:conditional, :conditional_function), 
                            (:exp, :exp), (:sin, :sin), (:cos, :cos), (:tan, :tan), 
                            (:sinh, :sinh), (:cosh, :cosh), (:tanh, :tanh), (:abs, :abs)]);
        if ex in specials
            newex = replacements[ex];
        elseif occursin("CALLBACK_", string(ex))
            # It's a callback function. Figure out the index in callback_functions
            fun_name = string(ex)[10:end];
            for i=1:length(finch_state.callback_functions)
                if finch_state.callback_functions[i].name == fun_name
                    newex = :(Finch.finch_state.callback_functions[$i].func);
                end
            end
        else
            newex = ex;
        end
    else
        newex = ex;
    end
    
    return newex;
end

function conditional_function(a, b, c)
    if a
        return b;
    else
        return c;
    end
end

# Order the factors in each term so that the test function is first and the unknown(when present) is second
function order_expression_for_fem(var, test, ex)
    # Assuming the upper levels of the tree look like t1+t2-t3+...
    # also possible: t1+(t2-t3)+... or just t1
    # Handle recursively in case there is something like t1+(t2+(t3+(t4...)))
    if typeof(ex) == Expr
        newex = copy(ex);
        if newex.head === :call
            if newex.args[1] === :+ || newex.args[1] === :.+ || newex.args[1] === :- || newex.args[1] === :.-
                # It is a + or - operation. Recursively order each term
                for i=2:length(newex.args)
                    newex.args[i] = order_expression_for_fem(var, test, newex.args[i])
                end
            elseif newex.args[1] === :* || newex.args[1] === :.*
                # search for the test function
                test_part = 0;
                var_part = 0;
                for i=2:length(newex.args)
                    if has_specific_part(newex.args[i], test)
                        test_part = i;
                    elseif has_specific_part(newex.args[i], var)
                        var_part = i;
                    end
                end
                
                # Place test first, var second
                if test_part > 2
                    tmp = newex.args[2];
                    newex.args[2] = newex.args[test_part];
                    newex.args[test_part] = tmp;
                    if var_part == 2
                        var_part = test_part;
                    end
                    if var_part > 3
                        tmp = newex.args[3];
                        newex.args[3] = newex.args[var_part];
                        newex.args[var_part] = tmp;
                    end
                    
                elseif var_part > 2
                    tmp = newex.args[2];
                    newex.args[2] = newex.args[var_part];
                    newex.args[var_part] = tmp;
                end
            end
            
        end
        
        return newex;
        
    else # ex is not an Expr
        # It is a single symentity. Do nothing.
        return ex;
    end
    
end

# Returns true if a.name matches any of the parts
function has_specific_part(a, parts)
    if typeof(a) == SymEntity
        if typeof(parts) <: Array
            for i=1:length(parts)
                if has_specific_part(a, parts[i])
                    return true;
                end
            end
        elseif typeof(parts) <: Variable
            return a.name == string(parts.symbol);
        elseif typeof(parts) == Coefficient
            return a.name == string(parts.symbol);
        end
    end
    return false;
end