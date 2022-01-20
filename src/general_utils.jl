#=
A place for general utilities that don't fit in elsewhere.

=#

# Traverses an expression, replacing variable symbols with the value arrays
function replace_var_symbols_with_values(ex)
    if typeof(ex) == Expr
        if ex.head === :call # function symbols (like :sin in sin(x)) don't need to be processed
            for i=2:length(ex.args)
                ex.args[i] = replace_var_symbols_with_values(ex.args[i]);
            end
        else
            for i=1:length(ex.args)
                ex.args[i] = replace_var_symbols_with_values(ex.args[i]);
            end
        end
        
    elseif typeof(ex) == Symbol
        # Check if it is a variable and substitute if needed
        var_index = 0;
        for v in variables
            if ex === v.symbol
                var_index = v.index;
            end
        end
        if var_index > 0
            if variables[var_index].type == SCALAR
                newex = :(Finch.variables[$var_index].values[node_index]);
            else
                newex = :(Finch.variables[$var_index].values[:,node_index]);
            end
            ex = newex;
        end
    end
    
    return ex;
end
