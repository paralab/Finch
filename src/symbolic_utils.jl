#=
Utilities for manipulating symbolic expressions that may also
be useful outside of SymbolicParser.
=#
export apply_negative, apply_flag_to_all_symbols, apply_old_to_symbol, apply_delta_to_symbol,
        get_root_symbol, apply_delta_to_symbols_in_expression, create_deriv_func

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

function get_root_symbol(s)
    str = string(s);
    # Find _u_ in FLAG_FLAG__u_n
    # Does it have any flags?
    st = 1;
    if str[1] == '_' # no flags
        st = 1;
        
    else # skip until "__"
        for i=2:length(str)
            tmp = st;
            if str[i-1] == '_' && str[i] == '_'
                st = i;
                break;
            end
        end
    end
    # find the end index
    en = st+1;
    for i=(st+1):length(str)
        tmp = st;
        if str[i] == '_'
            en = i;
            break;
        end
    end
    return Symbol(str[st:en]); # should be something like _u_
end

# Adds "OLD" to the variable symbol
function apply_old_to_symbol(ex)
    str = string(ex);
    # Find u in FLAG_FLAG__u_n
    # Does it have any flags?
    st = 1;
    if str[1] == '_' # no flags
        st = 2;
        
    else # skip until "__"
        for i=2:length(str)
            tmp = st;
            if str[i-1] == '_' && str[i] == '_'
                st = i+1;
                break;
            end
        end
    end
    
    if st > 1
        newstr = str[1:(st-1)] * "OLD" * str[st:end];
    else
        newstr = str;
        printerr("Applying OLD to a misformed symbol? " * str)
    end
    
    return symbols(newstr);
end
# Adds "DELTA" to the variable symbol
function apply_delta_to_symbol(ex)
    str = string(ex);
    # Find u in FLAG_FLAG__u_n
    # Does it have any flags?
    st = 1;
    if str[1] == '_' # no flags
        st = 2;
        
    else # skip until "__"
        for i=2:length(str)
            tmp = st;
            if str[i-1] == '_' && str[i] == '_'
                st = i+1;
                break;
            end
        end
    end
    
    if st > 1
        newstr = str[1:(st-1)] * "DELTA" * str[st:end];
    else
        newstr = str;
        printerr("Applying DELTA to a misformed symbol? " * str)
    end
    
    return symbols(newstr);
end
function apply_delta_to_symbols_in_expression(ex, var_symbol)
    root_symbol = get_root_symbol(var_symbol);
    
    new_ex = ex;
    if typeof(ex) == Basic
        # Extract a list of all symbols from ex
        allsymbols = SymEngine.free_symbols(ex);
        new_ex = copy(ex);
        
        # look for symbols containing the root symbol
        for s in allsymbols
            if get_root_symbol(s) === root_symbol
                new_ex = subs(new_ex, (s, apply_delta_to_symbol(s)))
            end
        end
    end
    return new_ex;
end

# Given a Basic term and variable, generate a function for the term.
# If method is "AD": 
#   Use AD to generate a derivative with respect to var.
#   Use this derivative to create a callback function.
#   First strip off the test function factor, add at end.
# If method is "symbolic":
#   Use symengine to make symbolic derivative with respect to var.
# 
function create_deriv_func(term, var, method)
    if method == "AD"
        # Figure out the present test function symbol
        allsymbols = SymEngine.free_symbols(term);
        test_symbol = nothing;
        for c in finch_state.test_functions
            for s in allsymbols
                if occursin("_"*string(c.symbol)*"_", string(s)) && test_symbol === nothing
                    test_symbol = s;
                    break;
                end
            end
        end
        
        # If present, trim off test function factor by dividing
        if !(test_symbol === nothing)
            new_term = term / test_symbol;
        else
            new_term = term;
        end
        
        # create a string for the needed args
        allsymbols = SymEngine.free_symbols(new_term);
        first_arg = string(var);
        other_args = "";
        arg_list = [string(var)];
        for s in allsymbols
            if !(s == var)
                other_args *= ", " * string(s);
                push!(arg_list, string(s));
            end
        end
        arg_str = first_arg * other_args;
        
        # Create a string for the function
        func_str = "("*arg_str*")->("*string(new_term)*")";
        
        # Make a derivative
        dfunc_str = "function df("*arg_str*") return Zygote.gradient("*func_str*", "*arg_str*")[1]; end"
        dfunc = eval(Meta.parse(dfunc_str));
        
        # Put it in a callback function
        i = length(finch_state.callback_functions);
        dfunc_name = "ADFUNCTION"*string(i+1);
        # callbackFunction(dfunc, name=dfunc_name);
        push!(finch_state.callback_functions, CallbackFunction(dfunc_name, arg_list, "", dfunc))
        log_entry("Added callback function for AD: "*dfunc_name, 2);
        
        # This will be placed in the term like ADFUNCTIONi(args)
        str = "CALLBACK_"*dfunc_name*"("*arg_str*")";
        
        # Create a symengine fun
        sfun = SymEngine.SymFunction("CALLBACK_"*dfunc_name);
        
        # If there was a test function, put it back now
        if !(test_symbol === nothing)
            new_term = Basic(str) * test_symbol;
        else
            new_term = Basic(str);
        end
        
    else # symbolic
        new_term = SymEngine.diff(term, var);
    end
    
    return new_term;
end