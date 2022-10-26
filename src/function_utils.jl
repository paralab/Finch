#=
# function utils
=#
export GenFunction, CallbackFunction, stringToFunction, add_genfunction, makeFunction, makeFunctions, makeCompleteFunction
# export @stringToFunction, @makeFunction, @makeFunctions

# Stores everything we could need about a generated function.
struct GenFunction      # example:
    name::String        # "genfunction_7"
    args::String        # "x"
    str::String         # "sin(x)"
    expr                # expression: sin(x) NOTE: could be an Expr, Symbol, Number,
    func                # handle for: genfunction_7(x) = sin(x)
end

struct CallbackFunction
    name::String        # "myFunction"
    args::Array{String} # ["x", "u", "c"] list of strings for args should match defined symbols
    body::String        # String of the function body: "a=u+2; return a;" 
    func                # The function
end

# Generates a function from a string
function stringToFunction(name, args, fun)
    return eval(Meta.parse(name*"("*args*")="*fun));
end

# Adds the GenFunction to a global array of generated functions
function add_genfunction(genfun)
    global genfunc_count += 1;
    global genfunctions = [genfunctions ; genfun];
    log_entry("Generated function: "*genfun.name);
end

# Makes either: a constant number, a genfunction, or an array of genfunctions
function makeFunctions(ex; args="x=0,y=0,z=0,t=0,node_index=1,face_index=1;indices=nothing")
    nfuns = 0;
    if typeof(ex) <: Array
        for i=1:length(ex)
            if typeof(ex[i]) == String
                tmp = makeFunction(args, ex[i]);
                nfuns += 1;
            # else # It could be a number or function handle
            end
        end
    else
        if typeof(ex) == String
            makeFunction(args, ex);
            nfuns = 1;
        # else # It could be a number or function handle
        end
    end
    return nfuns;
end

# Makes a GenFunction and adds it to the
# args is a string like "x,y,z"
# fun is a string like "sin(x)*y + 3*z" OR an Expr
function makeFunction(args, fun)
    name = "genfunction_"*string(Finch.genfunc_count+1);
    if typeof(fun) == Expr
        ex = fun;
        strfun = string(fun);
        nf = GenFunction(name, args, strfun, ex, stringToFunction(name, args, strfun));
    else
        ex = Meta.parse(fun);
        nf = GenFunction(name, args, fun, ex, stringToFunction(name, args, fun));
    end
    add_genfunction(nf);
end

# Takes a string or Expr for a complete function
# Extracts the name and args
# Creates a Genfunction and adds it to the genfunction array
function makeCompleteFunction(fun)
    if typeof(fun) == Expr
        ex = fun;
        strfun = string(fun);
        (name, args) = extract_function_name_and_args(fun);
        nf = GenFunction(name, args, strfun, ex, eval(fun));
    else
        printerr("makeCompleteFunction(fun) expects fun to be an Expr. Instead got\n  " * 
            string(stypeof(fun)) * ": " * string(fun), fatal=true);
    end
    add_genfunction(nf);
end

# Takes an Expr representing a function and returns strings for the 
# function name and arguments.
function extract_function_name_and_args(fun)
    # Find the piece of the Expr that contains the function def
    if fun.head === :function
        fun_args = fun.args[1].args;
    elseif fun.head === :block
        # the function was put into a block and may have a LineNumberNode
        for i=1:length(fun.args)
            if typeof(fun.args[i]) == Expr && fun.args[i].head === :function
                fun_args = fun.args[i].args[1].args;
                break;
            end
        end
    end
    # The name is a symbol in fun.args[1].args[1]
    name = string(fun_args[1]);
    # The args are in fun.args[1].args[2:end]
    nargs = length(fun_args) - 1; # This groups all keyword args into one item
    nkwargs = 0;
    args = "";
    kwargs = "";
    for i = 1:nargs
        ai = fun_args[1+i];
        if typeof(ai) == Symbol # like (a, b, c)
            if length(args) > 0
                args *= ", ";
            end
            args *= string(ai);
            
        elseif typeof(ai) == Expr 
            if ai.head === :kw # like (a=1, b=2, c=3)
                if length(args) > 0
                    args *= ", ";
                end
                args *= string(ai.args[1]) * "=" * string(ai.args[2]);
                
            elseif ai.head === :parameters # keyword args like (; kw1=4, kw2=5, kw3...)
                nkwargs = length(ai.args);
                for j = 1:nkwargs
                    kwi = ai.args[j];
                    if typeof(kwi) == Expr
                        if kwi.head === :kw
                            if length(kwargs) > 0
                                kwargs *= ", ";
                            end
                            kwargs *= string(kwi.args[1]) * "=" * string(kwi.args[2]);
                        elseif kwi.head === :(...)
                            if length(kwargs) > 0
                                kwargs *= ", ";
                            end
                            kwargs *= string(kwi.args[1]) * "...";
                        end
                    end
                end # kwargs
            end
        end
    end # args
    
    args = args * "; " * kwargs;
    
    return (name, args);
end

# Checks a string representing a function for time dependence.
# Assumes that a "t" that is not part of a word is time.
# ex could be a string or array containing strings or the result will be false.
# If any piece is true, all is true.
function check_time_dependence(ex)
    result = false;
    if typeof(ex) <: Array
        for i=1:length(ex)
            result = result || check_time_dependence(ex[i]);
        end
    elseif typeof(ex) == String
        # For coefficients, only x,y,z,t,node_index,face_index are allowed variables, so any t is time
        if !(findfirst(isequal('t'), ex) === nothing)
            result = true;
        end
    end
    return result;
end

# Interpret symbols in the arguments to a function for boundary or initial conditions.
# Replace them with something meaningful.
# Symbols could be:
# numbers/math ops -> no change
# x,y,z,t -> no change (provided as input arguments)
# variable -> Finch.variables[1].values[node_index]
# coefficient -> evaluate_coefficient(Finch.coefficients[1], x,y,z,t)
# indexer -> Finch.indexers[1].value
# parameter -> TODO replaced with its expression
# normal -> Finch.grid_data.faceNormals[:,face_index]
# callbackFunction -> callback_functions[1].func
function replace_symbols_in_conditions(ex)
    if typeof(ex) <: Array
        for i=1:length(ex)
            ex[i] = replace_symbols_in_conditions(ex[i]);
        end
    elseif typeof(ex) == Expr
        for i=1:length(ex.args)
            ex.args[i] = replace_symbols_in_conditions(ex.args[i]);
        end
        
    elseif typeof(ex) == Symbol
        # Check against the possibilities and replace as needed
        if !(ex in [:x, :y, :z, :t, :node_index, :face_index]) # don't change these
            v_index = variable_index_from_symbol(ex);
            c_index = coefficient_index_from_symbol(ex);
            i_index = indexer_index_from_symbol(ex);
            cb_index = callback_index_from_symbol(ex);
            if v_index > 0
                if variables[v_index].type == SCALAR
                    newex = :(Finch.variables[$v_index].values[node_index]);
                else
                    newex = :(Finch.variables[$v_index].values[:,node_index]);
                end
                ex = newex;
            elseif c_index > 0
                if coefficients[c_index].type == SCALAR
                    newex = :(evaluate_coefficient(Finch.coefficients[$c_index], 1, x, y, z, t, node_index, face_index));
                else
                    num_comps = length(Finch.coefficients[c_index].value);
                    newex = :([evaluate_coefficient(Finch.coefficients[$c_index], comp, x, y, z, t, node_index, face_index) for comp in 1:$num_comps]);
                end
                ex = newex;
            elseif i_index > 0
                index_tag = indexers[i_index].tag;
                ex = :(indices[$index_tag]);
            elseif cb_index > 0
                ex = :(Finch.callback_functions[$cb_index].func);
            elseif ex === :normal
                ex = :(Finch.grid_data.facenormals[:,face_index]);
            end
            
        end
    else
        # other types will not be changed
    end
    return ex;
end

function variable_index_from_symbol(s)
    for v in variables
        if s === v.symbol
            return v.index;
        end
    end
    return 0;
end

function coefficient_index_from_symbol(s)
    for c in coefficients
        if s === c.symbol
            return c.index;
        end
    end
    return 0;
end

function indexer_index_from_symbol(s)
    for i=1:length(indexers)
        if s === indexers[i].symbol
            return i;
        end
    end
    return 0;
end

function callback_index_from_symbol(s)
    for i=1:length(callback_functions)
        if string(s) == callback_functions[i].name
            return i;
        end
    end
    return 0;
end