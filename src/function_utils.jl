#=
# function utils
=#
export GenFunction, CallbackFunction, stringToFunction, add_genfunction, makeFunction, makeFunctions
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

# Makes a GenFunction and adds it to the
# args is a string like "x,y,z"
# fun is a string like "sin(x)*y + 3*z"
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

# Makes either: a constant number, a genfunction, or an array of genfunctions
function makeFunctions(ex; args="x=0,y=0,z=0,t=0,node_index=1,face_index=1")
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

# Interpret symbols in the arguments to a function for boundary or initial conditions.
# Replace them with something meaningful.
# Everything must be in terms of x, y, z, t, node_index, face_index
# Symbols could be:
# numbers/math ops -> no change
# x,y,z,t -> no change (provided as input arguments)
# variable -> Finch.variables[1].values[node_index]
# coefficient -> evaluate_coefficient(Finch.coefficients[1], x,y,z,t)
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
                    newex = :(evaluate_coefficient(Finch.coefficients[$c_index], 1, [x,y,z], t, node_index, face_index));
                else
                    num_comps = length(Finch.coefficients[c_index].value);
                    newex = :([evaluate_coefficient(Finch.coefficients[$c_index], comp, [x,y,z], t, node_index, face_index) for comp in 1:$num_comps]);
                end
                ex = newex;
            elseif cb_index > 0
                ex = :(callback_functions[$cb_index].func);
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

function callback_index_from_symbol(s)
    for i=1:length(callback_functions)
        if string(s) == callback_functions[i].name
            return i;
        end
    end
    return 0;
end