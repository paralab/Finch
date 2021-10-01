#=
# function utils
=#
export GenFunction, stringToFunction, add_genfunction, makeFunction, makeFunctions
# export @stringToFunction, @makeFunction, @makeFunctions

# Stores everything we could need about a generated function.
struct GenFunction      # example:
    name::String        # "genfunction_7"
    args::String        # "x"
    str::String         # "sin(x)"
    expr                # expression: sin(x) NOTE: could be an Expr, Symbol, Number,
    func                # handle for: genfunction_7(x) = sin(x)
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
    name = "genfunction_"*string(Finch.genfunc_count);
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
function makeFunctions(ex; args="x=0,y=0,z=0,t=0,node_index=1")
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

### Old macro versions below. keep for now.

# # Generates a function
# macro stringToFunction(name, args, fun)
#     return esc(quote
#         eval(Meta.parse($name*"("*$args*")="*$fun));
#     end)
# end

# # Makes a GenFunction and adds it to the
# # args is a string like "x,y,z"
# # fun is a string like "sin(x)*y + 3*z"
# macro makeFunction(args, fun)
#     return esc(quote
#         let
#             name = "genfunction_"*string(Finch.genfunc_count);
#             if typeof($fun) == Expr
#                 ex = $fun;
#                 strfun = string($fun);
#                 nf = GenFunction(name, $args, strfun, ex, @stringToFunction(name, $args, strfun));
#             else
#                 ex = Meta.parse($fun);
#                 nf = GenFunction(name, $args, $fun, ex, @stringToFunction(name, $args, $fun));
#             end
#             add_genfunction(nf);
#         end
#     end)
# end

# # Makes either: a constant number, a genfunction, or an array of genfunctions
# macro makeFunctions(ex)
#     return esc(quote
#         nfuns = 0;
#         args = "x=0,y=0,z=0,t=0";
#         if typeof($ex) <: Array
#             ### If it's an array, it can't be mixed type! At least for now.
#             # There's a problem with nfuns disappearing in the loop.
#             if typeof($ex[1]) <: Number
#                 nfuns = 0;
#             else
#                 for i=1:length($ex)
#                     if typeof($ex[i]) == String
#                         tmp = @makeFunction(args, $ex[i]);
#                     end
#                 end
#                 nfuns = length($ex)
#             end
#         else
#             if typeof($ex) == String
#                 @makeFunction(args, $ex);
#                 nfuns = 1;
#             elseif typeof($ex) <: Number
#                 nfuns = 0;
#             end
#         end
#         nfuns
#     end)
# end