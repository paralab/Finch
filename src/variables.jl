#=
# Variable struct and related functions
=#

# """
#     Variable

# Represents a value that is not defined by an independent function of coordinates.
# Note that this does not have to be an unknown that is solved for. It can hold
# values that are dependent on other variables such as transformed variables, or
# any other value that the user has control over.

# This is also useful for setting up data arrays that match those of unknowns.
# The values for a variable v can be directly accessed as `v.values` and has 
# dimensions [C, N] where C is the number of components (SCALAR=1, etc.) and N is 
# the number of nodes or cells.

# This should be built with the `variable` function.
# """
# mutable struct Variable
#     symbol::Symbol          # symbol used in expressions referring to this variable
#     symvar::Array{Basic,1}  # Array of symbolic layer variable symbols(Basic symbols)
#     index::Int              # index in the Finch list of variables
#     type::String            # constants for SCALAR, VECTOR, etc.
#     location::String        # constant for NODAL or CELL
#     discretization::String  # constant for CG, DG, FV
#     values::Array           # a C x N array, C is number of components (SCALAR=1, etc.)
    
#     # These are mainly for ARRAY type vars. Should there be a separate struct for those?
#     indexer::Union{Vector, Indexer} # indexer object or array of them
#     total_components::Int   # Number of components
#     #transposed::Bool        # If true, the values are actually trasposed -> N x C array (better for large C?)
    
#     # not currently used, may be removed
#     dependson::Vector       # a list of variables that this depends on
#     ready::Bool             # Is this variable's values ready? Can dependent variables use it?
# end

# """
#     VariableTransform

# Two sets of variables and a function that transforms one into the other.
# This is built with the `variableTransform` function.
# """
# struct VariableTransform
#     from::Union{Variable, Array{Variable}}  # Transformed from this
#     to::Union{Variable, Array{Variable}}    # to this
#     func::Function                          # using this function
#     # This function takes one number(or array matching from) 
#     # and returns one number(or an array matching to)
#     # NOTE: if "from" is an array, "to" must also be an array (size can be different)
# end

# For printing, write the variable symbol only
Base.show(io::IO, x::Variable) = print(io, string(x.symbol));

Base.:+(x::Finch.Variable, y::Finch.Variable) = x.values + y.values;
Base.:-(x::Finch.Variable, y::Finch.Variable) = x.values - y.values;
Base.:*(x::Finch.Variable, y::Finch.Variable) = x.values .* y.values;
Base.:/(x::Finch.Variable, y::Finch.Variable) = x.values ./ y.values;
Base.:^(x::Finch.Variable, y::Number) = x.values .^ y;

# Performs the transform on all of the variable values.
# Writes the results directly into "to".
function transform_variable_values(xform)
    is_array = false;
    if typeof(xform.from) <: Array
        from_len = length(xform.from);
        N = length(xform.from[1].values);
        is_array = true;
    else
        from_len = 1;
        N = length(xform.from.values);
    end
    if typeof(xform.to) <: Array
        to_len = length(xform.to);
    else
        to_len = 1;
    end
    
    if is_array
        datatype = typeof(xform.from[1].values[1]);
        input = zeros(datatype, from_len); # This is just the number of components, so shouldn't be too big
        output = zeros(datatype, to_len);
        for i=1:N
            for j=1:from_len
                input[j] = xform.from[j].values[i];
            end
            output = xform.func(input);
            for j=1:to_len
                xform.to[j].values[i] = output[j];
            end
        end
    else
        for i=1:N
            xform.to.values[i] = xform.func(xform.from.values[i]);
        end
    end
    return nothing;
end