#=
# Variable struct and related functions
=#

mutable struct Variable
    symbol::Symbol          # symbol used in expressions referring to this variable
    symvar::Array{Basic,1}  # Array of symbolic layer variable symbols(Basic symbols)
    index::Int              # index in the Finch list of variables
    type::String            # constants for SCALAR, VECTOR, etc.
    location::String        # constant for NODAL or CELL
    discretization::String  # constant for CG, DG, FV
    values::Array{Float64}  # a C x N array, C is number of components (SCALAR=1, etc.)
    
    # These are mainly for ARRAY type vars. Should there be a separate struct for those?
    indexer                  # indexer object or array of them
    total_components::Int   # Number of components
    #transposed::Bool        # If true, the values are actually trasposed -> N x C array (better for large C?)
    
    # not currently used, may be removed
    dependson               # a list of variables that this depends on
    ready::Bool             # Is this variable's values ready? Can dependent variables use it?
end

# For printing, write the variable symbol only
Base.show(io::IO, x::Variable) = print(io, string(x.symbol));

Base.:+(x::Finch.Variable, y::Finch.Variable) = x.values + y.values;
Base.:-(x::Finch.Variable, y::Finch.Variable) = x.values - y.values;
Base.:*(x::Finch.Variable, y::Finch.Variable) = x.values .* y.values;
Base.:/(x::Finch.Variable, y::Finch.Variable) = x.values ./ y.values;
Base.:^(x::Finch.Variable, y::Number) = x.values .^ y;

struct VariableTransform
    from::Union{Variable, Array{Variable}}  # Transformed from this
    to::Union{Variable, Array{Variable}}    # to this
    func                                    # using this function
    # This function takes one number(or array matching from) 
    # and returns one number(or and array matching to)
    # NOTE: if "from" is an array, "to" must also be an array (size can be different)
end

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
        input = zeros(from_len);
        output = zeros(to_len);
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