#=
# A coefficient to be used in the weak form.
# Can be a number or a generated function.
# Values are in an array to allow vector/tensor valued coefficients
=#

struct Coefficient
    symbol::Symbol          # symbol used in expressions referring to this coefficient
    symvar::Array{Basic,1}  # Array of symbolic layer symbols(Basic symbols)
    index::Int              # index in the Finch list of coefficients
    type::String            # constants for SCALAR, VECTOR, etc.
    location::String        # constant for NODAL, MODAL, CELL
    value::Array            # An array of either constant values(numbers) or genfunctions
    
    is_element_array::Bool  # Are the values specified for each element in an array?
end

# Returns the value of the coefficient for the specified component, coordinates and time.
# There are optional node index and face index inputs
function evaluate_coefficient(c::Coefficient, comp, x, t, nodeind=1, faceind=1)
    if typeof(comp) <: Array
        # This is likely an array type coefficient: comp=[a,b,c] -> value[a,b,c]
        coef_size = size(c.value);
        comp_index = comp[1];
        for i=2:length(comp)
            comp_index += comp_size[i-1]*(comp[i]-1);
        end
        
    else
        comp_index = comp;
    end
    if typeof(c.value[comp_index]) <: Number
        return c.value[comp_index];
    end
    
    # If not a number, it should be a genfunction
    dim = length(x);
    if dim == 1
        return c.value[comp_index].func(x[1],0,0,t,nodeind, faceind);
    elseif dim == 2
        return c.value[comp_index].func(x[1],x[2],0,t,nodeind, faceind);
    else
        return c.value[comp_index].func(x[1],x[2],x[3],t,nodeind, faceind);
    end
end