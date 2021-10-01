#=
# Variable struct and related functions
=#

mutable struct Variable
    symbol::Symbol          # symbol used in expressions referring to this variable
    symvar::Array{Basic,1}  # Array of symbolic layer variable symbols(Basic symbols)
    index::Int              # index in the Finch list of variables
    type::String            # constants for SCALAR, VECTOR, etc.
    location::String        # constant for NODAL or CELL
    values::Array{Float64}  # a C x N array, C is number of components (SCALAR=1, etc.)
    
    # These are mainly for ARRAY type vars. Should there be a separate struct for those?
    indexer                  # indexer object or array of them
    total_components::Int   # Number of components
    #transposed::Bool        # If true, the values are actually trasposed -> N x C array (better for large C?)
    
    # not currently used, may be removed
    dependson               # a list of variables that this depends on
    ready::Bool             # Is this variable's values ready? Can dependent variables use it?
end
