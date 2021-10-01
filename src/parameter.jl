#=
# A parameter to be used in the weak form.
# It is an array of Expr.
# Each represents a function of x,y,z,t, coefficients, or known variables
# Values are in an array to allow vector/tensor valued parameters
=#

struct Parameter
    symbol::Symbol          # symbol used in expressions
    index::Int              # index in the Finch list of parameters
    type::String            # constants for SCALAR, VECTOR, etc.
    value::Array            # An array of either constant values(numbers) or Expr expressions
end