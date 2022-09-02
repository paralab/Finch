#=
A struct containing problem information
=#
if !@isdefined(IRREGULAR)
    include("finch_constants.jl");
end

mutable struct Finch_prob
    # Boundary condition info
    bc_type::Array{String,2}        # DIRICHLET, etc. for each variable and bid
    bid::Array{Int,2}               # BID for each variable and boundary section
    bc_func::Array{Any,2}           # GenFunctions or numbers for each variable and bid
    
    # Reference points for pinning solutions without Dirichlet conditions
    ref_point::Array{Any,2}         # For each variable, [bool has_ref_point, [x,y,z], [vals]]
    
    # Time dependence info
    time_dependent::Bool            # Is this problem time dependent
    end_time::Float64               # If so, what is the final time
    initial::Array{Any,1}           # An array of initial condition GenFunctions, one for each variable, vector vars can have an array of functions
    
    # Nonlinear iteration info
    nonlinear::Bool                 # Is there a nonlinear iteration
    max_iters::Int                  # Max iterations
    relative_tol::Float64           # Relative tolerance for change between iterations (u(i) - u(i-1))
    absolute_tol::Float64           # Absolute tolerance
    
    # Constructor builds an empty prob.
    Finch_prob() = new(
        Array{String,2}(undef,(0,0)), 
        Array{Int,2}(undef,(0,0)), 
        Array{Any,2}(undef,(0,0)),
        Array{Any,2}(undef,(0,0)),
        false, 
        0, 
        Array{Any,1}(undef,(0)),
        false,
        1,
        1.0,
        1.0
    );
end
