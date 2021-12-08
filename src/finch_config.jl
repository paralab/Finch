#=
# A struct containing configuration information
# This is the most general Finch info that should be (almost)problem indepenent.
# Premade configs can be loaded for different classes of problems.
# Values can be set through various macros or automatically by problem specification.
=#
if !@isdefined(IRREGULAR)
    include("finch_constants.jl");
end

mutable struct Finch_config
    # Domain
    dimension::Int          # 1,2,3
    geometry::String        # square, irregular
    mesh_type::String       # unstructured, tree, uniform grid

    # FEM details
    solver_type::String     # cg, dg, fv
    trial_function::String  # Legendre
    test_function::String   # same as above
    elemental_nodes::String # uniform, gauss, lobatto (higher order node distribution within elements)
    quadrature::String      # uniform, gauss, lobatto (similar to above)
    p_adaptive::Bool        # Do adaptive p-refinement?
    basis_order_min::Int    # minimum order to use in p-refinement, or if p_adaptive is false
    basis_order_max::Int    # maximum order
    
    # Other solver details
    linear::Bool            # Is the equation linear?
    t_adaptive::Bool        # Do adaptive t_refinement?
    stepper::String         # Euler-explicit/implicit, RK4, LSRK4, etc. Type of time stepper to use
    linalg_matrixfree::Bool # Use matrix free methods?
    linalg_matfree_max::Int # max iters for matrix free
    linalg_matfree_tol::Float64 # tolerance for matrix free
    linalg_backend::String  # default, petsc, ?? (What to use for linear algebra)
    
    # Output
    output_format::String   # VTK, raw, custom (format for storing solutions)
    
    # Parallel details
    use_mpi::Bool           # Is MPI available?
    num_procs::Int;         # number of processes
    proc_rank::Int;         # this proccess rank
    num_threads::Int;       # number of available threads
    num_partitions::Int;    # number of mesh partitions
    partition_index::Int;   # this process's partition
    
    # Constructor builds a default config.
    Finch_config() = new(
        1,
        SQUARE,
        UNIFORM_GRID,
        CG,
        LEGENDRE,
        LEGENDRE,
        LOBATTO,
        GAUSS,
        false,
        1,
        1,
        true,
        false,
        EULER_IMPLICIT,
        false,
        1,
        1.0,
        DEFAULT_SOLVER,
        VTK,
        false,
        1,
        0,
        1,
        1,
        0
    );
end
