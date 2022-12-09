#=
This file contains most of the structs used by Finch.
=#
export FinchState, FinchConfig, FinchProblem, 
        Variable, VariableTransform, Coefficient, Indexer, Parameter, 
        GenFunction, CallbackFunction, FVInfo, MeshData, Grid, Refel,
        GeometricFactors, Jacobian, ParentMaps, Stepper, ParallelBuffers,
        SymExpression, SymEntity, SymOperator, IR_part
#

#=
# A struct containing configuration information
# This is the most general Finch info that should be (almost)problem indepenent.
# Premade configs can be loaded for different classes of problems.
# Values can be set through various macros or automatically by problem specification.
=#
mutable struct FinchConfig
    # Domain
    dimension::Int          # 1,2,3
    geometry::String        # square, irregular
    mesh_type::String       # unstructured, tree, uniform grid
    
    # Discretization details
    solver_type::String     # cg, dg, fv
    
    # FEM
    trial_function::String  # Legendre
    test_function::String   # same as above
    elemental_nodes::String # uniform, gauss, lobatto (higher order node distribution within elements)
    quadrature::String      # uniform, gauss, lobatto (similar to above)
    p_adaptive::Bool        # Do adaptive p-refinement?
    basis_order_min::Int    # minimum order to use in p-refinement, or if p_adaptive is false
    basis_order_max::Int    # maximum order
    
    # FVM
    fv_order::Int           # Order of reconstruction at faces
    
    # Time stepping
    t_adaptive::Bool        # Do adaptive t_refinement?
    stepper::String         # Euler-explicit/implicit, RK4, LSRK4, etc. Type of time stepper to use
    
    # Other solver details
    linear::Bool            # Is the equation linear?
    linalg_matrixfree::Bool             # Use matrix free methods?
    linalg_iterative::Bool              # Use an iterative solver?
    linalg_iterative_method::String     # GMRES or CG
    linalg_iterative_pc::String         # AMG, ILU, NONE
    linalg_iterative_maxiter::Int       # max iters for iterative solver
    linalg_iterative_abstol::Float64    # absolute tolerance
    linalg_iterative_reltol::Float64    # relative tolerance
    linalg_iterative_gmresRestart::Int  # GMRES restart iterations
    linalg_iterative_verbose::Bool      # print convergence info?
    
    linalg_usePetsc::Bool   # use PETSc?
    
    # Output
    output_format::String   # VTK, raw, custom (format for storing solutions)
    
    # Parallel details
    use_mpi::Bool           # Is MPI available?
    num_procs::Int;         # number of processes
    proc_rank::Int;         # this proccess rank
    num_threads::Int;       # number of available threads
    num_partitions::Int;    # number of mesh partitions
    partition_index::Int;   # this process's partition
    
    # Data types
    index_type::Type
    float_type::Type
    
    # Cachesim
    use_cachesim::Bool
    
    # Constructor builds a default config.
    FinchConfig() = new(
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
        
        1,
        
        false,
        EULER_IMPLICIT,
        
        true,
        false,
        true,
        "GMRES",
        "ILU",
        0,
        0,
        1e-8,
        0,
        false,
        false,
        
        VTK,
        
        false,
        1,
        0,
        1,
        1,
        0,
        
        Int64,
        Float64,
        
        false
    );
end

# Stores everything we could need about a generated function.
struct GenFunction      # example:
    name::String        # "genfunction_7"
    args::String        # "x"
    str::String         # "sin(x)"
    expr::Union{Expr,Symbol,Number} # expression: sin(x) NOTE: could be an Expr, Symbol, Number,
    func::Function      # handle for: genfunction_7(x) = sin(x)
end

#=
A struct containing problem information
=#
mutable struct FinchProblem
    # Boundary condition info
    bc_type::Matrix{String}        # DIRICHLET, etc. for each variable and bid
    bid::Vector{Int}               # BID for each boundary section
    bid_def::Vector{String}        # Description of boundary region if possible (ALL, XMIN, XMAX, (x-1)<eps, CUSTOM, ...)
    bc_func::Matrix{Vector{Union{Float64,GenFunction}}} # GenFunctions or numbers for each variable and bid
    
    # Reference points for pinning solutions without Dirichlet conditions
    has_ref_point::Vector{Bool}
    ref_point::Matrix{Vector{Union{Int,Float64,GenFunction}}} # For each variable, [bool has_ref_point, [x,y,z], [vals]]
    
    # Time dependence info
    time_dependent::Bool            # Is this problem time dependent
    end_time::Float64               # If so, what is the final time
    initial::Vector{Vector{Union{Float64,GenFunction}}} # An array of initial condition GenFunctions or values, one for each variable, vector vars can have an array of functions
    
    # pre and post-step functions
    pre_step_function               # user provided function to be called before
    post_step_function              # or after eash time step (initially nothing)
    
    # Nonlinear iteration info
    nonlinear::Bool                 # Is there a nonlinear iteration
    derivative_type::String         # "AD" or "symbolic"
    max_iters::Int                  # Max iterations
    relative_tol::Float64           # Relative tolerance for change between iterations (u(i) - u(i-1))
    absolute_tol::Float64           # Absolute tolerance
    relaxation::Float64             # Relaxation parameter
    
    # Constructor builds an empty prob.
    FinchProblem() = new(
        Matrix{String}(undef,(0,0)), 
        Vector{Int}(undef,(0)), 
        Vector{String}(undef,(0)), 
        Matrix{Vector}(undef,(0,0)),
        Vector{Bool}(undef,0),
        Matrix{Vector}(undef,(0,0)),
        false, 
        0, 
        zeros(0),
        nothing,
        nothing,
        false,
        "AD",
        1,
        1.0,
        1.0,
        1.0
    );
end

"""
    Indexer

An entity representing an index to be applied to a variable or coefficient. It has a 
symbol that can be used in expressions, a range of integer values, and a current 
value that can be accessed with `I.value` for an index labeled I.

This should be built with the `index` function.
"""
mutable struct Indexer
    symbol::Symbol   # symbol used in expressions
    range::Vector{Int}     # Range of integers.
    value::Int       # Current value
    tag::Int         # This indexer's position in the array of indexers?
end

"""
    Variable{T<:AbstractFloat}

Represents a value that is not defined by an independent function of coordinates.
Note that this does not have to be an unknown that is solved for. It can hold
values that are dependent on other variables such as transformed variables, or
any other value that the user has control over.

This is also useful for setting up data arrays that match those of unknowns.
The values for a variable v can be directly accessed as `v.values` and has 
dimensions [C, N] where C is the number of components (SCALAR=1, etc.) and N is 
the number of nodes or cells.

This should be built with the `variable` function.
"""
mutable struct Variable{T<:AbstractFloat}
    symbol::Symbol          # symbol used in expressions referring to this variable
    symvar::Vector{Basic}   # Array of symbolic layer variable symbols(Basic symbols)
    index::Int              # index in the Finch list of variables
    type::String            # constants for SCALAR, VECTOR, etc.
    location::String        # constant for NODAL or CELL
    discretization::String  # constant for CG, DG, FV
    values::Matrix{T}       # a C x N array, C is number of components (SCALAR=1, etc.)
    
    # These are mainly for ARRAY type vars. Should there be a separate struct for those?
    indexer::Vector{Indexer} # indexer object or array of them
    total_components::Int   # Number of components
    #transposed::Bool        # If true, the values are actually trasposed -> N x C array (better for large C?)
    
    # not currently used, may be removed
    dependson::Vector{Variable{T}} # a list of variables that this depends on
    ready::Bool             # Is this variable's values ready? Can dependent variables use it?
end

"""
    VariableTransform

Two sets of variables and a function that transforms one into the other.
This is built with the `variableTransform` function.
"""
struct VariableTransform{T<:AbstractFloat}
    from::Union{Variable{T}, Array{Variable{T}}}  # Transformed from this
    to::Union{Variable{T}, Array{Variable{T}}}    # to this
    func::Function                          # using this function
    # This function takes one number(or array matching from) 
    # and returns one number(or an array matching to)
    # NOTE: if "from" is an array, "to" must also be an array (size can be different)
end

"""
    Coefficient{T<:AbstractFloat}

Represents an independent value that can be defined in terms of coordinates
or assigned numerical values. Unlike variables, their values are typically not
available to the user since they may be in the form of generated functions 
rather than numbers.

This should be built with the `coefficient` function.
"""
struct Coefficient
    symbol::Symbol          # symbol used in expressions referring to this coefficient
    symvar::Vector{Basic}  # Array of symbolic layer symbols(Basic symbols)
    index::Int              # index in the Finch list of coefficients
    type::String            # constants for SCALAR, VECTOR, etc.
    location::String        # constant for NODAL, MODAL, CELL
    value::Vector{Union{Float64, GenFunction}} # An array of either constant values(numbers) or genfunctions
    
    is_element_array::Bool  # Are the values specified for each element in an array?
    is_time_dependent::Bool # Is it time dependent?
end

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
    value::Vector{Union{Float64, Expr}} # An array of either constant values(numbers) or Expr expressions
end

struct CallbackFunction
    name::String        # "myFunction"
    args::Array{String} # ["x", "u", "c"] list of strings for args should match defined symbols
    body::String        # String of the function body: "a=u+2; return a;" 
    func::Function      # The function
end

struct SymExpression
    tree                # The expression tree
    entities::Array     # The SymEntitys corresponding to the leaf nodes of the tree
end

struct SymEntity
    name::Union{Float64, String}    # The string or number 
    index::Union{Int64, Array}     # vector component, 1 for scalars, -1 for numbers, "INDEXEDBY..." for indexed vars
    derivs::Array                   # any derivatives applied
    flags::Array                    # any non-derivative modifiers or flags
end

# An operator struct that ties a symbol to a function.
struct SymOperator
    symbol::Symbol      # The name used in the input expression
    op::Function        # Function handle for the operator
end

# Struct for keeping the parent info.
struct ParentMaps
    num_parents::Int                # Number of parents including ghosts
    children_per_parent::Int        # Number of children per parent assuming similar element types
    faces_per_parent::Int           # Number of child faces in a parent
    patch_size::Int                 # Number of children in a patch (parent plus neighboring parents)
    
    child2parent::Array{Int,2}      # size = (2, allChildren) 1 = index of parent, 2 = location within parent
    parent2child::Array{Int,2}      # size = (myChildren, allParents) global index of each child in each parent
    parent2face::Array{Int,2}       # size = (myFaces, allParents) global index of each face in each parent
    cface2pface::Array{Int,3}       # size = (elementFaces, myChildren, allParents) relative face index in parent
    parent2neighbor::Array{Int,2}   # size = (outerFaces, allParents) index of neighboring parents
    
    patches::Array{Int, 2}          # size = (outerfaces*neighborChildren) local patch around each parent
    leftCells::Vector{Vector{Int}}          # Patch indices for left and right cell groups for each face
    rightCells::Vector{Vector{Int}}         #
    
    face_neighborhoods::Matrix{Vector{Int}} # indices of a group of elements to the left and right of each face
    
    ParentMaps() = new(0,0,0,0,zeros(Int,0,0),zeros(Int,0,0),zeros(Int,0,0),zeros(Int,0,0,0),
    zeros(Int,0,0),zeros(Int,0,0),Vector{Vector{Int}}(undef,0),Vector{Vector{Int}}(undef,0),Matrix{Vector{Int}}(undef,0,0))
    
    ParentMaps(num_parents,children_per_parent,face_per_parent,patch_size,child2parent,parent2child,parent2face,
    cface2pface, parent2neighbor,patches,leftCells,rightCells,face_neighborhoods) = 
    new(num_parents,children_per_parent,face_per_parent,patch_size,child2parent,parent2child,parent2face,
    cface2pface, parent2neighbor,patches,leftCells,rightCells,face_neighborhoods)
end

struct FVInfo{T<:AbstractFloat}
    fluxOrder::Int                  # Order for flux reconstruction
    cellCenters::Matrix{T}             # Coordinates of cell centers
    faceCenters::Matrix{T}             # Coordinates of face centers
    
    # cell2node::Vector{Vector{Int}}       # Neighboring cell indices for each node
    # cell2nodeWeight::Vector{Vector{T}} # Cell to node interpolation weights
    
    parentMaps::ParentMaps # For when there is a parent-child mesh
end

#=
# Stores the elemental jacobian
# Used by the "geometric_factors" function
=#
struct Jacobian{T<:AbstractFloat}
    rx::Vector{T}
    ry::Vector{T}
    rz::Vector{T}
    sx::Vector{T}
    sy::Vector{T}
    sz::Vector{T}
    tx::Vector{T}
    ty::Vector{T}
    tz::Vector{T}
end

#=
Stores all of the geometric factors for the grid.
=#
struct GeometricFactors{T<:AbstractFloat}
    J::Vector{Jacobian{T}}        # Jacobian for each element
    detJ::Matrix{T}      # Determinant of Jacobian for each element
    
    # These below are only computed if needed, otherwise empty arrays
    volume::Vector{T}    # Volume of each element (used by FV)
    face_detJ::Vector{T}  # Determinant of Jacobian for each face (used by DG and FV)
    area::Vector{T}      # Area of each face (used by FV)
end

struct MeshData
    #### Minimal required information ####
    # Nodes
    nx::Int;                    # Number of vertices
    nodes::Array{Float64,2};    # vertex locations (array has size (dim,nx))
    indices::Array{Int,1};      # vertex indices may not be in order
    # Elements
    nel::Int;                   # Number of elements
    elements::Array{Int,2};     # Element vertex mapping (array has size (Np, nel))*assumes only one element type
    etypes::Array{Int,1};       # Element types as defined by GMSH
    nv::Array{Int,1};           # Number of vertices for each element. Only different if they are different types,
    
    #### Optional information that will be built if not provided ####
    invind::Array{Int,1}        # Inverse of indices, maps vertex index to position in nodes array (invind[indices[i]] = i)
    face2vertex::Array{Int,2}   # Vertices defining each face (array has size (Nfp, Nfaces))
    face2element::Array{Int,2}  # Indices of elements on each side of the face. If 0, it is a boundary face. (size is (2,Nfaces))
    element2face::Array{Int,2}  # Indices of faces on each side of the element. (size is (NfacesPerElement, nel))
    normals::Array{Float64,2}   # Normal vectors for each face pointing from first to second in face2element order (size is (dim, Nfaces))
    bdryID::Array{Int,1};       # Boundary ID for each face (0=interior face)
    
    # The minimal constructor needs to build the optional information.
    # Note: Must uncomment to build.
    MeshData(n, x, ind, ne, el, et, v) = (
        # inv = invert_index(ind);
        # face2v = Array{Int,2}(undef,0,0);
        # face2e = Array{Int,2}(undef,0,0);
        # e2face = Array{Int,2}(undef,0,0);
        # norms = Array{Float64,2}(undef,0,0);
        # bdry = Array{Int,1}(undef,0);
        
        # uncomment these to compute. WARNING: can be slow
        inv = invert_index(ind);
        (face2v, face2e, e2face) = build_faces(ne, el, et);
        norms = find_normals(face2v, x);
        bdry = find_boundaries(face2e);
        new(n, x, ind, ne, el, et, v, inv, face2v, face2e, e2face, norms, bdry);
    )
    # The complete constructor
    MeshData(n, x, ind, ne, el, et, v, inv, face2v, face2e, e2face, norms, bdry) = (
        new(n, x, ind, ne, el, et, v, inv, face2v, face2e, e2face, norms, bdry);
    )
    # An empty mesh
    MeshData() = new(
        0, zeros(0,0), zeros(Int,0), 0, zeros(0,0), zeros(Int,0), zeros(Int,0), zeros(Int,0),
        zeros(Int,0,0), zeros(Int,0,0), zeros(Int,0,0), zeros(0,0), zeros(Int,0)
    )
end

struct Grid{T<:AbstractFloat}
    # nodes
    allnodes::Matrix{T}         # All node coordinates size = (dim, nnodes)
    
    # boundaries
    bdry::Vector{Vector{Int}}        # Indices of boundary nodes for each BID (bdry[bid][nodes])*note:array of arrays
    bdryface::Vector{Vector{Int}}    # Indices of faces touching each BID (bdryface[bid][faces])*note:array of arrays
    bdrynorm::Vector{Matrix{T}}    # Normal vector for boundary nodes for each BID (bdrynorm[bid][dim, nodes])*note:array of arrays
    bids::Vector{Int}                # BID corresponding to rows of bdrynodes
    nodebid::Vector{Int}             # BID for every node in allnodes order(interior=0)
    
    # elements
    loc2glb::Matrix{Int}             # local to global map for each element's nodes (size is (Np, nel))
    glbvertex::Matrix{Int}           # global indices of each elements' vertices (size is (Nvertex, nel))
    
    # faces (For CG, G=1. For DG, G=2)
    face2glb::Array{Int}             # local to global map for faces (size is (Nfp, G, Nfaces))
    element2face::Matrix{Int}        # face indices for each element (size is (Nfaces, nel))
    face2element::Matrix{Int}        # elements on both sides of a face, 0=boundary (size is (2, Nfaces))
    facenormals::Matrix{T}         # normal vector for each face
    faceRefelInd::Matrix{Int}        # Index for face within the refel for each side
    facebid::Vector{Int}             # BID of each face (0=interior face)
    
    # When partitioning the grid, this stores the ghost info.
    # Items specifying (for solver type) will be empty/0 for other types.
    is_subgrid::Bool            # Is this a partition of a greater grid?
    elemental_order::Vector{Int}     # Order used in elemental loops
    nel_global::Int             # Number of global elements
    nel_owned::Int              # Number of elements owned by this partition
    nel_ghost::Int              # Number of ghost elements (for FV)
    nface_owned::Int            # Number of faces owned by this partition
    nface_ghost::Int            # Number of ghost faces that are not owned (for FV)
    nnodes_global::Int          # Number of global nodes
    nnodes_borrowed::Int        # Number of nodes borrowed from another partition (for CG)
    element_owner::Vector{Int}       # The rank of each element's owner or -1 if locally owned (for FV)
    node_owner::Vector{Int}          # The rank of each node's owner (for FE)
    partition2global_element::Vector{Int}           # Map from partition elements to global mesh element index
    partition2global::Vector{Int}    # Global index of nodes (for CG,DG)
    global_bdry_index::Vector{Int8}   # Index in bids for every global node, or 0 for interior (Only proc 0 holds, only for FE)
    
    num_neighbor_partitions::Int   # number of partitions that share ghosts with this.
    neighboring_partitions::Vector{Int} # IDs of neighboring partitions
    ghost_counts::Vector{Int}           # How many ghost elements for each neighbor (for FV)
    ghost_index::Vector{Matrix{Int}}     # Lists of ghost elements to send/recv for each neighbor (for FV)
    
    # constructors
    Grid(T::DataType, allnodes, bdry, bdryfc, bdrynorm, bids, nodebid, loc2glb, glbvertex, f2glb, element2face, 
         face2element, facenormals, faceRefelInd, facebid) = 
     new{T}(allnodes, bdry, bdryfc, bdrynorm, bids, nodebid, loc2glb, glbvertex, f2glb, element2face, 
         face2element, facenormals, faceRefelInd, facebid, 
         false, Array(1:size(loc2glb,2)), size(loc2glb,2), size(loc2glb,2), 0,size(face2element,2), 0, 0, 0, zeros(Int,0), zeros(Int,0), 
         zeros(Int,0), zeros(Int,0), zeros(Int8,0), 0, zeros(Int,0), zeros(Int,0), [zeros(Int,2,0)]); # up to facebid only
     
    Grid(T::DataType, allnodes, bdry, bdryfc, bdrynorm, bids, nodebid, loc2glb, glbvertex, f2glb, element2face, 
         face2element, facenormals, faceRefelInd, facebid, 
         ispartitioned, el_order, nel_global, nel_owned, nel_ghost, nface_owned, nface_ghost, nnodes_global, nnodes_borrowed, element_owners, 
         node_owner, partition2global_element, partition2global, glb_bid, num_neighbors, neighbor_ids, ghost_counts, ghost_ind) = 
     new{T}(allnodes, bdry, bdryfc, bdrynorm, bids, nodebid, loc2glb, glbvertex, f2glb, element2face, 
         face2element, facenormals, faceRefelInd, facebid, 
         ispartitioned, el_order, nel_global, nel_owned, nel_ghost, nface_owned, nface_ghost, nnodes_global, nnodes_borrowed, element_owners, 
         node_owner, partition2global_element, partition2global, glb_bid, num_neighbors, neighbor_ids, ghost_counts, ghost_ind); # subgrid parts included
         
    # An empty Grid
    Grid(T::DataType) = new{T}(
        zeros(T,0,0),[zeros(Int,0)],[zeros(Int,0)],[zeros(T,0,0)],zeros(Int,0),zeros(Int,0),
        zeros(Int,0,0),zeros(Int,0,0),
        zeros(Int,0,0,0),
        zeros(Int,0,0), zeros(Int,0,0), zeros(T,0,0), zeros(Int,0,0),
        zeros(Int,0),
        false,zeros(Int,0),0,0,0,0,0,0,0,zeros(Int,0),zeros(Int,0),zeros(Int,0),zeros(Int,0),zeros(Int8,0),
        0,zeros(Int,0),zeros(Int,0),[zeros(Int,0,0)]
    )
end

mutable struct Refel{T<:AbstractFloat}
    dim::Int                # Dimension
    N::Int                  # Order of polynomials
    Np::Int                 # Number of nodes
    Nqp::Int                # Number of quadrature points
    Nfaces::Int             # Number of faces
    Nfp::Vector{Int}       # Number of nodes for each face
    
    ######################################
    # Volume nodes and quadrature matrices
    ######################################
    r1d::Vector{T}     # Node coordinates in 1D
    r::Matrix{T}       # dim-dim Node coordinates
    
    wr1d::Vector{T}    # r1d gll Quadrature weights
    wr::Vector{T}      # r gll Quadrature weights
    
    g1d::Vector{T}     # 1D Gauss points
    wg1d::Vector{T}    # 1D Gauss weights
    
    g::Matrix{T}       # dim-dim Gauss points
    wg::Vector{T}      # dim-dim Gauss weights
    
    V::Matrix{T}       # basis at r
    gradV::Matrix{T}   # grad of basis at r
    invV::Matrix{T}    # Inverse V
    
    Vg::Matrix{T}      # basis at Gauss
    gradVg::Matrix{T}  # grad of basis at g
    invVg::Matrix{T}   # Inverse Vg
    
    Dr::Matrix{T}      # Differentiation matrix for r
    Ds::Matrix{T}      # Differentiation matrix for s
    Dt::Matrix{T}      # Differentiation matrix for t
    Dg::Matrix{T}      # Differentiation matrix for g
    
    # Useful quadrature matrices for the volume integrals
    Q1d::Matrix{T}     # 1D quadrature matrix: like Vg*invV
    Q::Matrix{T}       # dim-dim quadrature matrix
    Qr::Matrix{T}      # quad of derivative matrix: like gradVg*invV
    Qs::Matrix{T}      # 
    Qt::Matrix{T}      # 
    
    Ddr::Matrix{T}      # Derivatives at the elemental nodes, not quadrature nodes
    Dds::Matrix{T}      # 
    Ddt::Matrix{T}      #
    
    #######################################
    # Surface nodes and quadrature matrices
    #######################################
    face2local::Vector{Vector{Int}}       # maps face nodes to local indices
    
    surf_r::Vector{Matrix{T}}       # surface node coordinates
    surf_wr::Vector{Vector{T}}      # surface gll weights
    
    surf_g::Vector{Matrix{T}}       # surface Gauss points
    surf_wg::Vector{Vector{T}}      # surface Gauss weights
    
    surf_V::Vector{Matrix{T}}       # basis at surf_r
    surf_gradV::Vector{Matrix{T}}   # grad of basis at surf_r
    surf_invV::Vector{Matrix{T}}    # Inverse surf_V
    
    surf_Vg::Vector{Matrix{T}}      # basis at surf_g
    surf_gradVg::Vector{Matrix{T}}  # grad of basis at surf_g
    surf_invVg::Vector{Matrix{T}}   # Inverse surf_Vg
    
    surf_Dr::Vector{Matrix{T}}      # Differentiation matrix for surf_r
    surf_Ds::Vector{Matrix{T}}      # Differentiation matrix for surf_s
    surf_Dt::Vector{Matrix{T}}      # Differentiation matrix for surf_t
    surf_Dg::Vector{Matrix{T}}      # Differentiation matrix for surf_g
    
    surf_Q::Vector{Matrix{T}}       # quadrature matrix
    surf_Qr::Vector{Matrix{T}}      # derivative quadrature matrix
    surf_Qs::Vector{Matrix{T}}      # 
    surf_Qt::Vector{Matrix{T}}      # 
    
    surf_Ddr::Vector{Matrix{T}}     # Derivatives at the elemental nodes, not quadrature nodes
    surf_Dds::Vector{Matrix{T}}     # 
    surf_Ddt::Vector{Matrix{T}}     #
    
    # Constructor needs at least this information
    Refel(T::DataType, dim, order, nnodes, nfaces, nfp) = new{T}(
        dim,
        order,
        nnodes,
        -1,
        nfaces,
        nfp,
        zeros(T,0),zeros(T,0,0),
        zeros(T,0),zeros(T,0),
        zeros(T,0),zeros(T,0),
        zeros(T,0,0),zeros(T,0),
        zeros(T,0,0),zeros(T,0,0),zeros(T,0,0),
        zeros(T,0,0),zeros(T,0,0),zeros(T,0,0),
        zeros(T,0,0),zeros(T,0,0),zeros(T,0,0),zeros(T,0,0),
        zeros(T,0,0),zeros(T,0,0),zeros(T,0,0),zeros(T,0,0),zeros(T,0,0),
        zeros(T,0,0),zeros(T,0,0),zeros(T,0,0),
        [zeros(Int,0)],
        [zeros(T,0,0)],[zeros(T,0)],
        [zeros(T,0,0)],[zeros(T,0)],
        [zeros(T,0,0)],[zeros(T,0,0)],[zeros(T,0,0)],
        [zeros(T,0,0)],[zeros(T,0,0)],[zeros(T,0,0)],
        [zeros(T,0,0)],[zeros(T,0,0)],[zeros(T,0,0)],[zeros(T,0,0)],
        [zeros(T,0,0)],[zeros(T,0,0)],[zeros(T,0,0)],[zeros(T,0,0)],
        [zeros(T,0,0)],[zeros(T,0,0)],[zeros(T,0,0)]
    )
end

mutable struct Stepper
    type::String;       # The constant for the stepper type
    implicit::Bool;     # implicit or explicit
    Nsteps::Int;        # number of steps
    dt::Float64;        # step size
    cfl::Float64;       # CFL number
    stages::Int;        # how many stages
    a::Array{Float64};  # for RK steppers
    b::Array{Float64};  # for RK steppers
    c::Array{Float64};  # for RK steppers
    
    Stepper(t, c) = new(t, false, 0, 0, c, 0, [], [], []);
end

# Buffers for the global system
mutable struct ParallelBuffers{T<:AbstractFloat}
    full_AI::Vector{Int}
    full_AJ::Vector{Int}
    full_AV::Vector{T}
    full_b::Vector{T}
    b_order::Vector{Int}
    
    vec_b::Vector{T}
end

# All of the IR node types are subtypes of this.
abstract type IR_part end;
#

# The Finch state is completely maintained within this struct.
# This should keep all data that is in any way related to the 
# computation out of the global scope.
mutable struct FinchState{T<:AbstractFloat}
    # Configuration
    config::FinchConfig
    project_name::String
    output_dir::String
    
    # Log
    use_log::Bool
    log_level::Int
    log_file::String
    log_line_index::Int
    
    # Timer
    timer_output::TimerOutput
    
    # Code gen target
    external_target::Bool
    target_language::String
    target_framework::String
    target_parameters::Dict
    
    # Problem specification
    prob::FinchProblem
    
    # Mesh
    mesh_data::MeshData # The basic element information as read from a MSH file or generated here.
    grid_data::Grid{T} # The full collection of nodes(including internal nodes) and other mesh info in the actual DOF ordering.
    fv_grid::Grid{T} # This FV version is only made if using FV.
    dg_grid::Grid{T} # This DG version is only made if using mixed CG/DG. Otherwise it is in grid_data.
    needed_grid_types::Vector{Bool}
    geo_factors::GeometricFactors{T} # Precomputed geometric factors for each element
    refel::Refel{T} # Reference element
    fv_geo_factors::GeometricFactors{T} # Simplified version for FV
    fv_refel::Refel{T} # Simplified version for FV
    fv_info::FVInfo{T} # FV specific data
    
    # Entities
    variables::Vector{Variable{T}}
    coefficients::Vector{Coefficient}
    parameters::Vector{Parameter}
    test_functions::Vector{Coefficient}
    indexers::Vector{Indexer}
    ordered_indexers::Vector{Indexer}
    variable_transforms::Vector{VariableTransform}
    
    # Time stepper
    time_stepper::Stepper
    time_steppers::Vector{Stepper} # only used for mixed methods with separate steppers
    use_specified_steps::Bool
    specified_dt::Float64
    specified_Nsteps::Int
    
    # Generated functions and code
    genfunctions::Vector{GenFunction}
    callback_functions::Vector{CallbackFunction}
    solve_functions::Vector{Union{GenFunction, IR_part, Nothing}}
    symexpressions::Vector{Vector}
    code_strings::Vector{String}
    
    # Misc
    use_cachesim::Bool
    parallel_buffers::ParallelBuffers{T}
    ops::Vector{SymOperator}
    
    # Default constructor
    FinchState(T::DataType, name::String) = new{T}(
        FinchConfig(),
        name,
        pwd(),
        
        false,
        2,
        "",
        1,
        
        TimerOutput(),
        
        false,
        JULIA,
        FINCH,
        Dict{String,Any}(),
        
        FinchProblem(),
        
        MeshData(),
        Grid(T),
        Grid(T),
        Grid(T),
        [false,false,false],
        GeometricFactors{T}([],zeros(0,0),[],[],[]),
        Refel(T,1,1,0,0,[1,1]),
        GeometricFactors{T}([],zeros(0,0),[],[],[]),
        Refel(T,1,1,0,0,[1,1]),
        FVInfo{T}(1,zeros(0,0),zeros(0,0),ParentMaps()),
        
        [], [], [], [], [], [], [],
        
        Stepper(EULER_EXPLICIT, 0),
        [],
        false,
        0,
        0,
        
        [], [], [], [[],[],[],[]], [],
        
        false,
        ParallelBuffers{T}([],[],[],[],[],[]),
        []
    )
end