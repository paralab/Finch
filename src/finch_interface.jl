#=
This file contains all of the common interface functions.
=#
export initFinch, generateFor, useLog, indexDataType, floatDataType,
        domain, solverType, functionSpace, finiteVolumeOrder,
        nodeType, timeStepper, setSteps, linAlgOptions, usePetsc, customOperator, customOperatorFile,
        mesh, exportMesh, variable, coefficient, parameter, testSymbol, index, boundary, addBoundaryID,
        referencePoint, timeInterval, initial, preStepFunction, postStepFunction, callbackFunction,
        variableTransform, transformVariable,
        weakForm, conservationForm, assemblyLoops,
        exportCode, importCode, printLatex,
        evalInitialConditions, nonlinear, 
        solve, cachesimSolve, 
        finalizeFinch, cachesim, outputValues,
        mortonNodes, hilbertNodes, tiledNodes, mortonElements, hilbertElements, 
        tiledElements, elementFirstNodes, randomNodes, randomElements

# Begin configuration setting functions

"""
initFinch(name = "unnamedProject", floatType::DataType=Float64)

This initializes and returns the Finch state.
The name of the project can be set here.
T is the data type to be used for floating point data.
T must be a subtype of AbstractFloat.
Note that while this generally applies to data arrays relevant to
the computation, some places may still use Float64, so the 
corresponding conversions should be defined.
"""
function initFinch(name="unnamedProject", floatType::DataType=Float64)
    return init_finch(floatType, name);
end

"""
    generateFor(lang; filename=project_name, header="", params...)

Specify the generation target. Lang could be one of the included target constants: 
(MATLAB, DENDRO) or the filename where the target is defined. The keyword argument
filename refers to the name to be applied to the generated code. The header text 
will be placed at the top of each generated code file. If the target requires 
some extra parameters, those are included in params and available to the target 
as a `Dict{Symbol, Any}`.
"""
function generateFor(lang; filename=finch_state.project_name, header="", params...)
    global finch_state;
    outputDirPath = pwd()*"/"*filename;
    if finch_state.config.proc_rank == 0 && !isdir(outputDirPath)
        mkdir(outputDirPath);
    end
    framew = 0;
    if !in(lang, [CPP,MATLAB,DENDRO, "Dendrite"])
        # lang should be a filename for a custom target
        # This file must include these three functions:
        # 1. get_external_language_elements() - file extensions, comment chars etc.
        # 3. generate_external_files(var, IR) - Writes all files based on generated code
        finch_state.target_framework = CUSTOM_GEN_TARGET;
        finch_state.target_language = CUSTOM_GEN_TARGET;
        include(lang);
        set_custom_gen_target(finch_state, get_external_language_elements, generate_external_files, outputDirPath, filename, head=header);
        
    else # Use an included target
        target_dir = @__DIR__
        if lang == DENDRO
            printerr("Sorry, this target is not ready for this version of Finch.", fatal=true)
            finch_state.target_framework = DENDRO;
            finch_state.target_language = CPP;
            target_file = "/targets/target_dendro.jl";
        elseif lang == MATLAB
            printerr("Sorry, this target is not ready for this version of Finch.", fatal=true)
            finch_state.target_framework = "";
            finch_state.target_language = MATLAB;
            target_file = "/targets/target_matlab.jl";
        elseif lang == "Dendrite"
            finch_state.target_framework = "Dendrite";
            finch_state.target_language = CPP;
            target_file = "/targets/target_dendrite.jl";
        else #CPP
            printerr("Sorry, this target is not ready for this version of Finch.", fatal=true)
            finch_state.target_framework = "";
            finch_state.target_language = CPP;
            target_file = "/targets/target_cpp.jl";
        end
        include(target_dir * target_file);
        set_custom_gen_target(finch_state, get_external_language_elements, generate_external_files, outputDirPath, filename, head=header);
    end
    
    set_codegen_parameters(finch_state, params);
end

"""
    useLog(name=project_name; dir=output_dir, level=2)

Turn on logging with the given file name and optional directory.
The verbosity level can be 1(basic progress info), 2(More details about
each step), or 3(everything).
"""
function useLog(name=finch_state.project_name; dir=finch_state.output_dir, level=2)
    init_log(finch_state, name, dir, level);
end

"""
    indexDataType(type<:Integer)

Set the data type to be used for indices.
The default is Int64.
"""
function indexDataType(type)
    printerr("indexDataType() is no longer available. Sorry")
end

"""
    floatDataType(type<:AbstractFloat)

This function has been removed. To set the float type, see 
initFinch(name, type)
"""
function floatDataType(type)
    printerr("floatDataType() is no longer available. To set float type, use initFinch(name, type).")
end

"""
    domain(dims; shape=SQUARE, grid=UNIFORM_GRID)

Set the dimensionality of the domain. The shape(SQUARE, IRREGULAR) and 
grid type(UNIFORM_GRID, UNSTRUCTURED, TREE) can be set, but may be changed
when building or importing the mesh.
"""
function domain(dims; shape=SQUARE, grid=UNIFORM_GRID)
    finch_state.config.dimension = dims;
    finch_state.config.geometry = shape;
    finch_state.config.mesh_type = grid;
end

"""
    solverType(method)

Select between CG, DG, or FV methods.
"""
function solverType(method)
    set_solver(finch_state, method);
end

"""
    functionSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)

Set the polynomial order and type of polynomials for FEM.
Some of these are placeholders, so only use order at this point.
"""
function functionSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
    finch_state.config.trial_function = space;
    finch_state.config.test_function = space;
    if orderMax > orderMin && orderMin >= 0
        finch_state.config.p_adaptive = true;
        finch_state.config.basis_order_min = orderMin;
        finch_state.config.basis_order_max = orderMax;
    else
        finch_state.config.basis_order_min = max(order, orderMin);
        finch_state.config.basis_order_max = max(order, orderMin);
    end
end

"""
    nodeType(type)

For FEM, set the nodal configuration within elements.
The default is LOBATTO. GAUSS and UNIFORM are available, but should be used with care.
"""
function nodeType(type)
    finch_state.config.elemental_nodes = type;
end

"""
    timeStepper(type; cfl=0)

Set the type of time stepping method and optionally the CFL number.
Options include `EULER_EXPLICIT`, `EULER_IMPLICIT`, `CRANK_NICHOLSON`, `RK4`, `LSRK4`.
There are some other options for specific targets, so see their documentation.
If no CFL number is provided, one will be chosen based on the mesh and stepper type.
"""
function timeStepper(type; cfl=0)
    set_stepper(finch_state, type, cfl);
end

"""
    timeInterval(T)

Set the ending time for time stepping. This is overridden if time steps
are manually specified.
"""
function timeInterval(T)
    finch_state.prob.time_dependent = true;
    if finch_state.time_stepper === nothing
        timeStepper(EULER_IMPLICIT);
    end
    finch_state.prob.end_time = T;
end

"""
    setSteps(dt, steps)

Manually set the time steps if desired.
"""
function setSteps(dt, steps)
    finch_state.prob.time_dependent = true;
    if finch_state.time_stepper === nothing
        timeStepper(EULER_IMPLICIT);
    end
    set_specified_steps(finch_state, dt, steps);
end

"""
    linAlgOptions(kwargs...)

Set options for solving the linear system.
Matrix-free must use an iterative method.
Iterative methods include GMRES or CG provided by IterativeSolvers.jl
Available preconditioners are AMG and ILU provided by
AlgebraicMultigrid.jl and IncompleteLU.jl.
Defaults maxiter=0, abstol=0, and gmresRestart=0 will use the 
corresponding defaults from IterativeSolvers.jl

**Keywords**
* `matrixFree=false`
* `iterative=false`
* `method="GMRES"` or `"CG"`
* `pc="ILU"` or `"AMG"` or `"NONE"`
* `maxiter::Int=0` 0 will result in `size(A, 2)`
* `abstol=0`
* `reltol=1e-8`
* `gmresRestart=0` 0 will result in `min(20, size(A, 2))`
* `verbose::Bool=false` Print convergence info for each iteration.

"""
function linAlgOptions(;matrixFree::Bool=false, iterative::Bool=false, method::String="GMRES", 
                    pc::String="ILU", maxiter::Int=0, abstol=0, reltol=1e-8, 
                    gmresRestart::Int=0, verbose::Bool=false)
    finch_state.config.linalg_matrixfree = matrixFree;
    finch_state.config.linalg_iterative = iterative;
    finch_state.config.linalg_iterative_method = method;
    finch_state.config.linalg_iterative_pc = pc;
    finch_state.config.linalg_iterative_maxiter = maxiter;
    finch_state.config.linalg_iterative_abstol = abstol;
    finch_state.config.linalg_iterative_reltol = reltol;
    finch_state.config.linalg_iterative_gmresRestart = gmresRestart;
    finch_state.config.linalg_iterative_verbose = verbose;
    
    if matrixFree
        log_entry("Using matrix-free with: "*method*"(pc="*pc*")", 2);
    elseif iterative
        log_entry("Set iterative solver to "*method*"(pc="*pc*")", 2);
    else
        log_entry("Using default A\\b", 2);
    end
end

"""
    usePetsc(useit::Bool=true)

If PETSc is available, use it.
If it is not available, this won't work.
"""
function usePetsc(useit::Bool=true)
    if @isdefined(PETSc)
        # Initialize using the first available library.
        # This should have been set up beforehand.
        if length(PETSc.petsclibs) == 0
            printerr(
"No PETSc library found. Please build PETSc first.
If you have a preferred PETSc library or the one supplied with PETSc.jl
is causing trouble, do this:
\$ export JULIA_PETSC_LIBRARY=/path/to/your/petsc_lib.so
julia> ]build PETSc", fatal=true);
        end
        petsclib = PETSc.petsclibs[1];
        PETSc.initialize(petsclib)
        
        finch_state.config.linalg.usePetsc = useit;
    else
        printerr("Cannot use PETSc. Set up PETSc.jl manually first. Proceeding with default.")
    end
end

"""
    customOperator(name, handle)

Define a new symbolic operator to be used in PDE expressions. The name is the
symbol that will be used in expressions. The function handle points to the 
operator's function.
"""
function customOperator(name, handle)
    s = Symbol(name);
    add_custom_op(s, handle);
end

"""
    customOperatorFile(filename)

Import a set of symbolic operators defined in a file. The file must contain
certain elements. See an example.
"""
function customOperatorFile(filename)
    log_entry("Adding custom operators from file: "*string(filename), 2);
    add_custom_op_file(filename);
end

# End configuration functions, begin problem definition functions

"""
    mesh(msh; elsperdim=5, bids=1, interval=[0,1], partitions=0)

Build or import a mesh. msh can be either a constant(LINEMESH, QUADMESH, HEXMESH)
to build a mesh with the built-in simple mesh generator, or a filename for a mesh
file. Currently GMSH files(.msh), either old or new versions, and MEDIT files(.mesh)
are supported.
"""
function mesh(msh; elsperdim=5, bids=1, interval=[0,1], partitions=0)
    
    @timeit finch_state.timer_output "Mesh" begin
    
    @timeit finch_state.timer_output "gen/read" begin
    if msh == LINEMESH
        log_entry("Building simple line mesh with nx elements, nx="*string(elsperdim));
        mshdat = simple_line_mesh(elsperdim.+1, bids, interval);
        
    elseif msh == QUADMESH
        if length(interval) == 2
            interval = [interval[1], interval[2], interval[1], interval[2]];
        end
        log_entry("Building simple quad mesh with nx*nx elements, nx="*string(elsperdim));
        mshdat = simple_quad_mesh(elsperdim.+1, bids, interval);
        
    elseif msh == TRIMESH
        if length(interval) == 2
            interval = [interval[1], interval[2], interval[1], interval[2]];
        end
        log_entry("Building simple triangle mesh with nx*nx*2 elements, nx="*string(elsperdim));
        mshdat = simple_tri_mesh(elsperdim.+1, bids, interval);
        finch_state.config.mesh_type = UNSTRUCTURED;
        
    elseif msh == HEXMESH
        if length(interval) == 2
            interval = [interval[1], interval[2], interval[1], interval[2], interval[1], interval[2]];
        end
        log_entry("Building simple hex mesh with nx*nx*nx elements, nx="*string(elsperdim));
        mshdat = simple_hex_mesh(elsperdim.+1, bids, interval);
        
    else # msh should be a mesh file name
        # open the file and read the mesh data
        mfile = open(msh, "r");
        log_entry("Reading mesh file: "*msh);
        mshdat=read_mesh(mfile);
        close(mfile);
        # Assume an irregular mesh
        finch_state.config.geometry = IRREGULAR;
        finch_state.config.mesh_type = UNSTRUCTURED;
    end
    
    end # gen/read timer
    
    # Set this in Finch and build the corresponding Grid
    add_mesh(finch_state, mshdat, partitions=partitions);
    
    # If bids>1 were specified for built meshes, add them here
    bid_defs = []
    if bids > 1
        @timeit finch_state.timer_output "bids" begin
        if msh == LINEMESH
            if bids == 2
                # already done in mesh_data
                # add_boundary_ID_to_grid(2, x -> (x >= interval[2]), finch_state.grid_data);
                add_boundary_ID_to_problem(1, "XMIN");
                add_boundary_ID_to_problem(2, "XMAX");
            end
            
        elseif msh == QUADMESH || msh == TRIMESH
            if bids == 2
                add_boundary_ID_to_grid(2, (x,y) -> (abs(y - interval[3]) <= eps()) || (abs(y - interval[4]) <= eps()), finch_state.grid_data);
                add_boundary_ID_to_problem(1, "XMIN || XMAX");
                add_boundary_ID_to_problem(2, "YMIN || YMAX");
            elseif bids == 3
                add_boundary_ID_to_grid(2, (x,y) -> (x >= interval[2]), finch_state.grid_data);
                add_boundary_ID_to_grid(3, (x,y) -> ((y <= interval[3]) || (y >= interval[4])) && (x > interval[1] && x < interval[2]), finch_state.grid_data);
                add_boundary_ID_to_problem(1, "XMIN");
                add_boundary_ID_to_problem(2, "XMAX");
                add_boundary_ID_to_problem(3, "YMIN || YMAX");
            elseif bids == 4
                add_boundary_ID_to_grid(2, (x,y) -> (x >= interval[2]), finch_state.grid_data);
                add_boundary_ID_to_grid(3, (x,y) -> (y <= interval[3] && (x > interval[1] && x < interval[2])), finch_state.grid_data);
                add_boundary_ID_to_grid(4, (x,y) -> (y >= interval[4] && (x > interval[1] && x < interval[2])), finch_state.grid_data);
                add_boundary_ID_to_problem(1, "XMIN");
                add_boundary_ID_to_problem(2, "XMAX");
                add_boundary_ID_to_problem(3, "YMIN");
                add_boundary_ID_to_problem(4, "YMAX");
            end
            
        elseif msh == HEXMESH
            tiny = 1e-13;
            if bids == 6
                # bids = [1,2,3,4,5,6]; # all separate
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny), finch_state.grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> (y <= interval[3]+tiny), finch_state.grid_data);
                add_boundary_ID_to_grid(4, (x,y,z) -> (y >= interval[4]-tiny), finch_state.grid_data);
                add_boundary_ID_to_grid(5, (x,y,z) -> (z <= interval[5]+tiny), finch_state.grid_data);
                add_boundary_ID_to_grid(6, (x,y,z) -> (z >= interval[6]-tiny), finch_state.grid_data);
                add_boundary_ID_to_problem(1, "XMIN");
                add_boundary_ID_to_problem(2, "XMAX");
                add_boundary_ID_to_problem(3, "YMIN");
                add_boundary_ID_to_problem(4, "YMAX");
                add_boundary_ID_to_problem(5, "ZMIN");
                add_boundary_ID_to_problem(6, "ZMAX");
            elseif bids == 5
                # bids = [1,2,3,4,5]; # combine z
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny), finch_state.grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> (y <= interval[3]+tiny), finch_state.grid_data);
                add_boundary_ID_to_grid(4, (x,y,z) -> (y >= interval[4]-tiny), finch_state.grid_data);
                add_boundary_ID_to_grid(5, (x,y,z) -> (z <= interval[5]+tiny) || (z >= interval[6]-tiny), finch_state.grid_data);
                add_boundary_ID_to_problem(1, "XMIN");
                add_boundary_ID_to_problem(2, "XMAX");
                add_boundary_ID_to_problem(3, "YMIN");
                add_boundary_ID_to_problem(4, "YMAX");
                add_boundary_ID_to_problem(5, "ZMIN || ZMAX");
            elseif bids == 4
                # bids = [1,2,3,4]; # combine y and z
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]), finch_state.grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> ((y <= interval[3]+tiny) || (y >= interval[4]-tiny)), finch_state.grid_data);
                add_boundary_ID_to_grid(4, (x,y,z) -> ((z <= interval[5]+tiny) || (z >= interval[6]-tiny)), finch_state.grid_data);
                add_boundary_ID_to_problem(1, "XMIN");
                add_boundary_ID_to_problem(2, "XMAX");
                add_boundary_ID_to_problem(3, "YMIN || YMAX");
                add_boundary_ID_to_problem(4, "ZMIN || ZMAX");
            elseif bids == 3
                # bids = [1,2,3]; # combine x,y,z
                add_boundary_ID_to_grid(2, (x,y,z) -> (y <= interval[3]+tiny) || (y >= interval[4]-tiny), finch_state.grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> (z <= interval[5]+tiny) || (z >= interval[6]-tiny), finch_state.grid_data);
                add_boundary_ID_to_problem(1, "XMIN || XMAX");
                add_boundary_ID_to_problem(2, "YMIN || YMAX");
                add_boundary_ID_to_problem(3, "ZMIN || ZMAX");
            elseif bids == 2
                # bids = [1,2]; # x=0, other
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny) || (y <= interval[3]+tiny) || (y >= interval[4]-tiny) || (z <= interval[5]+tiny) || (z >= interval[6]-tiny), finch_state.grid_data);
                add_boundary_ID_to_problem(1, "XMIN");
                add_boundary_ID_to_problem(2, "XMAX || YMIN || YMAX || ZMIN || ZMAX");
            end
            
        end
        end # bids timer
        
    else # only 1 bid
        dimension = finch_state.config.dimension;
        if dimension == 1
            add_boundary_ID_to_problem(1, "XMIN || XMAX");
        elseif dimension == 2
            add_boundary_ID_to_problem(1, "XMIN || XMAX || YMIN || YMAX");
        elseif dimension == 3
            add_boundary_ID_to_problem(1, "XMIN || XMAX || YMIN || YMAX || ZMIN || ZMAX");
        end
    end
    
    # If FV order was already set to > 1, set it here
    if finch_state.config.fv_order > 1
        @timeit finch_state.timer_output "parent-child maps" begin
        finiteVolumeOrder(finch_state.config.fv_order);
        end
    end
    
    end # timer block
end

"""
    exportMesh(filename, format=MSH_V2)

Export the mesh to a file with the given name. The format can be MSH_V2 or MSH_V4
for old and new style GMSH formats.
"""
function exportMesh(filename, format=MSH_V2)
    # open the file to write to
    mfile = open(filename, "w");
    log_entry("Writing mesh file: "*filename);
    output_mesh(finch_state, mfile, format);
    close(mfile);
end

"""
    finiteVolumeOrder(order)

Set the order of flux reconstruction for FVM.
For order > 1 this will cause the mesh to be subdivided into a parent/child mesh.
Take this into account when designing the mesh.
"""
function finiteVolumeOrder(order)
    if finch_state.config.dimension > 2 && order > 1
        printerr("Sorry, higher order FV is not ready for 3D (TODO: build parent/child mesh)\n Continuing with first order.");
        return;
    end
    
    finch_state.config.fv_order = order;
    if order > 1
        if length(finch_state.grid_data.allnodes) == 0
            # mesh hasn't been set yet
            set_parent_and_child(finch_state, nothing, nothing, order);
        else
            (parent, child) = divide_parent_grid(finch_state.grid_data, order);
            set_parent_and_child(finch_state, parent, child, order);
        end
    else # order == 1
        set_parent_and_child(finch_state, nothing, finch_state.fv_grid, order);
    end
end

"""
    variable(name; type=SCALAR, location=NODAL, method=CG, index=nothing)

Create a variable entity with name that will be used in expressions.
Type can be SCALAR, VECTOR, TENSOR, SYM_TENSOR, or VAR_ARRAY.
Location can be NODAL or CELL. Generally nodal will be used for FEM and 
cell for FVM. The method can be specified for this variable when using a
mixed solver, but will otherwise be the current solver type. It the type
is VAR_ARRAY, it should be indexed using a previously defined indexer or
an array of indexers.

Note that a variable does not have to represent an unknown. It can be used
to store any value, and can be used in expressions, but does not need to be 
solved for explicitly.
"""
function variable(name; type=SCALAR, location=NODAL, method=CG, index=nothing)
    varind = length(finch_state.variables) + 1;
    varsym = Symbol(name);
    # Not the default method?
    if !(finch_state.config.solver_type == MIXED) && !(finch_state.config.solver_type == method)
        method = finch_state.config.solver_type;
    end
    
    if index===nothing
        index = Vector{Indexer}(undef,0);
    elseif typeof(index) == Indexer
        index = [index];
    elseif !(typeof(index) <:Vector{Indexer})
        printerr("variable index must be an Indexer or vector of Indexers.", fatal=true);
    end
    # Just make an empty variable with the info known so far.
    ftype = finch_state.config.float_type;
    var = Variable(varsym, Vector{Basic}(undef,0), varind, type, location, method, zeros(ftype, 0,0), index, 
                    0, Vector{Variable{ftype}}(undef,0), false);
    add_variable(finch_state, var);
    return var;
end

"""
    coefficient(name, val; type=SCALAR, location=NODAL, element_array=false)

Create a coefficient entity with name that will be used in expressions.
The value can be a numerical constant, a string representing an expression of
coordinates(x, y, z, t), or an array of numbers corresponding to the location.
Type can be SCALAR, VECTOR, TENSOR, SYM_TENSOR, or VAR_ARRAY.
Location can be NODAL or CELL. Generally nodal will be used for FEM and 
cell for FVM. 
"""
function coefficient(name, val; type=SCALAR, location=NODAL, element_array=false, time_dependent=false)
    csym = Symbol(name);
    nfuns = makeFunctions(val); # if val is constant, nfuns will be 0
    time_dependent = time_dependent || check_time_dependence(val);
    return add_coefficient(finch_state, csym, type, location, val, nfuns, element_array, time_dependent);
end

"""
    parameter(name, val; type=SCALAR)

Create a parameter entity with name that will be used in expressions.
The value is a string expression that can include coordinates(x, y, z, t) and 
any variable and coefficient symbols. It is essentially a convenient object to
simplify more complicated expressions.
The type is not important as it will be determined by the symbolic expression it
represents.
"""
function parameter(name, val; type=SCALAR)
    if length(finch_state.parameters) == 0
        coefficient("parameterCoefficientForx", "x")
        coefficient("parameterCoefficientFory", "y")
        coefficient("parameterCoefficientForz", "z")
        coefficient("parameterCoefficientFort", "t")
    end
    if typeof(val) <: Number
        newval = [val];
    elseif typeof(val) == String
        # newval will be an array of expressions
        # search for x,y,z,t symbols and replace with special coefficients like parameterCoefficientForx
        newval = [swap_parameter_xyzt(Meta.parse(val))];
        
    elseif typeof(val) <: Array{String}
        newval = Array{Expr,1}(undef,length(val));
        for i=1:length(val)
            newval[i] = swap_parameter_xyzt(Meta.parse(val[i]));
        end
        newval = reshape(newval,size(val));
    else
        printerr("Error: use strings to define parameters", fatal=true);
        newval = 0;
    end
    
    return add_parameter(finch_state, Symbol(name), type, newval);
end

"""
    testSymbol(symbol; type=SCALAR)

Define a symbol for a test function when using FEM.
Type can be SCALAR, VECTOR, TENSOR, or SYM_TENSOR.
"""
function testSymbol(symbol; type=SCALAR)
    add_test_function(finch_state, Symbol(symbol), type);
end

"""
    index(name; range=[1])

Create an indexer entity to index variables and coefficients with a
VAR_ARRAY type. The range can be an array of integers or an array
with the min and max values. Though not strictly necessary, the range
should start at 1 and increase consecutively. Other configurations
may cause some issues in the generated code.

range=[1,2,3,4,5] is the same as range=[1,5]
"""
function index(name; range=[1])
    if length(range) == 2
        range = Vector(Int(range[1]):Int(range[2]));
    end
    idx = Indexer(Symbol(name), range, range[1], length(finch_state.indexers)+1)
    add_indexer(finch_state, idx);
    return idx;
end

"""
    variableTransform(var1, var2, func)

Define a function that will transform var1 into var2.
var1 and var2 can be variable entities or arrays of them.
The function should transform numerical values from one to the other.
"""
function variableTransform(var1, var2, func)
    # Make sure things are valid
    if typeof(var1) <: Array
        if !(typeof(var2) <: Array)
            return add_variable_transform(finch_state, var1, [var2], func);
        else
            return add_variable_transform(finch_state, var1, var2, func);
        end
    else
        if typeof(var2) <: Array
            return add_variable_transform(finch_state, [var1], var2, func);
        else
            return add_variable_transform(finch_state, var1, var2, func);
        end
    end
end

"""
    transformVariable(xform::VariableTransform)

Perform a variable transformation using a previously defined transform.
This can be done, for example, in a postStepFunction or in a callback
function for boundary conditions.
"""
function transformVariable(xform::VariableTransform)
    transform_variable_values(xform);
end
function transformVariable(var1, var2)
    # Look for a matching xform
    found = false;
    for i=1:length(finch_state.variable_transforms)
        if var1 == finch_state.variable_transforms[i].from || [var1] == finch_state.variable_transforms[i].from
            if var2 == finch_state.variable_transforms[i].to || [var2] == finch_state.variable_transforms[i].to
                found = true;
                transform_variable_values(variable_transforms[i]);
                break;
            end
        end
    end
end

"""
    boundary(var, bid, bc_type, bc_exp=0)

Set a boundary condition for a given variable on a boundary region with
this ID(bid). The type can be DIRICHLET, NEUMANN, NO_BC, FLUX.
bc_exp can be a constant number or a string expression. NO_BC type
does not need a bc_exp value.

Possible expressions can include coordinates(x, y, z, t), callback functions,
coefficients, variables, indexers, and certain other symbols such as "normal",
"node_index", "face_index".
"""
function boundary(var, bid, bc_type, bc_exp=0)
    # The expression may contain variable symbols.
    # Parse it, replace variable symbols with the appropriate parts
    # Type may change in bc_exp, so convert it to an array of Any
    newbc_exp = [];
    append!(newbc_exp, bc_exp);
    if typeof(bc_exp) <: Array
        nfuns = 0;
        for i=1:length(bc_exp)
            if typeof(bc_exp[i]) == String
                ex = Meta.parse(bc_exp[i]);
                ex = replace_symbols_in_conditions(ex);
                newbc_exp[i] = string(ex);
                nfuns += makeFunctions(newbc_exp[i]);
            elseif typeof(bc_exp[i]) <: Number
                newbc_exp[i] = Float64(bc_exp[i]);
            else
                # What else could it be?
            end
        end
    elseif typeof(bc_exp) == String
        ex = Meta.parse(bc_exp);
        ex = replace_symbols_in_conditions(ex);
        newbc_exp = string(ex);
        nfuns = makeFunctions(newbc_exp);
    elseif typeof(bc_exp) <: Number
        newbc_exp = Float64(bc_exp);
        nfuns = 0;
    else
        # ??
    end
    
    add_boundary_condition(finch_state, var, bid, bc_type, newbc_exp, nfuns);
end

"""
    addBoundaryID(bid::Int, trueOnBdry)

Create a new boundary region with the given ID number bid. It will be
assigned to all boundary faces where the center of the face satisfies
trueOnBdry. It will override any previously set ID for those faces.

trueOnBdry can be a function or a string expression of (x,y,z).
Note that it only applies to faces that are known boundary faces, not
interior faces.
"""
function addBoundaryID(bid::Int, trueOnBdry)
    # trueOnBdry(x, y, z) = something # points with x,y,z on this bdry segment evaluate true here
    if typeof(trueOnBdry) == String
        trueOnBdrystr = trueOnBdry;
        trueOnBdryfun = stringToFunction("trueOnBdry", "x,y=0,z=0", trueOnBdry);
    else
        trueOnBdryfun = trueOnBdry;
        trueOnBdrystr = "CUSTOM";
    end
    add_boundary_ID_to_grid(bid, trueOnBdryfun, finch_state.grid_data);
    add_boundary_ID_to_problem(bid, trueOnBdrystr);
end

"""
    referencePoint(var, pos, val)

Constrain a variable value at a single node to the given value.
This is needed for certain situations where boundary conditions do not
uniquely constrain a variable. The node closest to the position in pos
will be used and it can be a boundary or interior node. 
"""
function referencePoint(var, pos, val)
    add_reference_point(finch_state, var, pos, val);
end

"""
    initial(var, ics)

Set the initial condition for this variable. The value can be a constant
number or a string expression of coordinates(x, y, z).
This does not immediately set the variable values. To do so, use
evalInitialConditions(), or it will be done automatically before solving.
"""
function initial(var, value)
    nfuns = makeFunctions(value);
    add_initial_condition(finch_state, var.index, value, nfuns);
end

"""
    preStepFunction(fun)

Set a function to be called before each time step, or stage for multi-stage
steppers.
"""
function preStepFunction(fun)
    finch_state.prob.pre_step_function = fun;
end

"""
    postStepFunction(fun)

Set a function to be called after each time step, or stage for multi-stage
steppers.
"""
function postStepFunction(fun)
    finch_state.prob.post_step_function = fun;
end

"""
    callbackFunction(fun; name="", args=[], body="")

Include a callback function that can be included in expressions such as
the PDE or boundary conditions.
A better way to do this is with the macro @callbackFunction before the
function definition, which will extract the name, arguments, and body 
automatically.
"""
function callbackFunction(fun; name="", args=[], body="")
    if name==""
        name = string(fun);
    end
    # If args and body are not provided, this may still work internally
    # but it can't be generated for external targets.
    
    add_callback_function(finch_state, CallbackFunction(name, args, body, fun));
    
    log_entry("Added callback function: "*name, 2);
end

"""
    weakForm(var, wf)

Write the weak form of the PDE in residual form. This should be an expression
that is assumed to be equal to zero.
var can be a variable or an array of variables. When using arrays, wf must
also be an array of matching size.
wf is a string expression or array of them. It can include numbers, coefficients,
variables, parameters, indexers, and symbolic operators.
"""
function weakForm(var, wf)
    
    @timeit finch_state.timer_output "CodeGen" begin
    
    if !(typeof(var) <: Array)
        var = [var];
        wf = [wf];
    end
    wfvars = [];
    wfex = [];
    if !(length(var) == length(wf))
        printerr("Error in weak form: # of unknowns must equal # of equations. (example: weakform([a,b,c], [f1,f2,f3]))");
    end
    for vi=1:length(var)
        push!(wfvars, var[vi].symbol);
        push!(wfex, Meta.parse((wf)[vi]));
    end
    solver_type = var[1].discretization;
    
    # If nonlinear iteration is used, need to make copies of variables for DELTA versions
    # if prob.nonlinear
    #     deltavar = [];
    #     for i=1:length(var)
    #         deltaname = "DELTA"*string(var[i].symbol);
    #         push!(deltavar, variable(deltaname, type=var[i].type, location=var[i].location, method=var[i].discretization, index=var[i].indexer));
    #         # need to add zero dirichlet BC for every bid
    #         nbids = size(prob.bid, 2);
    #         for bid=1:nbids
    #             boundary(deltavar[end], prob.bid[1,bid], DIRICHLET, 0);
    #         end
    #     end
        
    # else
    #     deltavar = var;
    # end
    if finch_state.prob.nonlinear
        for i=1:length(var)
            oldname = "OLD"*string(var[i].symbol);
            variable(oldname, type=var[i].type, location=var[i].location, method=var[i].discretization, index=var[i].indexer);
        end
    end
    
    log_entry("Making weak form for variable(s): "*string(wfvars));
    log_entry("Weak form, input: "*string(wf));
    
    # This is the parsing step. It goes from an Expr to arrays of Basic
    result_exprs = sp_parse(wfex, wfvars);
    if length(result_exprs) == 4 # has surface terms
        (lhs_symexpr, rhs_symexpr, lhs_surf_symexpr, rhs_surf_symexpr) = result_exprs;
    else
        (lhs_symexpr, rhs_symexpr) = result_exprs;
        lhs_surf_symexpr = nothing;
        rhs_surf_symexpr = nothing;
    end
    
    # Here we set a SymExpression for each of the pieces. 
    # This is an Expr tree that is passed to the code generator.
    if length(result_exprs) == 4 # has surface terms
        set_symexpressions(finch_state, var, lhs_symexpr, LHS, "volume");
        set_symexpressions(finch_state, var, lhs_surf_symexpr, LHS, "surface");
        set_symexpressions(finch_state, var, rhs_symexpr, RHS, "volume");
        set_symexpressions(finch_state, var, rhs_surf_symexpr, RHS, "surface");
        
        log_entry("lhs volume symexpression:\n\t"*string(lhs_symexpr));
        log_entry("lhs surface symexpression:\n\t"*string(lhs_surf_symexpr));
        log_entry("rhs volume symexpression:\n\t"*string(rhs_symexpr));
        log_entry("rhs surface symexpression:\n\t"*string(rhs_surf_symexpr));
        
        log_entry("Latex equation:\n\t\t \\int_{K}"*symexpression_to_latex(lhs_symexpr)*
                    " dx + \\int_{\\partial K}"*symexpression_to_latex(lhs_surf_symexpr)*
                    " ds = \\int_{K}"*symexpression_to_latex(rhs_symexpr)*
                    " dx + \\int_{\\partial K}"*symexpression_to_latex(rhs_surf_symexpr)*" ds", 3);
    else
        set_symexpressions(finch_state, var, lhs_symexpr, LHS, "volume");
        set_symexpressions(finch_state, var, rhs_symexpr, RHS, "volume");
        
        log_entry("lhs symexpression:\n\t"*string(lhs_symexpr));
        log_entry("rhs symexpression:\n\t"*string(rhs_symexpr));
        
        log_entry("Latex equation:\n\t\t\$ \\int_{K}"*symexpression_to_latex(lhs_symexpr)*
                    " dx = \\int_{K}"*symexpression_to_latex(rhs_symexpr)*" dx");
    end
    
    # Init stepper here so it can be used in generation
    if finch_state.prob.time_dependent
        if length(finch_state.grid_data.allnodes) > 0
            # some measure of element size
            dim = finch_state.config.dimension;
            min_detj = minimum(finch_state.geo_factors.detJ);
            el_size = (2^dim * min_detj)^(1/dim);
            init_stepper(el_size, finch_state.time_stepper);
            if finch_state.use_specified_steps
                finch_state.time_stepper.dt = finch_state.specified_dt;
                finch_state.time_stepper.Nsteps = finch_state.specified_Nsteps;
            end
            share_time_step_info(finch_state.time_stepper, finch_state.config);
        else # no mesh was built for this target
            init_stepper(0.1, finch_state.time_stepper); # This should be manually set.
            if finch_state.use_specified_steps
                finch_state.time_stepper.dt = finch_state.specified_dt;
                finch_state.time_stepper.Nsteps = finch_state.specified_Nsteps;
            end
        end
    end
    
    # change symbolic layer into IR
    full_IR = build_IR_fem(lhs_symexpr, lhs_surf_symexpr, rhs_symexpr, rhs_surf_symexpr, 
                            var, finch_state.ordered_indexers, finch_state.config, finch_state.prob, finch_state.time_stepper);
    log_entry("Weak form IR: \n"*repr_IR(full_IR),3);
    
    # Generate code from IR
    code = generate_code_layer(var, full_IR, solver_type, finch_state.target_language, finch_state.target_framework);
    log_entry("Code layer: \n" * code, 3);
    set_code(finch_state, var, code, full_IR);
    
    end # timer block
end

"""
    conservationForm(var, cf)

Write the integral conservation form of the PDE. This should be an expression
that is assumed to be equal to the time derivative of the variable.
Surface integrals for the flux are wrapped in surface(), and all other terms are
assumed to be volume integrals for the source. Do not include the time derivative
as it is implicitly assumed.
var can be a variable or an array of variables. When using arrays, cf must
also be an array of matching size.
cf is a string expression or array of them. It can include numbers, coefficients,
variables, parameters, indexers, and symbolic operators.
"""
function conservationForm(var, cf)
    
    @timeit finch_state.timer_output "CodeGen" begin
    
    if !(typeof(var) <: Array)
        var = [var];
        cf = [cf];
    end
    cfvars = [];
    cfex = [];
    if !(length(var) == length(cf))
        printerr("Error in conservation form: # of unknowns must equal # of equations. (example: conservationform([a,b,c], [f1,f2,f3]))");
    end
    for vi=1:length(var)
        push!(cfvars, var[vi].symbol);
        push!(cfex, Meta.parse((cf)[vi]));
    end
    solver_type = var[1].discretization;
    
    log_entry("Making conservation form for variable(s): "*string(cfvars));
    log_entry("Conservation form, input: "*string(cf));
    
    # This is the parsing step. It goes from an Expr to arrays of Basic
    result_exprs = sp_parse(cfex, cfvars, is_FV=true);
    if length(result_exprs) == 4 # has surface terms
        (lhs_symexpr, rhs_symexpr, lhs_surf_symexpr, rhs_surf_symexpr) = result_exprs;
    else
        (lhs_symexpr, rhs_symexpr) = result_exprs;
    end
    
    # Here we set a SymExpression for each of the pieces. 
    # This is an Expr tree that is passed to the code generator.
    if length(result_exprs) == 4 # has surface terms
        set_symexpressions(finch_state, var, lhs_symexpr, LHS, "volume");
        set_symexpressions(finch_state, var, lhs_surf_symexpr, LHS, "surface");
        set_symexpressions(finch_state, var, rhs_symexpr, RHS, "volume");
        set_symexpressions(finch_state, var, rhs_surf_symexpr, RHS, "surface");
        
        log_entry("lhs volume symexpression:\n\t"*string(lhs_symexpr));
        log_entry("lhs surface symexpression:\n\t"*string(lhs_surf_symexpr));
        log_entry("rhs volume symexpression:\n\t"*string(rhs_symexpr));
        log_entry("rhs surface symexpression:\n\t"*string(rhs_surf_symexpr));
        
        log_entry("Latex equation:\n\t\t \\int_{K} -"*symexpression_to_latex(lhs_symexpr)*
                    " dx + \\int_{K}"*symexpression_to_latex(rhs_symexpr)*" dx"*
                    "\\int_{\\partial K}"*symexpression_to_latex(lhs_surf_symexpr)*
                    " ds - \\int_{\\partial K}"*symexpression_to_latex(rhs_surf_symexpr)*" ds", 3);
    else
        set_symexpressions(finch_state, var, lhs_symexpr, LHS, "volume");
        set_symexpressions(finch_state, var, rhs_symexpr, RHS, "volume");
        
        log_entry("lhs symexpression:\n\t"*string(lhs_symexpr));
        log_entry("rhs symexpression:\n\t"*string(rhs_symexpr));
        
        log_entry("Latex equation:\n\t\t\$ \\int_{K}"*symexpression_to_latex(lhs_symexpr)*
                    " dx = \\int_{K}"*symexpression_to_latex(rhs_symexpr)*" dx", 3);
    end
    
    # Init stepper here so it can be used in generation
    if finch_state.prob.time_dependent
        if length(finch_state.grid_data.allnodes) > 0
            # some measure of element size
            dim = finch_state.config.dimension;
            min_detj = minimum(finch_state.geo_factors.detJ);
            el_size = (2^dim * min_detj)^(1/dim);
            init_stepper(el_size, finch_state.time_stepper);
            if finch_state.use_specified_steps
                finch_state.time_stepper.dt = finch_state.specified_dt;
                finch_state.time_stepper.Nsteps = finch_state.specified_Nsteps;
            end
            share_time_step_info(finch_state.time_stepper, finch_state.config);
        else # no mesh was built for this target
            init_stepper(0.1, finch_state.time_stepper); # This should be manually set.
            if finch_state.use_specified_steps
                finch_state.time_stepper.dt = finch_state.specified_dt;
                finch_state.time_stepper.Nsteps = finch_state.specified_Nsteps;
            end
        end
    end
    
    # change symbolic layer into IR
    if length(result_exprs) == 4
        full_IR = build_IR_fvm(lhs_symexpr, lhs_surf_symexpr, rhs_symexpr, rhs_surf_symexpr, var, finch_state.ordered_indexers, finch_state.config, finch_state.prob, finch_state.time_stepper, finch_state.fv_info);
    else
        full_IR = build_IR_fvm(lhs_symexpr, nothing, rhs_symexpr, nothing, var, finch_state.ordered_indexers, finch_state.config, finch_state.prob, finch_state.time_stepper, finch_state.fv_info);
    end
    log_entry("Conservation form IR: "*repr_IR(full_IR),3);
    
    # Generate code from IR
    code = generate_code_layer(var, full_IR, solver_type, finch_state.target_language, finch_state.target_framework);
    log_entry("Code layer: \n" * code, 3);
    set_code(finch_state, var, code, full_IR);
    
    end # timer block
end

"""
    assemblyLoops(indices)

Specify the nesting order of loops for assembling the system.
This makes sense for problems with indexed variables.
Var is a variable of array of variables.
Indices is an array of indexer objects and a string "elements".
Loops will be nested in that order with outermost first.
"""
function assemblyLoops(indices)
    # find the associated indexer objects
    # and reorder ordered_indexers
    indexer_list = Vector{Indexer}(undef,0);
    for i=1:length(indices)
        if typeof(indices[i]) == String
            if indices[i] == "elements" || indices[i] == "cells"
                push!(indexer_list, Indexer(:FINCHELEMENTS, [0,0], 0, 0));
            else
                for j=1:length(finch_state.indexers)
                    if indices[i] == string(finch_state.indexers[j].symbol)
                        push!(indexer_list, finch_state.indexers[j]);
                        break;
                    end
                end
            end
            
        elseif typeof(indices[i]) == Indexer
            push!(indexer_list, indices[i]);
        end
    end
    set_ordered_indexers(finch_state, indexer_list);
    return nothing;
end

"""
    exportCode(filename)

Export all generated code including the elemental calculation and assembly loop
code for all variables if available.

"""
function exportCode(filename)
    # For now, only do this for Julia code because others are already output in code files.
    if finch_state.target_language == JULIA
        file = open(filename*".jl", "w");
        println(file, "#=\nGenerated functions for "*finch_state.project_name*"\n=#\n");
        
        for i=1:length(finch_state.variables)
            code = finch_state.code_strings[i];
            var = string(finch_state.variables[i].symbol);
            if length(code) > 1
                println(file, "# begin solve function for "*var*"\n");
                println(file, code);
                println(file, "# end solve function for "*var*"\n");
            else
                println(file, "# No code set for "*var*"\n");
            end
        end
        
        close(file);
        log_entry("Exported code to "*filename*".jl");
    else
        # Should we export for other targets?
    end
end

"""
    importCode(filename)

Import elemental calculation and assembly loop code for all variables if possible.
The "begin", "end", and "No code" comment lines must match a specific format to 
properly match them to the variables, so do not modify those lines from the
exported code.
"""
function importCode(filename)
    # For now, only do this for Julia code because others are already output in code files.
    if finch_state.target_language == JULIA
        file = open(filename*".jl", "r");
        lines = readlines(file, keep=true);
        
        # Loop over variables and check to see if a matching function is present.
        for i=1:length(finch_state.variables)
            # Scan the file for a pattern like
            #   # begin solve function for u
            #       ...
            #   # end solve function for u
            # OR
            #   # No code set for u
            var = string(finch_state.variables[i].symbol);
            func_begin_flag = "# begin solve function for "*var;
            func_end_flag = "# end solve function for "*var;
            func_none_flag = "# No code set for "*var;
            func_string = "";
            var_found = false;
            for st=1:length(lines)
                if occursin(func_begin_flag, lines[st])
                    # st+1 is the start of the function
                    for en=(st+1):length(lines)
                        if occursin(func_end_flag, lines[en])
                            # en-1 is the end of the function
                            st = en; # update st
                            break;
                        else
                            func_string *= lines[en];
                        end
                    end
                    var_found = true;
                elseif occursin(func_none_flag, lines[st])
                    # There is no function for this variable
                    log_entry("While importing, no solve function was found for "*var);
                    var_found = true;
                end
                if var_found
                    break;
                end
            end # lines loop
            
            # Generate the functions and set them in the right places
            if func_string == ""
                log_entry("While importing, nothing was found for "*var);
            else
                set_code(finch_state, finch_state.variables[i], func_string, IR_comment_node("Code imported from file, no IR available"));
            end
            
        end # vars loop
        
        close(file);
        log_entry("Imported code from "*filename*".jl");
        
    else
        # external code doesn't need to be imported
    end
end

"""
    printLatex(var)

Print a string of Latex formatted code for the symbolic layer form of the PDE.
This is somewhat limited and needs to be updated to a more useful output.
"""
function printLatex(var)
    if typeof(var) <: Array
        varname = "["*string(var[1].symbol);
        for i=2:length(var)
            varname *= ", "*string(var[i].symbol);
        end
        varname *= "]";
        var = var[1];
    else
        varname = string(var.symbol);
    end
    println("Latex equation for "*varname*":");
    
    symexpressions = finch_state.symexpressions;
    
    if finch_state.config.solver_type == FV
        result = raw"$" * " \\frac{d}{dt}\\bar{"*string(var.symbol)*"}";
        if !(symexpressions[1][var.index] === nothing)
            result *= " + \\int_{K}"*symexpression_to_latex(symexpressions[1][var.index])*" dx";
        end
        if !(symexpressions[2][var.index] === nothing)
            result *= " + \\int_{\\partial K}"*symexpression_to_latex(symexpressions[2][var.index])*" ds";
        end
        result *= " = ";
        if !(symexpressions[3][var.index] === nothing)
            result *= "\\int_{K}"*symexpression_to_latex(symexpressions[3][var.index])*" dx";
        end
        if !(symexpressions[4][var.index] === nothing)
            if !(symexpressions[3][var.index] === nothing)
                result *= " + ";
            end
            result *= "\\int_{\\partial K}"*symexpression_to_latex(symexpressions[4][var.index])*" ds";
        end
        if symexpressions[3][var.index] === nothing && symexpressions[4][var.index] === nothing
            result *= "0";
        end
        result *= raw"$";
        
    else # FEM
        if finch_state.prob.time_dependent
            result = raw"$ \left[";
            if !(symexpressions[1][var.index] === nothing)
                result *= "\\int_{K}"*symexpression_to_latex(symexpressions[1][var.index])*" dx";
            end
            if !(symexpressions[2][var.index] === nothing)
                if !(symexpressions[1][var.index] === nothing)
                    result *= " + ";
                end
                result *= "\\int_{\\partial K}"*symexpression_to_latex(symexpressions[2][var.index])*" ds";
            end
            if symexpressions[1][var.index] === nothing && symexpressions[2][var.index] === nothing
                result *= "0";
            end
            result *= "\\right]_{t+dt} = \\left[";
            if !(symexpressions[3][var.index] === nothing)
                result *= "\\int_{K}"*symexpression_to_latex(symexpressions[3][var.index])*" dx";
            end
            if !(symexpressions[4][var.index] === nothing)
                if !(symexpressions[3][var.index] === nothing)
                    result *= " + ";
                end
                result *= "\\int_{\\partial K}"*symexpression_to_latex(symexpressions[4][var.index])*" ds";
            end
            if symexpressions[3][var.index] === nothing && symexpressions[4][var.index] === nothing
                result *= "0";
            end
            result *= raw"\right]_{t} $";
            
        else # not time dependent
            result = raw"$";
            if !(symexpressions[1][var.index] === nothing)
                result *= "\\int_{K}"*symexpression_to_latex(symexpressions[1][var.index])*" dx";
            end
            if !(symexpressions[2][var.index] === nothing)
                if !(symexpressions[1][var.index] === nothing)
                    result *= " + ";
                end
                result *= "\\int_{\\partial K}"*symexpression_to_latex(symexpressions[2][var.index])*" ds";
            end
            if symexpressions[1][var.index] === nothing && symexpressions[2][var.index] === nothing
                result *= "0";
            end
            result *= " = ";
            if !(symexpressions[3][var.index] === nothing)
                result *= "\\int_{K}"*symexpression_to_latex(symexpressions[3][var.index])*" dx";
            end
            if !(symexpressions[4][var.index] === nothing)
                if !(symexpressions[3][var.index] === nothing)
                    result *= " + ";
                end
                result *= "\\int_{\\partial K}"*symexpression_to_latex(symexpressions[4][var.index])*" ds";
            end
            if symexpressions[3][var.index] === nothing && symexpressions[4][var.index] === nothing
                result *= "0";
            end
            result *= raw"$";
        end
    end
    println(result);
    println("");
    
    return result;
end

"""
    evalInitialConditions()

Evaluate initial conditions for all variables if possible.
This puts the initial values into each variable's values array.
This is called automatically by the solve step, but can be done manually here.
"""
function evalInitialConditions() 
    eval_initial_conditions(finch_state); 
end # Just for consistent style because this is also an internal function

"""
    nonlinear(;maxIters=100, relativeTol=1e-5, absoluteTol=1e-5)

Use to signal that this problem contains nonlinearity and needs an
iterative solution. Set the parameters for the iteration.
The tolerances are for the infinity norm of the change between 
iterations (u(i) - u(i-1)).
"""
function nonlinear(;maxIters=100, relativeTol=1e-5, absoluteTol=1e-5, relaxation=1, derivative="AD")
    finch_state.prob.nonlinear = true;
    finch_state.prob.derivative_type = derivative;
    finch_state.prob.max_iters = maxIters;
    finch_state.prob.relative_tol = relativeTol;
    finch_state.prob.absolute_tol = absoluteTol;
    finch_state.prob.relaxation = relaxation;
end

"""
    solve(var)

Either solve the problem using an internal target, or generate all code
files for an external target.
Var is a variable or array of variables to solve for.
The keyword arguments are for nonlinear equations which have very limited
support, so generally they won't be used.
"""
function solve(var)
    # if finch_state.use_cachesim
    #     return cachesimSolve(var);
    # end
    
    # Ensure an array of variables from finch_state
    
    if typeof(var) <: Variable
        var = [var];
    end
    nvars = length(var);
    vars = Vector{Variable{finch_state.config.float_type}}(undef,nvars);
    for i=1:nvars
        vars[i] = finch_state.variables[var[i].index];
    end
    
    dofs_per_node = 0;
    dofs_per_loop = 0;
    for vi=1:length(var)
        dofs_per_loop += length(var[vi].symvar);
        dofs_per_node += var[vi].total_components;
    end
    
    # Wrap this all in a timer
    @timeit finch_state.timer_output "Solve" begin
    
    # Generate files or solve directly
    if !(!finch_state.external_target && finch_state.target_language == JULIA) # if an external code gen target is ready
        varind = var[1].index;
        generate_all_files(var, finch_state.solve_functions[varind]);
        
    else
        varnames = "["*string(var[1].symbol);
        for vi=2:length(var)
            varnames = varnames*", "*string(var[vi].symbol);
        end
        varnames = varnames*"]";
        varind = var[1].index;
        
        # Evaluate initial conditions if not already done
        eval_initial_conditions(finch_state);
        
        # If any of the variables are indexed types, make sure ordered indexers has "elements"
        need_indexers = false;
        for vi=1:length(var)
            if !(var[vi].indexer === nothing)
                need_indexers = true;
            end
        end
        if need_indexers && length(finch_state.ordered_indexers) == length(finch_state.indexers)
            log_entry("Indexed variables detected, but no assembly loops specified. Using default.")
            assemblyLoops(["elements"; finch_state.ordered_indexers]);
        end
        
        if finch_state.prob.time_dependent
            if finch_state.use_specified_steps
                finch_state.time_stepper.dt = finch_state.specified_dt;
                finch_state.time_stepper.Nsteps = finch_state.specified_Nsteps;
            else
                # some measure of element size
                dim = finch_state.config.dimension;
                min_detj = minimum(finch_state.geo_factors.detJ);
                el_size = (2^dim * min_detj)^(1/dim);
                init_stepper(el_size, finch_state.time_stepper);
            end
            
            share_time_step_info(finch_state.time_stepper, finch_state.config);
        end
        
        # # If nonlinear, a matching set of variables should have been created. find them
        # if prob.nonlinear
        #     deltavar = Vector{Variable}(undef, 0);
        #     for i=1:length(var)
        #         deltaname = Symbol("DELTA"*string(var[i].symbol));
        #         for v in variables
        #             if deltaname === v.symbol
        #                 push!(deltavar, v);
        #                 v.values .= var[i].values; # copy initial conditions
        #                 break;
        #             end
        #         end
        #     end
        # end
        if finch_state.prob.nonlinear
            nl_var = Vector{Variable{finch_state.config.float_type}}(undef, 0);
            for i=1:length(var)
                oldname = Symbol("OLD"*string(var[i].symbol));
                for v in finch_state.variables
                    if oldname === v.symbol
                        push!(nl_var, v);
                        v.values .= var[i].values; # copy initial conditions
                        break;
                    end
                end
            end
        end
        
        # Use the appropriate solver
        if finch_state.config.solver_type == CG || finch_state.config.solver_type == DG
            func = finch_state.solve_functions[varind];
            
            if finch_state.prob.nonlinear
                @timeit finch_state.timer_output "FE_solve" func.func(vars, finch_state.grid_data, finch_state.refel, 
                                            finch_state.geo_factors, finch_state.config, finch_state.coefficients, 
                                            finch_state.variables, finch_state.test_functions, finch_state.ordered_indexers, 
                                            finch_state.prob, finch_state.time_stepper, finch_state.parallel_buffers,
                                            finch_state.timer_output, nl_var);
                
            else
                @timeit finch_state.timer_output "FE_solve" func.func(vars, finch_state.grid_data, finch_state.refel, 
                                            finch_state.geo_factors, finch_state.config, finch_state.coefficients, 
                                            finch_state.variables, finch_state.test_functions, finch_state.ordered_indexers, 
                                            finch_state.prob, finch_state.time_stepper, finch_state.parallel_buffers,
                                            finch_state.timer_output);
            end
            
            log_entry("Solved for "*varnames, 1);
            
        elseif finch_state.config.solver_type == FV
            func = finch_state.solve_functions[varind];
            
            if finch_state.prob.nonlinear
                @timeit finch_state.timer_output "FV_solve" func.func(vars, finch_state.fv_grid, finch_state.fv_refel, 
                                            finch_state.fv_geo_factors, finch_state.fv_info, finch_state.config, 
                                            finch_state.coefficients, finch_state.variables, finch_state.test_functions, 
                                            finch_state.ordered_indexers, finch_state.prob, finch_state.time_stepper, 
                                            finch_state.parallel_buffers, finch_state.timer_output, nl_var);
            else
                @timeit finch_state.timer_output "FV_solve" func.func(vars, finch_state.fv_grid, finch_state.fv_refel, 
                                            finch_state.fv_geo_factors, finch_state.fv_info, finch_state.config, 
                                            finch_state.coefficients, finch_state.variables, finch_state.test_functions, 
                                            finch_state.ordered_indexers, finch_state.prob, finch_state.time_stepper, 
                                            finch_state.parallel_buffers, finch_state.timer_output);
            end
            
            log_entry("Solved for "*varnames, 1);
            
        elseif finch_state.config.solver_type == MIXED
            println("Mixed solver is not ready. Please wait.");
            return nothing
            
            # # Need to determine the variable index for fe and fv
            # fe_var_index = 0;
            # fv_var_index = 0;
            # if typeof(var) <: Array
            #     for vi=1:length(var)
            #         if var[vi].discretization == FV && fv_var_index == 0
            #             fv_var_index += var[vi].index;
            #         elseif fe_var_index == 0;
            #             fe_var_index += var[vi].index;
            #         end
            #     end
            #     loop_func = [assembly_loops[fe_var_index], assembly_loops[fv_var_index]];
            # else
            #     # This shouldn't happen?
            #     printerr("mixed solver types specified, but only one variable being solved for?")
            #     if var.discretization == FV
            #         fv_var_index = var.index;
            #         loop_func = assembly_loops[fv_var_index];
            #     else
            #         fe_var_index = var.index;
            #         loop_func = assembly_loops[fe_var_index];
            #     end
            # end
            
        end
    end
    
    # At this point all of the final values for the variables in var should be placed in var.values.
    end # timer block
    
end

"""
    cachesimSolve(var)

When using the cache simulator target, this is used instead of solve().
"""
function cachesimSolve(var)
    if !(!finch_state.external_target && finch_state.target_language == JULIA)
        printerr("Cachesim solve is only ready for Julia target");
        
    else
        if typeof(var) <: Array
            varnames = "["*string(var[1].symbol);
            for vi=2:length(var)
                varnames = varnames*", "*string(var[vi].symbol);
            end
            varnames = varnames*"]";
            varind = var[1].index;
        else
            varnames = string(var.symbol);
            varind = var.index;
        end
        
        #TODO
        
        log_entry("Generated cachesim ouput for "*varnames*".(took "*string(t)*" seconds)", 1);
    end
end

"""
    outputValues(vars, filename; format="vtk", ascii=false)

Output variable values to a file in a spicified format.
vars can be a variable or array of variables.
Possible formats are "vtk", "csv", or "raw"
Set ascii to true to make ascii type vtk files instead of binary.
"""
function outputValues(vars, filename; format="vtk", ascii=false)
    available_formats = ["raw", "csv", "vtk", "try"];
    if format == "vtk"
        log_entry("Writing values to file: "*filename*".vtu");
        output_values_vtk(vars, filename, ascii)
        
    elseif format == "try"
        log_entry("Writing values to file: "*filename*".vtu");
        filename *= ".vtu";
        file = open(filename, "w");
        log_entry("Writing values to file: "*filename);
        output_values_myvtu(vars, file, ascii);
        
        close(file);
        
    else
        filename *= "."*format;
        file = open(filename, "w");
        log_entry("Writing values to file: "*filename);
        
        if format == "raw"
            output_values_raw(vars, file)
        elseif format == "csv"
            output_values_csv(vars, file)
        else
            println("Unknown output file format("*string(format)*"). Choose from: "*string(available_formats));
        end
        log_entry("Finished writing to "*filename);
        
        close(file);
    end
end

"""
finalizeFinch()

This closes all files and does any other finalization steps.
It does not deallocate any data, so further processing can be done after 
calling this. However, using Finch functions may cause issues because
files have been closed.
"""
function finalizeFinch()
    # Finalize generation
    finalize_code_generator();
    if finch_state.use_cachesim
        CachesimOut.finalize();
    end
    
    # mpi
    if finch_state.config.use_mpi
        MPI.Finalize(); # If MPI is finalized, it cannot be reinitialized without restarting Julia
    end
    
    if finch_state.config.proc_rank == 0
        # Print timer results to log and output
        buf = IOBuffer();
        print_timer(buf, finch_state.timer_output);
        log_entry(String(take!(buf)), 1);
        
        println("");
        show(finch_state.timer_output);
        println("");
        
        close_log(finch_state);
        println("Finch has completed.");
    end
end

### Other specialized functions ###

"""
    cachesim(use)

Toggle the cache simulator target. Set use=true to use the cachesim target.
If not using cachesim, don't use this function.
"""
function cachesim(use)
    log_entry("Using cachesim - Only cachesim output will be generated.", 1);
    finch_state.use_cachesim = use;
end

"""
    mortonNodes(griddim)

Reorder the nodes in memory to a Morton ordering in 2D or 3D.
This currently only works for a uniform grid such as the one generated
with Finch's internal utility.
griddim is an array representing the nodal grid size: like [n,n] for 2D or
[n,n,n] for 3D.
"""
function mortonNodes(griddim)
    t = @elapsed(finch_state.grid_data = reorder_grid_recursive!(finch_state.grid_data, griddim, "morton"));
    log_entry("Reordered nodes to Morton. Took "*string(t)*" sec.", 2);
end

"""
    mortonElements(griddim)

Reorder the elemental loop order to a spacial Morton ordering in 2D or 3D.
This currently only works for a uniform grid such as the one generated
with Finch's internal utility.
griddim is an array representing the elemental grid size: like [n,n] for 2D or
[n,n,n] for 3D.
"""
function mortonElements(griddim)
    finch_state.grid_data.elemental_order = get_recursive_order("morton", finch_state.config.dimension, griddim);
    log_entry("Reordered elements to Morton.", 2);
    ef_nodes();
end

"""
    hilbertNodes(griddim)

Reorder the nodes in memory to a Hilbert ordering in 2D or 3D.
This currently only works for a uniform grid such as the one generated
with Finch's internal utility.
griddim is an array representing the nodal grid size: like [n,n] for 2D or
[n,n,n] for 3D.
"""
function hilbertNodes(griddim)
    t = @elapsed(finch_state.grid_data = reorder_grid_recursive!(finch_state.grid_data, griddim, "hilbert"));
    log_entry("Reordered nodes to Hilbert. Took "*string(t)*" sec.", 2);
end

"""
    hilbertElements(griddim)

Reorder the elemental loop order to a spacial Hilbert ordering in 2D or 3D.
This currently only works for a uniform grid such as the one generated
with Finch's internal utility.
griddim is an array representing the elemental grid size: like [n,n] for 2D or
[n,n,n] for 3D.
"""
function hilbertElements(griddim)
    finch_state.grid_data.elemental_order = get_recursive_order("hilbert", finch_state.config.dimension, griddim);
    log_entry("Reordered elements to Hilbert.", 2);
    ef_nodes();
end

"""
    tiledNodes(griddim, tiledim)

Reorder the nodes in memory to a tiled ordering in 2D or 3D.
This currently only works for a uniform grid such as the one generated
with Finch's internal utility.
griddim is an array representing the nodal grid size: like [n,n] for 2D or
[n,n,n] for 3D.
tiledim is the desired tile dimensions such as [4,4] for a 4x4 tile in 2D.
"""
function tiledNodes(griddim, tiledim)
    t = @elapsed(finch_state.grid_data = reorder_grid_tiled(finch_state.grid_data, griddim, tiledim));
    log_entry("Reordered nodes to tiled. Took "*string(t)*" sec.", 2);
end

"""
    tiledElements(griddim, tiledim)

Reorder the elemental loop order to a spacial tiled ordering in 2D or 3D.
This currently only works for a uniform grid such as the one generated
with Finch's internal utility.
griddim is an array representing the elemental grid size: like [n,n] for 2D or
[n,n,n] for 3D.
tiledim is the desired tile dimensions such as [4,4] for a 4x4 tile in 2D.
"""
function tiledElements(griddim, tiledim)
    finch_state.grid_data.elemental_order = get_tiled_order(finch_state.config.dimension, griddim, tiledim, true);
    log_entry("Reordered elements to tiled("*string(tiledim)*").", 2);
    ef_nodes();
end

"""
    elementFirstNodes()

This is the default node ordering. Element first means the elements are given
some order and the nodes are added elementwise according to that. An element's
nodes are ordered according to the reference element.
"""
function elementFirstNodes()
    t = @elapsed(finch_state.grid_data = reorder_grid_element_first!(finch_state.grid_data, finch_state.config.basis_order_min));
    log_entry("Reordered nodes to EF. Took "*string(t)*" sec.", 2);
end

"""
    randomNodes(seed = 17)

Randomize nodes in memory for testing a worst-case arrangement.
The seed is for making results reproducible.
"""
function randomNodes(seed = 17)
    t = @elapsed(finch_state.grid_data = reorder_grid_random!(finch_state.grid_data, seed));
    log_entry("Reordered nodes to random. Took "*string(t)*" sec.", 2);
end

"""
    randomElements(seed = 17)

Randomize the order of the elemental loop for testing a worst-case arrangement.
The seed is for making results reproducible.
"""
function randomElements(seed = 17)
    finch_state.grid_data.elemental_order = random_order(size(finch_state.grid_data.loc2glb,2), seed);
    log_entry("Reordered elements to random.", 2);
    random_nodes(seed);
end
