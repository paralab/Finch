#=
This file contains all of the common interface functions.
Many of them simply call corresponding functions in jl.
=#
export generateFor, useLog, domain, solverType, functionSpace, trialSpace, testSpace, finiteVolumeOrder,
        nodeType, timeStepper, setSteps, matrixFree, customOperator, customOperatorFile,
        mesh, exportMesh, variable, coefficient, parameter, testSymbol, index, boundary, addBoundaryID,
        referencePoint, timeInterval, initial, preStepFunction, postStepFunction, callbackFunction,
        variableTransform, transformVariable,
        weakForm, conservationForm, flux, source, assemblyLoops,
        exportCode, importCode, printLatex,
        evalInitialConditions, solve, cachesimSolve, finalizeFinch, cachesim, outputValues,
        mortonNodes, hilbertNodes, tiledNodes, mortonElements, hilbertElements, 
        tiledElements, elementFirstNodes, randomNodes, randomElements,
        # These do not match the interface style, but are kept for legacy support. May be removed.
        finalize_finch, morton_nodes, hilbert_nodes, tiled_nodes, morton_elements, hilbert_elements, 
        tiled_elements, ef_nodes, random_nodes, random_elements

# Begin configuration setting functions

"""
    generateFor(lang; filename=project_name, header="", params=nothing)

Specify the generation target. Lang could be one of the included target constants: 
(MATLAB, DENDRO) or the filename where the target is defined. The keyword argument
filename refers to the name to be applied to the generated code. The header text 
will be placed at the top of each generated code file. If the target requires 
some extra parameters, those are included in params.
"""
function generateFor(lang; filename=project_name, header="", params=nothing)
    outputDirPath = pwd()*"/"*filename;
    if config.proc_rank == 0 && !isdir(outputDirPath)
        mkdir(outputDirPath);
    end
    framew = 0;
    if !in(lang, [CPP,MATLAB,DENDRO,HOMG])
        # lang should be a filename for a custom target
        # This file must include these three functions:
        # 1. get_external_language_elements() - file extensions, comment chars etc.
        # 2. generate_external_code_layer(var, entities, terms, lorr, vors) - Turns symbolic expressions into code
        # 3. generate_external_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf) - Writes all files based on generated code
        include(lang);
        set_custom_gen_target(get_external_language_elements, generate_external_code_layer, generate_external_files, outputDirPath, filename, head=header);
    else
        # Use an included target
        target_dir = @__DIR__
        if lang == DENDRO
            framew = DENDRO;
            lang = CPP;
            target_file = "/targets/target_dendro_cg.jl";
        elseif lang == MATLAB
            framew = MATLAB;
            lang = MATLAB;
            target_file = "/targets/target_matlab_cg.jl";
        else
            framew = 0;
            target_file = "/targets/target_matlab_cg.jl";
        end
        include(target_dir * target_file);
        set_custom_gen_target(get_external_language_elements, generate_external_code_layer, generate_external_files, outputDirPath, filename, head=header);
    end
    if !(params === nothing)
        set_codegen_parameters(params);
    end
end

"""
    useLog(name=project_name; dir=output_dir, level=2)

Turn on logging with the given file name and optional directory.
The verbosity level can be 1(basic progress info), 2(More details about
each step), or 3(everything).
"""
function useLog(name=project_name; dir=output_dir, level=2)
    init_log(name, dir, level);
end

"""
    domain(dims; shape=SQUARE, grid=UNIFORM_GRID)

Set the dimensionality of the domain. The shape(SQUARE, IRREGULAR) and 
grid type(UNIFORM_GRID, UNSTRUCTURED, TREE) can be set, but may be changed
when building or importing the mesh.
"""
function domain(dims; shape=SQUARE, grid=UNIFORM_GRID)
    config.dimension = dims;
    config.geometry = shape;
    config.mesh_type = grid;
end

"""
    solverType(method, backend=DEFAULT_SOLVER)

Select between CG, DG, FV, or MIXED methods. The backend refers to the tools used
for solving linear systems: DEFAULT_SOLVER, PETSC_SOLVER, CUDA_SOLVER. 
Not all combinations of options are available yet.
"""
function solverType(method, backend=DEFAULT_SOLVER)
    set_solver(method, backend);
end

"""
    functionSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)

Set the polynomial order and type of polynomials for FEM.
Some of these are placeholders, so only use order at this point.
"""
function functionSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
    config.trial_function = space;
    config.test_function = space;
    if orderMax > orderMin && orderMin >= 0
        config.p_adaptive = true;
        config.basis_order_min = orderMin;
        config.basis_order_max = orderMax;
    else
        config.basis_order_min = max(order, orderMin);
        config.basis_order_max = max(order, orderMin);
    end
end

function trialSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
    #TODO
    functionSpace(space=space, order=order, orderMin=orderMin, orderMax=orderMax);
end

function testSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
    #TODO
    functionSpace(space=space, order=order, orderMin=orderMin, orderMax=orderMax);
end

"""
    nodeType(type)

For FEM, set the nodal configuration within elements.
The default is LOBATTO. GAUSS and UNIFORM are available, but should be used with care.
"""
function nodeType(type)
    config.elemental_nodes = type;
end

"""
    timeStepper(type; cfl=0)

Set the type of time stepping method and optionally the CFL number.
Options include EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICHOLSON, RK4, LSRK4, PECE.
If no CFL number is provided, one will be chosen based on the mesh and stepper type.
"""
function timeStepper(type; cfl=0)
    set_stepper(type, cfl);
end

"""
    timeInterval(T)

Set the ending time for time stepping. This is overridden if time steps
are manually specified.
"""
function timeInterval(T)
    prob.time_dependent = true;
    if time_stepper === nothing
        timeStepper(EULER_IMPLICIT);
    end
    prob.end_time = T;
end

"""
    setSteps(dt, steps)

Manually set the time steps if desired.
"""
function setSteps(dt, steps)
    prob.time_dependent = true;
    if time_stepper === nothing
        timeStepper(EULER_IMPLICIT);
    end
    set_specified_steps(dt, steps);
end

"""
    matrixFree(shallwe=true; maxiters=100, tol=1e-6)

Select a matrix free method for FEM with the given max iterations and tolerance.
This will use a basic conjugate gradient method.
"""
function matrixFree(shallwe=true; maxiters=100, tol=1e-6)
    config.linalg_matrixfree = shallwe;
    config.linalg_matfree_max = maxiters;
    config.linalg_matfree_tol = tol;
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
    
    @timeit timer_output "Mesh" begin
    
    if msh == LINEMESH
        log_entry("Building simple line mesh with nx elements, nx="*string(elsperdim));
        meshtime = @elapsed(mshdat = simple_line_mesh(elsperdim.+1, bids, interval));
        log_entry("Mesh building took "*string(meshtime)*" seconds");
        
    elseif msh == QUADMESH
        if length(interval) == 2
            interval = [interval[1], interval[2], interval[1], interval[2]];
        end
        log_entry("Building simple quad mesh with nx*nx elements, nx="*string(elsperdim));
        meshtime = @elapsed(mshdat = simple_quad_mesh(elsperdim.+1, bids, interval));
        log_entry("Mesh building took "*string(meshtime)*" seconds");
        
    elseif msh == TRIMESH
        if length(interval) == 2
            interval = [interval[1], interval[2], interval[1], interval[2]];
        end
        log_entry("Building simple triangle mesh with nx*nx*2 elements, nx="*string(elsperdim));
        meshtime = @elapsed(mshdat = simple_tri_mesh(elsperdim.+1, bids, interval));
        log_entry("Mesh building took "*string(meshtime)*" seconds");
        config.mesh_type = UNSTRUCTURED;
        
    elseif msh == HEXMESH
        if length(interval) == 2
            interval = [interval[1], interval[2], interval[1], interval[2], interval[1], interval[2]];
        end
        log_entry("Building simple hex mesh with nx*nx*nx elements, nx="*string(elsperdim));
        meshtime = @elapsed(mshdat = simple_hex_mesh(elsperdim.+1, bids, interval));
        log_entry("Mesh building took "*string(meshtime)*" seconds");
        
    else # msh should be a mesh file name
        # open the file and read the mesh data
        mfile = open(msh, "r");
        log_entry("Reading mesh file: "*msh);
        meshtime = @elapsed(mshdat=read_mesh(mfile));
        log_entry("Mesh reading took "*string(meshtime)*" seconds");
        close(mfile);
        # Assume an irregular mesh
        config.geometry = IRREGULAR;
        config.mesh_type = UNSTRUCTURED;
    end
    
    # Set this in Finch and build the corresponding Grid
    add_mesh(mshdat, partitions=partitions);
    
    # If bids>1 were specified for built meshes, add them here
    if bids > 1
        if msh == LINEMESH
            # already done in mesh_data
            # if bn == 2
            #     add_boundary_ID_to_grid(2, x -> (x >= interval[2]), grid_data);
            # end
            
        elseif msh == QUADMESH
            if bids == 2
                add_boundary_ID_to_grid(2, (x,y) -> (y <= interval[3]) || (y >= interval[4]), grid_data);
            elseif bids == 3
                add_boundary_ID_to_grid(2, (x,y) -> (x >= interval[2]), grid_data);
                add_boundary_ID_to_grid(3, (x,y) -> ((y <= interval[3]) || (y >= interval[4])) && (x > interval[1] && x < interval[2]), grid_data);
            elseif bids == 4
                add_boundary_ID_to_grid(2, (x,y) -> (x >= interval[2]), grid_data);
                add_boundary_ID_to_grid(3, (x,y) -> (y <= interval[3] && (x > interval[1] && x < interval[2])), grid_data);
                add_boundary_ID_to_grid(4, (x,y) -> (y >= interval[4] && (x > interval[1] && x < interval[2])), grid_data);
            end
            
        elseif msh == HEXMESH
            tiny = 1e-13;
            if bids == 6
                # bids = [1,2,3,4,5,6]; # all separate
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny), grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> (y <= interval[3]+tiny), grid_data);
                add_boundary_ID_to_grid(4, (x,y,z) -> (y >= interval[4]-tiny), grid_data);
                add_boundary_ID_to_grid(5, (x,y,z) -> (z <= interval[5]+tiny), grid_data);
                add_boundary_ID_to_grid(6, (x,y,z) -> (z >= interval[6]-tiny), grid_data);
            elseif bids == 5
                # bids = [1,2,3,4,5]; # combine z
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny), grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> (y <= interval[3]+tiny), grid_data);
                add_boundary_ID_to_grid(4, (x,y,z) -> (y >= interval[4]-tiny), grid_data);
                add_boundary_ID_to_grid(5, (x,y,z) -> (z <= interval[5]+tiny) || (z >= interval[6]-tiny), grid_data);
            elseif bids == 4
                # bids = [1,2,3,4]; # combine y and z
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]), grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> ((y <= interval[3]+tiny) || (y >= interval[4]-tiny)), grid_data);
                add_boundary_ID_to_grid(4, (x,y,z) -> ((z <= interval[5]+tiny) || (z >= interval[6]-tiny)), grid_data);
            elseif bids == 3
                # bids = [1,2,3]; # combine x,y,z
                add_boundary_ID_to_grid(2, (x,y,z) -> (y <= interval[3]+tiny) || (y >= interval[4]-tiny), grid_data);
                add_boundary_ID_to_grid(3, (x,y,z) -> (z <= interval[5]+tiny) || (z >= interval[6]-tiny), grid_data);
            elseif bids == 2
                # bids = [1,2]; # x=0, other
                add_boundary_ID_to_grid(2, (x,y,z) -> (x >= interval[2]-tiny) || (y <= interval[3]+tiny) || (y >= interval[4]-tiny) || (z <= interval[5]+tiny) || (z >= interval[6]-tiny), grid_data);
            end
            
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
    output_mesh(mfile, format);
    close(mfile);
end

"""
    finiteVolumeOrder(order)

Set the order of flux reconstruction for FVM.
For order > 1 this will cause the mesh to be subdivided into a parent/child mesh.
Take this into account when designing the mesh.
"""
function finiteVolumeOrder(order)
    if config.dimension > 2
        printerr("Sorry, higher order FV is not ready for 3D (TODO: build parent/child grid)\n Continuing with first order.");
        return;
    end
    (parent, child) = divide_parent_grid(grid_data, order);
    set_parent_and_child(parent, child, order);
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
    varind = var_count + 1;
    varsym = Symbol(name);
    # Not the default method?
    if !(config.solver_type == MIXED) && !(config.solver_type == method)
        method = config.solver_type;
    end
    # Just make an empty variable with the info known so far.
    var = Variable(varsym, [], varind, type, location, method, [], index, 0, [], false);
    add_variable(var);
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
    return add_coefficient(csym, type, location, val, nfuns, element_array, time_dependent);
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
    if length(parameters) == 0
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
        
    elseif typeof(val) <: Array
        newval = Array{Expr,1}(undef,length(val));
        for i=1:length(val)
            newval[i] = swap_parameter_xyzt(Meta.parse(val[i]));
        end
        newval = reshape(newval,size(val));
    else
        println("Error: use strings to define parameters");
        newval = 0;
    end
    
    return add_parameter(Symbol(name), type, newval);
end

"""
    testSymbol(symbol; type=SCALAR)

Define a symbol for a test function when using FEM.
Type can be SCALAR, VECTOR, TENSOR, or SYM_TENSOR.
"""
function testSymbol(symbol; type=SCALAR)
    add_test_function(Symbol(symbol), type);
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
        range = Array(range[1]:range[2]);
    end
    idx = Indexer(Symbol(name), range, range[1])
    add_indexer(idx);
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
            return add_variable_transform(var1, [var2], func);
        else
            return add_variable_transform(var1, var2, func);
        end
    else
        if typeof(var2) <: Array
            return add_variable_transform([var1], var2, func);
        else
            return add_variable_transform(var1, var2, func);
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
    for i=1:length(variable_transforms)
        if var1 == variable_transforms[i].from || [var1] == variable_transforms[i].from
            if var2 == variable_transforms[i].to || [var2] == variable_transforms[i].to
                found = true;
                transform_variable_values(variable_transforms[i]);
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
                # do nothing
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
        nfuns = 0;
    else
        # ??
    end
    
    add_boundary_condition(var, bid, bc_type, newbc_exp, nfuns);
end

"""
    addBoundaryID(bid, trueOnBdry)

Create a new boundary region with the given ID number bid. It will be
assigned to all boundary faces where the center of the face satisfies
trueOnBdry. It will override any previously set ID for those faces.

trueOnBdry can be a function or a string expression of (x,y,z).
Note that it only applies to faces that are known boundary faces, not
interior faces.
"""
function addBoundaryID(bid, trueOnBdry)
    # trueOnBdry(x, y, z) = something # points with x,y,z on this bdry segment evaluate true here
    if typeof(trueOnBdry) == String
        trueOnBdry = stringToFunction("trueOnBdry", "x,y=0,z=0", trueOnBdry);
    end
    add_boundary_ID_to_grid(bid, trueOnBdry, grid_data);
end

"""
    referencePoint(var, pos, val)

Constrain a variable value at a single node to the given value.
This is needed for certain situations where boundary conditions do not
uniquely constrain a variable. The node closest to the position in pos
will be used and it can be a boundary or interior node. 
"""
function referencePoint(var, pos, val)
    add_reference_point(var, pos, val);
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
    add_initial_condition(var.index, value, nfuns);
end

"""
    preStepFunction(fun)

Set a function to be called before each time step, or stage for multi-stage
steppers.
"""
function preStepFunction(fun)
    solver.set_pre_step(fun);
end

"""
    postStepFunction(fun)

Set a function to be called after each time step, or stage for multi-stage
steppers.
"""
function postStepFunction(fun)
    solver.set_post_step(fun);
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
    
    add_callback_function(CallbackFunction(name, args, body, fun));
    
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
    
    @timeit timer_output "CodeGen-weakform" begin
    
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
    
    log_entry("Making weak form for variable(s): "*string(wfvars));
    log_entry("Weak form, input: "*string(wf));
    
    # This is the parsing step. It goes from an Expr to arrays of Basic
    result_exprs = sp_parse(wfex, wfvars);
    if length(result_exprs) == 4 # has surface terms
        (lhs_symexpr, rhs_symexpr, lhs_surf_symexpr, rhs_surf_symexpr) = result_exprs;
    else
        (lhs_symexpr, rhs_symexpr) = result_exprs;
    end
    
    # Here we set a SymExpression for each of the pieces. 
    # This is an Expr tree that is passed to the code generator.
    if length(result_exprs) == 4 # has surface terms
        set_symexpressions(var, lhs_symexpr, LHS, "volume");
        set_symexpressions(var, lhs_surf_symexpr, LHS, "surface");
        set_symexpressions(var, rhs_symexpr, RHS, "volume");
        set_symexpressions(var, rhs_surf_symexpr, RHS, "surface");
        
        log_entry("lhs volume symexpression:\n\t"*string(lhs_symexpr));
        log_entry("lhs surface symexpression:\n\t"*string(lhs_surf_symexpr));
        log_entry("rhs volume symexpression:\n\t"*string(rhs_symexpr));
        log_entry("rhs surface symexpression:\n\t"*string(rhs_surf_symexpr));
        
        log_entry("Latex equation:\n\t\t \\int_{K}"*symexpression_to_latex(lhs_symexpr)*
                    " dx + \\int_{\\partial K}"*symexpression_to_latex(lhs_surf_symexpr)*
                    " ds = \\int_{K}"*symexpression_to_latex(rhs_symexpr)*
                    " dx + \\int_{\\partial K}"*symexpression_to_latex(rhs_surf_symexpr)*" ds", 3);
    else
        set_symexpressions(var, lhs_symexpr, LHS, "volume");
        set_symexpressions(var, rhs_symexpr, RHS, "volume");
        
        log_entry("lhs symexpression:\n\t"*string(lhs_symexpr));
        log_entry("rhs symexpression:\n\t"*string(rhs_symexpr));
        
        log_entry("Latex equation:\n\t\t\$ \\int_{K}"*symexpression_to_latex(lhs_symexpr)*
                    " dx = \\int_{K}"*symexpression_to_latex(rhs_symexpr)*" dx");
    end
    
    # Init stepper here so it can be used in generation
    if prob.time_dependent
        init_stepper(grid_data.allnodes, time_stepper);
        if use_specified_steps
            time_stepper.dt = specified_dt;
            time_stepper.Nsteps = specified_Nsteps;
        end
    end
    
    # change symbolic layer into IR
    if length(result_exprs) == 4
        full_IR = build_IR_fem(lhs_symexpr, lhs_surf_symexpr, rhs_symexpr, rhs_surf_symexpr, var, config, prob, time_stepper);
    else
        full_IR = build_IR_fem(lhs_symexpr, nothing, rhs_symexpr, nothing, var, config, prob, time_stepper);
    end
    log_entry("Weak form IR: "*string(full_IR),3);
    
    # Generate code from IR
    code = generate_code_layer(var, full_IR, solver_type, language, gen_framework);
    log_entry("Code layer: \n" * code, 3);
    set_code(var, code);
    
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
    
    @timeit timer_output "CodeGen-conservationform" begin
    
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
        set_symexpressions(var, lhs_symexpr, LHS, "volume");
        set_symexpressions(var, lhs_surf_symexpr, LHS, "surface");
        set_symexpressions(var, rhs_symexpr, RHS, "volume");
        set_symexpressions(var, rhs_surf_symexpr, RHS, "surface");
        
        log_entry("lhs volume symexpression:\n\t"*string(lhs_symexpr));
        log_entry("lhs surface symexpression:\n\t"*string(lhs_surf_symexpr));
        log_entry("rhs volume symexpression:\n\t"*string(rhs_symexpr));
        log_entry("rhs surface symexpression:\n\t"*string(rhs_surf_symexpr));
        
        log_entry("Latex equation:\n\t\t \\int_{K} -"*symexpression_to_latex(lhs_symexpr)*
                    " dx + \\int_{K}"*symexpression_to_latex(rhs_symexpr)*" dx"*
                    "\\int_{\\partial K}"*symexpression_to_latex(lhs_surf_symexpr)*
                    " ds - \\int_{\\partial K}"*symexpression_to_latex(rhs_surf_symexpr)*" ds", 3);
    else
        set_symexpressions(var, lhs_symexpr, LHS, "volume");
        set_symexpressions(var, rhs_symexpr, RHS, "volume");
        
        log_entry("lhs symexpression:\n\t"*string(lhs_symexpr));
        log_entry("rhs symexpression:\n\t"*string(rhs_symexpr));
        
        log_entry("Latex equation:\n\t\t\$ \\int_{K}"*symexpression_to_latex(lhs_symexpr)*
                    " dx = \\int_{K}"*symexpression_to_latex(rhs_symexpr)*" dx", 3);
    end
    
    # Init stepper here so it can be used in generation
    if prob.time_dependent
        init_stepper(grid_data.allnodes, time_stepper);
        if use_specified_steps
            time_stepper.dt = specified_dt;
            time_stepper.Nsteps = specified_Nsteps;
        end
    end
    
    # change symbolic layer into IR
    if length(result_exprs) == 4
        full_IR = build_IR_fvm(lhs_symexpr, lhs_surf_symexpr, rhs_symexpr, rhs_surf_symexpr, var, config, prob, time_stepper, fv_info);
    else
        full_IR = build_IR_fvm(lhs_symexpr, nothing, rhs_symexpr, nothing, var, config, prob, time_stepper, fv_info);
    end
    log_entry("Conservation form IR: "*string(full_IR),3);
    
    # Generate code from IR
    code = generate_code_layer(var, full_IR, solver_type, language, gen_framework);
    log_entry("Code layer: \n" * code, 3);
    set_code(var, code);
    
    end # timer block
end

"""
    flux(var, fex)

Write the surface integral term of the PDE. This is the flux after transforming 
using the divergence theorem into a surface integral.
var can be a variable or an array of variables. When using arrays, fex must
also be an array of matching size.
fex is a string expression or array of them. It can include numbers, coefficients,
variables, parameters, indexers, and symbolic operators.
"""
function flux(var, fex)
    if typeof(var) <: Array
        # multiple simultaneous variables
        symvars = [];
        symfex = [];
        if !(length(var) == length(fex))
            printerr("Error in flux function: # of unknowns must equal # of equations. (example: flux([a,b,c], [f1,f2,f3]))");
        end
        for vi=1:length(var)
            push!(symvars, var[vi].symbol);
            push!(symfex, Meta.parse((fex)[vi]));
        end
    else
        symfex = Meta.parse(fex);
        symvars = var.symbol;
    end
    
    log_entry("Making flux for variable(s): "*string(symvars));
    log_entry("flux, input: "*string(fex));
    
    # The parsing step
    (lhs_symexpr, rhs_symexpr) = sp_parse(symfex, symvars, is_FV=true, is_flux=true);
    
    set_symexpressions(var, lhs_symexpr, LHS, "surface");
    set_symexpressions(var, rhs_symexpr, RHS, "surface");
    
    log_entry("flux lhs symexpression:\n\t"*string(lhs_symexpr));
    log_entry("flux rhs symexpression:\n\t"*string(rhs_symexpr));
    log_entry("Latex flux equation:\n\t\t \\int_{\\partial K}"*symexpression_to_latex(lhs_symexpr)*
                    " ds - \\int_{\\partial K}"*symexpression_to_latex(rhs_symexpr)*" ds");
    
    # change symbolic layer into code layer
    (lhs_string, lhs_code) = generate_code_layer(lhs_symexpr, var, LHS, "surface", FV, language, gen_framework);
    (rhs_string, rhs_code) = generate_code_layer(rhs_symexpr, var, RHS, "surface", FV, language, gen_framework);
    log_entry("flux, code layer: \n  LHS = "*string(lhs_string)*" \n  RHS = "*string(rhs_string));
    
    if language == JULIA || language == 0
        args = "args; kwargs...";
        makeFunction(args, string(lhs_code));
        set_lhs_surface(var, lhs_string);
        
        makeFunction(args, string(rhs_code));
        set_rhs_surface(var, rhs_string);
        
    else
        set_lhs_surface(var, lhs_code);
        set_rhs_surface(var, rhs_code);
    end
end

"""
    source(var, sex)

Write the volume integral term of the PDE.
var can be a variable or an array of variables. When using arrays, fex must
also be an array of matching size.
sex is a string expression or array of them. It can include numbers, coefficients,
variables, parameters, indexers, and symbolic operators.
"""
function source(var, sex)
    if typeof(var) <: Array
        # multiple simultaneous variables
        symvars = [];
        symsex = [];
        if !(length(var) == length(sex))
            printerr("Error in source function: # of unknowns must equal # of equations. (example: source([a,b,c], [s1,s2,s3]))");
        end
        for vi=1:length(var)
            push!(symvars, var[vi].symbol);
            push!(symsex, Meta.parse((sex)[vi]));
        end
    else
        symsex = Meta.parse(sex);
        symvars = var.symbol;
    end
    
    log_entry("Making source for variable(s): "*string(symvars));
    log_entry("source, input: "*string(sex));
    
    # The parsing step
    (lhs_symexpr, rhs_symexpr) = sp_parse(symsex, symvars, is_FV=true);
    
    set_symexpressions(var, lhs_symexpr, LHS, "volume");
    set_symexpressions(var, rhs_symexpr, RHS, "volume");
    
    log_entry("source lhs symexpression:\n\t"*string(lhs_symexpr));
    log_entry("source rhs symexpression:\n\t"*string(rhs_symexpr));
    log_entry("Latex source equation:\n\t\t \\int_{K} -"*symexpression_to_latex(lhs_symexpr)*
                    " dx + \\int_{K}"*symexpression_to_latex(rhs_symexpr)*" dx");
    
    # change symbolic layer into code code_layer
    (lhs_string, lhs_code) = generate_code_layer(lhs_symexpr, var, LHS, "volume", FV, language, gen_framework);
    (rhs_string, rhs_code) = generate_code_layer(rhs_symexpr, var, RHS, "volume", FV, language, gen_framework);
    log_entry("source, code layer: \n  LHS = "*string(lhs_string)*" \n  RHS = "*string(rhs_string));
    
    if language == JULIA || language == 0
        args = "args; kwargs...";
        makeFunction(args, string(lhs_code));
        set_lhs(var, lhs_string);
        
        makeFunction(args, string(rhs_code));
        set_rhs(var, rhs_string);
        
    else
        set_lhs(var, lhs_code);
        set_rhs(var, rhs_code);
    end
end

"""
    assemblyLoops(var, indices, parallel_type=[])

Specify the nesting order of loops for assembling the system.
This makes sense for problems with indexed variables.
Var is a variable of array of variables.
Indices is an array of indexer objects and a string "elements".
Loops will be nested in that order with outermost first.
If parallel_type is specified for the loops, they will be generated with 
that type of parallel strategy. Possibilities are "none", "mpi", "threads".
"""
function assemblyLoops(var, indices, parallel_type=[])
    if length(parallel_type) < length(indices)
        parallel_type = fill("none", length(indices));
    end
    # find the associated indexer objects
    indexer_list = [];
    for i=1:length(indices)
        if typeof(indices[i]) == String
            if indices[i] == "elements" || indices[i] == "cells"
                push!(indexer_list, "elements");
            else
                for j=1:length(indexers)
                    if indices[i] == string(indexers[j].symbol)
                        push!(indexer_list, indexers[j]);
                    end
                end
            end
            
        elseif typeof(indices[i]) == Indexer
            push!(indexer_list, indices[i]);
        end
        
    end
    
    (loop_string, loop_code) = generate_assembly_loops(var, indexer_list, config.solver_type, language, parallel_type);
    log_entry("assembly loop function: \n\t"*string(loop_string));
    
    if language == JULIA || language == 0
        if config.solver_type == FV
            args = "var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0; is_explicit=true";
        elseif config.solver_type == CG
            args = "var, bilinear, linear, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0; rhs_only=false";
        end
        makeFunction(args, string(loop_code));
        set_assembly_loops(var, loop_string);
    else
        set_assembly_loops(var, loop_code);
    end
end

"""
    exportCode(filename)

Export all generated code including the elemental calculation and assembly loop
code for all variables if available.

"""
function exportCode(filename)
    # For now, only do this for Julia code because others are already output in code files.
    if language == JULIA || language == 0
        file = open(filename*".jl", "w");
        println(file, "#=\nGenerated functions for "*project_name*"\n=#\n");
        
        for i=1:length(variables)
            code = code_strings[1][i];
            var = string(variables[i].symbol);
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
    if language == JULIA || language == 0
        file = open(filename*".jl", "r");
        lines = readlines(file, keep=true);
        
        # Loop over variables and check to see if a matching function is present.
        for i=1:length(variables)
            # Scan the file for a pattern like
            #   # begin solve function for u
            #       ...
            #   # end solve function for u
            # OR
            #   # No code set for u
            var = string(variables[i].symbol);
            func_begin_flag = "# begin solve function for "*var;
            func_end_flag = "# end solve function for "*var;
            func_none_flag = "# No code set for "*var;
            func_string = "";
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
                elseif occursin(func_none_flag, lines[st])
                    # There is no function for this variable
                    log_entry("While importing, no solve function was found for "*var);
                end
            end # lines loop
            
            # Generate the functions and set them in the right places
            if func_string == ""
                log_entry("While importing, nothing was found for "*var);
            else
                set_code(variables[i], func_string);
            end
            
        end # vars loop
        
        close(file);
        log_entry("Imported code from "*filename*".jl");
        
    else
        # TODO non-julia
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
    
    if config.solver_type == FV
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
    else
        if prob.time_dependent
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
            
        else
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
    eval_initial_conditions(); 
end # Just for consistent style because this is also an internal function

"""
    solve(var, nlvar=nothing; nonlinear=false)

Either solve the problem using an internal target, or generate all code
files for an external target.
Var is a variable or array of variables to solve for.
The keyword arguments are for nonlinear equations which have very limited
support, so generally they won't be used.
"""
function solve(var, nlvar=nothing; nonlinear=false)
    if use_cachesim
        return cachesimSolve(var);
    end
    
    global time_stepper; # This should not be necessary. It will go away eventually
    
    # Wrap this all in a timer
    @timeit timer_output "Solve" begin
    
    # Generate files or solve directly
    if !(!generate_external && (language == JULIA || language == 0)) # if an external code gen target is ready
        if typeof(var) <: Array
            varind = var[1].index;
        else
            varind = var.index;
        end
        generate_all_files(var, bilinears[varind], face_bilinears[varind], linears[varind], face_linears[varind]);
        
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
        
        # Evaluate initial conditions if not already done
        eval_initial_conditions();
        
        # If any of the variables are indexed types, generate assembly loops
        var_indexers = [];
        if typeof(var) <: Array
            for vi=1:length(var)
                if !(var[vi].indexer === nothing)
                    push!(var_indexers, var[vi].indexer);
                end
            end
        else
            if !(var.indexer === nothing)
                var_indexers = [var.indexer];
            end
        end
        if length(var_indexers)>0 && assembly_loops[varind] === nothing
            log_entry("Indexed variables detected, but no assembly loops specified. Using default.")
            assemblyLoops(var, ["elements"; var_indexers]);
        end
        
        # Use the appropriate solver
        if config.solver_type == CG
            lhs = bilinears[varind];
            rhs = linears[varind];
            func = solve_function[varind];
            loop_func = assembly_loops[varind];
            
            if prob.time_dependent
                # time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                # if use_specified_steps
                #     time_stepper.dt = specified_dt;
				#     time_stepper.Nsteps = specified_Nsteps;
                # end
                if (nonlinear)
                    if time_stepper.type == EULER_EXPLICIT || time_stepper.type == LSRK4
                        printerr("Warning: Use implicit stepper for nonlinear problem. cancelling solve.");
                        return;
                    end
                	t = @elapsed(result = CGSolver.nonlinear_solve(var, nlvar, lhs, rhs, time_stepper, assemble_func=loop_func));
				else
                	# t = @elapsed(result = CGSolver.linear_solve(var, lhs, rhs, time_stepper, assemble_func=loop_func));
                    t = @elapsed(result = CGSolver.linear_solve(var, func, time_stepper));
				end
                # result is already stored in variables
            else
                # solve it!
				if (nonlinear)
                	t = @elapsed(result = CGSolver.nonlinear_solve(var, nlvar, lhs, rhs, assemble_func=loop_func));
                else
                    t = @elapsed(result = CGSolver.linear_solve(var, func));
				end
                
                # values should already be placed in variable arrays
            end
            
            log_entry("Solved for "*varnames*".(took "*string(t)*" seconds)", 1);
            
        elseif config.solver_type == DG
            lhs = bilinears[varind];
            rhs = linears[varind];
            slhs = face_bilinears[varind];
            srhs = face_linears[varind];
            
            if prob.time_dependent
                time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                if use_specified_steps
                    time_stepper.dt = specified_dt;
				    time_stepper.Nsteps = specified_Nsteps;
                end
				if (nonlinear)
                	# t = @elapsed(result = DGSolver.nonlinear_solve(var, nlvar, lhs, rhs, slhs, srhs, time_stepper));
                    printerr("Nonlinear solver not ready for DG");
                    return;
				else
                	t = @elapsed(result = DGSolver.linear_solve(var, lhs, rhs, slhs, srhs, time_stepper));
				end
                # result is already stored in variables
            else
                # solve it!
				if (nonlinear)
                	# t = @elapsed(result = DGSolver.nonlinear_solve(var, nlvar, lhs, rhs, slhs, srhs));
                    printerr("Nonlinear solver not ready for DG");
                else
                    t = @elapsed(result = DGSolver.linear_solve(var, lhs, rhs, slhs, srhs));
				end
                
                # place the values in the variable value arrays
                if typeof(var) <: Array && length(result) > 1
                    tmp = 0;
                    totalcomponents = 0;
                    for vi=1:length(var)
                        totalcomponents = totalcomponents + length(var[vi].symvar);
                    end
                    for vi=1:length(var)
                        components = length(var[vi].symvar);
                        for compi=1:components
                            var[vi].values[compi,:] = result[(compi+tmp):totalcomponents:end];
                            tmp = tmp + 1;
                        end
                    end
                elseif length(result) > 1
                    components = length(var.symvar);
                    for compi=1:components
                        var.values[compi,:] = result[compi:components:end];
                    end
                end
            end
            
            log_entry("Solved for "*varnames*".(took "*string(t)*" seconds)", 1);
            
        elseif config.solver_type == FV
            slhs = bilinears[varind];
            srhs = linears[varind];
            flhs = face_bilinears[varind];
            frhs = face_linears[varind];
            loop_func = assembly_loops[varind];
            
            if prob.time_dependent
                if fv_grid === nothing
                    time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                else
                    time_stepper = init_stepper(fv_grid.allnodes, time_stepper);
                end
                if use_specified_steps
                    time_stepper.dt = specified_dt;
				    time_stepper.Nsteps = specified_Nsteps;
                end
				if (nonlinear)
                	printerr("Nonlinear solver not ready for FV");
                    return;
				else
                	t = @elapsed(result = solver.linear_solve(var, slhs, srhs, flhs, frhs, time_stepper, loop_func));
				end
                # result is already stored in variables
                log_entry("Solved for "*varnames*".(took "*string(t)*" seconds)", 1);
            else
                # does this make sense?
                printerr("FV assumes time dependence. Set initial conditions etc.");
            end
            
        elseif config.solver_type == MIXED
            # Need to determine the variable index for fe and fv
            fe_var_index = 0;
            fv_var_index = 0;
            if typeof(var) <: Array
                for vi=1:length(var)
                    if var[vi].discretization == FV && fv_var_index == 0
                        fv_var_index += var[vi].index;
                    elseif fe_var_index == 0;
                        fe_var_index += var[vi].index;
                    end
                end
                loop_func = [assembly_loops[fe_var_index], assembly_loops[fv_var_index]];
            else
                # This shouldn't happen?
                printerr("mixed solver types specified, but only one variable being solved for?")
                if var.discretization == FV
                    fv_var_index = var.index;
                    loop_func = assembly_loops[fv_var_index];
                else
                    fe_var_index = var.index;
                    loop_func = assembly_loops[fe_var_index];
                end
            end
            
            vol_lhs = [];
            vol_rhs = [];
            surf_lhs = [];
            surf_rhs = [];
            if fe_var_index > 0
                push!(vol_lhs, bilinears[fe_var_index]);
                push!(vol_rhs, linears[fe_var_index]);
                push!(surf_lhs, face_bilinears[fe_var_index]);
                push!(surf_rhs, face_linears[fe_var_index]);
            else
                push!(vol_lhs, nothing);
                push!(vol_rhs, nothing);
                push!(surf_lhs, nothing);
                push!(surf_rhs, nothing);
            end
            
            if fv_var_index > 0
                push!(vol_lhs, bilinears[fv_var_index]);
                push!(vol_rhs, linears[fv_var_index]);
                push!(surf_lhs, face_bilinears[fv_var_index]);
                push!(surf_rhs, face_linears[fv_var_index]);
            else
                push!(vol_lhs, nothing);
                push!(vol_rhs, nothing);
                push!(surf_lhs, nothing);
                push!(surf_rhs, nothing);
            end
            
            if prob.time_dependent
                if typeof(time_stepper) <: Array
                    time_stepper[1] = init_stepper(grid_data.allnodes, time_stepper[1]);
                    time_stepper[2] = init_stepper(grid_data.allnodes, time_stepper[2]);
                    
                    if use_specified_steps
                        time_stepper[1].dt = specified_dt;
                        time_stepper[1].Nsteps = specified_Nsteps;
                        time_stepper[2].dt = specified_dt;
                        time_stepper[2].Nsteps = specified_Nsteps;
                    else
                        min_dt = min(time_stepper[1].dt, time_stepper[2].dt);
                        max_Nsteps = max(time_stepper[1].Nsteps, time_stepper[2].Nsteps);
                        time_stepper[1].dt = min_dt;
                        time_stepper[1].Nsteps = max_Nsteps;
                        time_stepper[2].dt = min_dt;
                        time_stepper[2].Nsteps = max_Nsteps;
                    end
                else
                    time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                    if use_specified_steps
                        time_stepper.dt = specified_dt;
                        time_stepper.Nsteps = specified_Nsteps;
                    end
                end
                
				if (nonlinear)
                	printerr("Nonlinear solver not ready for mixed solver");
                    return;
				else
                	t = @elapsed(result = solver.linear_solve(var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, time_stepper, loop_func));
				end
                # result is already stored in variables
                log_entry("Solved for "*varnames*".(took "*string(t)*" seconds)", 1);
            else
                # does this make sense?
                printerr("FV assumes time dependence. Set initial conditions etc.");
            end
            
        end
    end
    
    # At this point all of the final values for the variables in var should be placed in var.values.
    end # timer block
    
end

"""
    cachesimSolve(var, nlvar=nothing; nonlinear=false)

When using the cache simulator target, this is used instead of solve().
"""
function cachesimSolve(var, nlvar=nothing; nonlinear=false)
    if !(!generate_external && (language == JULIA || language == 0))
        printerr("Cachesim solve is only ready for Julia direct solve");
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

        lhs = bilinears[varind];
        rhs = linears[varind];
        
        t = @elapsed(result = linear_solve_cachesim(var, lhs, rhs));
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
function outputValues(vars, filename; format="vtk", ascii=false) output_values(vars, filename, format=format, ascii=ascii); end
function output_values(vars, filename; format="vtk", ascii=false)
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
function finalizeFinch() finalize_finch() end
function finalize_finch()
    # Finalize generation
    finalize_code_generator();
    if use_cachesim
        CachesimOut.finalize();
    end
    # timeroutput
    show(timer_output)
    # log
    close_log();
    # # mpi
    # if config.use_mpi
    #     MPI.Finalize(); # If MPI is finalized, it cannot be reinitialized without restarting Julia
    # end
    if config.proc_rank == 0
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
    global use_cachesim = use;
end

"""
    mortonNodes(griddim)

Reorder the nodes in memory to a Morton ordering in 2D or 3D.
This currently only works for a uniform grid such as the one generated
with Finch's internal utility.
griddim is an array representing the nodal grid size: like [n,n] for 2D or
[n,n,n] for 3D.
"""
function mortonNodes(griddim) morton_nodes(griddim) end
function morton_nodes(griddim)
    t = @elapsed(global grid_data = reorder_grid_recursive!(grid_data, griddim, MORTON_ORDERING));
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
function mortonElements(griddim) morton_elements(griddim) end
function morton_elements(griddim)
    global elemental_order = get_recursive_order(MORTON_ORDERING, config.dimension, griddim);
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
function hilbertNodes(griddim) hilbert_nodes(griddim) end
function hilbert_nodes(griddim)
    t = @elapsed(global grid_data = reorder_grid_recursive!(grid_data, griddim, HILBERT_ORDERING));
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
function hilbertElements(griddim) hilbert_elements(griddim) end
function hilbert_elements(griddim)
    global elemental_order = get_recursive_order(HILBERT_ORDERING, config.dimension, griddim);
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
function tiledNodes(griddim, tiledim) tiled_nodes(griddim, tiledim) end
function tiled_nodes(griddim, tiledim)
    t = @elapsed(global grid_data = reorder_grid_tiled(grid_data, griddim, tiledim));
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
function tiledElements(griddim, tiledim) tiled_elements(griddim, tiledim) end
function tiled_elements(griddim, tiledim)
    global elemental_order = get_tiled_order(config.dimension, griddim, tiledim, true);
    log_entry("Reordered elements to tiled("*string(tiledim)*").", 2);
    ef_nodes();
end

"""
    elementFirstNodes()

This is the default node ordering. Element first means the elements are given
some order and the nodes are added elementwise according to that. An element's
nodes are ordered according to the reference element.
"""
function elementFirstNodes() ef_nodes() end
function ef_nodes()
    t = @elapsed(global grid_data = reorder_grid_element_first!(grid_data, config.basis_order_min, elemental_order));
    log_entry("Reordered nodes to EF. Took "*string(t)*" sec.", 2);
end

"""
    randomNodes(seed = 17)

Randomize nodes in memory for testing a worst-case arrangement.
The seed is for making results reproducible.
"""
function randomNodes(seed = 17) random_nodes(seed) end
function random_nodes(seed = 17)
    t = @elapsed(global grid_data = reorder_grid_random!(grid_data, seed));
    log_entry("Reordered nodes to random. Took "*string(t)*" sec.", 2);
end

"""
    randomElements(seed = 17)

Randomize the order of the elemental loop for testing a worst-case arrangement.
The seed is for making results reproducible.
"""
function randomElements(seed = 17) random_elements(seed) end
function random_elements(seed = 17)
    global elemental_order = random_order(size(grid_data.loc2glb,2), seed);
    log_entry("Reordered elements to random.", 2);
    random_nodes(seed);
end