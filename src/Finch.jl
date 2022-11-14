#=
The main module for Finch.
=#
module Finch

export init_finch, set_language, set_custom_gen_target, dendro, set_solver, set_stepper, set_specified_steps, set_matrix_free,
        reformat_for_stepper, reformat_for_stepper_fv,
        add_mesh, output_mesh, add_boundary_ID, add_test_function, 
        add_initial_condition, add_boundary_condition, add_reference_point, set_rhs, set_lhs, set_lhs_surface, set_rhs_surface,
        set_symexpressions, set_assembly_loops
export build_cache_level, build_cache, build_cache_auto
export sp_parse
export generate_code_layer, generate_code_layer_surface, generate_code_layer_fv
export Variable, add_variable, VariableTransform
export Coefficient, add_coefficient
export Parameter, add_parameter
export Indexer
export FVInfo

### Module's global variables ###
# configuration, output files, etc.
config = nothing;
prob = nothing;
project_name = "unnamedProject";
output_dir = pwd();
language = 0;
gen_framework = 0;
generate_external = false;
codegen_params = nothing;
#log
use_log = false;
log_file = "";
log_line_index = 1;
#mesh
mesh_data = nothing;    # The basic element information as read from a MSH file or generated here.
needed_grid_types = [false, false, false]; # [CG, DG, FV] true if that type is needed
# FEM specific
grid_data = nothing;    # The full collection of nodes(including internal nodes) and other mesh info in the actual DOF ordering.
dg_grid_data = nothing; # This DG version is only made if using mixed CG/DG. Otherwise it is in grid_data.
geo_factors = nothing;  # Geometric factors
refel = nothing;        # Reference element(s?)
# FVM specific
fv_order = 1;           # Order of reconstruction for FV
fv_info = nothing;      # Finite volume info
fv_refel = nothing;     # Reference element for FV
fv_grid = nothing;      # Similar to grid_data, but only first order CG elements, and includes ghost info when partitioned.
fv_geo_factors = nothing;# Geometric factors for FV

#problem variables
var_count = 0;
variables = [];
coefficients = [];
parameters = [];
test_functions = [];
indexers = [];
ordered_indexers = [];
variable_transforms = [];
#generated functions
genfunc_count = 0;
genfunctions = [];
callback_functions = [];
# solve functions for each variable
solve_function = [];
#assembly loop functions
assembly_loops = [];
#symbolic layer
symexpressions = [[],[],[],[]];
# string versions of generated code for printing
code_strings = [[],[],[],[],[]];
#time stepper
time_stepper = nothing;
specified_dt = 0;
specified_Nsteps = 0;
use_specified_steps = false;

use_cachesim = false;

#handles for custom code gen functions
custom_gen_funcs = [];

timer_output = nothing;

###############################################################
include("finch_includes.jl");
include("macros.jl");
include("finch_interface.jl"); # included here after globals are defined
###############################################################

config = Finch_config(); # These need to be initialized here
prob = Finch_prob();

# This has already been initialized on loading, so this is more of a reset function
function init_finch(name="unnamedProject")
    global config = Finch_config();
    global prob = Finch_prob();
    global project_name = name;
    global language = JULIA;
    global framework = 0;
    global generate_external = false;
    global codegen_params = nothing;
    global log_file = "";
    global use_log = false;
    global log_line_index = 1;
    global mesh_data = nothing;
    global needed_grid_types = [false, false, false];
    global grid_data = nothing;
    global dg_grid_data = nothing;
    global geo_factors = nothing;
    global fv_order = 1;
    global fv_info = nothing;
    global refel = nothing;
    global fv_refel = nothing;
    global fv_grid = nothing;
    global fv_geo_factors = nothing;
    global var_count = 0;
    global variables = [];
    global coefficients = [];
    global parameters = [];
    global test_functions = [];
    global indexers = [];
    global ordered_indexers = [];
    global variable_transforms = [];
    global genfunc_count = 0;
    global genfunctions = [];
    global callback_functions = [];
    global symexpressions = [[],[],[],[]];
    global code_strings = [[],[],[],[],[]];
    global time_stepper = nothing;
    global specified_dt = 0;
    global specified_Nsteps = 0;
    global use_specified_steps = false;
    global use_cachesim = false;
    
    global timer_output = TimerOutput();
    
    # check for MPI and initialize
    if @isdefined(MPI)
        MPI.Init();
        config.num_procs = MPI.Comm_size(MPI.COMM_WORLD);
        if config.num_procs < 2
            # If only one process, ignore MPI.
            # MPI.Finalize();
        else
            config.use_mpi = true;
            config.proc_rank = MPI.Comm_rank(MPI.COMM_WORLD);
        end
    end
    
    # check for thread availability
    config.num_threads = Threads.nthreads();
    if (config.num_threads > 1 || config.num_procs > 1) && config.proc_rank == 0
        println("Initialized with "*string(config.num_procs)*" processes and "*string(config.num_threads)*" threads per proc.");
    end
end

# Sets the code generation target
function set_included_gen_target(lang, framework, dirpath, name; head="")
    #TODO
    println("Premade generation targets currently under construction. Sorry. Reverting to Julia target.");
    return;
    
    global language = lang;
    global gen_framework = framework;
    global output_dir = dirpath;
    global project_name = name;
    global generate_external = true;
    init_codegenerator(dirpath, name, head);
    
    # Need to set these three functions
    
    set_generation_target(get_external_language_elements, generate_external_files);
end

# Setting a custom target requires three functions
# 1. get_external_language_elements() - file extensions, comment chars etc.
# 3. generate_external_files(var, IR) - Writes all files based on generated code
function set_custom_gen_target(lang_elements, file_maker, dirpath, name; head="", params=nothing)
    global language = CUSTOM_GEN_TARGET;
    global gen_framework = CUSTOM_GEN_TARGET;
    global output_dir = dirpath;
    global project_name = name;
    global generate_external = true;
    init_code_generator(dirpath, name, head);
    set_generation_target(lang_elements, file_maker);
end

# Set parameters used by code generation. This is specific to the target.
function set_codegen_parameters(params)
    global codegen_params = params;
end

# Set the needed discretizations and backend
function set_solver(stype, backend)
    config.linalg_backend = backend;
    if typeof(stype) <: Array
        config.solver_type = MIXED;
        for s in stype
            if s == CG
                global needed_grid_types[1] = true;
            elseif s == DG
                global needed_grid_types[2] = true;
            elseif s == FV
                global needed_grid_types[3] = true;
            end
        end
    elseif stype == DG
        config.solver_type = stype;
        global needed_grid_types[2] = true;
    elseif stype == CG
        config.solver_type = stype;
        global needed_grid_types[1] = true;
    elseif stype == FV
        config.solver_type = stype;
        global needed_grid_types[3] = true;
    end
    
    if backend == PETSC_SOLVER && !(PETSc===nothing)
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
    end
end

# Sets the time stepper
function set_stepper(type, cfl)
    if typeof(type) <: Array
        first_stepper = Stepper(type[1], cfl);
        second_stepper = Stepper(type[2], cfl);
        
        global time_stepper = [first_stepper, second_stepper];
        global config.stepper = MIXED_STEPPER;
        log_entry("Set mixed time stepper to ["*type[1]*", "*type[2]*"]", 1);
    else
        global time_stepper = Stepper(type, cfl);
        global config.stepper = type;
        log_entry("Set time stepper to "*type, 1);
    end
end

# Time steps can be chosen automatically or specified. This is for manual specifying.
function set_specified_steps(dt, steps)
    global specified_dt = dt;
    global specified_Nsteps = steps;
    global use_specified_steps = true;
    prob.time_dependent = true;
    prob.end_time = dt*steps;
    log_entry("Set time stepper values to dt="*string(dt)*", Nsteps="*string(steps), 1);
end

# Uses matrix free iterative solver with the given max iterations and error tolerance.
function set_matrix_free(max, tol)
    config.linalg_matrixfree = true;
    config.linalg_matfree_max = max;
    config.linalg_matfree_tol = tol;
end

# Adds a mesh and builds the full grid and reference elements.
function add_mesh(mesh; partitions=0)
    global mesh_data = mesh;
    global refel;
    global fv_refel;
    global grid_data;
    global dg_grid_data;
    global fv_grid;
    global geo_factors;
    global fv_geo_factors;
    global fv_info;
    
    # If no method has been specified, assume CG
    if !(needed_grid_types[1] || needed_grid_types[2] || needed_grid_types[3])
        set_solver(CG, DEFAULT_SOLVER);
    end
    
    if partitions==0
        np = config.num_procs; # By default each process gets a partition
        config.num_partitions = np;
        config.partition_index = config.proc_rank;
    else
        # If other partitioning strategies are desired.
        # More than one proc may be assigned to each mesh partition.
        # They are grouped by partition number: 0,0,0,1,1,1,...,p,p,p
        np = partitions;
        config.num_partitions = np;
        config.partition_index = Int(floor(((config.proc_rank+0.5) * np) / config.num_procs));
    end
    if config.use_mpi && np > 1
        if config.proc_rank == 0
            # Partition the mesh and send the partition info to each proc to build their subgrids.
            # For now, partition into the number of procs. This will be modified in the future.
            epart = get_element_partitions(mesh, np);
            
        else # other ranks recieve the partition info
            epart = fill(Cint(-1), mesh.nel);
            
        end
        # Broadcast partition info to all
        MPI.Bcast!(epart, 0, MPI.COMM_WORLD);
        
        # Each proc will build a subgrid containing only their elements and ghost neighbors.
        # Make a Grid struct for each needed type
        if needed_grid_types[1] # CG
            (refel, grid_data) = partitioned_grid_from_mesh(mesh_data, epart, grid_type=CG, order=config.basis_order_min);
        end
        if needed_grid_types[2] # DG
            if needed_grid_types[1] # CG also needed
                (refel, dg_grid_data) = partitioned_grid_from_mesh(mesh_data, epart, grid_type=DG, order=config.basis_order_min);
            else
                (refel, grid_data) = partitioned_grid_from_mesh(mesh_data, epart, grid_type=DG, order=config.basis_order_min);
            end
        end
        if needed_grid_types[3] # FV
            (fv_refel, fv_grid) = partitioned_grid_from_mesh(mesh_data, epart, grid_type=FV, order=1);
            # If only FV is used, also set grid_data and refel to this. Just in case.
            if !(needed_grid_types[1] || needed_grid_types[2])
                grid_data = fv_grid;
                refel = fv_refel;
            end
        end
        
    else
        # Make a Grid struct for each needed type
        if needed_grid_types[1] # CG
            (refel, grid_data) = grid_from_mesh(mesh_data, grid_type=CG, order=config.basis_order_min);
        end
        if needed_grid_types[2] # DG
            if needed_grid_types[1] # CG also needed
                (refel, dg_grid_data) = grid_from_mesh(mesh_data, grid_type=DG, order=config.basis_order_min);
            else
                (refel, grid_data) = grid_from_mesh(mesh_data, grid_type=DG, order=config.basis_order_min);
            end
        end
        if needed_grid_types[3] # FV
            (fv_refel, fv_grid) = grid_from_mesh(mesh_data, grid_type=FV, order=1);
            # If only FV is used, also set grid_data and refel to this. Just in case.
            if !(needed_grid_types[1] || needed_grid_types[2])
                grid_data = fv_grid;
                refel = fv_refel;
            end
        end
    end
    
    # regular parallel sided elements or simplexes have constant Jacobians, so only store one value per element.
    if ((config.geometry == SQUARE && config.mesh_type == UNIFORM_GRID) ||
        config.dimension == 1 ||
        (config.dimension == 2 && refel.Nfaces == 3) ||
        (config.dimension == 3 && refel.Nfaces == 4) )
        constantJ = true;
    else
        constantJ = false;
    end
    do_faces = false;
    do_vol = false;
    if config.solver_type == DG || config.solver_type == FV || config.solver_type == MIXED
        do_faces = true;
    end
    if config.solver_type == FV || config.solver_type == MIXED
        do_vol = true;
    end
    
    # FE version is always made?
    geo_factors = build_geometric_factors(refel, grid_data, do_face_detj=do_faces, do_vol_area=do_vol, constant_jacobian=constantJ);
    if needed_grid_types[3] # FV
        fv_geo_factors = build_geometric_factors(fv_refel, fv_grid, do_face_detj=do_faces, do_vol_area=do_vol, constant_jacobian=constantJ);
        fv_info = build_FV_info(fv_grid);
    end
    
    log_entry("Added mesh with "*string(mesh_data.nx)*" vertices and "*string(mesh_data.nel)*" elements.", 1);
    if np > 1
        e_count = zeros(Int, np);
        for ei=1:length(epart)
            e_count[epart[ei]+1] += 1;
        end
        log_entry("Number of elements in each partition: "*string(e_count), 2);
        log_entry("Full grid has "*string(length(grid_data.global_bdry_index))*" nodes.", 2);
        
    else
        log_entry("Full grid has "*string(size(grid_data.allnodes,2))*" nodes.", 2);
    end
end

# Write the mesh to a MSH file
function output_mesh(file, format)
    if config.proc_rank == 0
        write_mesh(file, format, mesh_data);
        log_entry("Wrote mesh data to file.", 1);
    end
end

# For higher order FV, a finer child map is used for calculation, but the 
# coarse parent map is partitioned.
function set_parent_and_child(p_maps, c_grid, order)
    global fv_order = order;
    # If the mesh has not meen made yet, this will be done after that.
    if c_grid === nothing
        return;
    end
    
    global fv_grid = c_grid;
    dim = config.dimension;
    nfaces = size(fv_grid.element2face,1);
    global fv_refel = build_refel(dim, 1, nfaces, config.elemental_nodes);
    global fv_geo_factors = build_geometric_factors(fv_refel, fv_grid, do_face_detj=true, do_vol_area=true, constant_jacobian=true);
    global fv_info = build_FV_info(fv_grid, order, p_maps);
    log_entry("Set child grid with "*string(size(c_grid.allnodes,2))*" nodes and "*string(size(c_grid.loc2glb,2))*" elements.", 2);
    
    # If CELL variables exist, resize their values
    N = size(fv_grid.loc2glb, 2);
    for i=1:length(variables)
        if variables[i].location == CELL
            if variables[i].type == SCALAR
                val_size = (1, N);
            elseif variables[i].type == VECTOR
                val_size = (config.dimension, N);
            elseif variables[i].type == TENSOR
                val_size = (config.dimension*config.dimension, N);
            elseif variables[i].type == SYM_TENSOR
                val_size = (Int((config.dimension*(config.dimension+1))/2), N);
            elseif variables[i].type == VAR_ARRAY
                if typeof(variables[i].indexer) <: Array
                    comps = 1;
                    for j=1:length(variables[i].indexer)
                        comps *= length(variables[i].indexer[j].range);
                    end
                elseif typeof(variables[i].indexer) == Indexer
                    comps = length(variables[i].indexer.range);
                else
                    comps = 1;
                end
                val_size = (comps, N);
            end
            
            variables[i].values = zeros(config.float_type, val_size);
        end
    end
end

# Maybe remove. This is in grid.jl
function add_boundary_ID(bid, on_bdry)
    add_boundary_ID_to_grid(bid, on_bdry, grid_data);
end

# Defines a test function symbol by creating a special coefficient object.
function add_test_function(v, type)
    varind = length(test_functions) + 1;
    # make SymType
    dim = config.dimension;
    components = 1;
    if type == SCALAR
        components = 1;
    elseif type == VECTOR
        components = dim;
    elseif type == TENSOR
        components = dim*dim;
    elseif type == SYM_TENSOR
        if dim == 1
            components = 1;
        elseif dim == 2
            components = 3;
        elseif dim == 3
            components = 6;;
        elseif dim == 4
            components = 10;
        end
    elseif type == VAR_ARRAY
        components = 1;
    end
    symvar = sym_var(string(v), type, components);

    push!(test_functions, Finch.Coefficient(v, symvar, varind, type, NODAL, [], false, false););
    log_entry("Set test function symbol: "*string(v)*" of type: "*type, 2);
end

# Adds a variable and allocates everything associated with it.
function add_variable(var)
    global var_count += 1;

    # adjust values arrays
    if grid_data === nothing
        N = 1;
    else
        if var.location == CELL
            if fv_grid === nothing
                N = size(grid_data.loc2glb, 2);
            else
                N = size(fv_grid.loc2glb, 2);
            end
        else
            N = size(grid_data.allnodes,2);
        end
    end
    
    if var.type == SCALAR
        val_size = (1, N);
        comps = 1;
    elseif var.type == VECTOR
        val_size = (config.dimension, N);
        comps = config.dimension;
    elseif var.type == TENSOR
        val_size = (config.dimension*config.dimension, N);
        comps = config.dimension*config.dimension;
    elseif var.type == SYM_TENSOR
        val_size = (Int((config.dimension*(config.dimension+1))/2), N);
        comps = val_size[1];
    elseif var.type == VAR_ARRAY
        if typeof(var.indexer) <: Array
            comps = 1;
            for i=1:length(var.indexer)
                comps *= length(var.indexer[i].range);
            end
        elseif typeof(var.indexer) == Indexer
            comps = length(var.indexer.range);
        else
            comps = 1;
        end
        val_size = (comps, N);
    end
    var.total_components = comps;
    
    if language == JULIA || language == 0
        var.values = zeros(config.float_type, val_size);
    else
        var.values = zeros(config.float_type, 1,0);
    end
    
    # make symbolic layer variable symbols
    var.symvar = sym_var(string(var.symbol), var.type, var.total_components);

    global variables = [variables; var];
    
    global solve_function = [solve_function; nothing];
    global symexpressions[1] = [symexpressions[1]; nothing];
    global symexpressions[2] = [symexpressions[2]; nothing];
    global symexpressions[3] = [symexpressions[3]; nothing];
    global symexpressions[4] = [symexpressions[4]; nothing];
    global code_strings[1] = [code_strings[1]; ""];
    global code_strings[2] = [code_strings[2]; ""];
    global code_strings[3] = [code_strings[3]; ""];
    global code_strings[4] = [code_strings[4]; ""];
    global code_strings[5] = [code_strings[5]; ""];

    log_entry("Added variable: "*string(var.symbol)*" of type: "*var.type*", location: "*var.location, 2);
end

# Adds a coefficient with either constant value or some generated function of (x,y,z,t)
function add_coefficient(c, type, location, val, nfuns, element_array=false, time_dependent=false)
    global coefficients;
    # The values of c will have the same array structure as val
    if typeof(val) <: Array
        vals = Array{Any}(undef,size(val));
    else
        vals = [];
    end
    if nfuns == 0 # constant values
        vals = val;
        if length(vals) == 1 && !(typeof(vals) <: Array)
            vals = [val];
        end

    else # genfunction values
        if typeof(val) <: Array
            ind = length(genfunctions) - nfuns + 1;
            for i=1:length(val)
                if typeof(val[i]) == String
                    vals[i] = genfunctions[ind];
                    ind += 1;
                else
                    vals[i] = val[i];
                end
            end
        else
            push!(vals, genfunctions[end]);
        end
    end
    # If element_array, the last index in size(vals) should be per element.
    # c[n1, nel] for n1 components by nel elements.
    if element_array
        sz = size(vals);
        components = 0;
        for ind=1:length(sz)-1
            components += sz[ind];
        end
        if components == 0
            components = 1;
        end
    else
        components = length(vals);
    end
    
    symvar = sym_var(string(c), type, components);

    index = length(coefficients) + 1;
    push!(coefficients, Coefficient(c, symvar, index, type, location, vals, element_array, time_dependent));
    
    if element_array
        log_entry("Added coefficient "*string(c)*" : (array of elemental values)", 2);
    else
        log_entry("Added coefficient "*string(c)*" : "*string(val), 2);
    end
    
    return coefficients[end];
end

# Adds a parameter entity.
function add_parameter(p, type, val)
    index = length(parameters);
    push!(parameters, Parameter(p, index, type, val));

    log_entry("Added parameter "*string(p)*" : "*string(val), 2);

    return parameters[end];
end

# Adds an Indexer entity
function add_indexer(indexer)
    push!(indexers, indexer);
    push!(ordered_indexers, indexer);
    if length(indexer.range) < 3
        log_entry("Added indexer "*string(indexer.symbol)*" : "*string(indexer.range));
    else
        log_entry("Added indexer "*string(indexer.symbol)*" : ["*string(indexer.range[1])*", ... "*string(indexer.range[end])*"]");
    end
end

# Sets the ordered indexers
function set_ordered_indexers(ind)
    global ordered_indexers = ind;
end

# Defines a variable transform
function add_variable_transform(var1, var2, func)
    push!(variable_transforms, VariableTransform(var1, var2, func));
    log_entry("Added transform "*string(func)*": "*string(var1)*" -> "*string(var2));
    return variable_transforms[end];
end

# Needed to insert parameter expressions into the weak form expressions.
function swap_parameter_xyzt(ex)
    if typeof(ex) == Symbol
        if ex === :x
            return :parameterCoefficientForx ;
        elseif ex === :y
            return :parameterCoefficientFory ;
        elseif ex === :z
            return :parameterCoefficientForz ;
        elseif ex === :t
            return :parameterCoefficientFort ;
        end
    elseif typeof(ex) == Expr
        for i=1:length(ex.args)
            ex.args[i] = swap_parameter_xyzt(ex.args[i]); # Recursively swap
        end
    end
    return ex;
end

# Sets initial condition for the given variable, but does not initialize.
function add_initial_condition(varindex, ex, nfuns)
    global prob;
    while length(prob.initial) < varindex
        prob.initial = [prob.initial; nothing];
    end
    if typeof(ex) <: Array
        vals = [];
        ind = length(genfunctions) - nfuns + 1;
        for i=1:length(ex)
            if typeof(ex[i]) == String
                push!(vals, genfunctions[ind]);
                ind += 1;
            else
                push!(vals, ex[i]);
            end
        end
        prob.initial[varindex] = vals;
    else
        var_components = size(variables[varindex].values,1);
        if typeof(ex) == String
            prob.initial[varindex] = fill(genfunctions[end], var_components);
        else
            prob.initial[varindex] = fill(ex, var_components);
        end
    end
    
    if length(prob.initial[varindex]) < 10
        log_entry("Initial condition for "*string(variables[varindex].symbol)*" : "*string(prob.initial[varindex]), 2);
    else
        log_entry("Initial condition for "*string(variables[varindex].symbol)*" : "*string(prob.initial[varindex][1:6])*" (truncated for printing)", 2);
    end
    # hold off on initializing till solve or generate is determined.
end

# Evaluates all available intitial conditions setting the values for variables.
function eval_initial_conditions()
    dim = config.dimension;

    # build initial conditions
    for vind=1:length(variables)
        if !(variables[vind].ready) && vind <= length(prob.initial)
            if !(prob.initial[vind] === nothing)
                if variables[vind].location == CELL && !(fv_grid === nothing)
                    # Need to use the fv_grid instead of grid_data
                    this_grid_data = fv_grid;
                    this_geo_factors = fv_geo_factors;
                    this_refel = fv_refel;
                else
                    this_grid_data = grid_data;
                    this_geo_factors = geo_factors;
                    this_refel = refel;
                end
                
                # Evaluate at nodes
                nodal_values = zeros(config.float_type, length(prob.initial[vind]), size(this_grid_data.allnodes,2));
                for ci=1:length(prob.initial[vind])
                    for ni=1:size(this_grid_data.allnodes,2)
                        if typeof(prob.initial[vind][ci]) <: Number
                            nodal_values[ci,ni] = prob.initial[vind][ci];
                        elseif dim == 1
                            nodal_values[ci,ni] = prob.initial[vind][ci].func(this_grid_data.allnodes[ni],0.0,0.0,0);
                        elseif dim == 2
                            nodal_values[ci,ni] = prob.initial[vind][ci].func(this_grid_data.allnodes[1,ni],this_grid_data.allnodes[2,ni],0.0,0);
                        elseif dim == 3
                            nodal_values[ci,ni] = prob.initial[vind][ci].func(this_grid_data.allnodes[1,ni],this_grid_data.allnodes[2,ni],this_grid_data.allnodes[3,ni],0);
                        end
                    end
                end
                
                # compute cell averages using nodal values if needed
                if variables[vind].location == CELL
                    nel = size(this_grid_data.loc2glb, 2);
                    for ei=1:nel
                        e = ei;
                        glb = this_grid_data.loc2glb[:,e];
                        vol = this_geo_factors.volume[e];
                        detj = this_geo_factors.detJ[e];
                        
                        for ci=1:length(prob.initial[vind])
                            variables[vind].values[ci,e] = detj / vol * (this_refel.wg' * this_refel.Q * (nodal_values[ci,glb][:]))[1];
                        end
                    end
                else
                    variables[vind].values = nodal_values;
                end
                
                variables[vind].ready = true;
                log_entry("Built initial conditions for: "*string(variables[vind].symbol));
            end
        end
    end
end

# Sets boundary condition for a variable and BID.
function add_boundary_condition(var, bid, type, ex, nfuns)
    global prob;
    # make sure the arrays are big enough
    if size(prob.bc_func, 1) < var_count || size(prob.bc_func, 2) < bid
        if !(grid_data===nothing)
            nbid = length(grid_data.bids);
        else
            # For some targets no grid is created
            nbid = 1;
        end
        
        tmp1 = fill(NO_BC, var_count, nbid);
        tmp2 = Array{Any,2}(undef, (var_count, nbid)); # need to keep this array open to any type
        fill!(tmp2, 0);
        tmp3 = zeros(Int, (var_count, nbid));
        
        for i=1:size(prob.bc_func,1)
            for j=1:size(prob.bc_func,2)
                tmp1[i,j] = prob.bc_type[i,j];
                tmp2[i,j] = prob.bc_func[i,j];
                tmp3[i,j] = prob.bid[i,j];
            end
        end
        prob.bc_type = tmp1;
        prob.bc_func = tmp2;
        prob.bid = tmp3;
    end
    
    # Add this boundary condition to the struct
    valstr = "";
    if typeof(ex) <: Array
        vals = [];
        valstr = "[";
        ind = length(genfunctions) - nfuns + 1;
        for i=1:length(ex)
            if typeof(ex[i]) <: Number
                push!(vals, ex[i]);
                valstr *= string(ex[i]);
            elseif typeof(ex[i]) == String
                push!(vals, genfunctions[ind]);
                valstr *= genfunctions[ind].name;
                ind += 1;
            else # callback
                push!(vals, ex[i]);
                valstr *= ex[i].name;
            end
            if i < length(ex)
                valstr *= ", ";
            end
        end
        
        # If the variable has mor components than specified, just expand vals with its last element
        var_components = size(var.values,1);
        if var_components > length(ex)
            for i=(length(ex)+1):var_components
                push!(vals, vals[i-1]);
            end
        end
        
        valstr *= "]";
        prob.bc_func[var.index, bid] = vals;
        
    else
        var_components = size(var.values,1);
        if typeof(ex) <: Number
            prob.bc_func[var.index, bid] = fill(ex, var_components); 
            valstr = string(ex);
        elseif typeof(ex) == String
            prob.bc_func[var.index, bid] = fill(genfunctions[end], var_components); 
            valstr = genfunctions[end].name;
        else # callback
            prob.bc_func[var.index, bid] = fill(ex, var_components); 
            valstr = ex.name;
        end
    end
    prob.bc_type[var.index, bid] = type;
    prob.bid[var.index, bid] = bid;

    log_entry("Boundary condition: var="*string(var.symbol)*" bid="*string(bid)*" type="*type*" val="*valstr, 2);
end

# A reference point is a single point with a defined value.
# It's treated like a Dirichlet boundary with one boundary node.
function add_reference_point(var, pos, val)
    global prob;
    # make sure the array is big enough
    if size(prob.ref_point,1) < var_count
        tmp = Array{Any,2}(undef, (var_count, 3));
        for i=1:size(prob.ref_point,1)
            tmp[i,:] = prob.ref_point[i,:];
        end
        for i=(size(prob.ref_point,1)+1):var_count
            tmp[i,1] = false;
            tmp[i,2] = [0,0];
            tmp[i,3] = [0];
        end
        prob.ref_point = tmp;
    end
    if typeof(pos) <: Number
        pos = pos*ones(config.dimension);
    end
    if typeof(val) <: Number
        val = val*ones(length(var.symvar));
    end
    
    # Find the closest vertex to pos
    # The stored pos is actually the index into glbvertex pointing to the closest vertex
    ind = [1,1];
    mindist = 12345;
    nel = size(grid_data.loc2glb,2);
    for ei=1:nel
        for i=1:size(grid_data.glbvertex,1)
            d = 0;
            for comp=1:length(pos)
                d = d + abs(grid_data.allnodes[comp, grid_data.glbvertex[i, ei]] - pos[comp]);
            end
            if d<mindist
                ind = [i,ei];
                mindist = d;
            end
        end
    end
    
    prob.ref_point[var.index, 1] = true;
    prob.ref_point[var.index, 2] = ind;
    prob.ref_point[var.index, 3] = val;
    
    log_entry("Reference point: var="*string(var.symbol)*" position="*string(pos)*" value="*string(val), 2);
end

# adds a callback function to the 
function add_callback_function(f)
    push!(callback_functions, f);
end

# Generates a function from a code string and sets that as the code for the variable(s).
function set_code(var, code, IR)
    if language == JULIA || language == 0
        code_expr = CodeGenerator.code_string_to_expr(code);
        # args = "args; kwargs...";
        # makeFunction(args, code_expr);
        makeCompleteFunction(code_expr);
        if typeof(var) <:Array
            for i=1:length(var)
                solve_function[var[i].index] = genfunctions[end];
                code_strings[1][var[i].index] = code;
            end
        else
            solve_function[var.index] = genfunctions[end];
            code_strings[1][var.index] = code;
        end
    else
        if typeof(var) <:Array
            for i=1:length(var)
                solve_function[var[i].index] = IR;
                code_strings[1][var[i].index] = code;
            end
        else
            solve_function[var.index] = IR;
            code_strings[1][var.index] = code;
        end
    end
end

function set_symexpressions(var, ex, lorr, vors)
    # If the ex is empty, don't set anything
    if ex == [] || ex == [[]]
        return;
    end
    
    if typeof(var) <:Array
        for i=1:length(var)
            set_symexpressions(var[i], ex, lorr, vors);
        end
        
    else
        global symexpressions;
        if lorr == LHS
            if vors == "volume"
                ind = 1;
            else
                ind = 2;
            end
        else # rhs
            if vors == "volume"
                ind = 3;
            else
                ind = 4;
            end
        end
        symexpressions[ind][var.index] = ex;
    end
end

end # module
