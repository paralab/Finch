#=
The main module for Finch.
=#
module Finch

include("finch_imports.jl");

# Need to define constants and all other structs before FinchState
include("finch_constants.jl");
include("finch_structs.jl");

# A global FinchState is maintained
# finch_state = FinchState("");
finch_state = FinchState(Float64, "");

###############################################################
include("finch_includes.jl");
###############################################################

# Create a new FinchState and initialize MPI and threads
function init_finch(T::DataType, name="unnamedProject")
    state = FinchState(T, name);
    state.config.float_type = T;
    state.config.index_type = Int;
    
    state.ops = load_basic_ops();
    
    # check for MPI and initialize
    if @isdefined(MPI)
        MPI.Init();
        state.config.num_procs = MPI.Comm_size(MPI.COMM_WORLD);
        if state.config.num_procs < 2
            # If only one process, ignore MPI.
            # MPI.Finalize();
        else
            state.config.use_mpi = true;
            state.config.proc_rank = MPI.Comm_rank(MPI.COMM_WORLD);
        end
    end
    
    # check for thread availability
    state.config.num_threads = Threads.nthreads();
    if (state.config.num_threads > 1 || state.config.num_procs > 1) && state.config.proc_rank == 0
        println("Initialized with "*string(state.config.num_procs)*" processes and "*
                    string(state.config.num_threads)*" threads per proc.");
    end
    
    global finch_state = state;
    return state;
end

# Sets the code generation target
function set_included_gen_target(state::FinchState, lang::String, framework::String, dirpath::String, name::String; head::String="")
    #TODO
    println("Premade generation targets currently under construction. Sorry. Reverting to Julia target.");
    return;
    
    # state.target_language = lang;
    # state.gen_framework = framework;
    # state.output_dir = dirpath;
    # state.project_name = name;
    # state.generate_external = true;
    # init_codegenerator(dirpath, name, head);
    
    # # Need to set these three functions
    
    # set_generation_target(get_external_language_elements, generate_external_files);
end

# Setting a custom target requires three functions
# 1. get_external_language_elements() - file extensions, comment chars etc.
# 3. generate_external_files(var, IR) - Writes all files based on generated code
function set_custom_gen_target(state::FinchState, lang_elements, file_maker, dirpath::String, name::String; head::String="", params=nothing)
    state.output_dir = dirpath;
    state.project_name = name;
    state.external_target = true;
    init_code_generator(dirpath, name, head);
    set_generation_target(lang_elements, file_maker);
end

# Set parameters used by code generation. This is specific to the target.
# params is the named pairs passed in the interface
function set_codegen_parameters(state::FinchState, params)
    state.target_parameters = Dict(params);
end

# Set the needed discretizations and backend
function set_solver(state::FinchState, stype)
    if typeof(stype) <: Array
        printerr("Mixed discretizations are not ready in this version of Finch.", fatal=true);
        state.config.solver_type = MIXED;
        for s in stype
            if s == CG
                state.needed_grid_types[1] = true;
            elseif s == DG
                state.needed_grid_types[2] = true;
            elseif s == FV
                state.needed_grid_types[3] = true;
            end
        end
    elseif stype == DG
        state.config.solver_type = stype;
        state.needed_grid_types[2] = true;
    elseif stype == CG
        state.config.solver_type = stype;
        state.needed_grid_types[1] = true;
    elseif stype == FV
        state.config.solver_type = stype;
        state.needed_grid_types[3] = true;
    end
end

# Sets the time stepper
function set_stepper(state::FinchState, type, cfl)
    if typeof(type) <: Array
        first_stepper = Stepper(type[1], cfl);
        second_stepper = Stepper(type[2], cfl);
        
        state.time_stepper = [first_stepper, second_stepper];
        state.config.stepper = MIXED_STEPPER;
        log_entry("Set mixed time stepper to ["*type[1]*", "*type[2]*"]", 1);
    else
        state.time_stepper = Stepper(type, cfl);
        state.config.stepper = type;
        log_entry("Set time stepper to "*type, 1);
    end
end

# Time steps can be chosen automatically or specified. This is for manual specifying.
function set_specified_steps(state::FinchState, dt, steps)
    state.specified_dt = dt;
    state.specified_Nsteps = steps;
    state.use_specified_steps = true;
    state.prob.time_dependent = true;
    state.prob.end_time = dt*steps;
    log_entry("Set time stepper values to dt="*string(dt)*", Nsteps="*string(steps), 1);
end

# Adds a mesh and builds the full grid and reference elements.
function add_mesh(state::FinchState, mesh; partitions=0)
    state.mesh_data = mesh;
    
    @timeit finch_state.timer_output "mesh2grid" begin
    # If no method has been specified, assume CG
    if !(state.needed_grid_types[1] || state.needed_grid_types[2] || state.needed_grid_types[3])
        set_solver(state, CG);
    end
    
    if partitions==0
        np = state.config.num_procs; # By default each process gets a partition
        state.config.num_partitions = np;
        state.config.partition_index = state.config.proc_rank;
    else
        # If other partitioning strategies are desired.
        # More than one proc may be assigned to each mesh partition.
        # They are grouped by partition number: 0,0,0,1,1,1,...,p,p,p
        np = partitions;
        state.config.num_partitions = np;
        state.config.partition_index = Int(floor(((state.config.proc_rank+0.5) * np) / state.config.num_procs));
    end
    if state.config.use_mpi && np > 1
        if state.config.proc_rank == 0
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
        if state.needed_grid_types[1] # CG
            (state.refel, state.grid_data) = partitioned_grid_from_mesh(state.mesh_data, epart, grid_type=CG, order=state.config.basis_order_min);
        end
        if state.needed_grid_types[2] # DG
            if state.needed_grid_types[1] # CG also needed
                (state.refel, state.dg_grid) = partitioned_grid_from_mesh(state.mesh_data, epart, grid_type=DG, order=state.config.basis_order_min);
            else
                (state.refel, state.grid_data) = partitioned_grid_from_mesh(state.mesh_data, epart, grid_type=DG, order=state.config.basis_order_min);
            end
        end
        if state.needed_grid_types[3] # FV
            (state.fv_refel, state.fv_grid) = partitioned_grid_from_mesh(state.mesh_data, epart, grid_type=FV, order=1);
            # If only FV is used, also set grid_data and refel to this. Just in case.
            if !(state.needed_grid_types[1] || state.needed_grid_types[2])
                state.grid_data = state.fv_grid;
                state.refel = state.fv_refel;
            end
        end
        
    else
        # Make a Grid struct for each needed type
        if state.needed_grid_types[1] # CG
            (state.refel, state.grid_data) = grid_from_mesh(state.mesh_data, grid_type=CG, order=state.config.basis_order_min);
        end
        if state.needed_grid_types[2] # DG
            if state.needed_grid_types[1] # CG also needed
                (state.refel, state.dg_grid) = grid_from_mesh(state.mesh_data, grid_type=DG, order=state.config.basis_order_min);
            else
                (state.refel, state.grid_data) = grid_from_mesh(state.mesh_data, grid_type=DG, order=state.config.basis_order_min);
            end
        end
        if state.needed_grid_types[3] # FV
            (state.fv_refel, state.fv_grid) = grid_from_mesh(state.mesh_data, grid_type=FV, order=1);
            # If only FV is used, also set grid_data and refel to this. Just in case.
            if !(state.needed_grid_types[1] || state.needed_grid_types[2])
                state.grid_data = state.fv_grid;
                state.refel = state.fv_refel;
            end
        end
    end
    end#timer
    
    # regular parallel sided elements or simplexes have constant Jacobians, so only store one value per element.
    if ((state.config.geometry == SQUARE && state.config.mesh_type == UNIFORM_GRID) ||
        state.config.dimension == 1 ||
        (state.config.dimension == 2 && state.refel.Nfaces == 3) ||
        (state.config.dimension == 3 && state.refel.Nfaces == 4) )
        constantJ = true;
    else
        constantJ = false;
    end
    do_faces = false;
    do_vol = false;
    if state.config.solver_type == DG || state.config.solver_type == FV || state.config.solver_type == MIXED
        do_faces = true;
    end
    if state.config.solver_type == FV || state.config.solver_type == MIXED
        do_vol = true;
    end
    
    # FE version is always made?
    @timeit finch_state.timer_output "geo factors" begin
    state.geo_factors = build_geometric_factors(state.refel, state.grid_data, do_face_detj=do_faces, do_vol_area=do_vol, constant_jacobian=constantJ);
    if state.needed_grid_types[3] # FV
        if state.refel.Np == state.fv_refel.Np
            state.fv_geo_factors = state.geo_factors;
        else
            state.fv_geo_factors = build_geometric_factors(state.fv_refel, state.fv_grid, do_face_detj=do_faces, do_vol_area=do_vol, constant_jacobian=constantJ);
        end
        @timeit finch_state.timer_output "fv_info" begin
        state.fv_info = build_FV_info(state.fv_grid);
        end#timer
    end
    end#timer
    
    log_entry("Added mesh with "*string(state.mesh_data.nx)*" vertices and "*string(state.mesh_data.nel)*" elements.", 1);
    if np > 1
        e_count = zeros(Int, np);
        for ei=1:length(epart)
            e_count[epart[ei]+1] += 1;
        end
        log_entry("Number of elements in each partition: "*string(e_count), 2);
        log_entry("Full grid has "*string(length(state.grid_data.global_bdry_index))*" nodes.", 2);
        
    else
        log_entry("Full grid has "*string(size(state.grid_data.allnodes,2))*" nodes.", 2);
    end
end

# Write the mesh to a MSH file
function output_mesh(state::FinchState, file::String, format::String)
    if state.config.proc_rank == 0
        write_mesh(file, format, state.mesh_data);
        log_entry("Wrote mesh data to file.", 1);
    end
end

# For higher order FV, a finer child map is used for calculation, but the 
# coarse parent map is partitioned.
function set_parent_and_child(state::FinchState, p_maps, c_grid, order)
    # If the mesh has not meen made yet, this will be done after that.
    if c_grid === nothing
        return;
    end
    
    state.fv_grid = c_grid;
    dim = state.config.dimension;
    nfaces = size(state.fv_grid.element2face,1);
    state.fv_refel = build_refel(dim, 1, nfaces, state.config.elemental_nodes);
    state.fv_geo_factors = build_geometric_factors(state.fv_refel, state.fv_grid, do_face_detj=true, do_vol_area=true, constant_jacobian=true);
    state.fv_info = build_FV_info(state.fv_grid, order, p_maps);
    log_entry("Set child grid with "*string(size(c_grid.allnodes,2))*" nodes and "*string(size(c_grid.loc2glb,2))*" elements.", 2);
    
    # If CELL variables exist, resize their values
    N = size(state.fv_grid.loc2glb, 2);
    for i=1:length(state.variables)
        if state.variables[i].location == CELL
            if state.variables[i].type == SCALAR
                val_size = (1, N);
            elseif state.variables[i].type == VECTOR
                val_size = (dim, N);
            elseif state.variables[i].type == TENSOR
                val_size = (dim * dim, N);
            elseif state.variables[i].type == SYM_TENSOR
                val_size = (Int((dim * (dim+1))/2), N);
            elseif state.variables[i].type == VAR_ARRAY
                if typeof(state.variables[i].indexer) <: Array
                    comps = 1;
                    for j=1:length(state.variables[i].indexer)
                        comps *= length(state.variables[i].indexer[j].range);
                    end
                elseif typeof(state.variables[i].indexer) == Indexer
                    comps = length(state.variables[i].indexer.range);
                else
                    comps = 1;
                end
                val_size = (comps, N);
            end
            
            state.variables[i].values = zeros(state.config.float_type, val_size);
        end
    end
end

# Maybe remove. This is in grid.jl
function add_boundary_ID(state::FinchState, bid::Int, on_bdry)
    add_boundary_ID_to_grid(bid, on_bdry, state.grid_data);
end

# Defines a test function symbol by creating a special coefficient object.
function add_test_function(state::FinchState, v, type::String)
    varind = length(state.test_functions) + 1;
    # make SymType
    dim = state.config.dimension;
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

    push!(state.test_functions, Coefficient(v, symvar, varind, type, NODAL, Vector{Union{Float64,GenFunction}}(undef,0), false, false););
    log_entry("Set test function symbol: "*string(v)*" of type: "*type, 2);
end

# Adds a variable and allocates everything associated with it.
function add_variable(state::FinchState, var)
    # adjust values arrays
    if length(state.grid_data.allnodes) == 0
        N = 1;
    else
        if var.location == CELL
            if length(state.fv_grid.allnodes) == 0
                N = size(state.grid_data.loc2glb, 2);
            else
                N = size(state.fv_grid.loc2glb, 2);
            end
        else
            N = size(state.grid_data.allnodes,2);
        end
    end
    
    dim = state.config.dimension;
    if var.type == SCALAR
        val_size = (1, N);
        comps = 1;
    elseif var.type == VECTOR
        val_size = (dim, N);
        comps = dim;
    elseif var.type == TENSOR
        val_size = (dim * dim, N);
        comps = dim * dim;
    elseif var.type == SYM_TENSOR
        val_size = (Int((dim * (dim+1))/2), N);
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
    
    if state.target_language == JULIA
        var.values = zeros(state.config.float_type, val_size);
    else
        var.values = zeros(state.config.float_type, 1,0);
    end
    
    # make symbolic layer variable symbols
    var.symvar = sym_var(string(var.symbol), var.type, var.total_components);

    push!(state.variables, var);
    
    push!(state.solve_functions, nothing);
    push!(state.symexpressions[1], nothing);
    push!(state.symexpressions[2], nothing);
    push!(state.symexpressions[3], nothing);
    push!(state.symexpressions[4], nothing);
    push!(state.code_strings, "");

    log_entry("Added variable: "*string(var.symbol)*" of type: "*var.type*", location: "*var.location, 2);
end

# Adds a coefficient with either constant value or some generated function of (x,y,z,t)
function add_coefficient(state::FinchState, c, type, location, val, nfuns, element_array=false, time_dependent=false)
    # The values of c will have the same array structure as val
    ftype = finch_state.config.float_type;
    if typeof(val) <: Array
        vals = Vector{Union{Float64,GenFunction}}(undef,length(val));
    else
        vals = Vector{Union{Float64,GenFunction}}(undef,1);
        val = [val];
    end
    
    ind = length(state.genfunctions) - nfuns + 1;
    for i=1:length(val)
        if typeof(val[i]) == String
            vals[i] = state.genfunctions[ind];
            ind += 1;
        else
            vals[i] = Float64(val[i]);
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

    index = length(state.coefficients) + 1;
    push!(state.coefficients, Coefficient(c, symvar, index, type, location, vals, element_array, time_dependent));
    
    if element_array
        log_entry("Added coefficient "*string(c)*" : (array of elemental values)", 2);
    else
        log_entry("Added coefficient "*string(c)*" : "*string(val), 2);
    end
    
    return state.coefficients[end];
end

# Adds a parameter entity.
function add_parameter(state::FinchState, p, type, val)
    index = length(state.parameters);
    push!(state.parameters, Parameter(p, index, type, val));

    log_entry("Added parameter "*string(p)*" : "*string(val), 2);

    return state.parameters[end];
end

# Adds an Indexer entity
function add_indexer(state::FinchState, indexer)
    push!(state.indexers, indexer);
    push!(state.ordered_indexers, indexer);
    if length(indexer.range) < 3
        log_entry("Added indexer "*string(indexer.symbol)*" : "*string(indexer.range));
    else
        log_entry("Added indexer "*string(indexer.symbol)*" : ["*string(indexer.range[1])*", ... "*string(indexer.range[end])*"]");
    end
end

# Sets the ordered indexers
function set_ordered_indexers(state::FinchState, ind::Vector{Indexer})
    state.ordered_indexers = ind;
end

# Defines a variable transform
function add_variable_transform(state::FinchState, var1, var2, func)
    push!(state.variable_transforms, VariableTransform(var1, var2, func));
    log_entry("Added transform "*string(func)*": "*string(var1)*" -> "*string(var2));
    return state.variable_transforms[end];
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
function add_initial_condition(state::FinchState, varindex, ex, nfuns)
    while length(state.prob.initial) < varindex
        push!(state.prob.initial, [0.0]);
    end
    if typeof(ex) <: Array
        vals = [];
        ind = length(state.genfunctions) - nfuns + 1;
        for i=1:length(ex)
            if typeof(ex[i]) == String
                push!(vals, state.genfunctions[ind]);
                ind += 1;
            else
                push!(vals, Float64(ex[i]));
            end
        end
        state.prob.initial[varindex] = vals;
    else
        var_components = size(state.variables[varindex].values,1);
        if typeof(ex) == String
            state.prob.initial[varindex] = fill(state.genfunctions[end], var_components);
        else
            state.prob.initial[varindex] = fill(Float64(ex), var_components);
        end
    end
    
    if length(state.prob.initial[varindex]) < 10
        log_entry("Initial condition for "*string(state.variables[varindex].symbol)*" : "*string(state.prob.initial[varindex]), 2);
    else
        log_entry("Initial condition for "*string(state.variables[varindex].symbol)*" : "*string(state.prob.initial[varindex][1:6])*" (truncated for printing)", 2);
    end
    # hold off on initializing till solve or generate is determined.
end

# Evaluates all available intitial conditions setting the values for variables.
function eval_initial_conditions(state::FinchState)
    dim = state.config.dimension;
    # build initial conditions
    for vind=1:length(state.variables)
        if !(state.variables[vind].ready) && vind <= length(state.prob.initial)
            if length(state.prob.initial[vind]) == state.variables[vind].total_components
                if state.variables[vind].location == CELL && !(length(state.fv_grid.allnodes) == 0)
                    # Need to use the fv_grid instead of grid_data
                    this_grid_data = state.fv_grid;
                    this_geo_factors = state.fv_geo_factors;
                    this_refel = state.fv_refel;
                else
                    this_grid_data = state.grid_data;
                    this_geo_factors = state.geo_factors;
                    this_refel = state.refel;
                end
                
                # Evaluate at nodes
                nodal_values = zeros(state.config.float_type, length(state.prob.initial[vind]), size(this_grid_data.allnodes,2));
                for ci=1:length(state.prob.initial[vind])
                    for ni=1:size(this_grid_data.allnodes,2)
                        if typeof(state.prob.initial[vind][ci]) <: Number
                            nodal_values[ci,ni] = state.prob.initial[vind][ci];
                        elseif dim == 1
                            nodal_values[ci,ni] = state.prob.initial[vind][ci].func(this_grid_data.allnodes[ni],0.0,0.0,0.0, 0,0,zeros(Int,0));
                        elseif dim == 2
                            nodal_values[ci,ni] = state.prob.initial[vind][ci].func(this_grid_data.allnodes[1,ni],this_grid_data.allnodes[2,ni],0.0,0.0, 0,0,zeros(Int,0));
                        elseif dim == 3
                            nodal_values[ci,ni] = state.prob.initial[vind][ci].func(this_grid_data.allnodes[1,ni],this_grid_data.allnodes[2,ni],this_grid_data.allnodes[3,ni],0.0, 0,0,zeros(Int,0));
                        end
                    end
                end
                
                # compute cell averages using nodal values if needed
                if state.variables[vind].location == CELL
                    nel = size(this_grid_data.loc2glb, 2);
                    for ei=1:nel
                        e = ei;
                        glb = this_grid_data.loc2glb[:,e];
                        vol = this_geo_factors.volume[e];
                        detj = this_geo_factors.detJ[e];
                        
                        for ci=1:length(state.prob.initial[vind])
                            state.variables[vind].values[ci,e] = detj / vol * (this_refel.wg' * this_refel.Q * (nodal_values[ci,glb][:]))[1];
                        end
                    end
                else
                    state.variables[vind].values = nodal_values;
                end
                
                state.variables[vind].ready = true;
                log_entry("Built initial conditions for: "*string(state.variables[vind].symbol));
            end
        end
    end
end

# Sets boundary condition for a variable and BID.
function add_boundary_condition(state::FinchState, var, bid, type, ex, nfuns)
    var_count = length(state.variables);
    # make sure the arrays are big enough
    if size(state.prob.bc_func, 1) < var_count || size(state.prob.bc_func, 2) < bid
        if length(state.grid_data.bids) > 0
            nbid = length(state.grid_data.bids);
        else
            # For some targets no grid is created
            nbid = bid;
        end
        
        tmp1 = fill(NO_BC, var_count, nbid);
        tmp2 = Matrix{Vector{Union{Float64,GenFunction}}}(undef, (var_count, nbid));
        fill!(tmp2, [0.0]);
        tmp3 = zeros(Int, nbid);
        
        for i=1:size(state.prob.bc_func,1)
            for j=1:size(state.prob.bc_func,2)
                tmp1[i,j] = state.prob.bc_type[i,j];
                tmp2[i,j] = state.prob.bc_func[i,j];
                tmp3[j] = state.prob.bid[j];
            end
        end
        state.prob.bc_type = tmp1;
        state.prob.bc_func = tmp2;
        state.prob.bid = tmp3;
    end
    
    # Add this boundary condition to the struct
    valstr = "";
    if typeof(ex) <: Array
        vals = [];
        valstr = "[";
        ind = length(state.genfunctions) - nfuns + 1;
        for i=1:length(ex)
            if typeof(ex[i]) <: Number
                push!(vals, Float64(ex[i]));
                valstr *= string(ex[i]);
            elseif typeof(ex[i]) == String
                push!(vals, state.genfunctions[ind]);
                valstr *= state.genfunctions[ind].name;
                ind += 1;
            # else # callback
            #     push!(vals, ex[i]);
            #     valstr *= ex[i].name;
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
        state.prob.bc_func[var.index, bid] = vals;
        
    else
        var_components = size(var.values,1);
        if typeof(ex) <: Number
            state.prob.bc_func[var.index, bid] = fill(Float64(ex), var_components); 
            valstr = string(ex);
        elseif typeof(ex) == String
            state.prob.bc_func[var.index, bid] = fill(state.genfunctions[end], var_components); 
            valstr = state.genfunctions[end].name;
        # else # callback
        #     state.prob.bc_func[var.index, bid] = fill(ex, var_components); 
        #     valstr = ex.name;
        end
    end
    state.prob.bc_type[var.index, bid] = type;
    state.prob.bid[bid] = bid;

    log_entry("Boundary condition: var="*string(var.symbol)*" bid="*string(bid)*" type="*type*" val="*valstr, 2);
end

# Add the boundary region to the FinchProblem struct
function add_boundary_ID_to_problem(bid, trueOnBdry)
    prob = finch_state.prob;
    nbid = length(prob.bid) + 1;
    var_count = length(finch_state.variables);
    # make sure the arrays are big enough
    if size(prob.bc_func, 1) < var_count || size(prob.bc_func, 2) < nbid
        tmp1 = fill(NO_BC, var_count, nbid);
        tmp2 = Matrix{Vector{Union{Float64,GenFunction}}}(undef, (var_count, nbid));
        fill!(tmp2, [0.0]);
        tmp3 = zeros(Int, nbid);
        tmp4 = Vector{String}(undef, nbid);
        
        for i=1:size(prob.bc_func,1)
            for j=1:size(prob.bc_func,2)
                tmp1[i,j] = prob.bc_type[i,j];
                tmp2[i,j] = prob.bc_func[i,j];
                tmp3[j] = prob.bid[j];
                tmp4[j] = prob.bid_def[j];
            end
        end
        prob.bc_type = tmp1;
        prob.bc_func = tmp2;
        prob.bid = tmp3;
        prob.bid_def = tmp4;
    end
    prob.bid[nbid] = bid;
    prob.bid_def[nbid] = trueOnBdry;
end

# A reference point is a single point with a defined value.
# It's treated like a Dirichlet boundary with one boundary node.
function add_reference_point(state::FinchState, var, pos, val)
    var_count = length(state.variables)
    # make sure the array is big enough
    if size(state.prob.ref_point,1) < var_count
        tmp = Matrix{Vector{Union{Int,Float64,GenFunction}}}(undef, (var_count, 2));
        tmp2 = Vector{Bool}(undef, var_count);
        for i=1:size(state.prob.ref_point,1)
            tmp[i,:] = state.prob.ref_point[i,:];
            tmp2[i] = state.prob.has_ref_point[i];
        end
        for i=(size(state.prob.ref_point,1)+1):var_count
            tmp2[i] = false;
            tmp[i,1] = [0,0];
            tmp[i,2] = [0.0];
        end
        state.prob.ref_point = tmp;
        state.prob.has_ref_point = tmp2;
    end
    if typeof(pos) <: Number
        pos = pos*ones(Int,state.config.dimension);
    end
    if typeof(val) <: Number
        val = val*ones(length(var.symvar));
    end
    
    # Find the closest vertex to pos
    # The stored pos is actually the index into glbvertex pointing to the closest vertex
    ind = [1,1];
    mindist = 12345;
    nel = size(state.grid_data.loc2glb,2);
    for ei=1:nel
        for i=1:size(state.grid_data.glbvertex,1)
            d = 0;
            for comp=1:length(pos)
                d = d + abs(state.grid_data.allnodes[comp, state.grid_data.glbvertex[i, ei]] - pos[comp]);
            end
            if d<mindist
                ind = [i,ei];
                mindist = d;
            end
        end
    end
    
    state.prob.has_ref_point[var.index] = true;
    state.prob.ref_point[var.index, 1] = ind;
    state.prob.ref_point[var.index, 2] = val;
    
    log_entry("Reference point: var="*string(var.symbol)*" position="*string(pos)*" value="*string(val), 2);
end

# adds a callback function to the 
function add_callback_function(state::FinchState, f)
    push!(state.callback_functions, f);
end

# Generates a function from a code string and sets that as the code for the variable(s).
function set_code(state::FinchState, var, code, IR)
    if state.target_language == JULIA
        code_expr = CodeGenerator.code_string_to_expr(code);
        # args = "args; kwargs...";
        # makeFunction(args, code_expr);
        makeCompleteFunction(code_expr);
        if typeof(var) <:Array
            for i=1:length(var)
                state.solve_functions[var[i].index] = state.genfunctions[end];
                state.code_strings[var[i].index] = code;
            end
        else
            state.solve_functions[var.index] = state.genfunctions[end];
            state.code_strings[var.index] = code;
        end
    else
        if typeof(var) <:Array
            for i=1:length(var)
                state.solve_functions[var[i].index] = IR;
                state.code_strings[var[i].index] = code;
            end
        else
            state.solve_functions[var.index] = IR;
            state.code_strings[var.index] = code;
        end
    end
end

function set_symexpressions(state::FinchState, var, ex, lorr, vors)
    # If the ex is empty, don't set anything
    if ex == [] || ex == [[]]
        return;
    end
    
    if typeof(var) <:Array
        for i=1:length(var)
            set_symexpressions(state, var[i], ex, lorr, vors);
        end
        
    else
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
        state.symexpressions[ind][var.index] = ex;
    end
end

end # module
