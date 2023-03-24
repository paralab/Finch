#=
This version is specific to GPU accelerated code.
Take an IR and generate julia code.
Two main tasks are:
- directly translate the IR into code
- create any target specfic support code to surround it

The GPU version creates:
1. a "generated_solve_function_for_" function that does all setup, allocation, time stepping
2. a "gpu_kernel_for_" function that is the gpu kernel for system assembly
=#

# This creates the full code string including support code.
# If desired this can make the code as a complete, stand-alone function
# or as just the body of a function that will be generated(default).
function generate_code_layer_julia_gpu(var::Vector{Variable{FT}}, IR::IR_part, solver, wrap_in_function=true) where FT<:AbstractFloat
    # This will hold the code string to be returned
    codes =fill("", length(var));
    code ="";
    aux_code ="";
    
    # Set up useful numbers
    # Count variables, dofs, and store offsets
    varcount = length(var);
    dofs_per_node = var[1].total_components;
    dofs_per_loop = length(var[1].symvar);
    offset_ind = [0];
    if varcount > 1
        offset_ind = zeros(Int, varcount);
        for i=2:varcount
            offset_ind[i] = dofs_per_node;
            dofs_per_node += var[i].total_components;
            dofs_per_loop += length(var[i].symvar);
        end
    end
    
    # 
    if solver == CG || solver == DG
        args = "(var::Vector{Variable{FT}}, mesh::Grid, refel::Refel, geometric_factors::GeometricFactors, "*
                "config::FinchConfig, coefficients::Vector{Coefficient}, variables::Vector{Variable{FT}}, "*
                "test_functions::Vector{Coefficient}, ordered_indexers::Vector{Indexer}, prob::FinchProblem, "*
                "time_stepper::Stepper, buffers::ParallelBuffers, timer_output::TimerOutput, nl_var=nothing) where FT<:AbstractFloat";
        disc_specific = "
# FEM specific pieces
Q = refel.Q;
wg = refel.wg;
gness = size(mesh.face2glb,2); # CG=1, DG=2

# For partitioned meshes
(partitioned_order, partitioned_sizes) = get_partitioned_ordering(dofs_per_node, mesh, config);
"
        
    elseif solver == FV
        args = "(var::Vector{Variable{FT}}, mesh::Grid, refel::Refel, geometric_factors::GeometricFactors, "*
                "fv_info::FVInfo, config::FinchConfig, coefficients::Vector{Coefficient}, variables::Vector{Variable{FT}}, "*
                "test_functions::Vector{Coefficient}, ordered_indexers::Vector{Indexer}, prob::FinchProblem, "*
                "time_stepper::Stepper, buffers::ParallelBuffers, timer_output::TimerOutput, nl_var=nothing) where FT<:AbstractFloat";
        disc_specific = "
# FVM specific pieces
dofs_global = fv_dofs_global;
# boundary values for flux on each boundary face
boundary_flux = zeros(Float64, num_bdry_faces * dofs_per_node);
boundary_dof_index = zeros(Int64, num_bdry_faces * dofs_per_node);

"
        
    else
        printerr("This solver type is not ready: "*solver);
        args = "()";
        disc_specific = "";
    end
    
    if wrap_in_function
        for vi=1:varcount
            codes[vi] = "function generated_solve_function_for_"*string(var[vi].symbol) * args * "\n";
        end
    end
    
    # Generated functions will be placed inline like:
    # coeffficient_1_value = 2 * x
    genfunction_names = Vector{String}(undef,0);
    genfunction_defs = Vector{String}(undef,0);
    for c in finch_state.coefficients
        nc = size(c.value,1);
        for i=1:nc
            if typeof(c.value[i]) == GenFunction
                already_included = false;
                for j=1:length(genfunction_list)
                    if c.value[i].name == genfunction_list[j]
                        already_included = true;
                        break;
                    end
                end
                
                if !already_included
                    push!(genfunction_names, c.value[i].name);
                    push!(genfunction_defs, c.value[i].str);
                end
            end
        end
    end
    
    # data types
    itype = string(finch_state.config.index_type);
    ftype = string(finch_state.config.float_type);
    if !(ftype=="Float64"||ftype=="Float32"||ftype==Float16)
        ftype *= "\n    $ftype = FT"
    end
    
    # index bounds for any indexers
    index_bounds = "";
    index_bound_array = "[";
    for ind in finch_state.ordered_indexers
        if !(ind.symbol === :FINCHELEMENTS)
            index_bounds *= "num_"*string(ind.symbol)*"_indices = $(length(ind.range));\n";
            index_bound_array *= "num_"*string(ind.symbol)*"_indices, ";
        else
            index_bounds *= "num_elements = mesh.nel_owned;\n";
            index_bound_array *= "num_elements, ";
        end
    end
    index_bound_array *= "]";
    if length(index_bounds) > 1
        index_bounds *= "tmp_index_ranges = $index_bound_array;\n";
    end
    
    code *="
# User specified data types for int and float
# int type is $itype
# float type is $ftype

# pre/post step functions if defined
pre_step_function = prob.pre_step_function;
post_step_function = prob.post_step_function;

$index_bounds

# Prepare some useful numbers
# dofs_per_node = "*string(dofs_per_node)*";
# dofs_per_loop = "*string(dofs_per_loop)*";
# dof_offsets = "*string(offset_ind)*";

varcount = length(var);
dofs_per_node = var[1].total_components;
dofs_per_loop = length(var[1].symvar);
dof_offsets = zeros(Int, varcount);
for i=2:varcount
dof_offsets[i] = dofs_per_node;
dofs_per_node += var[i].total_components;
dofs_per_loop += length(var[i].symvar);
end

nnodes_partition = size(mesh.allnodes,2);
nnodes_global = nnodes_partition;
num_elements = mesh.nel_owned;
num_elements_global = mesh.nel_global;
num_elements_ghost = mesh.nel_ghost;
num_faces = mesh.nface_owned + mesh.nface_ghost;

dofs_global = dofs_per_node * nnodes_global;
fv_dofs_global = dofs_per_node * num_elements_global;
dofs_partition = dofs_per_node * nnodes_partition;
fv_dofs_partition = dofs_per_node * (num_elements + num_elements_ghost);
num_partitions = config.num_partitions;
proc_rank = config.proc_rank;

nodes_per_element = refel.Np;
qnodes_per_element = refel.Nqp;
faces_per_element = refel.Nfaces;
nodes_per_face = refel.Nfp[1];
dofs_per_element = dofs_per_node * nodes_per_element;
local_system_size = dofs_per_loop * nodes_per_element;

num_bdry_faces = 0;
nbids = length(mesh.bids);
for bi=1:nbids
num_bdry_faces += length(mesh.bdryface[bi]);
end

"
    code *= disc_specific;
    
    # Allocate GPU storage based on needs of assembly block
    assembly_block = extract_specified_block(IR, "assembly");
    gpu_allocate = get_gpu_allocations(assembly_block, IR_entry_types());
    gpu_allocate *= "global_vector_gpu = CuArray(global_vector);\n"
    
    # Remove duplicates
    gpu_allocate = remove_duplicate_lines(gpu_allocate);
    # Add all allocated things to the list of args
    kernel_args = get_args_from_allocations(gpu_allocate);
    push!(kernel_args, "dofs_global");
    push!(kernel_args, "faces_per_element");
    push!(kernel_args, "index_ranges");
    # push!(kernel_args, "global_vector_gpu");
    
    # Generate GPU kernel code
    aux_code = generate_gpu_kernel(var, kernel_args, IR);
    
    # Do CPU allocations first
    allocation_block = IR.parts[1];
    code *= generate_from_IR_julia_gpu(allocation_block, []);
    
    # then GPU allocations
    code *= "@timeit timer_output \"GPU alloc\" begin\n" * 
            "# Allocate and transfer things to the GPU\n" * gpu_allocate * 
            "index_ranges = CuArray(tmp_index_ranges);\nend\n";
    
    # The rest of the CPU code
    for i=2:length(IR.parts)
        code *= generate_from_IR_julia_gpu(IR.parts[i], kernel_args);
    end
    
    # finish
    code *="\nreturn nothing;\n"
    
    if wrap_in_function
        for vi=1:varcount
            codes[vi] *= code * "\nend # function\n";
            codes[vi] = set_indent(codes[vi]);
        end
        
    else
        for vi=1:varcount
            codes[vi] = set_indent(code);
        end
    end
    
    aux_code = set_indent(aux_code);
    
    return (codes, aux_code);
end

# Scan the IR and extract all needed arrays on the gpu as well as args for kernel call
function get_gpu_allocations(IR, IRtypes)
    alloc = "";
    node_type = typeof(IR);
    if node_type == IR_data_node # data nodes will look like "x" or "x[i,2]" or "x[a[1],17]"
        # skip certain ones
        skipit = false;
        to_skip = [:source, :flux, :flux_tmp];
        for s in to_skip
            if IR.label === s
                skipit = true;
            end
        end
        
        if !skipit
            var_name = string(IR.label);
            
            if length(IR.size) > 0 && length(IR.index) > 0
                alloc *= var_name * "_gpu = CuArray("*var_name*");\n"
            end
        end
        
    elseif node_type == IR_operation_node
        if IR.type == IRtypes.assign_op
            if typeof(IR.args[2]) <: IR_part
                alloc *= get_gpu_allocations(IR.args[2], IRtypes);
            end
            
        elseif IR.type == IRtypes.math_assign_op # x += ...
            if typeof(IR.args[3]) <: IR_part
                alloc *= get_gpu_allocations(IR.args[3], IRtypes);
            end
            
        elseif IR.type == IRtypes.function_op # f(...)
            # any of the args could need this
            for i=2:length(IR.args)
                if typeof(IR.args[i]) <: IR_part
                    alloc *= get_gpu_allocations(IR.args[i], IRtypes);
                end
            end
            
        elseif IR.type == IRtypes.math_op # (a * b * c) or sin(a)
            for i=2:length(IR.args)
                if typeof(IR.args[i]) <: IR_part
                    alloc *= get_gpu_allocations(IR.args[i], IRtypes);
                end
            end
            
        elseif IR.type == IRtypes.member_op # refel.Q
            # something like refel.Q needs to change to refel_Q
            # so member_op{refel, Q} -> data_node{refel_Q}
            # and a conversion needs to be made like: refel_Q = CuArray(refel.Q)
            code = string(IR.args[1]) * "." * generate_from_IR_julia_gpu(IR.args[2], [], IRtypes);
            original = string(IR.args[1]) * ".";
            converted = string(IR.args[1]) * "_";
            if typeof(IR.args[2]) == IR_data_node
                data_part = string(IR.args[2].label);
                if length(IR.args[2].size) > 0 # array data
                    alloc *= converted * data_part * "_gpu = CuArray("*original*data_part*");\n"
                else
                    alloc *= converted * data_part * "_gpu = "*original*data_part*";\n"
                end
                
            elseif typeof(IR.args[2]) == Symbol
                data_part = string(IR.args[2]);
                alloc *= converted * data_part * "_gpu = "*original*data_part*";\n"
                
            elseif typeof(IR.args[2]) == IR_operation_node && IR.args[2].type == IRtypes.member_op
                converted *= string(IR.args[2].args[1]) * "_";
                original *= string(IR.args[2].args[1]) * ".";
                data_part = string(IR.args[2].args[2].label);
                
                if length(IR.args[2].args[2].size) > 0 # array data
                    alloc *= converted * data_part * "_gpu = CuArray("*original*data_part*");\n"
                else
                    alloc *= converted * data_part * "_gpu = "*original*data_part*";\n"
                end
            else
                printerr("unexpected IR in conversion to CuArray: " * string(IR));
                return "ERROR";
            end
            
        elseif IR.type == IRtypes.named_op # handled case by case
            op = IR.args[1];
            if op === :COEF_EVAL
                # If the coefficient value is a number or array of numbers
                # coefficients_1_value_gpu = coefficients[1].value;
                # coefficients_1_value_gpu = CuArray(coefficients[1].value);
                coef = finch_state.coefficients[IR.args[2]];
                ncomponents = length(coef.value);
                vals = zeros(ncomponents);
                has_func = false;
                for i=1:ncomponents
                    if typeof(coef.value[i]) <:Number
                        vals[i] = coef.value[i];
                    else
                        has_func = true;
                    end
                end
                
                gpu_name = "coefficients_$(IR.args[2])_value_gpu";
                if has_func
                    alloc *= gpu_name * " = CuArray(" * string(vals) * ");\n";
                else
                    alloc *= gpu_name * " = CuArray(Array{Float64}(coefficients[$(IR.args[2])].value));\n";
                end
                
            elseif op === :KNOWN_VAR
                gpu_name = "variables_$(IR.args[2])_values_gpu";
                alloc *= gpu_name * " = CuArray(variables[$(IR.args[2])].values);\n";
                
            elseif op === :TIMER
                alloc *= get_gpu_allocations(IR.args[3], IRtypes);
                
            elseif op === :GLOBAL_FORM_MATRIX
                alloc *= "global_matrix_gpu = CuArray(global_matrix);\n"
                
            elseif op === :GLOBAL_GATHER_VECTOR
                alloc *= "global_vector_gpu = CuArray(global_vector);\n"
                
            elseif op === :GLOBAL_GATHER_SYSTEM
                # ??
                # There any many other possibilities. Figure them out as needed.
            end
        end
        
    elseif node_type == IR_block_node # A collection of statements. Do them one line at a time
        for i=1:length(IR.parts)
            alloc *= get_gpu_allocations(IR.parts[i], IRtypes);
        end
        
    elseif node_type == IR_loop_node
        alloc *= get_gpu_allocations(IR.body, IRtypes);
        
    elseif node_type == IR_conditional_node
        alloc *= get_gpu_allocations(IR.condition, IRtypes);
        alloc *= get_gpu_allocations(IR.body, IRtypes);
        if !(IR.elsepart===nothing)
            alloc *= get_gpu_allocations(IR.elsepart, IRtypes);
        end
        
    end
    
    return alloc;
end

# Extract args for the gpu kernel from the allocations
function get_args_from_allocations(code)
    lines = split(code, "\n");
    nlines = length(lines);
    args = fill("", 0);
    
    for i=1:nlines
        if occursin(" = ", lines[i])
            arg = split(lines[i], " = ")[1];
            push!(args, arg);
        end
    end
    
    return args;
end

# Generate the GPU kernel from the IR
# This only generates the part for the assembly loops for now
function generate_gpu_kernel(var, args, IR)
    IRtypes = IR_entry_types();
    
    # The function start
    nargs = length(args);
    args_per_line = 5;
    code = "function gpu_assembly_kernel(";
    for i=1:nargs
        code *= args[i];
        if i < nargs
            code *= ", ";
            if mod(i,args_per_line) == 0
                code *= "\n            ";
            end
        end
    end
    code *= ")\n";
    
    varcount = length(var);
    dofs_per_node = var[1].total_components;
    dofs_per_loop = length(var[1].symvar);
    offset_ind = [0];
    if varcount > 1
        offset_ind = zeros(Int, varcount);
        for i=2:varcount
            offset_ind[i] = dofs_per_node;
            dofs_per_node += var[i].total_components;
            dofs_per_loop += length(var[i].symvar);
        end
    end
    
    # Extract index counts and values from current dof
    index_counts = "";
    index_values = "";
    multiplier = "";
    num_indices = length(finch_state.ordered_indexers);
    if dofs_per_loop > 1
        # TODO
        # VAR_INDEX = Int(mod(current_dof-1, dofs_per_loop) + 1);
        # multiplier, index_values, index_counts needed
    end
    if num_indices > 1
        for indi = num_indices:-1:1 # work in reverse
            ind = finch_state.ordered_indexers[indi];
            if indi == num_indices && dofs_per_loop == 1
                if !(ind.symbol === :FINCHELEMENTS)
                    index_symbol = string(ind.symbol);
                    index_values = "INDEX_VAL_$index_symbol = Int(mod(current_dof-1, num_$index_symbol) + 1);\n";
                    multiplier = "num_$index_symbol";
                    index_counts *= "num_$index_symbol = index_ranges[$indi];\n";
                else
                    index_values *= "eid = Int(mod(current_dof-1, num_elements) + 1);\n";
                    multiplier = "num_elements";
                    index_counts *= "num_elements = index_ranges[$indi];\n";
                end
                
            else
                if !(ind.symbol === :FINCHELEMENTS)
                    index_symbol = string(ind.symbol);
                    index_values *= "INDEX_VAL_$index_symbol = Int(floor(mod(current_dof-1, num_$index_symbol * $multiplier) / ($multiplier)) + 1);\n";
                    multiplier *= " * num_$index_symbol";
                    index_counts *= "num_$index_symbol = index_ranges[$indi];\n";
                else
                    index_values *= "eid = Int(floor(mod(current_dof-1, num_elements * $multiplier) / ($multiplier)) + 1);\n";
                    multiplier *= " * num_elements";
                    index_counts *= "num_elements = index_ranges[$indi];\n";
                end
            end
        end
        
    else
        index_counts = "num_elements = index_ranges[1];\n";
        index_values = "eid = current_dof;\n";
    end
    
    # Set up the dof loop
    code *="
# The dofs handled by this thread are determined by the global thread ID
# max_threads = blockDim().x * gridDim().x * blockDim().y * gridDim().y * blockDim().z * gridDim().z;
# block_id = blockIdx().x + (blockIdx().y - 1) * gridDim().x + (blockIdx().z - 1) * gridDim().x * gridDim().y;
# thread_id = threadIdx().x + (block_id - 1) * blockDim().x * blockDim().y * blockDim().z;

# simplified version with 1D grid and blocks
max_threads = blockDim().x * gridDim().x;
block_id = blockIdx().x;
thread_id = threadIdx().x + (block_id - 1) * blockDim().x;

# Index counts
$index_counts

# Loop over all dofs assigned to this thread
# strided by max_threads
current_dof = thread_id;
while current_dof <= dofs_global

# extract index values
$index_values

"
    
    # Add the body of the dof loop
    assembly_block = extract_specified_block(IR, "assembly");
    if assembly_block === nothing
        assenbly_block = IR_block_node([IR_comment_node("No assembly block found!")]);
    end
    # Skip the index values part until the comment "Begin assembly code"
    assembly_offset = 1;
    
    for i=1:length(assembly_block.parts)
        if typeof(assembly_block.parts[i]) == IR_comment_node && assembly_block.parts[i].string == "Begin assembly code"
            assembly_offset = i;
            break;
        end
    end
    assembly_block = IR_block_node(assembly_block.parts[assembly_offset:end]);
    # Parts with these symbols need to be changed into scalars
    for s in [:source, :flux, :flux_tmp]
        scalarize_data_nodes(assembly_block, s);
    end
    code *= generate_from_IR_gpu_assembly(assembly_block);
    
    # end of the kernel
    code *= "

# go to the next assigned dof
current_dof = current_dof + max_threads;

end # dof loop

return nothing;
end # GPU kernel
";
    
    return code;
end

# Extract a block with the specified label from the IR.
# Returns nothing if block is not present.
function extract_specified_block(IR, label)
    specified_block = nothing;
    irtype = typeof(IR);
    # traverse the IR looking for a block labeled "assembly"
    if irtype == IR_block_node && IR.label == label
        specified_block = IR;
    else
        # Look deeper into blocks, loops, etc
        if irtype == IR_block_node
            for part in IR.parts
                specified_block = extract_specified_block(part, label);
                if !(specified_block === nothing)
                    break;
                end
            end
        elseif irtype == IR_loop_node
            specified_block = extract_specified_block(IR.body, label);
        elseif irtype == IR_conditional_node
            specified_block = extract_specified_block(IR.body, label);
            if specified_block === nothing
                specified_block = extract_specified_block(IR.elsepart, label);
            end
        elseif irtype == IR_operation_node
            for part in IR.args
                specified_block = extract_specified_block(part, label);
                if !(specified_block === nothing)
                    break;
                end
            end
        end
    end
    
    return specified_block;
end

# Generate the assembly loop code for the GPU kernel
function generate_from_IR_gpu_assembly(IR, IRtypes::Union{IR_entry_types, Nothing} = nothing)
    if IRtypes === nothing
        IRtypes = IR_entry_types();
    end
    code = "";
    node_type = typeof(IR);
    
    if node_type == IR_data_node # data nodes will look like "x" or "x[i,2]" or "x[a[1],17]"
        var_name = string(IR.label) * "_gpu";
        
        if length(IR.size) == 0
            code = var_name;
        elseif length(IR.size) > 0 && length(IR.index) > 0
            code = var_name*"[";
            for i=1:length(IR.index)
                if typeof(IR.index[i]) <: IR_part
                    code *= generate_from_IR_gpu_assembly(IR.index[i], IRtypes);
                else
                    code *= string(IR.index[i]);
                end
                if i<length(IR.index)
                    code *= ", ";
                end
            end
            code *= "]";
        else
            code = var_name;
        end
        
    elseif node_type == IR_operation_node
        if IR.type == IRtypes.allocate_op
            # No allocation should happen here
            code = "# Invalid allocation\n";
            
        elseif IR.type == IRtypes.assign_op # x = ...
            if ((typeof(IR.args[2]) == IR_operation_node) && (IR.args[2].type == IRtypes.allocate_op))
                code = "# Invalid allocation";
                
            else
                if typeof(IR.args[1]) <: IR_part
                    code = generate_from_IR_gpu_assembly(IR.args[1], IRtypes);
                else
                    code = string(IR.args[1]);
                end
                code *= " = ";
                if typeof(IR.args[2]) <: IR_part
                    code *= generate_from_IR_gpu_assembly(IR.args[2], IRtypes);
                else
                    code *= string(IR.args[2]);
                end
            end
            
        elseif IR.type == IRtypes.math_assign_op # x += ...
            # What is the math op?
            math_op = IR.args[1];
            if math_op in [:+, :-, :/, :*, :(.+), :(.-), :(.*), :(./)]
                math_op = string(math_op);
            else
                printerr("Unsupported math op in math-assignment: "*string(math_op)*"=");
                math_op = "";
            end
            
            # Note that if the symbol (IR.args[1]) is not yet on the declared list, 
            # it will be here.
            if typeof(IR.args[2]) <: IR_part
                code = generate_from_IR_gpu_assembly(IR.args[2], IRtypes);
            else
                code = string(IR.args[2]);
            end
            
            code *= " "*math_op*"= ";
            if typeof(IR.args[3]) <: IR_part
                code *= generate_from_IR_gpu_assembly(IR.args[3], IRtypes);
            else
                code *= string(IR.args[3]);
            end
            
        elseif IR.type == IRtypes.function_op # f(...)
            if typeof(IR.args[1]) <: IR_part
                code = generate_from_IR_gpu_assembly(IR.args[1], IRtypes);
            else
                code = string(IR.args[1]);
            end
            
            # Conditional function
            if code == "conditional_function"
                # three args
                condition = generate_from_IR_gpu_assembly(IR.args[2], IRtypes);
                case1 = generate_from_IR_gpu_assembly(IR.args[3], IRtypes);
                case2 = generate_from_IR_gpu_assembly(IR.args[4], IRtypes);
                code = "(($condition) ? ($case1) : ($case2))"
                
            else
                code *= "(";
                for i=2:length(IR.args)
                    if typeof(IR.args[i]) <: IR_part
                        code *= generate_from_IR_gpu_assembly(IR.args[i], IRtypes);
                    else
                        code *= string(IR.args[i]);
                    end
                    if i < length(IR.args)
                        code *= ", ";
                    end
                end
                code *= ")";
            end
            
        elseif IR.type == IRtypes.math_op # (a * b * c) or sin(a)
            if IR.args[1] in [:+, :-, :*, :/, :&&, :||, :<, :>, :(==), :(>=), :(<=)]
                if length(IR.args) < 3 # -a, +a
                    code = "(" * string(IR.args[1]);
                    if typeof(IR.args[2]) <: IR_part
                        code *= generate_from_IR_gpu_assembly(IR.args[2], IRtypes);
                    else
                        code *= string(IR.args[2]);
                    end
                    code *= ")";
                else # a + b + c + d
                    code = "(";
                    for i=2:length(IR.args)
                        if typeof(IR.args[i]) <: IR_part
                            code *= generate_from_IR_gpu_assembly(IR.args[i], IRtypes);
                        else
                            code *= string(IR.args[i]);
                        end
                        if i < length(IR.args)
                            code *= " " * string(IR.args[1]) * " ";
                        end
                    end
                    code *= ")";
                end
                
            else # sin(a)
                
                if string(IR.args[1]) == "conditional_function"
                    # three args
                    condition = generate_from_IR_gpu_assembly(IR.args[2], IRtypes);
                    case1 = generate_from_IR_gpu_assembly(IR.args[3], IRtypes);
                    case2 = generate_from_IR_gpu_assembly(IR.args[4], IRtypes);
                    code = "(($condition) ? ($case1) : ($case2))"
                    
                else
                    code = string(IR.args[1]) * "(";
                    for i=2:length(IR.args)
                        if typeof(IR.args[i]) <: IR_part
                            code *= generate_from_IR_gpu_assembly(IR.args[i], IRtypes);
                        else
                            code *= string(IR.args[i]);
                        end
                        if i < length(IR.args)
                            code *= ", ";
                        end
                    end
                    code *= ")";
                end
                
            end
            
        elseif IR.type == IRtypes.member_op # refel.Q -> refel_Q
            code = string(IR.args[1]) * "_" * generate_from_IR_gpu_assembly(IR.args[2], IRtypes);
            
        elseif IR.type == IRtypes.named_op # handled case by case
            code = generate_named_op_gpu_kernel(IR, IRtypes);
        end
        
    elseif node_type == IR_block_node # A collection of statements. Do them one line at a time
        code = "";
        # If it is the boundary conditions block, skip this on the GPU
        if IR.label == "boundary conditions"
            code = "# boundary conditions handled on cpu side\n"
            # create check for flux BCs
            bc_types = finch_state.prob.bc_type;
            bids = finch_state.prob.bid;
            nbids = length(bids);
            bid_check = "";
            for i=1:nbids
                if !(bc_types[1][i] == NO_BC) # TODO this assumes first variable only
                    if length(bid_check) > 0
                        bid_check *= " || ";
                    end
                    bid_check *= "fbid_gpu==$(bids[i])";
                end
            end
            if length(bid_check) == 0
                bid_check = "false";
            end
            if finch_state.config.solver_type == FV
                code *= "flux_tmp_gpu = ($(bid_check)) ? 0.0 : flux_tmp_gpu\n\n";
            else
                
            end
            
        elseif IR.label == "local to global"
            code = "# Row index is current_dof\n";
            # Note: this is only for explicit fvm. TODO loc2glb system
            code *= "global_vector_gpu[current_dof] = source_gpu + flux_gpu\n";
            
        else
            for i=1:length(IR.parts)
                if typeof(IR.parts[i]) == IR_operation_node
                    code *= generate_from_IR_gpu_assembly(IR.parts[i], IRtypes);
                    code *= "\n";
                else
                    code *= generate_from_IR_gpu_assembly(IR.parts[i], IRtypes);
                end
                
            end
        end
        
    elseif node_type == IR_loop_node
        # while and for loops
        if IR.type == IRtypes.while_loop
            condition = generate_from_IR_gpu_assembly(IR.last, IRtypes);
            code = string(IR.iterator) * " = " * string(IR.first) * ";\n";
            code *= "while (" * condition * ")\n";
            code *= string(IR.iterator) * " += 1;\n";
            code *= generate_from_IR_gpu_assembly(IR.body, IRtypes);
            code *= "end\n";
            
        else # for loop
            if (IR.first == IR.last) && IR.first==0
                range = string(IR.collection);
            else
                range = string(IR.first) * ":" * string(IR.last);
            end
            
            code = "for "* string(IR.iterator) * " = " * range * "\n";
            code *= generate_from_IR_gpu_assembly(IR.body, IRtypes);
            code *= "end\n";
        end
        
        
    elseif node_type == IR_conditional_node
        code = "if " * generate_from_IR_gpu_assembly(IR.condition, IRtypes) * "\n";
        code *= generate_from_IR_gpu_assembly(IR.body, IRtypes);
        if !(IR.elsepart===nothing)
            code *= "else\n" * generate_from_IR_gpu_assembly(IR.elsepart, IRtypes);
        end
        code *= "end\n";
        
    elseif node_type == IR_comment_node
        code = "#= " * IR.string * " =#\n";
        
    elseif node_type <: IR_part
        code = IR_string(IR);
        
    else
        code = string(IR);
    end
    
    return code;
end

# Directly translate IR into julia code
# Named operators are treated specially
function generate_from_IR_julia_gpu(IR, kernel_args, IRtypes::Union{IR_entry_types, Nothing} = nothing)
    if IRtypes === nothing
        IRtypes = IR_entry_types();
    end
    code = "";
    indent = "";
    node_type = typeof(IR);
    
    if node_type == IR_data_node # data nodes will look like "x" or "x[i,2]" or "x[a[1],17]"
        var_name = string(IR.label);
        
        if length(IR.size) == 0
            code = var_name;
        elseif length(IR.size) > 0 && length(IR.index) > 0
            code = var_name*"[";
            for i=1:length(IR.index)
                if typeof(IR.index[i]) <: IR_part
                    code *= generate_from_IR_julia_gpu(IR.index[i], kernel_args, IRtypes);
                else
                    code *= string(IR.index[i]);
                end
                if i<length(IR.index)
                    code *= ", ";
                end
            end
            code *= "]";
        else
            code = var_name;
        end
        
    elseif node_type == IR_operation_node
        if IR.type == IRtypes.allocate_op # zeros(Float64, 2,5)
            type_name = IRtypes.name[IR.args[1]];
            if type_name == "CustomInt"
                type_name = string(finch_state.config.index_type);
            elseif type_name == "CustomFloat"
                type_name = string(finch_state.config.float_type);
            end
            code = "zeros("*type_name;
            for i=2:length(IR.args)
                code *= ", ";
                if typeof(IR.args[i]) <: IR_part
                    code *= generate_from_IR_julia_gpu(IR.args[i], kernel_args, IRtypes);
                else
                    code *= string(IR.args[i]);
                end
            end
            code *= ")";
            
        elseif IR.type == IRtypes.assign_op # x = ...
            if ((typeof(IR.args[2]) == IR_operation_node) && (IR.args[2].type == IRtypes.allocate_op))
                type_name = IRtypes.name[IR.args[2].args[1]];
                if type_name == "CustomInt"
                    type_name = string(finch_state.config.index_type);
                elseif type_name == "CustomFloat"
                    type_name = string(finch_state.config.float_type);
                end
                adims = length(IR.args[2].args) - 1;
                if adims == 1
                    type_name = "::Vector{"*type_name*"}";
                elseif adims == 2
                    type_name = "::Matrix{"*type_name*"}";
                elseif adims > 2
                    type_name = "::Array{"*type_name*", "*string(adims)*"}";
                else
                    type_name = "::"*type_name;
                end
                
            elseif ((typeof(IR.args[2]) == IR_operation_node) && (IR.args[2].type == IRtypes.named_op) && (IR.args[2].args[1] == :COEF_EVAL))
                type_name = "::" * string(finch_state.config.float_type);
            else
                type_name = "";
            end
            if typeof(IR.args[1]) <: IR_part
                code = generate_from_IR_julia_gpu(IR.args[1], kernel_args, IRtypes);
            else
                code = string(IR.args[1]);
            end
            code *= type_name * " = ";
            if typeof(IR.args[2]) <: IR_part
                code *= generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes);
            else
                code *= string(IR.args[2]);
            end
            
        elseif IR.type == IRtypes.math_assign_op # x += ...
            # What is the math op?
            math_op = IR.args[1];
            if math_op in [:+, :-, :/, :*, :(.+), :(.-), :(.*), :(./)]
                math_op = string(math_op);
            else
                printerr("Unsupported math op in math-assignment: "*string(math_op)*"=");
                math_op = "";
            end
            
            # Note that if the symbol (IR.args[1]) is not yet on the declared list, 
            # it will be here.
            if typeof(IR.args[2]) <: IR_part
                code = generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes);
            else
                code = string(IR.args[2]);
            end
            
            code *= " "*math_op*"= ";
            if typeof(IR.args[3]) <: IR_part
                code *= generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes);
            else
                code *= string(IR.args[3]);
            end
            
        elseif IR.type == IRtypes.function_op # f(...)
            if typeof(IR.args[1]) <: IR_part
                code = generate_from_IR_julia_gpu(IR.args[1], kernel_args, IRtypes);
            else
                code = string(IR.args[1]);
            end
            code *= "(";
            for i=2:length(IR.args)
                if typeof(IR.args[i]) <: IR_part
                    code *= generate_from_IR_julia_gpu(IR.args[i], kernel_args, IRtypes);
                else
                    code *= string(IR.args[i]);
                end
                if i < length(IR.args)
                    code *= ", ";
                end
            end
            code *= ")";
            
        elseif IR.type == IRtypes.math_op # (a * b * c) or sin(a)
            if IR.args[1] in [:+, :-, :*, :/, :&&, :||, :<, :>, :(==), :(>=), :(<=)]
                if length(IR.args) < 3 # -a, +a
                    code = "(" * string(IR.args[1]);
                    if typeof(IR.args[2]) <: IR_part
                        code *= generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes);
                    else
                        code *= string(IR.args[2]);
                    end
                    code *= ")";
                else # a + b + c + d
                    code = "(";
                    for i=2:length(IR.args)
                        if typeof(IR.args[i]) <: IR_part
                            code *= generate_from_IR_julia_gpu(IR.args[i], kernel_args, IRtypes);
                        else
                            code *= string(IR.args[i]);
                        end
                        if i < length(IR.args)
                            code *= " " * string(IR.args[1]) * " ";
                        end
                    end
                    code *= ")";
                end
                
            else # sin(a)
                code = string(IR.args[1]) * "(";
                for i=2:length(IR.args)
                    if typeof(IR.args[i]) <: IR_part
                        code *= generate_from_IR_julia_gpu(IR.args[i], kernel_args, IRtypes);
                    else
                        code *= string(IR.args[i]);
                    end
                    if i < length(IR.args)
                        code *= ", ";
                    end
                end
                code *= ")";
            end
            
        elseif IR.type == IRtypes.member_op # refel.Q
            code = string(IR.args[1]) * "." * generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes);
            
        elseif IR.type == IRtypes.named_op # handled case by case
            code = generate_named_op_gpu(IR, kernel_args, IRtypes);
        end
        
    elseif node_type == IR_block_node # A collection of statements. Do them one line at a time
        code = "";
        
        ### For GPU target, the assembly block is split into assembly and boundary conditions.
        ### Assembly is done in the GPU kernel. Boundary conditions are done here.
        ### 
        if IR.label == "assembly loop"
            # This needs to be figured out TODO
            # varcount = length(var);
            # dofs_per_loop = length(var[1].symvar);
            # if varcount > 1
            #     for i=2:varcount
            #         dofs_per_loop += length(var[i].symvar);
            #     end
            # end
            
            dofs_per_loop = 1;
            
            # Create the index loops
            index_loops = "";
            index_loop_ends = "";
            index_offset = "";
            index_values = "";
            multiplier = "";
            num_indices = length(finch_state.ordered_indexers);
            if dofs_per_loop > 1
                # TODO
                printerr("TODO: more dofs per loop");
                # VAR_INDEX = Int(mod(current_dof-1, dofs_per_loop) + 1);
                # multiplier, index_values, index_counts needed
            end
            if num_indices > 1
                for indi = num_indices:-1:1 # work in reverse
                    ind = finch_state.ordered_indexers[indi];
                    if indi == num_indices
                        if !(ind.symbol === :FINCHELEMENTS)
                            index_symbol = string(ind.symbol);
                            index_index = ind.tag;
                            index_loops = "for INDEX_VAL_$index_symbol = 1:num_$(index_symbol)_indices\n";
                            index_loop_ends = "end\n";
                            index_offset = "INDEX_VAL_$index_symbol";
                            index_values = "index_values[$(index_index)] = INDEX_VAL_$index_symbol\n";
                            multiplier = "num_$(index_symbol)_indices";
                        end
                        
                    else
                        if !(ind.symbol === :FINCHELEMENTS)
                            index_symbol = string(ind.symbol);
                            index_index = ind.tag;
                            index_loops *= "for INDEX_VAL_$index_symbol = 1:num_$(index_symbol)_indices\n";
                            index_loop_ends *= "end\n";
                            index_offset *= " + $(multiplier) * (INDEX_VAL_$index_symbol - 1)";
                            index_values *= "index_values[$(index_index)] = INDEX_VAL_$index_symbol\n";
                            multiplier *= " * num_$index_symbol";
                        end
                    end
                end
                index_offset = index_offset * " - 1";
                
            else
                index_offset = "0";
            end
            
            # Args for the kernel
            nargs = length(kernel_args);
            args_per_line = 5;
            variable_args = [];
            args_str = "(";
            for i=1:nargs
                args_str *= kernel_args[i];
                if i < nargs
                    args_str *= ", ";
                    if mod(i,args_per_line) == 0
                        args_str *= "\n            ";
                    end
                end
                if occursin("variables_", kernel_args[i])
                    push!(variable_args, kernel_args[i]);
                end
            end
            args_str *= ")\n";
            
            # update the variables on the gpu
            var_update = "# Send needed values back to gpu\n";
            for i=1:length(variable_args)
                vind = parse(Int, variable_args[i][11]);
                var_update *= "copyto!($(variable_args[i]), variables[$(vind)].values);\n"
            end
            var_update *= "CUDA.synchronize();\n"

            code *="

$(var_update)

# This is done on gpu
@cuda threads=256 blocks=min(4096,ceil(Int, dofs_global/256)) gpu_assembly_kernel$(args_str)

# Asynchronously compute boundary values on cpu
@timeit timer_output \"bdry_vals\" begin
next_bdry_index = 1;
for bi=1:nbids
nfaces = length(mesh.bdryface[bi]);
for fi=1:nfaces
fid = mesh.bdryface[bi][fi];
eid = mesh.face2element[1,fid];
fbid = mesh.bids[bi];
volume = geometric_factors.volume[eid]
area = geometric_factors.area[fid]
area_over_volume = (area / volume)

$(index_loops)
$(index_values)
index_offset = $(index_offset)

row_index = index_offset + 1 + dofs_per_node * (eid - 1);

apply_boundary_conditions_face_rhs(var, eid, fid, fbid, mesh, refel, geometric_factors, fv_info, prob, 
                                    t, dt, flux_tmp, bdry_done, index_offset, index_values)
#
# store it
boundary_flux[next_bdry_index] = flux_tmp[1] * area_over_volume;
boundary_dof_index[next_bdry_index] = row_index;
next_bdry_index += 1;
$(index_loop_ends)
end
end
end # timer bdry_vals

# Then get global_vector from gpu
CUDA.synchronize()
copyto!(global_vector, global_vector_gpu)
CUDA.synchronize()

# And add BCs to global vector
for update_i = 1:(num_bdry_faces * dofs_per_node)
row_index = boundary_dof_index[update_i]
global_vector[row_index] = global_vector[row_index] + boundary_flux[update_i];
end

"
            
        else
            for i=1:length(IR.parts)
                if typeof(IR.parts[i]) == IR_operation_node
                    code *= generate_from_IR_julia_gpu(IR.parts[i], kernel_args, IRtypes);
                    code *= "\n";
                else
                    code *= generate_from_IR_julia_gpu(IR.parts[i], kernel_args, IRtypes);
                end
            end
        end
        
    elseif node_type == IR_loop_node
        # while and for loops
        if IR.type == IRtypes.while_loop
            condition = generate_from_IR_julia_gpu(IR.last, kernel_args, IRtypes);
            code = string(IR.iterator) * " = " * string(IR.first) * ";\n";
            code *= "while (" * condition * ")\n";
            code *= string(IR.iterator) * " += 1;\n";
            code *= generate_from_IR_julia_gpu(IR.body, kernel_args, IRtypes);
            code *= "end\n";
            
        else # for loop
            if (IR.first == IR.last) && IR.first==0
                range = string(IR.collection);
            else
                range = string(IR.first) * ":" * string(IR.last);
            end
            
            code = "for "* string(IR.iterator) * " = " * range * "\n";
            code *= generate_from_IR_julia_gpu(IR.body, kernel_args, IRtypes);
            code *= "end\n";
        end
        
        
    elseif node_type == IR_conditional_node
        code = "if " * generate_from_IR_julia_gpu(IR.condition, kernel_args, IRtypes) * "\n";
        code *= generate_from_IR_julia_gpu(IR.body, kernel_args, IRtypes);
        if !(IR.elsepart===nothing)
            code *= "else\n" * generate_from_IR_julia_gpu(IR.elsepart, kernel_args, IRtypes);
        end
        code *= "end\n";
        
    elseif node_type == IR_comment_node
        code = "#= " * IR.string * " =#\n";
        
    elseif node_type <: IR_part
        code = IR_string(IR);
        
    else
        code = string(IR);
    end
    
    return code;
end

#=
A list of named ops that really needs to be shortened:
COEF_EVAL
KNOWN_VAR
ROWCOL_TO_INDEX
TIMER
FILL_ARRAY
INIT_MATRIX_IJ_FV
GLOBAL_FORM_MATRIX
GLOBAL_GATHER_VECTOR
GLOBAL_GATHER_SYSTEM
GHOST_EXCHANGE_FV
GLOBAL_SOLVE
GLOBAL_DISTRIBUTE_VECTOR
GATHER_VARS
BDRY_TO_VECTOR
BDRY_TO_VAR
SCATTER_VARS
LOCAL2GLOBAL
LOCAL2GLOBAL_VEC
ADD_GLOBAL_VECTOR_AND_NORM
UPDATE_GLOBAL_VECTOR_AND_NORM
LINALG_VEC_BLOCKS
LINALG_MAT_BLOCKS
LINALG_MATMAT_BLOCKS
LINALG_MATVEC_BLOCKS
LINALG_TDM
LINALG_Tv
=#
function generate_named_op_gpu(IR::IR_operation_node, kernel_args, IRtypes::Union{IR_entry_types, Nothing} = nothing)
    code = "";
    
    op = IR.args[1];
    if op === :PRINT_STRING
        if typeof(IR.args[2]) == String
            code = "print(\""*IR.args[2]*"\")"
        else
            code = "print(string("*generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes)*"))"
        end
        
    elseif op === :COEF_EVAL
        # Rather than use evaluate_coefficient, try calling one of the nested functions
        type_name = string(finch_state.config.float_type);
        use_eval = true;
        if typeof(IR.args[2]) == Int
            coef = finch_state.coefficients[IR.args[2]];
            if typeof(IR.args[3]) == Int
                val = coef.value[IR.args[3]];
                if typeof(val) <:Number
                    code = type_name*"("*string(val)*")";
                    use_eval = false;
                elseif typeof(val) == GenFunction
                    code = type_name*"("*val.name * "(x,y,z,t,"*string(IR.args[8])*", "*string(IR.args[9])*", index_values))";
                    use_eval = false;
                end
            else
                # component is unknown.
            end
        end
        
        if use_eval
            if type_name == "Float64"
                type_name = ""
            end
            code = type_name*"(evaluate_coefficient(coefficients[" * string(IR.args[2]) * "], " * string(IR.args[3]) * 
                ", " * string(IR.args[4]) * ", " * string(IR.args[5]) * ", " * string(IR.args[6]) * 
                ", " * string(IR.args[7]) * ", " * string(IR.args[8]) * ", "*string(IR.args[9])*", index_values))";
        end
        
    elseif op === :KNOWN_VAR
        if length(IR.args) > 3
            code = "variables[" * string(IR.args[2]) * "].values["*generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes)*", "*string(IR.args[4])*"]";
        elseif length(IR.args) == 2
            code = "variables[" * string(IR.args[2]) * "].values";
        end
        
    elseif op === :ROWCOL_TO_INDEX
        # code = string(IR.args[2]) * " + (" * string(IR.args[3]) * "-1)*" * string(IR.args[4]);
        code = string(IR.args[2]) * ", " * string(IR.args[3]);
        
    elseif op === :TIMER
        # A timer has two more args, the label and the content
        code = "@timeit timer_output \""*string(IR.args[2])*"\" begin\n";
        code *= generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes);
        code *= "\n" * "end # timer:"*string(IR.args[2])*"\n";
        
    elseif op === :FILL_ARRAY
        # args[2] is the array, args[3] is the value
        code = generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) * " .= " * string(IR.args[3]);
        
    elseif op === :INIT_MATRIX_IJ_FV
        # args[2] is the I matrix, args[3] is the J matrix
        code = "set_matrix_indices!(" * generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) * ", " * 
                    generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes) * ", dofs_per_node, mesh)";
        
    elseif op === :GLOBAL_FORM_MATRIX
        # This will eventually be more complicated, but let's do the default for now
        code = "global_matrix = sparse(@view(global_matrix_I[1:(next_nonzero_index-1)]), @view(global_matrix_J[1:(next_nonzero_index-1)]), @view(global_matrix_V[1:(next_nonzero_index-1)]));\n"
        # code = "global_matrix = sparse(global_matrix_I, global_matrix_J, global_matrix_V);\n"
        
    elseif op === :GLOBAL_GATHER_VECTOR
        if length(IR.args) > 1 && IR.args[2] === :FV
            # FV has elemental dofs
            code *= "global_vector = gather_system_FV(nothing, nothing, nothing, global_vector, mesh.nel_owned, dofs_per_node, config, buffers);"
        else
            # FE has nodal dofs
            code = "global_vector = gather_system(nothing, nothing, nothing, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config, buffers);"
        end
        
    elseif op === :GLOBAL_GATHER_SYSTEM
        if length(IR.args) > 1 && IR.args[2] === :FV
            # FV has elemental dofs
            # code *= "(global_matrix_I, global_matrix_J, global_matrix_V, global_vector) = gather_system_FV(global_matrix_I, global_matrix_J, global_matrix_V, global_vector, mesh.nel_owned, dofs_per_node, config, buffers);"
            code *= "(full_global_matrix, full_global_vector) = gather_system_FV(global_matrix_I, global_matrix_J, global_matrix_V, global_vector, num_elements, dofs_per_node, config, buffers);"
        else
            # FE has nodal dofs
            # code *= "(global_matrix_I, global_matrix_J, global_matrix_V, global_vector) = gather_system(global_matrix_I, global_matrix_J, global_matrix_V, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config, buffers);"
            code *= "(global_matrix, global_vector) = gather_system(global_matrix_I, global_matrix_J, global_matrix_V, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config, buffers);\n"
            code *= "if length(global_vector) > length(global_solution)\n"
            code *= "global_solution = zeros(length(global_vector));\n"
            code *= "end"
        end
        
    elseif op === :GHOST_EXCHANGE_FV
        if length(IR.args) < 2
            code *= "exchange_ghosts_fv(var, mesh, dofs_per_node, ti, config);"
        else
            code *= "exchange_ghosts_fv(var, mesh, dofs_per_node, "*generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes)*", config);"
        end
        
    elseif op === :GLOBAL_SOLVE
        # not in-place
        # code =  generate_from_IR_julia_gpu(IR.args[2], IRtypes) * " .= linear_system_solve("* generate_from_IR_julia_gpu(IR.args[3], IRtypes) *", "* generate_from_IR_julia_gpu(IR.args[4], IRtypes) *", config);"
        # in-place
        code = "linear_system_solve!("* generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes) *", "* generate_from_IR_julia_gpu(IR.args[4], kernel_args, IRtypes) *", "* generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) * ", config);"
        
    elseif op === :GLOBAL_DISTRIBUTE_VECTOR
        if length(IR.args) > 3 && IR.args[4] === :FV
            # FV has elemental dofs
            code = "distribute_solution_FV!("*generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) * ", "* generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes) * ", mesh.nel_owned, dofs_per_node, config, buffers);"
        else
            # FE has nodal dofs
            code = "distribute_solution!("*generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) * ", "* generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes) * ", nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config, buffers);"
        end
        
    elseif op === :GATHER_VARS
        # place variable arrays in global vector
        if length(IR.args) < 3
            code = generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) * " = get_var_vals(var, "* generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) * ", );";
        else
            code = generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) * 
                " = get_var_vals("* generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes) * ", "* generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) * ");";
        end
        
    elseif op === :BDRY_TO_VECTOR
        # FV_copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node, prob);
        if length(IR.args) < 3
            code = "copy_bdry_vals_to_vector(var, "* generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) *", mesh, dofs_per_node, prob);";
        else
            code = "copy_bdry_vals_to_vector("* generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes) *", "* 
                            generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) *", mesh, dofs_per_node, prob);";
        end
        
    elseif op === :BDRY_TO_VAR
        # copy_bdry_vals_to_variables(var, solution, mesh, dofs_per_node, prob, true)
        if length(IR.args) < 4
            code = "copy_bdry_vals_to_variables(var, "* generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) *
                            ", mesh, dofs_per_node, prob, "* generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes) *");";
        else
            code = "copy_bdry_vals_to_variables("* generate_from_IR_julia_gpu(IR.args[4], kernel_args, IRtypes) *", "* 
                            generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) *
                            ", mesh, dofs_per_node, prob, "* generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes) *");";
        end
        
    elseif op === :SCATTER_VARS
        # place global vector in variable arrays
        if length(IR.args) < 3
            code = "place_vector_in_vars(var, "* generate_from_IR_julia_gpu(IR.args[2], kernel_args, IRtypes) *");";
        else
            code = "place_vector_in_vars("* generate_from_IR_julia_gpu(IR.args[3], kernel_args, IRtypes) *", "* 
                            generate_from_IR_julia_gpu(IR.args[2], IRtypes) *");";
        end
        
    elseif op === :LOCAL2GLOBAL
        # put elemental matrix and vector in global system
        code = generate_from_IR_julia_gpu(generate_local_to_global_fem(IR.args[2]), kernel_args, IRtypes);
        
    elseif op === :LOCAL2GLOBAL_VEC
        # put elemental vector in global system
        code = generate_from_IR_julia_gpu(generate_local_to_global_fem(IR.args[2], vec_only=true), kernel_args, IRtypes);
        
    elseif op === :ADD_GLOBAL_VECTOR_AND_NORM
        # u = u + du 
        # and find absolute and relative norms of du
        # args are: u, du, abs_residual, rel_residual
        code = generate_from_IR_julia_gpu(generate_residual_norms_and_update(
            IR.args[2], IR.args[3], IR.args[4], IR.args[5]), kernel_args, IRtypes);
        
    elseif op === :UPDATE_GLOBAL_VECTOR_AND_NORM
        # b = a
        # and find absolute and relative norms of (a-b)
        # args are: a, b, abs_residual, rel_residual
        code = generate_from_IR_julia_gpu(generate_difference_norms_and_update(
            IR.args[2], IR.args[3], IR.args[4], IR.args[5]), kernel_args, IRtypes);
        
    elseif op === :LINALG_VEC_BLOCKS
        # This sets blocks of a vector
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_julia_gpu(IR.args[4], kernel_args, IRtypes);
        
        init_lines = "";
        compute_lines = "";
        code = "";
        if blocksize > 1
            # loop over blocksize TODO
        end
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*2 + 1];
            comp  = IR.args[4 + (blk-1)*2 + 2];
            
            if blocksize > 1
                # TODO block loop : row_index = r_ind > 1 ? (string(r_ind-1)*"*"*string(blocksize)*" + row") : "row";
                row_index = string(r_ind);
            else # blocksize == 1
                row_index = string(r_ind);
            end
            
            content = generate_from_IR_julia_gpu(comp, kernel_args, IRtypes);
            code *= "$vecname[$row_index] = $content;\n";
        end
        
        if blocksize > 1
            # end # loop over blocksize
        end
    
    elseif op === :LINALG_MAT_BLOCKS
        # This sets blocks of a matrix
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_julia_gpu(IR.args[4], kernel_args, IRtypes);
        
        init_lines = "";
        compute_lines = "";
        code = "";
        if blocksize > 1
            # loop over blocksize TODO
        end
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*3 + 1];
            c_ind = IR.args[4 + (blk-1)*3 + 2];
            comp  = IR.args[4 + (blk-1)*3 + 3];
            
            if blocksize > 1
                # TODO block loop : row_index = r_ind > 1 ? (string(r_ind-1)*"*"*string(blocksize)*" + row") : "row";
                row_index = generate_from_IR_julia_gpu(r_ind, kernel_args, IRtypes);
                col_index = generate_from_IR_julia_gpu(c_ind, kernel_args, IRtypes);
            else # blocksize == 1
                row_index = generate_from_IR_julia_gpu(r_ind, kernel_args, IRtypes);
                col_index = generate_from_IR_julia_gpu(c_ind, kernel_args, IRtypes);
            end
            
            content = generate_from_IR_julia_gpu(comp, kernel_args, IRtypes);
            code *= "$vecname[$row_index, $col_index] = $content;\n";
        end
        
        if blocksize > 1
            # end # loop over blocksize
        end
    
    elseif op === :LINALG_MATMAT_BLOCKS
        # This is not just a matmat, 
        # it is the loop structure for a block matmat that
        # computes blocks that are computed like small 
        # mat*mat loops. A_ij = sum_k( something(i,j,k) )
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        matname = generate_from_IR_julia_gpu(IR.args[4], kernel_args, IRtypes);
        
        init_lines = "";
        compute_lines = "";
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*3 + 1];
            c_ind = IR.args[4 + (blk-1)*3 + 2];
            comp  = IR.args[4 + (blk-1)*3 + 3];
            row_index = r_ind > 1 ? (string(r_ind-1)*"*nodes_per_element + row") : "row";
            col_index = c_ind > 1 ? (string(c_ind-1)*"*nodes_per_element + col") : "col";
            content = generate_from_IR_julia_gpu(comp, kernel_args, IRtypes);
            
            init_lines *=    "$matname[$row_index, $col_index] = 0;\n";
            compute_lines *= "$matname[$row_index, $col_index] += $content;\n";
        end
        
        # content = generate_from_IR_julia_gpu(IR.args[6], IRtypes);
        code = "
@inbounds begin
for row=1:nodes_per_element
for col=1:nodes_per_element
$init_lines
for i=1:qnodes_per_element
$compute_lines
end
end
end
end";
    
    elseif op === :LINALG_MATVEC_BLOCKS
        # This is not just a matvec, 
        # it is the loop structure for a block matvec that
        # computes blocks that are computed like small 
        # mat*vec loops. b_i = sum_j( something(i,j) )
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_julia_gpu(IR.args[4], kernel_args, IRtypes);
        
        init_lines = "";
        compute_lines = "";
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*2 + 1];
            comp  = IR.args[4 + (blk-1)*2 + 2];
            row_index = r_ind > 1 ? (string(r_ind-1)*"*nodes_per_element + row") : "row";
            content = generate_from_IR_julia_gpu(comp, kernel_args, IRtypes);
            
            init_lines *=    "$vecname[$row_index] = 0;\n";
            compute_lines *= "$vecname[$row_index] += $content;\n";
        end
        
        # content = generate_from_IR_julia_gpu(IR.args[5], IRtypes);
        code = "
@inbounds begin
for row=1:nodes_per_element
$init_lines
for col=1:qnodes_per_element
$compute_lines
end
end
end";
    
    elseif op === :LINALG_TDM
        # Tcode = generate_from_IR_julia_gpu(IR.args[2], IRtypes) * "[i + (row-1)*qnodes_per_element]";
        # Mcode = generate_from_IR_julia_gpu(IR.args[4], IRtypes) * "[i + (col-1)*qnodes_per_element]";
        # Dcode = generate_from_IR_julia_gpu(apply_indexed_access(IR.args[3], [:i], IRtypes), IRtypes);
        # code = "*(" * Tcode *", "* Dcode *", "* Mcode * ")";
        
        i_index = :row;
        j_index = :col;
        k_index = :i;
        
        code = generate_from_IR_julia_gpu(generate_linalg_TDM_product(IR.args[2], IR.args[3], IR.args[4], i_index, j_index, k_index), kernel_args, IRtypes);
            
    elseif op === :LINALG_Tv
        # Tcode = generate_from_IR_julia_gpu(IR.args[2], IRtypes) * "[col + (row-1)*qnodes_per_element]";
        # vcode = generate_from_IR_julia_gpu(apply_indexed_access(IR.args[3], [:col], IRtypes), IRtypes);
        # code = "*(" * Tcode *", "* vcode * ")";
        
        i_index = :row;
        j_index = :col;
        
        code = generate_from_IR_julia_gpu(generate_linalg_Tv_product(IR.args[2], IR.args[3], i_index, j_index), kernel_args, IRtypes);
    end
    
    return code;
end

###############################################################################################################################
###############################################################################################################################
function generate_named_op_gpu_kernel(IR::IR_operation_node, IRtypes::Union{IR_entry_types, Nothing} = nothing)
    code = "";
    
    op = IR.args[1];
    if op === :PRINT_STRING
        # don't print from GPU
        
    elseif op === :COEF_EVAL
        # Rather than use evaluate_coefficient, try calling one of the nested functions
        type_name = string(finch_state.config.float_type);
        use_eval = true;
        if typeof(IR.args[2]) == Int
            coef = finch_state.coefficients[IR.args[2]];
            if typeof(IR.args[3]) == Int
                val = coef.value[IR.args[3]];
                if typeof(val) <:Number
                    code = type_name*"("*string(val)*")";
                    use_eval = false;
                elseif typeof(val) == GenFunction
                    code = type_name*"("*val.name * "(x,y,z,t,"*string(IR.args[8])*", "*string(IR.args[9])*", index_values))";
                    use_eval = false;
                end
            else
                # component is index
                code = "coefficients_$(IR.args[2])_value_gpu[" * string(IR.args[3]) * "]"
                use_eval = false;
            end
        end
        
        if use_eval
            if type_name == "Float64"
                type_name = ""
            end
            code = type_name*"(evaluate_coefficient(coefficients[" * string(IR.args[2]) * "], " * string(IR.args[3]) * 
                ", " * string(IR.args[4]) * ", " * string(IR.args[5]) * ", " * string(IR.args[6]) * 
                ", " * string(IR.args[7]) * ", " * string(IR.args[8]) * ", "*string(IR.args[9])*", index_values))";
        end
        
    elseif op === :KNOWN_VAR
        if length(IR.args) > 3
            code = "variables_" * string(IR.args[2]) * "_values_gpu["*generate_from_IR_gpu_assembly(IR.args[3], IRtypes)*", "*string(IR.args[4])*"]";
        elseif length(IR.args) == 2
            code = "variables_" * string(IR.args[2]) * "_values_gpu";
        end
        
    elseif op === :ROWCOL_TO_INDEX
        # code = string(IR.args[2]) * " + (" * string(IR.args[3]) * "-1)*" * string(IR.args[4]);
        code = string(IR.args[2]) * ", " * string(IR.args[3]);
        
    elseif op === :TIMER
        # Don't do on GPU
        code *= generate_from_IR_gpu_assembly(IR.args[3], IRtypes);
        
    elseif op === :FILL_ARRAY
        # args[2] is the array, args[3] is the value
        # code = generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * " .= " * string(IR.args[3]);
        if typeof(IR.args[2]) == IR_data_node && length(IR.args[2].size) > 0
            code = generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * " .= " * string(IR.args[3]);
        else
            code = generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * " = " * string(IR.args[3]);
        end
        
    elseif op === :INIT_MATRIX_IJ_FV
        # Don't do on GPU
        
    elseif op === :GLOBAL_FORM_MATRIX
        # This will eventually be more complicated, but let's do the default for now
        code = "global_matrix = sparse(@view(global_matrix_I[1:(next_nonzero_index-1)]), @view(global_matrix_J[1:(next_nonzero_index-1)]), @view(global_matrix_V[1:(next_nonzero_index-1)]));\n"
        # code = "global_matrix = sparse(global_matrix_I, global_matrix_J, global_matrix_V);\n"
        
    elseif op === :GLOBAL_GATHER_VECTOR
        if length(IR.args) > 1 && IR.args[2] === :FV
            # FV has elemental dofs
            code *= "global_vector = gather_system_FV(nothing, nothing, nothing, global_vector, mesh.nel_owned, dofs_per_node, config, buffers);"
        else
            # FE has nodal dofs
            code = "global_vector = gather_system(nothing, nothing, nothing, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config, buffers);"
        end
        
    elseif op === :GLOBAL_GATHER_SYSTEM
        if length(IR.args) > 1 && IR.args[2] === :FV
            # FV has elemental dofs
            # code *= "(global_matrix_I, global_matrix_J, global_matrix_V, global_vector) = gather_system_FV(global_matrix_I, global_matrix_J, global_matrix_V, global_vector, mesh.nel_owned, dofs_per_node, config, buffers);"
            code *= "(full_global_matrix, full_global_vector) = gather_system_FV(global_matrix_I, global_matrix_J, global_matrix_V, global_vector, num_elements, dofs_per_node, config, buffers);"
        else
            # FE has nodal dofs
            # code *= "(global_matrix_I, global_matrix_J, global_matrix_V, global_vector) = gather_system(global_matrix_I, global_matrix_J, global_matrix_V, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config, buffers);"
            code *= "(global_matrix, global_vector) = gather_system(global_matrix_I, global_matrix_J, global_matrix_V, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config, buffers);\n"
            code *= "if length(global_vector) > length(global_solution)\n"
            code *= "global_solution = zeros(length(global_vector));\n"
            code *= "end"
        end
        
    elseif op === :GHOST_EXCHANGE_FV
        if length(IR.args) < 2
            code *= "exchange_ghosts_fv(var, mesh, dofs_per_node, ti, config);"
        else
            code *= "exchange_ghosts_fv(var, mesh, dofs_per_node, "*generate_from_IR_gpu_assembly(IR.args[2], IRtypes)*", config);"
        end
        
    elseif op === :GLOBAL_SOLVE
        # not in-place
        # code =  generate_from_IR_julia_gpu(IR.args[2], IRtypes) * " .= linear_system_solve("* generate_from_IR_julia_gpu(IR.args[3], IRtypes) *", "* generate_from_IR_julia_gpu(IR.args[4], IRtypes) *", config);"
        # in-place
        code = "linear_system_solve!("* generate_from_IR_gpu_assembly(IR.args[3], IRtypes) *", "* generate_from_IR_gpu_assembly(IR.args[4], IRtypes) *", "* generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * ", config);"
        
    elseif op === :GLOBAL_DISTRIBUTE_VECTOR
        if length(IR.args) > 3 && IR.args[4] === :FV
            # FV has elemental dofs
            code = "distribute_solution_FV!("*generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * ", "* generate_from_IR_gpu_assembly(IR.args[3], IRtypes) * ", mesh.nel_owned, dofs_per_node, config, buffers);"
        else
            # FE has nodal dofs
            code = "distribute_solution!("*generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * ", "* generate_from_IR_gpu_assembly(IR.args[3], IRtypes) * ", nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config, buffers);"
        end
        
    elseif op === :GATHER_VARS
        # place variable arrays in global vector
        if length(IR.args) < 3
            code = generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * " = get_var_vals(var, "* generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * ", );";
        else
            code = generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * 
                " = get_var_vals("* generate_from_IR_gpu_assembly(IR.args[3], IRtypes) * ", "* generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * ");";
        end
        
    elseif op === :BDRY_TO_VECTOR
        # FV_copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node, prob);
        if length(IR.args) < 3
            code = "copy_bdry_vals_to_vector(var, "* generate_from_IR_gpu_assembly(IR.args[2], IRtypes) *", mesh, dofs_per_node, prob);";
        else
            code = "copy_bdry_vals_to_vector("* generate_from_IR_gpu_assembly(IR.args[3], IRtypes) *", "* 
                            generate_from_IR_gpu_assembly(IR.args[2], IRtypes) *", mesh, dofs_per_node, prob);";
        end
        
    elseif op === :BDRY_TO_VAR
        # copy_bdry_vals_to_variables(var, solution, mesh, dofs_per_node, prob, true)
        if length(IR.args) < 4
            code = "copy_bdry_vals_to_variables(var, "* generate_from_IR_gpu_assembly(IR.args[2], IRtypes) *
                            ", mesh, dofs_per_node, prob, "* generate_from_IR_gpu_assembly(IR.args[3], IRtypes) *");";
        else
            code = "copy_bdry_vals_to_variables("* generate_from_IR_gpu_assembly(IR.args[4], IRtypes) *", "* 
                            generate_from_IR_gpu_assembly(IR.args[2], IRtypes) *
                            ", mesh, dofs_per_node, prob, "* generate_from_IR_gpu_assembly(IR.args[3], IRtypes) *");";
        end
        
    elseif op === :SCATTER_VARS
        # place global vector in variable arrays
        if length(IR.args) < 3
            code = "place_vector_in_vars(var, "* generate_from_IR_gpu_assembly(IR.args[2], IRtypes) *");";
        else
            code = "place_vector_in_vars("* generate_from_IR_gpu_assembly(IR.args[3], IRtypes) *", "* 
                            generate_from_IR_gpu_assembly(IR.args[2], IRtypes) *");";
        end
        
    elseif op === :LOCAL2GLOBAL
        # put elemental matrix and vector in global system
        # code = generate_from_IR_gpu_assembly(generate_local_to_global_fem(IR.args[2]), IRtypes);
        code = "# TODO: local 2 global system\n"
        
    elseif op === :LOCAL2GLOBAL_VEC
        # put elemental vector in global system
        # code = generate_from_IR_gpu_assembly(generate_local_to_global_fem(IR.args[2], vec_only=true), IRtypes);
        code = "global_vector_gpu[current_dof] = source_gpu + flux_gpu\n"
        
    elseif op === :ADD_GLOBAL_VECTOR_AND_NORM
        # u = u + du 
        # and find absolute and relative norms of du
        # args are: u, du, abs_residual, rel_residual
        code = generate_from_IR_gpu_assembly(generate_residual_norms_and_update(
            IR.args[2], IR.args[3], IR.args[4], IR.args[5]), IRtypes);
        
    elseif op === :UPDATE_GLOBAL_VECTOR_AND_NORM
        # b = a
        # and find absolute and relative norms of (a-b)
        # args are: a, b, abs_residual, rel_residual
        code = generate_from_IR_gpu_assembly(generate_difference_norms_and_update(
            IR.args[2], IR.args[3], IR.args[4], IR.args[5]), IRtypes);
        
    elseif op === :LINALG_VEC_BLOCKS
        # This sets blocks of a vector
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_gpu_assembly(IR.args[4], IRtypes);
        
        init_lines = "";
        compute_lines = "";
        code = "";
        if blocksize > 1
            # loop over blocksize TODO
        end
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*2 + 1];
            comp  = IR.args[4 + (blk-1)*2 + 2];
            
            if blocksize > 1
                # TODO block loop : row_index = r_ind > 1 ? (string(r_ind-1)*"*"*string(blocksize)*" + row") : "row";
                row_index = string(r_ind);
            else # blocksize == 1
                row_index = string(r_ind);
            end
            
            content = generate_from_IR_gpu_assembly(comp, IRtypes);
            
            # This is where multiple dofs per loop will need to be figured out
            # for now, assume 1
            #code *= "$vecname[$row_index] = $content;\n";
            code *= "$vecname = $content;\n";
        end
        
        if blocksize > 1
            # end # loop over blocksize
        end
    
    elseif op === :LINALG_MAT_BLOCKS
        # This sets blocks of a matrix
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_gpu_assembly(IR.args[4], IRtypes);
        
        init_lines = "";
        compute_lines = "";
        code = "";
        if blocksize > 1
            # loop over blocksize TODO
        end
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*3 + 1];
            c_ind = IR.args[4 + (blk-1)*3 + 2];
            comp  = IR.args[4 + (blk-1)*3 + 3];
            
            if blocksize > 1
                # TODO block loop : row_index = r_ind > 1 ? (string(r_ind-1)*"*"*string(blocksize)*" + row") : "row";
                row_index = generate_from_IR_gpu_assembly(r_ind, IRtypes);
                col_index = generate_from_IR_gpu_assembly(c_ind, IRtypes);
            else # blocksize == 1
                row_index = generate_from_IR_gpu_assembly(r_ind, IRtypes);
                col_index = generate_from_IR_gpu_assembly(c_ind, IRtypes);
            end
            
            content = generate_from_IR_gpu_assembly(comp, IRtypes);
            code *= "$vecname[$row_index, $col_index] = $content;\n";
        end
        
        if blocksize > 1
            # end # loop over blocksize
        end
    
    elseif op === :LINALG_MATMAT_BLOCKS
        # This is not just a matmat, 
        # it is the loop structure for a block matmat that
        # computes blocks that are computed like small 
        # mat*mat loops. A_ij = sum_k( something(i,j,k) )
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        matname = generate_from_IR_gpu_assembly(IR.args[4], IRtypes);
        
        init_lines = "";
        compute_lines = "";
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*3 + 1];
            c_ind = IR.args[4 + (blk-1)*3 + 2];
            comp  = IR.args[4 + (blk-1)*3 + 3];
            row_index = r_ind > 1 ? (string(r_ind-1)*"*nodes_per_element + row") : "row";
            col_index = c_ind > 1 ? (string(c_ind-1)*"*nodes_per_element + col") : "col";
            content = generate_from_IR_gpu_assembly(comp, IRtypes);
            
            init_lines *=    "$matname[$row_index, $col_index] = 0;\n";
            compute_lines *= "$matname[$row_index, $col_index] += $content;\n";
        end
        
        # content = generate_from_IR_gpu_assembly(IR.args[6], IRtypes);
        code = "
@inbounds begin
for row=1:nodes_per_element
for col=1:nodes_per_element
$init_lines
for i=1:qnodes_per_element
$compute_lines
end
end
end
end";
    
    elseif op === :LINALG_MATVEC_BLOCKS
        # This is not just a matvec, 
        # it is the loop structure for a block matvec that
        # computes blocks that are computed like small 
        # mat*vec loops. b_i = sum_j( something(i,j) )
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_gpu_assembly(IR.args[4], IRtypes);
        
        init_lines = "";
        compute_lines = "";
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*2 + 1];
            comp  = IR.args[4 + (blk-1)*2 + 2];
            row_index = r_ind > 1 ? (string(r_ind-1)*"*nodes_per_element + row") : "row";
            content = generate_from_IR_gpu_assembly(comp, IRtypes);
            
            init_lines *=    "$vecname[$row_index] = 0;\n";
            compute_lines *= "$vecname[$row_index] += $content;\n";
        end
        
        # content = generate_from_IR_gpu_assembly(IR.args[5], IRtypes);
        code = "
@inbounds begin
for row=1:nodes_per_element
$init_lines
for col=1:qnodes_per_element
$compute_lines
end
end
end";
    
    elseif op === :LINALG_TDM
        # Tcode = generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * "[i + (row-1)*qnodes_per_element]";
        # Mcode = generate_from_IR_gpu_assembly(IR.args[4], IRtypes) * "[i + (col-1)*qnodes_per_element]";
        # Dcode = generate_from_IR_gpu_assembly(apply_indexed_access(IR.args[3], [:i], IRtypes), IRtypes);
        # code = "*(" * Tcode *", "* Dcode *", "* Mcode * ")";
        
        i_index = :row;
        j_index = :col;
        k_index = :i;
        
        code = generate_from_IR_gpu_assembly(generate_linalg_TDM_product(IR.args[2], IR.args[3], IR.args[4], i_index, j_index, k_index), IRtypes);
            
    elseif op === :LINALG_Tv
        # Tcode = generate_from_IR_gpu_assembly(IR.args[2], IRtypes) * "[col + (row-1)*qnodes_per_element]";
        # vcode = generate_from_IR_gpu_assembly(apply_indexed_access(IR.args[3], [:col], IRtypes), IRtypes);
        # code = "*(" * Tcode *", "* vcode * ")";
        
        i_index = :row;
        j_index = :col;
        
        code = generate_from_IR_gpu_assembly(generate_linalg_Tv_product(IR.args[2], IR.args[3], i_index, j_index), IRtypes);
    end
    
    return code;
end
