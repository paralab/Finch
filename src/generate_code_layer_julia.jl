#=
Take an IR and generate julia code.
Two main tasks are:
- directly translate the IR into code
- create any target specfic support code to surround it

The support code should define sevaral numbers and data structures
with specific names that are used in the IR, such as mesh, refel, num_elements, etc.
(TODO create a concrete list of these required parts)
=#

# This creates the full code string including support code.
# If desired this can make the code as a complete, stand-alone function
# or as just the body of a function that will be generated(default).
function generate_code_layer_julia(var::Vector{Variable}, IR::IR_part, solver, wrap_in_function=true)
    # This will hold the code string to be returned
    code ="";
    
    # Set up useful numbers
    # Count variables, dofs, and store offsets
    varcount = 1;
    dofs_per_node = var[1].total_components;
    dofs_per_loop = length(var[1].symvar);
    offset_ind = [0];
    if length(var) > 1
        varcount = length(var);
        offset_ind = zeros(Int, varcount);
        for i=2:length(var)
            offset_ind[i] = dofs_per_node;
            dofs_per_node += var[i].total_components;
            dofs_per_loop += length(var[i].symvar);
        end
    end
    
    # 
    if solver == CG || solver == DG
        args = "(var, mesh, refel, geometric_factors, config, coefficients, variables, test_functions, indexers, prob, time_stepper, nl_var=nothing)";
        disc_specific = "
    # FEM specific pieces
    Q = refel.Q;
    wg = refel.wg;
    
    # For partitioned meshes
    (partitioned_order, partitioned_sizes) = get_partitioned_ordering(dofs_per_node, mesh, config);
    "
        
    elseif solver == FV
        args = "(var, mesh, refel, geometric_factors, fv_info, config, coefficients, variables, test_functions, indexers, prob, time_stepper, nl_var=nothing)";
        disc_specific = "
    # FVM specific pieces
    
    "
        
    else
        printerr("This colver type is not ready: "*solver);
        args = "()";
        disc_specific = "";
    end
    
    if wrap_in_function
        code = "function generated_solve_function_for_"*string(var[1].symbol) * args * "\n";
    end
    
    code *="
    # pre/post step functions if defined
    pre_step_function = prob.pre_step_function;
    post_step_function = prob.post_step_function;
    
    # Prepare some useful numbers
    dofs_per_node = "*string(dofs_per_node)*";
    dofs_per_loop = "*string(dofs_per_loop)*";
    dof_offsets = "*string(offset_ind)*";
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
    dofs_per_element = dofs_per_node * nodes_per_element;
    local_system_size = dofs_per_loop * nodes_per_element;
    
    "
    code *= disc_specific;
    
    # The assembly is taken from the IR
    code *= generate_from_IR_julia(IR);
    
    # Post assembly part
    code *="
    
    return nothing;\n"
    
    if wrap_in_function
        code *= "\nend # function\n";
    end
    
    return code;
end

# Directly translate IR into julia code
# Named operators are treated specially
function generate_from_IR_julia(IR, IRtypes::Union{IR_entry_types, Nothing} = nothing, indent="")
    if IRtypes === nothing
        IRtypes = IR_entry_types();
    end
    code = "";
    node_type = typeof(IR);
    
    if node_type == IR_data_node # data nodes will look like "x" or "x[i,2]" or "x[a[1],17]"
        var_name = string(IR.label);
        
        if length(IR.size) == 0
            code = var_name;
        elseif length(IR.size) > 0 && length(IR.index) > 0
            code = var_name*"[";
            for i=1:length(IR.index)
                if typeof(IR.index[i]) <: IR_part
                    code *= generate_from_IR_julia(IR.index[i], IRtypes);
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
            code = indent * "zeros("*IRtypes.name[IR.args[1]];
            for i=2:length(IR.args)
                code *= ", ";
                if typeof(IR.args[i]) <: IR_part
                    code *= generate_from_IR_julia(IR.args[i], IRtypes);
                else
                    code *= string(IR.args[i]);
                end
            end
            code *= ")";
            
        elseif IR.type == IRtypes.assign_op # x = ...
            if typeof(IR.args[1]) <: IR_part
                code = indent * generate_from_IR_julia(IR.args[1], IRtypes);
            else
                code = indent * string(IR.args[1]);
            end
            code *= " = ";
            if typeof(IR.args[2]) <: IR_part
                code *= generate_from_IR_julia(IR.args[2], IRtypes);
            else
                code *= string(IR.args[2]);
            end
            
        elseif IR.type == IRtypes.function_op # f(...)
            if typeof(IR.args[1]) <: IR_part
                code = generate_from_IR_julia(IR.args[1], IRtypes);
            else
                code = string(IR.args[1]);
            end
            code *= "(";
            for i=2:length(IR.args)
                if typeof(IR.args[i]) <: IR_part
                    code *= generate_from_IR_julia(IR.args[i], IRtypes);
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
                        code *= generate_from_IR_julia(IR.args[2], IRtypes);
                    else
                        code *= string(IR.args[2]);
                    end
                    code *= ")";
                else # a + b + c + d
                    code = "(";
                    for i=2:length(IR.args)
                        if typeof(IR.args[i]) <: IR_part
                            code *= generate_from_IR_julia(IR.args[i], IRtypes);
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
                        code *= generate_from_IR_julia(IR.args[i], IRtypes);
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
            code = string(IR.args[1]) * "." * generate_from_IR_julia(IR.args[2], IRtypes);
            
        elseif IR.type == IRtypes.named_op # handled case by case
            code = generate_named_op(IR, IRtypes, indent);
        end
        
    elseif node_type == IR_block_node # A collection of statements. Do them one line at a time
        code = "";
        for i=1:length(IR.parts)
            if typeof(IR.parts[i]) == IR_operation_node
                code *= generate_from_IR_julia(IR.parts[i], IRtypes, indent * "    ");
                code *= "\n";
            elseif typeof(IR.parts[i]) == IR_block_node
                code *= generate_from_IR_julia(IR.parts[i], IRtypes, indent);
            else
                code *= generate_from_IR_julia(IR.parts[i], IRtypes, indent*"    ");
            end
            
        end
        
    elseif node_type == IR_loop_node
        # while and for loops
        if IR.type == IRtypes.while_loop
            condition = generate_from_IR_julia(IR.last, IRtypes);
            code = indent * string(IR.iterator) * " = " * string(IR.first) * ";\n";
            code *= indent * "while (" * condition * ")\n";
            code *= indent * "    " * string(IR.iterator) * " += 1;\n";
            code *= generate_from_IR_julia(IR.body, IRtypes, indent);
            code *= indent * "end\n";
            
        else # for loop
            if (IR.first == IR.last) && IR.first==0
                range = string(IR.collection);
            else
                range = string(IR.first) * ":" * string(IR.last);
            end
            
            code = indent * "for "* string(IR.iterator) * " = " * range * "\n";
            code *= generate_from_IR_julia(IR.body, IRtypes, indent);
            code *= indent * "end\n";
        end
        
        
    elseif node_type == IR_conditional_node
        code = indent * "if " * generate_from_IR_julia(IR.condition, IRtypes) * "\n";
        code *= generate_from_IR_julia(IR.body, IRtypes, indent);
        if !(IR.elsepart===nothing)
            code *= indent * "else\n" * generate_from_IR_julia(IR.elsepart, IRtypes, indent);
        end
        code *= indent * "end\n";
        
    elseif node_type == IR_comment_node
        code = indent * "#= " * IR.string * " =#\n";
        
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
function generate_named_op(IR::IR_operation_node, IRtypes::Union{IR_entry_types, Nothing} = nothing, indent="")
    code = "";
    
    op = IR.args[1];
    if op === :COEF_EVAL
        code = "evaluate_coefficient(coefficients[" * string(IR.args[2]) * "], " * string(IR.args[3]) * 
            ", " * string(IR.args[4]) * ", " * string(IR.args[5]) * ", " * string(IR.args[6]) * 
            ", " * string(IR.args[7]) * ", " * string(IR.args[8]) * ")";
        
    elseif op === :KNOWN_VAR
        code = "variables[" * string(IR.args[2]) * "].values["*generate_from_IR_julia(IR.args[3], IRtypes)*", "*string(IR.args[4])*"]";
        
    elseif op === :ROWCOL_TO_INDEX
        # code = string(IR.args[2]) * " + (" * string(IR.args[3]) * "-1)*" * string(IR.args[4]);
        code = string(IR.args[2]) * ", " * string(IR.args[3]);
        
    elseif op === :TIMER
        # A timer has two more args, the label and the content
        code = indent * "@timeit timer_output \""*string(IR.args[2])*"\" begin\n";
        code *= generate_from_IR_julia(IR.args[3], IRtypes, indent);
        code *= "\n" * indent * "end # timer:"*string(IR.args[2])*"\n";
        
    elseif op === :FILL_ARRAY
        # args[2] is the array, args[3] is the value
        code = indent * generate_from_IR_julia(IR.args[2], IRtypes) * " .= " * string(IR.args[3]);
        
    elseif op === :INIT_MATRIX_IJ_FV
        # args[2] is the I matrix, args[3] is the J matrix
        code = indent * "set_matrix_indices!(" * generate_from_IR_julia(IR.args[2], IRtypes) * ", " * 
                    generate_from_IR_julia(IR.args[3], IRtypes) * ", dofs_per_node, mesh)";
        
    elseif op === :GLOBAL_FORM_MATRIX
        # This will eventually be more complicated, but let's do the default for now
        code = indent * "global_matrix = sparse(global_matrix_I, global_matrix_J, global_matrix_V);\n"
        
    elseif op === :GLOBAL_GATHER_VECTOR
        code = indent * "global_vector = gather_system(nothing, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config);"
        
    elseif op === :GLOBAL_GATHER_SYSTEM
        code *= indent * "(global_matrix, global_vector) = gather_system(global_matrix, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config);"
        
    elseif op === :GHOST_EXCHANGE_FV
        code *= indent * "exchange_ghosts_fv(var, mesh, dofs_per_node, ti);"
        
    elseif op === :GLOBAL_SOLVE
        code = indent * generate_from_IR_julia(IR.args[2], IRtypes) * " = linear_system_solve("* generate_from_IR_julia(IR.args[3], IRtypes) *", "* generate_from_IR_julia(IR.args[4], IRtypes) *", config);"
        
    elseif op === :GLOBAL_DISTRIBUTE_VECTOR
        code = indent * generate_from_IR_julia(IR.args[2], IRtypes) * " = distribute_solution("*generate_from_IR_julia(IR.args[2], IRtypes) * ", nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config);"
        
    elseif op === :GATHER_VARS
        # place variable arrays in global vector
        if length(IR.args) < 3
            code = indent * generate_from_IR_julia(IR.args[2], IRtypes) * " = get_var_vals(var, "* generate_from_IR_julia(IR.args[2], IRtypes) * ");";
        else
            code = indent * generate_from_IR_julia(IR.args[2], IRtypes) * 
                " = get_var_vals("* generate_from_IR_julia(IR.args[3], IRtypes) * ", "* generate_from_IR_julia(IR.args[2], IRtypes) * ");";
        end
        
    elseif op === :BDRY_TO_VECTOR
        # FV_copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node);
        if length(IR.args) < 3
            code = indent * "copy_bdry_vals_to_vector(var, "* generate_from_IR_julia(IR.args[2], IRtypes) *", mesh, dofs_per_node);";
        else
            code = indent * "copy_bdry_vals_to_vector("* generate_from_IR_julia(IR.args[3], IRtypes) *", "* 
                            generate_from_IR_julia(IR.args[2], IRtypes) *", mesh, dofs_per_node);";
        end
        
    elseif op === :BDRY_TO_VAR
        # copy_bdry_vals_to_variables(var, solution, mesh, dofs_per_node, true)
        if length(IR.args) < 4
            code = indent * "copy_bdry_vals_to_variables(var, "* generate_from_IR_julia(IR.args[2], IRtypes) *
                            ", mesh, dofs_per_node, "* generate_from_IR_julia(IR.args[3], IRtypes) *");";
        else
            code = indent * "copy_bdry_vals_to_variables("* generate_from_IR_julia(IR.args[4], IRtypes) *", "* 
                            generate_from_IR_julia(IR.args[2], IRtypes) *
                            ", mesh, dofs_per_node, "* generate_from_IR_julia(IR.args[3], IRtypes) *");";
        end
        
    elseif op === :SCATTER_VARS
        # place global vector in variable arrays
        # place_sol_in_vars(var, sol, stepper);
        if length(IR.args) < 3
            code = indent * "place_sol_in_vars(var, "* generate_from_IR_julia(IR.args[2], IRtypes) *");";
        else
            code = indent * "place_sol_in_vars("* generate_from_IR_julia(IR.args[3], IRtypes) *", "* 
                            generate_from_IR_julia(IR.args[2], IRtypes) *");";
        end
        
    elseif op === :LOCAL2GLOBAL
        # put elemental matrix and vector in global system
        code = generate_from_IR_julia(IntermediateRepresentation.generate_local_to_global_fem(IR.args[2]), IRtypes, indent);
        
    elseif op === :LOCAL2GLOBAL_VEC
        # put elemental vector in global system
        code = generate_from_IR_julia(IntermediateRepresentation.generate_local_to_global_fem(IR.args[2], vec_only=true), IRtypes, indent);
        
    elseif op === :ADD_GLOBAL_VECTOR_AND_NORM
        # u = u + du 
        # and find absolute and relative norms of du
        # args are: u, du, abs_residual, rel_residual
        code = generate_from_IR_julia(IntermediateRepresentation.generate_residual_norms_and_update(
            IR.args[2], IR.args[3], IR.args[4], IR.args[5]), IRtypes, indent);
        
    elseif op === :UPDATE_GLOBAL_VECTOR_AND_NORM
        # b = a
        # and find absolute and relative norms of (a-b)
        # args are: a, b, abs_residual, rel_residual
        code = generate_from_IR_julia(IntermediateRepresentation.generate_difference_norms_and_update(
            IR.args[2], IR.args[3], IR.args[4], IR.args[5]), IRtypes, indent);
        
    elseif op === :LINALG_VEC_BLOCKS
        # This sets blocks of a vector
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_julia(IR.args[4], IRtypes);
        
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
            
            content = generate_from_IR_julia(comp, IRtypes);
            code *= indent * "$vecname[$row_index] = $content;\n";
        end
        
        if blocksize > 1
            # end # loop over blocksize
        end
    
    elseif op === :LINALG_MAT_BLOCKS
        # This sets blocks of a matrix
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_julia(IR.args[4], IRtypes);
        
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
                row_index = generate_from_IR_julia(r_ind, IRtypes);
                col_index = generate_from_IR_julia(c_ind, IRtypes);
            else # blocksize == 1
                row_index = generate_from_IR_julia(r_ind, IRtypes);
                col_index = generate_from_IR_julia(c_ind, IRtypes);
            end
            
            content = generate_from_IR_julia(comp, IRtypes);
            code *= indent * "$vecname[$row_index, $col_index] = $content;\n";
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
        matname = generate_from_IR_julia(IR.args[4], IRtypes);
        
        init_lines = "";
        compute_lines = "";
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*3 + 1];
            c_ind = IR.args[4 + (blk-1)*3 + 2];
            comp  = IR.args[4 + (blk-1)*3 + 3];
            row_index = r_ind > 1 ? (string(r_ind-1)*"*nodes_per_element + row") : "row";
            col_index = c_ind > 1 ? (string(c_ind-1)*"*nodes_per_element + col") : "col";
            content = generate_from_IR_julia(comp, IRtypes);
            
            init_lines *=    "                $matname[$row_index, $col_index] = 0;\n";
            compute_lines *= "                    $matname[$row_index, $col_index] += $content;\n";
        end
        
        # content = generate_from_IR_julia(IR.args[6], IRtypes);
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
        vecname = generate_from_IR_julia(IR.args[4], IRtypes);
        
        init_lines = "";
        compute_lines = "";
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*2 + 1];
            comp  = IR.args[4 + (blk-1)*2 + 2];
            row_index = r_ind > 1 ? (string(r_ind-1)*"*nodes_per_element + row") : "row";
            content = generate_from_IR_julia(comp, IRtypes);
            
            init_lines *=    "            $vecname[$row_index] = 0;\n";
            compute_lines *= "                $vecname[$row_index] += $content;\n";
        end
        
        # content = generate_from_IR_julia(IR.args[5], IRtypes);
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
        # Tcode = generate_from_IR_julia(IR.args[2], IRtypes) * "[i + (row-1)*qnodes_per_element]";
        # Mcode = generate_from_IR_julia(IR.args[4], IRtypes) * "[i + (col-1)*qnodes_per_element]";
        # Dcode = generate_from_IR_julia(IntermediateRepresentation.apply_indexed_access(IR.args[3], [:i], IRtypes), IRtypes);
        # code = "*(" * Tcode *", "* Dcode *", "* Mcode * ")";
        
        i_index = :row;
        j_index = :col;
        k_index = :i;
        
        code = generate_from_IR_julia(IntermediateRepresentation.generate_linalg_TDM_product(IR.args[2], IR.args[3], IR.args[4], i_index, j_index, k_index), IRtypes);
            
    elseif op === :LINALG_Tv
        # Tcode = generate_from_IR_julia(IR.args[2], IRtypes) * "[col + (row-1)*qnodes_per_element]";
        # vcode = generate_from_IR_julia(IntermediateRepresentation.apply_indexed_access(IR.args[3], [:col], IRtypes), IRtypes);
        # code = "*(" * Tcode *", "* vcode * ")";
        
        i_index = :row;
        j_index = :col;
        
        code = generate_from_IR_julia(IntermediateRepresentation.generate_linalg_Tv_product(IR.args[2], IR.args[3], i_index, j_index), IRtypes);
    end
    
    return code;
end
