#=
This generates a set of C++ files from a given IR

The following two functions are required for a target:
1. get_external_language_elements() - file extensions, comment chars etc.
2. generate_external_files(var, IR) - Writes all files based on the configuration and IR

You will have access to anything in the Finch module scope because this file will be included there.
=#

function get_external_language_elements()
    file_extension = ".cpp";
    comment_char = "//";
    block_comment = ["/*"; "*/"];
    
    return (file_extension, comment_char, block_comment);
end

# Create all of the code files.
# The input are code strings generated above. They will be inserted into their appropriate files.
function generate_external_files(var, IR)
    config = finch_state.config;
    
    # If multiple procs, only rank 0 does this
    if config.num_procs == 1 || config.proc_rank == 0
        # Write the static files (see the end of this file)
        dendrite_static_files();
        # Build and readme files
        dendrite_build_files();
        
        # The generated code
        dendrite_main_file(var);
        dendrite_equation_file(var, IR);
        dendrite_boundary_file(var);
        dendrite_nodedata_file(var);
        # dendrite_inputdata_file(var);
        dendrite_refine_file(var);
    end
    
    if config.num_procs > 1
        MPI.Barrier(MPI.COMM_WORLD);
    end
end

#=
Translate the IR into code directly.
IR - The IR_part holding the full IR

Note that this should be independent of problem setup, discretization, etc.
It may depend on target-specific parameters.
This will probably only be called by one of the file generation functions.
=#
function generate_from_IR_external(IR, IRtypes::Union{IR_entry_types, Nothing} = nothing; indent="", declared=[])
    if IRtypes === nothing
        IRtypes = IR_entry_types();
    end
    code = "";
    node_type = typeof(IR);
    
    # Handle different IR part types
    if node_type == IR_data_node # data nodes will look like "x" or "x[i,2]" or "x[a[1],17]"
        var_name = string(IR.label);
        
        if length(IR.size) == 0 # scalars
            code = var_name;
        elseif length(IR.index) > 0 # indexed element in an array
            code = var_name*"[";
            for i=1:length(IR.index)
                if typeof(IR.index[i]) <: IR_part
                    code *= generate_from_IR_external(IR.index[i], IRtypes);
                elseif typeof(IR.index[i]) <: Int
                    code *= string(IR.index[i]-1);
                else
                    code *= string(IR.index[i]);
                end
                if i<length(IR.index)
                    code *= ", ";
                end
            end
            code *= "]";
        else # arrays without index
            code = var_name;
        end
        
    elseif node_type == IR_operation_node
        if IR.type == IRtypes.allocate_op # new <type>[dim1 * dim2...]
            typestring = cpp_type_name(IRtypes.name[IR.args[1]]);
            code = indent * "new " * typestring * "[";
            for i=2:length(IR.args)
                if typeof(IR.args[i]) <: IR_part
                    code *= generate_from_IR_external(IR.args[i], IRtypes);
                else
                    code *= string(IR.args[i]);
                end
                if i<length(IR.args)
                    code *= " * "
                end
            end
            code *= "]";
            
        elseif IR.type == IRtypes.assign_op # x = ...
            ## NOTE: Since this uses PETSC for global vectors and matrix
            ## They need special creation/assign steps. For now don't allocate 
            ## or directly assign to them here.
            if typeof(IR.args[1]) <: IR_data_node && (
                IR.args[1].label in [:global_matrix_I, :global_matrix_J, :global_matrix_V, :global_vector, :solution] ||
                (length(IR.args[1].size)>0 && IR.args[1].size[1] in [:allocated_nonzeros, :dofs_global]))
                
                return "GLOBAL_ARRAY_ASSIGNMENT";
            end
            
            # Note that if the symbol (IR.args[1]) is not yet on the declared list, 
            # it will be here.
            if typeof(IR.args[1]) <: IR_data_node
                typestring = cpp_type_name(IRtypes.name[IR.args[1].type]);
                name = IR.args[1].label;
                # check if it has been declared
                foundit = false;
                for d in declared
                    if name === d
                        foundit = true;
                        break;
                    end
                end
                if !foundit
                    code = indent * typestring * " " * generate_from_IR_external(IR.args[1], IRtypes);
                    push!(declared, name);
                else
                    code = indent * generate_from_IR_external(IR.args[1], IRtypes);
                end
                
            elseif typeof(IR.args[1]) <: IR_part
                code = indent * generate_from_IR_external(IR.args[1], IRtypes);
            else
                
                code = indent * string(IR.args[1]);
            end
            code *= " = ";
            if typeof(IR.args[2]) <: IR_part
                code *= generate_from_IR_external(IR.args[2], IRtypes);
            else
                code *= string(IR.args[2]);
            end
            
        elseif IR.type == IRtypes.math_assign_op # x += ...
            # What is the math op?
            math_op = IR.args[1];
            if math_op in [:+, :-, :/, :*]
                math_op = string(math_op);
            else
                printerr("Unsupported math op in math-assignment: "*string(math_op)*"=");
                math_op = "";
            end
            
            ## NOTE: Since this uses PETSC for global vectors and matrix
            ## They need special creation/assign steps. For now don't allocate 
            ## or directly assign to them here.
            if typeof(IR.args[2]) <: IR_data_node && (
                IR.args[2].label in [:global_matrix_I, :global_matrix_J, :global_matrix_V, :global_vector, :solution] ||
                (length(IR.args[2].size)>0 && IR.args[2].size[1] in [:allocated_nonzeros, :dofs_global]))
                
                return "GLOBAL_ARRAY_ASSIGNMENT";
            end
            
            if typeof(IR.args[2]) <: IR_part
                code = indent * generate_from_IR_external(IR.args[2], IRtypes);
            else
                code = indent * string(IR.args[2]);
            end
            
            code *= " "*math_op*"= ";
            if typeof(IR.args[3]) <: IR_part
                code *= generate_from_IR_external(IR.args[3], IRtypes);
            else
                code *= string(IR.args[3]);
            end
            
        elseif IR.type == IRtypes.function_op # f(...)
            if typeof(IR.args[1]) <: IR_part
                code = generate_from_IR_external(IR.args[1], IRtypes);
            else
                code = string(IR.args[1]);
            end
            code *= "(";
            for i=2:length(IR.args)
                if typeof(IR.args[i]) <: IR_part
                    code *= generate_from_IR_external(IR.args[i], IRtypes);
                else
                    code *= string(IR.args[i]);
                end
                if i < length(IR.args)
                    code *= ", ";
                end
            end
            code *= ")";
            
        elseif IR.type == IRtypes.math_op # (a * b * 2) or sin(a)
            # is it one of +-*/
            if IR.args[1] in [:+, :-, :*, :/, :&&, :||, :<, :>, :(==), :(>=), :(<=)]
                if length(IR.args) < 3 # -a, +a
                    code = "(" * string(IR.args[1]);
                    if typeof(IR.args[2]) <: IR_part
                        code *= generate_from_IR_external(IR.args[2], IRtypes);
                    else
                        code *= string(IR.args[2]);
                    end
                    code *= ")";
                else # a + b + c + d
                    code = "(";
                    for i=2:length(IR.args)
                        if typeof(IR.args[i]) <: IR_part
                            code *= generate_from_IR_external(IR.args[i], IRtypes);
                        else
                            code *= string(IR.args[i]);
                        end
                        if i < length(IR.args)
                            code *= " " * string(IR.args[1]) * " ";
                        end
                    end
                    code *= ")";
                end
                
            elseif IR.args[1] === :^ # power is special
                code = "pow(";
                for i=2:length(IR.args)
                    if typeof(IR.args[i]) <: IR_part
                        code *= generate_from_IR_external(IR.args[i], IRtypes);
                    else
                        code *= string(IR.args[i]);
                    end
                    if i < length(IR.args)
                        code *= ", ";
                    end
                end
                code *= ")";
                
            else # sin(a)
                code = string(IR.args[1]) * "(";
                for i=2:length(IR.args)
                    if typeof(IR.args[i]) <: IR_part
                        code *= generate_from_IR_external(IR.args[i], IRtypes);
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
            code = string(IR.args[1]) * "." * generate_from_IR_external(IR.args[2], IRtypes);
            
        elseif IR.type == IRtypes.named_op # handled case by case
            code = generate_named_op(IR, IRtypes, indent);
        end
        
    elseif node_type == IR_block_node # A collection of statements. Do them one line at a time
        code = "";
        for i=1:length(IR.parts)
            if typeof(IR.parts[i]) == IR_operation_node
                code *= generate_from_IR_external(IR.parts[i], IRtypes, indent = indent);
                code *= ";\n";
            elseif typeof(IR.parts[i]) == IR_block_node
                code *= generate_from_IR_external(IR.parts[i], IRtypes, indent = indent);
            else
                code *= generate_from_IR_external(IR.parts[i], IRtypes, indent = indent);
            end
            
        end
        
    elseif node_type == IR_loop_node
        # while and for loops
        if IR.type == IRtypes.while_loop
            condition = generate_from_IR_external(IR.last, IRtypes);\
            # int iterator = first; while(condition){iterator++; body}
            code = indent * "int " * string(IR.iterator) * " = " * string(IR.first) * ";\n";
            code *= indent * "while(" * condition * "){\n";
            code *= indent * "    " * string(IR.iterator) * "++;\n";
            code *= generate_from_IR_external(IR.body, IRtypes, indent = indent*"    ");
            code *= indent * "}\n";
            
        else # for loop
            # for(unsigned int i = 1; i <= n; i++){}
            code = indent * "for(unsigned int "* string(IR.iterator) * " = " * string(IR.first) * 
                "; " * string(IR.iterator) * "<=" * string(IR.last) * "; " * string(IR.iterator) * "++){\n";
            code *= generate_from_IR_external(IR.body, IRtypes, indent = indent*"    ");
            code *= indent * "}\n";
        end
        
        
    elseif node_type == IR_conditional_node
        code = indent * "if(" * generate_from_IR_external(IR.condition, IRtypes) * "){\n";
        code *= generate_from_IR_external(IR.body, IRtypes, indent = indent*"    ");
        if !(IR.elsepart===nothing)
            code *= indent * "}else{\n" * generate_from_IR_external(IR.elsepart, IRtypes, indent = indent*"    ");
        end
        code *= indent * "}\n";
    elseif node_type == IR_comment_node
        code = indent * "// " * IR.string * " //\n";
    elseif node_type <: IR_part
        code = IR_string(IR);
        
    else
        code = string(IR);
    end
    
    return code;
end

#=
Some special operations are given a name.
What that means depends on the target.
=#
function generate_named_op(IR::IR_operation_node, IRtypes::Union{IR_entry_types, Nothing} = nothing, indent="")
    code = "";
    
    op = IR.args[1];
    if op === :COEF_EVAL
        # "genfunction_"*string(cval)*"(p)"
        gen_func_name = finch_state.coefficients[IR.args[2]].value[IR.args[3]].name;
        
        code = gen_func_name * "(p, currentT)";
        
    elseif op === :KNOWN_VAR
        # code = "p_data_->valueFEM(fe, dofind)";
        if IR.args[4] == 1
            code = "p_data_->valueFEM(fe, " * string(IR.args[2]-1) * ")";
        else
            code = "p_data_->valueFEM(fe, " * string(IR.args[2]-1) * " + NUM_VARS)";
        end
        
    elseif op === :ROWCOL_TO_INDEX
        code = string(IR.args[2]) * " + (" * string(IR.args[3]) * "-1)*" * string(IR.args[4]);
        
    elseif op === :TIMER
        # A timer has two more args, the label and the content
        code = indent * "double "*string(IR.args[2])*"_timer = MPI_Wtime();\n";
        code *= generate_from_IR_external(IR.args[3], IRtypes, indent=indent);
        code *= "\n" * indent * string(IR.args[2])*"_timer = MPI_Wtime() - " * string(IR.args[2])*"_timer;\n"
        
    elseif op === :FILL_ARRAY
        ## NOTE: Since this uses PETSC for global vectors and matrix
            ## They need special creation/assign steps. For now don't allocate 
            ## or directly assign to them here.
            if typeof(IR.args[2]) <: IR_data_node && (
                IR.args[2].label in [:global_matrix_I, :global_matrix_J, :global_matrix_V, :global_vector, :solution] ||
                (length(IR.args[2].size)>0 && IR.args[2].size[1] in [:allocated_nonzeros, :dofs_global]))
                
                return "";
            end
        # args[2] is the array, args[3] is the value, args[4] is the length
        # for(unsigned int i=1; i<= IR4; i++){
        #    IR2[i] = IR3
        # }
        limit = generate_from_IR_external(IR.args[4], IRtypes);
        code = indent * "for(unsigned int fill_i=1; fill_i<="*limit*"; fill_i++){\n"
        code *= indent * "    " * generate_from_IR_external(IR.args[2], IRtypes) * "[fill_i] = " * string(IR.args[3]) * ";\n";
        code *= indent * "}\n";
        
    elseif op === :INIT_MATRIX_IJ_FV
        code = "";
        
    elseif op === :GLOBAL_FORM_MATRIX
        # // Petsc assembles the system
        code = "";
        
    elseif op === :GLOBAL_GATHER_VECTOR
        # code = indent * "global_vector = gather_system(nothing, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config);"
        code = indent * "// GLOBAL_GATHER_VECTOR named op not ready";
        
    elseif op === :GLOBAL_GATHER_SYSTEM
        # code *= indent * "(global_matrix, global_vector) = gather_system(global_matrix, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config);"
        code = indent * "// GLOBAL_GATHER_SYSTEM named op not ready";
        
    elseif op === :GHOST_EXCHANGE_FV
        # code *= indent * "exchange_ghosts_fv(var, mesh, dofs_per_node, ti);"
        code = indent * "// GHOST_EXCHANGE_FV named op not ready";
        
    elseif op === :GLOBAL_SOLVE
        code = "";
        
    elseif op === :GLOBAL_DISTRIBUTE_VECTOR
        code = "";
        
    elseif op === :GATHER_VARS
        # place variable arrays in global vector
        code = "";
        
    elseif op === :BDRY_TO_VECTOR
        # FV_copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node);
        code = "";
        
    elseif op === :BDRY_TO_VAR
        # copy_bdry_vals_to_variables(var, solution, mesh, dofs_per_node, true)
        code = "";
        
    elseif op === :SCATTER_VARS
        # place solution in variable arrays
        code = "";
        
    elseif op === :LOCAL2GLOBAL
        # put elemental matrix and vector in global system
        code = "";
        
    elseif op === :LOCAL2GLOBAL_VEC
        # put elemental vector in global system
        code = "";
        
    elseif op === :ADD_GLOBAL_VECTOR_AND_NORM
        # u = u + du 
        # and find absolute and relative norms of du
        # args are: u, du, abs_residual, rel_residual
        code = "";
        
    elseif op === :UPDATE_GLOBAL_VECTOR_AND_NORM
        # b = a
        # and find absolute and relative norms of (a-b)
        # args are: a, b, abs_residual, rel_residual
        code = "";
        
    elseif op === :LINALG_VEC_BLOCKS
        # This sets blocks of a vector
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_external(IR.args[4], IRtypes);
        
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
            
            content = generate_from_IR_external(comp, IRtypes);
            code *= indent * "$vecname[$row_index] = $content;\n";
        end
        
        if blocksize > 1
            # end # loop over blocksize
        end
    
    elseif op === :LINALG_MAT_BLOCKS
        # This sets blocks of a matrix
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_external(IR.args[4], IRtypes);
        
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
                row_index = generate_from_IR_external(r_ind, IRtypes);
                col_index = generate_from_IR_external(c_ind, IRtypes);
            else # blocksize == 1
                row_index = generate_from_IR_external(r_ind, IRtypes);
                col_index = generate_from_IR_external(c_ind, IRtypes);
            end
            
            content = generate_from_IR_external(comp, IRtypes);
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
        matname = string(IR.args[4]);
        
        init_lines = "";
        compute_lines = "";
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*3 + 1];
            c_ind = IR.args[4 + (blk-1)*3 + 2];
            comp  = IR.args[4 + (blk-1)*3 + 3];
            row_index = r_ind > 1 ? (string(r_ind-1)*"*nodes_per_element + row") : "row";
            col_index = c_ind > 1 ? (string(c_ind-1)*"*nodes_per_element + col") : "col";
            content = generate_from_IR_external(comp, IRtypes);
            
            init_lines *=    "                $matname[($row_index - 1)*dofs_per_element + $col_index - 1] = 0;\n";
            compute_lines *= "                    $matname[($row_index - 1)*dofs_per_element + $col_index - 1] += $content;\n";
        end
        
        # content = generate_from_IR_external(IR.args[6], IRtypes);
        code = "
        for(unsigned int row=1; row<=nodes_per_element; row++){
            for(unsigned int col=1; col<=nodes_per_element; col++){
$init_lines
                for(unsigned int i=1; i<=qnodes_per_element; i++){
$compute_lines
                }
            }
        }\n";
    
    elseif op === :LINALG_MATVEC_BLOCKS
        # This is not just a matvec, 
        # it is the loop structure for a block matvec that
        # computes blocks that are computed like small 
        # mat*vec loops. b_i = sum_j( something(i,j) )
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = generate_from_IR_external(IR.args[4], IRtypes);
        
        init_lines = "";
        compute_lines = "";
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*2 + 1];
            comp  = IR.args[4 + (blk-1)*2 + 2];
            row_index = r_ind > 1 ? (string(r_ind-1)*"*nodes_per_element + row") : "row";
            content = generate_from_IR_external(comp, IRtypes);
            
            init_lines *=    "            $vecname[$row_index] = 0;\n";
            compute_lines *= "                $vecname[$row_index] += $content;\n";
        end
        
        # content = generate_from_IR_external(IR.args[5], IRtypes);
        code = "
        @inbounds begin
        for row=1:nodes_per_element
$init_lines
            for col=1:qnodes_per_element
$compute_lines
            end
        end
        end";
    
    elseif op === :LINALG_VECTOR_BLOCK
        n_blocks = IR.args[2];
        blocksize = IR.args[3];
        vecname = string(IR.args[4]);
        
        init_lines = "";
        compute_lines = "";
        
        for blk = 1:n_blocks
            r_ind = IR.args[4 + (blk-1)*2 + 1];
            comp  = IR.args[4 + (blk-1)*2 + 2];
            row_index = r_ind > 1 ? (string(r_ind-1)*"*nodes_per_element + row") : "row";
            content = generate_from_IR_external(comp, IRtypes);
            
            init_lines *=    "            $vecname[$row_index] = 0;\n";
            compute_lines *= "                $vecname[$row_index] += $content;\n";
        end
        
        # content = generate_from_IR_external(IR.args[5], IRtypes);
        code = "
        for(unsigned int row=1; row<=nodes_per_element; row++){
$init_lines
            for(unsigned int row=1; row<=qnodes_per_element; row++){
$compute_lines
            }
        }\n";
    
    elseif op === :LINALG_TDM
        # Tcode = generate_from_IR_external(IR.args[2], IRtypes) * "[i + (row-1)*qnodes_per_element]";
        # Mcode = generate_from_IR_external(IR.args[4], IRtypes) * "[i + (col-1)*qnodes_per_element]";
        # Dcode = generate_from_IR_external(IntermediateRepresentation.apply_indexed_access(IR.args[3], [:i], IRtypes), IRtypes);
        # code = "*(" * Tcode *", "* Dcode *", "* Mcode * ")";
        
        i_index = :row;
        j_index = :col;
        k_index = :i;
        
        code = generate_from_IR_external(IntermediateRepresentation.generate_linalg_TDM_product(IR.args[2], IR.args[3], IR.args[4], i_index, j_index, k_index), IRtypes);
            
    elseif op === :LINALG_Tv
        # Tcode = generate_from_IR_external(IR.args[2], IRtypes) * "[col + (row-1)*qnodes_per_element]";
        # vcode = generate_from_IR_external(IntermediateRepresentation.apply_indexed_access(IR.args[3], [:col], IRtypes), IRtypes);
        # code = "*(" * Tcode *", "* vcode * ")";
        
        i_index = :row;
        j_index = :col;
        
        code = generate_from_IR_external(IntermediateRepresentation.generate_linalg_Tv_product(IR.args[2], IR.args[3], i_index, j_index), IRtypes);
    
    else # An unknown named op
        code = indent * "// UNKNOWN NAMED OP: "*string(op);
    end
    
    return code;
end

#########################################################
# code writing utilities
#########################################################

# returns the string for the C++ equivalent type
function cpp_type_name(T::String)
    if T == "Int32"
        return "int32_t";
    elseif T == "Int64" || T == "CustomInt"
        return "int64_t";
    elseif T == "Float32"
        return "float";
    elseif T == "Float64" || T == "CustomFloat"
        return "double";
    elseif T == "Bool"
        return "bool";
    end
    # What else could it be?
    printerr("Unknown type encountered in C++ generation: "*T*", check generated code for UNKNOWNTYPE");
    return "UNKNOWNTYPE";
end

# numbers: 2 -> "2"
# strings: "thing" -> "\"thing\""
# arrays: [1 2; 3 4] -> "{{1,2},{3,4}}"
function cpp_gen_string(v)
    if typeof(v) == String
        return "\""*v*"\"";
        
    elseif typeof(v) <: Number
        return string(v);
        
    elseif typeof(v) <: Array
        if ndims(v) == 1 # "{a, b, c}"
            n = length(v);
            str = "{";
            for i=1:n
                str = str*cpp_gen_string(v[i]);
                if i < n
                    str = str*", ";
                end
            end
            str = str*"}";
        elseif ndims(v) == 2 # "{{a, b}, {c, d}}
            (n,m) = size(v);
            str = "{";
            for i=1:n
                for j=1:m
                    str = str*cpp_gen_string(v[i,j]);
                    if i*j < n*m
                        str = str*", ";
                    end
                end
            end
            str = str*"}";
        end
        return str;
        
    elseif typeof(v) == GenFunction
        return "finch::" * v.name;
    else
        return string(v);
    end
end

#= Generates:
std::function<returntype(argtypes)> name = [captures](args){
    (content)
};
# Note: content should be an array of lines to allow indentation
=#
function cpp_functional(indent, name, args, argtypes, ret, rettype, captures, content)
    lines = fill("", length(content)+2);
    inner_indent = indent*"    ";
    
    arg = string(argtypes[1])*" "*string(args[1]);
    argtype = string(argtypes[1]);
    for i=2:length(args)
        arg = arg*", "*string(argtypes[i])*" "*string(args[i]);
        argtype = argtype*", "*string(argtypes[i]);
    end
    ret = string(ret);
    
    lines[1] = indent*"std::function<"*rettype*"("*argtype*")> "*name*" = ["*captures*"]("*arg*"){";
    for i=1:length(content)
        lines[i+1] = inner_indent*content[i];
    end
    lines[end] = indent*"};";
    
    return lines;
end

#= Generates:
returntype name(args){
    (content)
}
# Note: content should be an array of lines to allow indentation
=#
function cpp_function_def(indent, name, args, argtypes, rettype, content)
    lines = fill("", length(content)+2);
    inner_indent = indent*"    ";
    
    arg = string(argtypes[1])*" "*string(args[1]);
    argtype = string(argtypes[1]);
    for i=2:length(args)
        arg = arg*", "*string(argtypes[i])*" "*string(args[i]);
        argtype = argtype*", "*string(argtypes[i]);
    end
    
    lines[1] = indent*rettype*" "*name*"("*arg*"){";
    for i=1:length(content)
        lines[i+1] = inner_indent*content[i];
    end
    lines[end] = indent*"}";
    
    return lines;
end

#= Generates:
for(iterator=range[1]; iterator<range[2]; iterator+=step){
    (content)
}
# Note: content should be an array of lines to allow indentation
=#
function cpp_for_loop(indent, iterator, range, step, content)
    lines = [];
    inner_indent = indent*"    ";
    push!(lines, indent*"for("*iterator*" = "*string(range[1])*";"*iterator*" < "*string(range[2])*";"*iterator*" += "*string(step)*"){");
    for i=1:length(content)
        push!(lines, inner_indent*content[i]);
    end
    push!(lines, indent*"}")
    
    return lines;
end

# changes symbol "a" to symbol "b" in expression ex
function cpp_swap_symbol(a, b, ex)
    if typeof(ex) == Symbol
        if ex === a
            return b;
        else
            return ex;
        end
    elseif typeof(ex) <: Number
        return ex;
    elseif typeof(ex) == Expr && length(ex.args) > 1
        swapped = copy(ex);
        for i=1:length(ex.args)
            swapped.args[i] = cpp_swap_symbol(a,b,ex.args[i]);
        end
        return swapped;
    else
        return ex;
    end
end

# changes the math operators in the expr
function cpp_change_math_ops(ex)
    if typeof(ex) == Expr
        if ex.head === :.
            # a broadcast operator, change to call
            ex.head = :call
            ex.args = [ex.args[1]; ex.args[2].args];
        end
        for i=1:length(ex.args)
            ex.args[i] = cpp_change_math_ops(ex.args[i]);
        end
    elseif typeof(ex) == Symbol
        if ex === :.+   ex = :+; end
        if ex === :.-   ex = :-; end
        if ex === :.*   ex = :*; end
        if ex === :./   ex = :/; end
        if (ex === :.^ || ex === :^)   ex = :pow; end
        if ex === :abs   ex = :fabs; end
    end
    return ex;
end

# builds a c++ functional for a constant
# When would this ever be useful?
function cpp_number_to_function(name, val)
    indent = "";
    args = ["x"; "y"; "z"; "var"];
    argtypes = ["double"; "double"; "double"; "double*"];
    ret = "";
    rettype = "void";
    captures = "gridX_to_X,gridY_to_Y,gridZ_to_Z";
    content = ["var[0] = "*string(val)*";"];
    return cpp_functional(indent, name, args, argtypes, ret, rettype, captures, content);
end

function cpp_conditional_to_oneline_string(ex)
    if typeof(ex) == Expr
        if ex.head === :if
            # :(a ? b : c)
            # get strings for a, b, c
            a_str = string(ex.args[1]);
            b_str = string(ex.args[2]);
            c_str = string(ex.args[3]); # This could potentially be an error is there is no else, but there should be an else.
            # Now ex will just be a string expression
            ex = "($(a_str) ? $(b_str) : $(c_str))";
        else
            for i=1:length(ex.args)
                ex.args[i] = cpp_conditional_to_oneline_string(ex.args[i]);
            end
        end
    end
    return ex;
end

function cpp_genfunction_to_string(genfun)
    newex = cpp_change_math_ops(genfun.expr); # change operators to match C++
    newex = cpp_swap_symbol(:pi, :M_PI, newex); # swap pi for M_PI
    newex = cpp_swap_symbol(:x, :(pt.x()), newex); # swap x for pt.x()
    newex = cpp_swap_symbol(:y, :(pt.y()), newex); # swap x for pt.y()
    newex = cpp_swap_symbol(:z, :(pt.z()), newex); # swap x for pt.z()
    newex = cpp_conditional_to_oneline_string(newex); # make conditionals like (a ? b : c)
    
    s = string(newex);
    ns = replace(s, r"([\d)])([(A-Za-z])" => s"\1*\2"); # explicitly multiply with "*" (2*x not 2x)
    return ns;
end

# Returns the C++ string corresponding to a bid_def string.
function cpp_bid_def(str::String)
    result = str;
    # Check for CUSTOM
    if occursin("CUSTOM", str) || length(str)<1
        result = "false; // Must manually specify custom boundaries";
        printerr("Custom boundary regions must be specified manually. See XXXBoundaryConditions.h");
    else
        # Build by replacing keywords with expressions
        # Julia 1.7+ has a more elegant way to do this, but for now...
        result = replace(str, "abs(" => "fabs(")
        result = replace(result, "XMIN" => "(fabs(x - domainMin.x(0)) < eps)")
        result = replace(result, "XMAX" => "(fabs(x - domainMax.x(0)) < eps)")
        result = replace(result, "YMIN" => "(fabs(y - domainMin.x(1)) < eps)")
        result = replace(result, "YMAX" => "(fabs(y - domainMax.x(1)) < eps)")
        result = replace(result, "ZMIN" => "(fabs(z - domainMin.x(2)) < eps)")
        result = replace(result, "ZMAX" => "(fabs(z - domainMax.x(2)) < eps)")
        result = replace(result, "eps()" => "eps")
    end
    
    return result;
end

#######################################################
# Write code files

#######  ######  ##        #######   #####
##         ##    ##        ##       ###   #
#######    ##    ##        ######     ### 
##         ##    ##        ##       #   ###
##       ######  ########  #######   #####

# dendrite_main_file(var, config)
# dendrite_equation_file(var, config, IR)
# dendrite_boundary_file(var)
# dendrite_nodedata_file(var, config)
# dendrite_inputdata_file(var)
#######################################################

#=
The main function sets everything up and includes the solve or time stepping code.
=#
function dendrite_main_file(var)
    config = finch_state.config;
    prob = finch_state.prob;
    project_name = finch_state.project_name;
    file = add_generated_file(project_name*".cpp", dir="src");
    
    if !(typeof(var) <:Array)
        var = [var];
    end
    
    # gather important numbers
    dofs_per_node = 0;
    dofs_per_loop = 0;
    dof_names = [];
    varcount = length(var);
    offset_ind = zeros(Int, varcount);
    for i=1:length(var)
        offset_ind[i] = dofs_per_node;
        dofs_per_node += var[i].total_components;
        dofs_per_loop += length(var[i].symvar);
        for j=1:length(var[i].total_components)
            push!(dof_names, string(var[i].symbol)*"_"*string(j));
        end
    end
    
    # Things to delete/destroy at the end
    delete_part = "";
    
    # Make the solve or time stepping part
    solution_step = ""
    if prob.time_dependent
        stepper = finch_state.time_stepper;
        # Make a function for the initial conditions
        # std::function<void(const double *, double *)> initial_condition = [](const double *x, double *var) {
        #     pt[0] = oct_pt[0];
        #     pt[1] = oct_pt[1];
        #     octToPhys.convertCoordsToPhys(nonconst_pt, 1);
        #     var[0] = sin(M_PI * x[0]) * sin(M_PI * x[1]) * sin(M_PI * x[2]);
        #     var[1] = ...
        # };
        solution_step *= "    std::function<void(const double *, double *)> initial_condition = [pt,octToPhys](const double *oct_pt, double *var) {\n";
        solution_step *= "        pt[0] = oct_pt[0];\n";
        solution_step *= "        pt[1] = oct_pt[1];\n";
        if config.dimension > 2
            solution_step *= "        pt[2] = oct_pt[2];\n";
        end
        solution_step *= "        octToPhys.convertCoordsToPhys(pt, 1);\n";
        
        dof_ind = 0;
        for vi=1:varcount
            for ci=1:var[vi].total_components
                ic = prob.initial[var[vi].index][ci];
                if typeof(ic) <: Number
                    solution_step *= "        var[$(dof_ind)] = $(ic);\n";
                else # genfunction
                    newex = cpp_change_math_ops(ic.expr); # change operators to match C++
                    newex = cpp_swap_symbol(:pi, :M_PI, newex); # swap pi for M_PI
                    newex = cpp_swap_symbol(:x, :(pt[0]), newex); # swap x for pt[0]
                    newex = cpp_swap_symbol(:y, :(pt[1]), newex); # swap x for pt[1]
                    newex = cpp_swap_symbol(:z, :(pt[2]), newex); # swap x for pt[2]
                    newex = cpp_conditional_to_oneline_string(newex); # make conditionals like (a ? b : c)
                    
                    s = string(newex);
                    ns = replace(s, r"([\d)])([(A-Za-z])" => s"\1*\2"); # explicitly multiply with "*" (2*x not 2x)
                    solution_step *= "        var[$(dof_ind)] = $(ns);\n";
                end
                dof_ind += 1;
            end
        end
        solution_step *= "    };\n";
        
        # A progress meter
        solution_step *= "    // Progress meter\n";
        solution_step *= "    int progress_step_size = 10;\n";
        solution_step *= "    int last_progress = 0;\n";
        solution_step *= "    TALYFEMLIB::PrintStatus(\"Time step progress (%) 0\");\n";
        progress_update = """
        if((currentStep * 100.0 / nSteps) >= (last_progress + progress_step_size)){
            last_progress += progress_step_size;
            TALYFEMLIB::PrintStatus(last_progress, "%");
        }
"""
        
        # Set up the time stepping loop
        # Prepare the required storage
        time_storage = ["prev_solution"];
        set_initial = [true];
        step_loop = "";
        if stepper.type == EULER_IMPLICIT || stepper.type == CRANK_NICHOLSON
            step_loop = """
    while (currentStep < nSteps) {
        currentT += dt;
        timeInfo.increment();
        $(project_name)Eq->equation()->advanceTime(dt);
        currentStep++;
        $(project_name)Solver->solve();
        VecCopy($(project_name)Solver->getCurrentSolution(), prev_solution);
        
$(progress_update)
    }
"""
        elseif stepper.type == EULER_EXPLICIT
            step_loop = """
    while (currentStep < nSteps) {
        $(project_name)Solver->solve();
        currentT += dt;
        timeInfo.increment();
        $(project_name)Eq->equation()->advanceTime(dt);
        currentStep++;
        VecCopy($(project_name)Solver->getCurrentSolution(), prev_solution);
        
$(progress_update)
    }
"""
        elseif stepper.type == LSRK4
            printerr("TIme stepper "*string(stepper.type)*" is not supported by this target yet.")
            step_loop = "UNSUPPORTED TIME STEPPING METHOD "*string(stepper.type);
            
        elseif stepper.type == RK4
            append!(time_storage, ["tmp_solution", "rk_k_1", "rk_k_2", "rk_k_3"]);
            append!(set_initial, [false, false, false, false]);
            
            step_loop = """
    PetscScalar dt_half = dt/2;
    PetscScalar dt_full = dt;
    PetscScalar* rk_b = {0.16666666666666667, 0.33333333333333333, 0.33333333333333333, 0.16666666666666667};
    Vec rk_vecs[4];
    rk_vecs[0] = rk_k_1;
    rk_vecs[1] = rk_k_2;
    rk_vecs[2] = rk_k_3;
    while (currentStep < nSteps) {
        currentStep++;
        VecCopy(prev_solution, tmp_solution);
        
        $(project_name)Solver->solve();
        VecCopy($(project_name)Solver->getCurrentSolution(), rk_k_1);
        
        currentT += dt/2;
        $(project_name)Eq->equation()->advanceTime(dt/2);
        VecWAXPY(prev_solution, dt_half, rk_k_1, tmp_solution);
        $(project_name)Solver->solve();
        VecCopy($(project_name)Solver->getCurrentSolution(), rk_k_2);
        
        VecWAXPY(prev_solution, dt_half, rk_k_2, tmp_solution);
        $(project_name)Solver->solve();
        VecCopy($(project_name)Solver->getCurrentSolution(), rk_k_3);
        
        currentT += dt/2;
        timeInfo.increment();
        $(project_name)Eq->equation()->advanceTime(dt/2);
        VecWAXPY(prev_solution, dt_full, rk_k_3, tmp_solution);
        $(project_name)Solver->solve();
        
        rk_vecs[3] = $(project_name)Solver->getCurrentSolution();
        VecMAXPY(tmp_solution, 4, rk_b, rk_vecs);
        VecCopy(tmp_solution, prev-solution);
        
$(progress_update)
    }
    
"""
            
        elseif stepper.type == "BDF2"
            push!(time_storage, "prev_prev_solution");
            push!(set_initial, true);
            step_loop = """
    while (currentStep < nSteps) {
        currentT += dt;
        timeInfo.increment();
        $(project_name)Eq->equation()->advanceTime(dt);
        currentStep++;
        $(project_name)Solver->solve();
        VecCopy(prev_solution, prev_prev_solution);
        VecCopy($(project_name)Solver->getCurrentSolution(), prev_solution);
        
$(progress_update)
    }
"""
        else
            printerr("TIme stepper "*string(stepper.type)*" is not supported by this target yet.")
            step_loop = "UNSUPPORTED TIME STEPPING METHOD "*string(stepper.type) * "\n";
        end
        
        # Set up the vectors
        solution_step *= "    // Set up storage used by time stepper.\n";
        nvecs = length(time_storage);
        solution_step *= "    Vec " * time_storage[1];
        delete_part *= "    VecDestroy(&" * time_storage[1] * ");\n";
        for i=2:nvecs
            solution_step *= ", "*time_storage[i];
            delete_part *= "    VecDestroy(&" * time_storage[i] * ");\n";
        end
        solution_step *= ";\n";
        
        for i=1:nvecs
            solution_step *= "    octDA->petscCreateVector("*time_storage[i]*", false, false, NUM_VARS);\n";
        end
        
        # set initial conditions where needed
        solution_step *= "    // Set initial condtions.\n";
        solution_step *= "    octDA->petscSetVectorByFunction(prev_solution, initial_condition, false, false, NUM_VARS);\n";
        for i=2:nvecs
            if set_initial[i]
                solution_step *= "    VecCopy(prev_solution, "*time_storage[i]*");\n";
            end
        end
        
        # Set vectors in Eq
        vecOfVecs = "{VecInfo("*time_storage[1]*", NUM_VARS, 0)";
        for i=2:nvecs
            j = i-1;
            vecOfVecs *= ", VecInfo("*time_storage[i]*", NUM_VARS, $(j))";
        end
        vecOfVecs *= "}";
        solution_step *= "    $(project_name)Eq->setVectors($(vecOfVecs), SYNC_TYPE::VECTOR_ONLY);\n";
        
        # Time stepping loop
        solution_step *= "    \n    timers.Stop(timer_tags[\"Setup\"]);\n\n";
        solution_step *= "    // end setup, start solve /////////////////////////////////////////////////////////////////////////////\n\n";
        solution_step *= "    timers.Start(timer_tags[\"Solve\"]);\n\n"
        solution_step *= "    // Beginning time steps\n";
        solution_step *= step_loop;
        solution_step *= "    \n    timers.Stop(timer_tags[\"Solve\"]);\n\n";
        
    else # not time dependent
        # A single solve step
        solution_step = "    \n    timers.Stop(timer_tags[\"Setup\"]);\n\n";
        solution_step *= "    // end setup, start solve /////////////////////////////////////////////////////////////////////////////\n\n";
        solution_step *= "    timers.Start(timer_tags[\"Solve\"]);\n\n"
        solution_step *= "    // Solve it\n    " * project_name * "Solver->solve();\n"
        solution_step *= "    \n    timers.Stop(timer_tags[\"Solve\"]);\n\n";
    end # solution_step
    
    # Output
    output_part = "    // Output to vtk\n";
    output_part *= "    timers.Start(timer_tags[\"FileIO\"]);\n"
    output_part *= "    static const char *varname[]{\"" * dof_names[1] * "\"";
    for i=2:length(dof_names)
        output_part *= ", \"" * dof_names[i] * "\"";
    end
    output_part *= "};\n";
    if prob.time_dependent
        #output_part *= "    petscVectopvtu(octDA, prev_solution, \"$(project_name)\", varname, physDomain, false, false, $(project_name)NodeData::NUM_VARS);\n";
        #            PETSc::petscVectopvtu(octDA, treePartition, vec, folder, fname, varName, subDomain.domainExtents(), false, false, ndof);
        output_part *= "    util_funcs::save_data(octDA, dTree.getTreePartFiltered(), prev_solution, inputData, timeInfo, subDomain, varname);"  
    else
        # output_part *= "    petscVectopvtu(octDA, $(project_name)Solver->getCurrentSolution(), \"$(project_name)\", varname, physDomain, false, false, $(project_name)NodeData::NUM_VARS);\n";
        output_part *= "    util_funcs::save_data(octDA, dTree.getTreePartFiltered(), $(project_name)Solver->getCurrentSolution(), inputData, timeInfo, subDomain, varname);" 
    end
    output_part *= "    timers.Stop(timer_tags[\"FileIO\"]);\n"
    
    content = """
/**
* Main file for $(project_name).
*/ 
// General
#include <iostream>
#include <point.h>
// Dendrite/Taly
#include <DendriteUtils.h>
#include <TalyEquation.h>
#include <Traversal/Analytic.h>
#include <IO/VTU.h>
#include <IMGA/IMGA.h>
#include <IMGA/Marker.h>
// Petsc
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <PETSc/IO/petscVTU.h>

// generated
#include <$(project_name)Equation.h>
#include <$(project_name)BoundaryConditions.h>
#include <$(project_name)NodeData.h>
#include <$(project_name)InputData.h>
#include "util.h"

using namespace PETSc;

int main(int argc, char* argv[]) {
    /// initialize
    dendrite_init(argc, argv);
    int rank = TALYFEMLIB::GetMPIRank();
    
    // timers //////////////////////////////////////////////////////////////////////////////
    TimerGroup<MPITimer> timers;

    std::vector<std::string>
        timer_labels = {"Total", "Setup", "Solve", "FileIO"};
    std::map<std::string, int> timer_tags;
    for (int i = 0; i < timer_labels.size(); i++){
        timer_tags.insert(std::pair<std::string, int>(timer_labels[i], i));
        timers.AddTimer(timer_labels[i]);
    }
    timers.Start(timer_tags["Total"]);
    timers.Start(timer_tags["Setup"]);
    //////////////////////////////////////////////////////////////////////////////////////////
    
    /// read parameters from config.txt //////////////////////////////////////////////////////
    $(project_name)InputData inputData;
    std::ifstream configFile("config.txt");
    if (configFile.good()) {
        if (!inputData.ReadFromFile()) {
            throw std::runtime_error("[ERR] Error reading input data, check the config file!");
        }
        if (!inputData.CheckInputData()) {
            throw std::runtime_error("[ERR] Problem with input data, check the config file!");
        }
    } else{
        TALYFEMLIB::PrintStatus("No config.txt file found. Copy the generated config.txt into the build directory");
        return -1;
    }
    
    // Print config
    inputData.PrintInputData();
    //////////////////////////////////////////////////////////////////////////////////////////
    
    // Configuration variables
    unsigned int eleOrder = inputData.basisFunction;
    unsigned int level = inputData.meshDef.refineLevel_base;
    bool mfree = inputData.ifMatrixFree;
    
    // build the mesh ////////////////////////////////////////////////////////////////////////
    bool resume_from_checkpoint = false;
    DomainExtents domainExtents(inputData.meshDef.fullDADomain, inputData.meshDef.physDomain);
    DA *octDA = nullptr;
    DistTREE dTree;
    SubDomain subDomain(domainExtents, resume_from_checkpoint);
    
    IMGA *imga = new IMGA(domainExtents, IBM_METHOD::SBM);

    /// Load stl
    std::vector<GEOMETRY::Geometry *> curve_geoms;
    std::vector<GEOMETRY::STL *> stls; // 3D
    std::vector<GEOMETRY::MSH *> mshs; // 2D
    std::array<DENDRITE_REAL, DIM> shift;
    std::vector<GeomRefinement> ibm_refinements;

    for (const auto &geom_def : inputData.carved_out_geoms_def){
#if (DIM == 2)
        if (geom_def.type == CarvedOutGeom::Type::MESHOBJECT_2D){
            mshs.push_back(new GEOMETRY::MSH(geom_def.mesh_path, GEOMETRY::InOutTest2D::RAY_TRACING_2D));
            shift[0] = geom_def.InitialDisplacement[0];
            shift[1] = geom_def.InitialDisplacement[1];
            auto geom_retain_side = RetainSide::OUT;
            if (geom_def.outer_boundary){
                geom_retain_side = RetainSide::IN;
            }
            curve_geoms.push_back(new GEOMETRY::Geometry(mshs.back(), Point<DIM>(shift), geom_retain_side));
            ibm_refinements.emplace_back(geom_def.geomRefine);
        }
#endif
#if (DIM == 3)
        if (geom_def.type == CarvedOutGeom::Type::MESHOBJECT){
            stls.push_back(new GEOMETRY::STL(geom_def.mesh_path, GEOMETRY::InOutTest::RAY_TRACING));
            shift[0] = geom_def.InitialDisplacement[0];
            shift[1] = geom_def.InitialDisplacement[1];
            shift[2] = geom_def.InitialDisplacement[2];
            auto geom_retain_side = RetainSide::OUT;
            if (geom_def.outer_boundary){
                geom_retain_side = RetainSide::IN;
            }
            curve_geoms.push_back(new GEOMETRY::Geometry(stls.back(), Point<DIM>(shift), geom_retain_side));
            ibm_refinements.emplace_back(geom_def.geomRefine);
        }
#endif
    }
    
    for (const auto &c : curve_geoms){
        subDomain.addObject(c);
    }
    for (int i = 0; i < curve_geoms.size(); i++){
        imga->addGeometry(curve_geoms.at(i), ibm_refinements.at(i));
    }
    
    /// carving subda
    std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *octCoords, double scale){
        return (subDomain.functionToRetain(octCoords, scale));
    };
    
    octDA = createSubDA(dTree, functionToRetain, level, eleOrder);
    subDomain.finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
    util_funcs::performRefinementSubDA(octDA, dTree.getTreePartFiltered(), domainExtents, dTree, inputData, &subDomain);
    
    TALYFEMLIB::PrintStatus("total No of nodes in the mesh = ", octDA->getGlobalNodeSz());
    
    subDomain.finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
    IO::writeBoundaryElements(octDA, dTree.getTreePartFiltered(), "boundary", "subDA", domainExtents);
    SubDomainBoundary boundary(&subDomain, octDA, domainExtents);
    
    //////////////////////////////////////////////////////////////////////////////////////////
    
    // Problem specification
    static const unsigned int NUM_VARS = $(dofs_per_node);
    DomainInfo physDomain;
    physDomain.min = inputData.meshDef.min;
    physDomain.max = inputData.meshDef.max;
    
    // Time steps (define even if no time stepping)
    unsigned int currentStep = 0;
    double currentT = 0.0;
    double dt = inputData.dt[0];
    unsigned int nSteps = inputData.nSteps[0];
    double totalT = nSteps * dt;
    TimeInfo timeInfo(0.0, inputData.dt, inputData.totalT);
    
    /// Boundary condition
    $(project_name)BoundaryConditions $(project_name)BC(&boundary, &inputData, &timeInfo);

    /// Gridfield setup
    TalyMesh<$(project_name)NodeData> talyMesh(octDA->getElementOrder());
    SubDomainBoundary *subDomainBoundary = &boundary;
    
    /// imga setup is not currently done. When is this needed?
    // imga->initIMGAComputation(octDA, dTree.getTreePartFiltered());
    
    /// Equation and solver setup
    Marker *elementMarker = new Marker(octDA, dTree.getTreePartFiltered(), domainExtents, imga, MarkerType::GAUSS_POINT);
    std::vector<PetscInt> dirichletNodes;
    
    auto $(project_name)Eq = new TalyEquation<$(project_name)Equation, $(project_name)NodeData>
                                                (&talyMesh, octDA, dTree.getTreePartFiltered(),
                                                subDomain.domainExtents(), NUM_VARS, &timeInfo, true, subDomainBoundary,
                                                &inputData, imga);
    
    $(project_name)Eq->equation()->setBoundaryCondition(&$(project_name)BC);
    
    // Setup PETSc solver
    LinearSolver *$(project_name)Solver = setLinearSolver($(project_name)Eq, octDA, NUM_VARS, mfree);
    /// apply solver parameters from config.txt
    inputData.solverOptions.apply_to_petsc_options("-$(project_name)_");
    KSP m_ksp = $(project_name)Solver->ksp();
    KSPSetOptionsPrefix(m_ksp, "$(project_name)_");
    KSPSetFromOptions(m_ksp);
    
    // Boundary conditions
    $(project_name)Solver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary {
        Boundary b;
        $(project_name)BC.setBoundaryConditions(b, pos);
        return b;
    });
    
    $(project_name)Solver->setIBMDirichletNodes(dirichletNodes);
    $(project_name)Eq->assignIBMConstructs(imga, elementMarker->getMarkers().data());
    
    // initial condition function needs physical coordinates
    OctToPhysical octToPhys(domainExtents);
    double *pt = new double[DIM];
    
$(solution_step)

$(output_part)
    
    timers.Stop(timer_tags["Total"]);
    timers.PrintTotalTimeSeconds();
    
$(delete_part)

    delete $(project_name)Eq;
    delete $(project_name)Solver;
    dendrite_finalize(octDA);
    return 0;
}   
"""
    println(file, content);
end

#=
This file contains the elemental matrix and vector computations.
=#
function dendrite_equation_file(var, IR)
    config = finch_state.config;
    project_name = finch_state.project_name;
    file = add_generated_file(project_name*"Equation.h", dir="include");
    
    dofs_per_node = 0;
    for i=1:length(var)
        dofs_per_node += var[i].total_components;
    end
    
    # Generate the elemental matrix and vector calculation
    matrix_part = "// Matrix computation not found in IR";
    vector_part = "// Vector computation not found in IR";
    matrix_coef_part = "// No coefficients to compute";
    vector_coef_part = "// No coefficients to compute";
    # Find the right blocks in the IR
    if typeof(IR) == IR_block_node
        for i=1:length(IR.parts)
            if typeof(IR.parts[i]) == IR_block_node
                if IR.parts[i].label == "elemental matrix"
                    matrix_part = generate_from_IR_external(IR.parts[i], indent="                ");
                elseif IR.parts[i].label == "elemental vector"
                    vector_part = generate_from_IR_external(IR.parts[i], indent="            ");
                elseif IR.parts[i].label == "prepare matrix"
                    matrix_coef_part = generate_from_IR_external(IR.parts[i], indent="        ");
                elseif IR.parts[i].label == "prepare vector"
                    vector_coef_part = generate_from_IR_external(IR.parts[i], indent="        ");
                end
            end
        end
    end
    
    # Any generated functions used for coefficients
    genfunctions = finch_state.genfunctions;
    function_defs = "";
    indent = "    ";
    args ="(const TALYFEMLIB::ZEROPTV &pt, const double t)";
    for i = 1:length(genfunctions)
        str = cpp_genfunction_to_string(genfunctions[i]);
        fun = "    return "*str*";";
        
        function_defs *= indent * "double " * genfunctions[i].name * args * "{\n";
        function_defs *= indent * fun * "\n"*indent*"}\n";
    end
    
    content = """
#pragma once

#include <talyfem/fem/cequation.h>
#include "$(project_name)NodeData.h"
#include "$(project_name)InputData.h"
#include <$(project_name)BoundaryConditions.h>
#include "util.h"
#include "SBMcalc.h"
#include <DataTypes.h>
#include <Basis/MatVec.h>
#include <Basis/Vec.h>
#include <Basis/Mat.h>
class $(project_name)Equation : public TALYFEMLIB::CEquation<$(project_name)NodeData> {
    
    public:
    explicit $(project_name)Equation(const $(project_name)InputData * idata, const IMGA *imga)
            : TALYFEMLIB::CEquation<$(project_name)NodeData>(false, TALYFEMLIB::kAssembleGaussPoints), imga_(imga){
        idata_ = idata;
        dt = idata->dt[0];
        currentT = 0.0;
        NUM_VARS = $(dofs_per_node);
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // These do nothing or redirect to other functions. Needed by talyfem? ////////////////////////
    void Solve(double dt, double t) override {
        assert(false);
    }
    void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                    TALYFEMLIB::ZEROARRAY<double> &be) override {
        assert(false);
    }
    void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, const double *h){
        Integrands_Ae(fe, Ae);
    }
    void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const double *h){
        Integrands_be(fe, be);
    }
    void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZeroMatrix<double> &Ae,
                          const double *h){
        Integrands4side_Ae(fe, side_idx, id, Ae);
    }
    void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZEROARRAY<double> &be,
                          const double *h){
        Integrands4side_be(fe, side_idx, id, be);
    }
    void ibm_Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROMATRIX<double> &Ae,
                                const NodeAndValues<DENDRITE_REAL> &gpinfo,
                                const TALYFEMLIB::ZEROPTV &position,
                                const TALYFEMLIB::ZEROPTV &h){
        return;
    }
    void ibm_Integrands4side_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be,
                                const NodeAndValues<DENDRITE_REAL> &gpinfo,
                                const TALYFEMLIB::ZEROPTV &position,
                                const TALYFEMLIB::ZEROPTV &h){
        return;
    }
    void ibm_Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROMATRIX<double> &Ae,
                                const NodeAndValues<DENDRITE_REAL> &gpinfo,
                                const TALYFEMLIB::ZEROPTV &position,
                                const TALYFEMLIB::ZEROPTV &h,
                                const std::vector<double> &surface_values){
        ibm_Integrands4side_Ae(fe, Ae, gpinfo, position, h);
    }
    void ibm_Integrands4side_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &Ae,
                                const NodeAndValues<DENDRITE_REAL> &gpinfo,
                                const TALYFEMLIB::ZEROPTV &position,
                                const TALYFEMLIB::ZEROPTV &h,
                                const std::vector<double> &surface_values){
        ibm_Integrands4side_be(fe, Ae, gpinfo, position, h);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    // Volume integrals //////////////////////////////////////////////////////////////////////
    /*
    This is called once for each quadrature point (1D???)
    Computes every element of Ae for one quadrature point.
    Will be inside this loop:
    while (fe.next_itg_pt()) {
        Integrands_Ae(fe, Ae);
    }
    
    Why do it this way?
    */
    void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {
        using namespace TALYFEMLIB;
        // # of dimensions: 1, 2, or 3
        // const int n_dimensions = fe.nsd(); // <- config.dimension
        const int n_dimensions = """*string(config.dimension)*""";
        // # of basis functions
        const int n_basis_functions = fe.nbf(); // <- refel.Np ?
        
        // (determinant of J) cross W
        const double wdetj = fe.detJxW();
        // coordinates of this point
        const ZEROPTV p = fe.position();
        
"""*matrix_coef_part*"""
        
        for (int row = 0; row < n_basis_functions; row++) {
            for (int col = 0; col < n_basis_functions; col++) {
                ////////////////////////////////////////////////////////////////////////
                // This is generated from the input expressions
"""*matrix_part*"""
                ////////////////////////////////////////////////////////////////////////
            }
        }
    }
    
    void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
        using namespace TALYFEMLIB;
        // # of basis functions
        const int n_basis_functions = fe.nbf();
        // (determinant of J) cross W
        const double wdetj = fe.detJxW();
        // coordinates of this point
        const ZEROPTV p = fe.position();
        
"""*vector_coef_part*"""
        
        for (int row = 0; row < n_basis_functions; row++) {
            ////////////////////////////////////////////////////////////////////////
            // This is generated from the input expressions
"""*vector_part*"""
            ////////////////////////////////////////////////////////////////////////
        }
        
    }
    
    // Surface integrals //////////////////////////////////////////////////////////////////////
    void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZeroMatrix<double> &Ae){
        
        GPpos_.push_back(fe.position());
        
        /// from integral by parts, there is no Ae for the boundary term,
        /// because -[w\\dot \\alpha \\grad T] are all known variable, and they go into the be
        /// for the weak BC
        const int nsd = DIM;
        const double Cb_e = idata_->Cb_e;
        const double detSideJxW = fe.detJxW();

        double h = util_funcs::ElementSize(fe);
        double alpha = Cb_e;
        double d[DIM];

        SBMCalc sbmCalc(fe, idata_, imga_);
        sbmCalc.Dist2Geo(d);
        
        DENDRITE_REAL secondOrderTerm_a(0), secondOrderTerm_b(0);

        for (int a = 0; a < fe.nbf(); a++){
            DENDRITE_REAL gradWdotn = 0;
            DENDRITE_REAL gradWdotd = 0;

#if (DIM == 2)
            if (idata_->elemOrder == 2 and idata_->SecondOrderTaylorQBF){
                secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1]) +
                                     d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1])) / 2;
            }else{
                secondOrderTerm_a = 0;
            }
#endif
#if (DIM == 3)
            if (idata_->elemOrder == 2 and idata_->SecondOrderTaylorQBF){
                secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1] + fe.d2N(a, 0, 2) * d[2]) + 
                                     d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1] + fe.d2N(a, 1, 2) * d[2]) + 
                                     d[2] * (fe.d2N(a, 2, 0) * d[0] + fe.d2N(a, 2, 1) * d[1] + fe.d2N(a, 2, 2) * d[2])) / 2;
            }else{
                secondOrderTerm_a = 0;
            }
#endif

            for (int k = 0; k < DIM; k++){
                gradWdotn += fe.dN(a, k) * (fe.surface()->normal().data()[k]);
                gradWdotd += fe.dN(a, k) * d[k];
            }

            for (int b = 0; b < fe.nbf(); b++){
#if (DIM == 2)
                if (idata_->elemOrder == 2 and idata_->SecondOrderTaylorQBF){
                    secondOrderTerm_b = (d[0] * (fe.d2N(b, 0, 0) * d[0] + fe.d2N(b, 0, 1) * d[1]) +
                                        d[1] * (fe.d2N(b, 1, 0) * d[0] + fe.d2N(b, 1, 1) * d[1])) / 2;
                }else{
                    secondOrderTerm_b = 0;
                }
#endif
#if (DIM == 3)
                if (idata_->elemOrder == 2 and idata_->SecondOrderTaylorQBF){
                    secondOrderTerm_b = (d[0] * (fe.d2N(b, 0, 0) * d[0] + fe.d2N(b, 0, 1) * d[1] + fe.d2N(b, 0, 2) * d[2]) + 
                                        d[1] * (fe.d2N(b, 1, 0) * d[0] + fe.d2N(b, 1, 1) * d[1] + fe.d2N(b, 1, 2) * d[2]) + 
                                        d[2] * (fe.d2N(b, 2, 0) * d[0] + fe.d2N(b, 2, 1) * d[1] + fe.d2N(b, 2, 2) * d[2])) / 2;
                }else{
                    secondOrderTerm_b = 0;
                }
#endif

                DENDRITE_REAL gradUdotd = 0;
                DENDRITE_REAL gradUdotn = 0;

                for (int k = 0; k < DIM; k++){
                    gradUdotn += fe.dN(b, k) * (fe.surface()->normal().data()[k]);
                    gradUdotd += fe.dN(b, k) * d[k];
                }

                Ae(a, b) += -fe.N(a) * gradUdotn * fe.detJxW();

                if (idata_->IfAdjointConsistency){
                    Ae(a, b) += -gradWdotn * (fe.N(b) + gradUdotd + secondOrderTerm_b) * fe.detJxW();
                }

                Ae(a, b) += alpha / h * (fe.N(a) + gradWdotd + secondOrderTerm_a) *
                            (fe.N(b) + gradUdotd + secondOrderTerm_b) * fe.detJxW();
            }
        }
    }
    
    void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZEROARRAY<double> &be){
        
        using namespace TALYFEMLIB;
        const int nsd = DIM;
        const double detSideJxW = fe.detJxW();
        const double Cb_e = idata_->Cb_e;
        
        double h = util_funcs::ElementSize(fe);

        double alpha = Cb_e;
        DENDRITE_REAL d[DIM];
        
        // Need to find the Dirichlet boundary value.
        DENDRITE_REAL boundary_value;
        SBMCalc sbmCalc(fe, idata_, imga_);
        sbmCalc.Dist2Geo(d);
        // The position of the GP
        const ZEROPTV pt = fe.position();
        // The corresponding position on the true boundary
        double x_true = pt.x() + d[0];
        double y_true = pt.y() + d[1];
        boundaryConditions->getBoundaryValue(ZEROPTV(x_true,y_true), &boundary_value);
        
        DENDRITE_REAL secondOrderTerm_a(0);

        for (int a = 0; a < fe.nbf(); a++){
            DENDRITE_REAL gradWdotn = 0;
            DENDRITE_REAL gradWdotd = 0;

#if (DIM == 2)
            if (idata_->elemOrder == 2 and idata_->SecondOrderTaylorQBF){
                secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1]) +
                                    d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1])) / 2;
            }else{
                secondOrderTerm_a = 0;
            }
#endif
#if (DIM == 3)
            if (idata_->elemOrder == 2 and idata_->SecondOrderTaylorQBF){
                secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1] + fe.d2N(a, 0, 2) * d[2]) + 
                                    d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1] + fe.d2N(a, 1, 2) * d[2]) + 
                                    d[2] * (fe.d2N(a, 2, 0) * d[0] + fe.d2N(a, 2, 1) * d[1] + fe.d2N(a, 2, 2) * d[2])) / 2;
            }else{
                secondOrderTerm_a = 0;
            }
#endif

            for (int k = 0; k < DIM; k++){
                gradWdotn += fe.dN(a, k) * (fe.surface()->normal().data()[k]);
                gradWdotd += fe.dN(a, k) * d[k];
            }

            if (idata_->IfAdjointConsistency){
                be(a) += -gradWdotn * boundary_value * fe.detJxW();
            }

            be(a) += alpha / h * (fe.N(a) + gradWdotd + secondOrderTerm_a) * boundary_value * fe.detJxW();
        }
    }
    
    void advanceTime(const double dti){
        currentT += dti;
    }
    
    void GetPos(std::vector<ZEROPTV> &GPpos){
        GPpos = GPpos_;
    }
    
    void setBoundaryCondition($(project_name)BoundaryConditions *bcs){
      boundaryConditions = bcs;
    }
    
    protected:
    const $(project_name)InputData *idata_;
    $(project_name)BoundaryConditions *boundaryConditions;
    double dt;
    double currentT;
    unsigned int NUM_VARS;
    
    const IMGA *imga_;
    std::vector<ZEROPTV> GPpos_;
    
    ////////////////////////////////////////////////////////////////////////
    // Coefficient functions
"""*function_defs*"""
    ////////////////////////////////////////////////////////////////////////
};
    
"""
    println(file, content);
end

#=
Sets boundary conditions.
=#
function dendrite_boundary_file(var)
    config = finch_state.config;
    prob = finch_state.prob;
    project_name = finch_state.project_name;
    file = add_generated_file(project_name*"BoundaryConditions.h", dir="include");
    
    dofs_per_node = 0;
    nvars = length(var);
    for i=1:nvars
        dofs_per_node += var[i].total_components;
    end
    
    if config.dimension == 1
        extract_coords = """
        double x = position.x();
        Point<1> domainMin(inputData_->meshDef.physDomain.min);
        Point<1> domainMax(inputData_->meshDef.physDomain.max);
"""
    elseif config.dimension == 2
        extract_coords = """
        double x = position.x();
        double y = position.y();
        Point<2> domainMin(inputData_->meshDef.physDomain.min);
        Point<2> domainMax(inputData_->meshDef.physDomain.max);
"""
    elseif config.dimension == 3
        extract_coords = """
        double x = position.x();
        double y = position.y();
        double z = position.z();
        Point<3> domainMin(inputData_->meshDef.physDomain.min);
        Point<3> domainMax(inputData_->meshDef.physDomain.max);
"""
    end
    
    # bool on_BID_1 = (fabs(y - domainMax.x(1)) < eps);
    bid_defs = "";
    bid_array = "{";
    nbids = length(prob.bid);
    for i=1:nbids
        bid_name = "on_bid_"*string(prob.bid[i]);
        bid_array *= bid_name;
        if i < nbids
            bid_array *= ", ";
        end
        bid_location = cpp_bid_def(prob.bid_def[i]);
        bid_defs *= "    bool " * bid_name * " = " * bid_location * ";\n";
    end
    bid_array *= "}";
    
    # genfunctions
    bdry_genfunctions = [];
    
    # if on_BID_1 {
    #     b.addDirichlet($(project_name)NodeData::U_1, 0.0);
    # } else if on_BID_2 {
    #     b.addDirichlet($(project_name)NodeData::U_1, 0.0);
    # }
    dirichlet_bc = "";
    dirichlet_bv = "";
    dirichlet_bc_for_bids = fill("", nbids);
    for i=1:nbids
        bid_name = "on_bid_"*string(prob.bid[i]);
        # One line for each var with a dirichlet bdry on this bid
        add_dirichlet = "";
        get_dirichlet = "";
        for vi = 1:length(var)
            var_ind = var[vi].index;
            # One line for each component
            for ci=1:var[vi].total_components
                if prob.bc_type[var_ind, i] == DIRICHLET
                    # Constant or genfunction values
                    if typeof(prob.bc_func[var_ind, i][ci]) <: Number
                        bc_val = string(prob.bc_func[var_ind, i][ci]);
                    elseif typeof(prob.bc_func[var_ind, i][ci]) == GenFunction
                        bc_val = prob.bc_func[var_ind, i][ci].name * "(position, t)";
                        push!(bdry_genfunctions, prob.bc_func[var_ind, i][ci]);
                    else # uh oh
                        printerr("Boundary values that are not constant or [x,y,z,t] gen functions must be entered manually.");
                        bc_val = "0.0";
                    end
                    
                    dof_name = string(var[vi].symbol) * "_" * string(ci) * "_dofind"; # u_1_dofind
                    add_dirichlet *= "b.addDirichlet($(project_name)NodeData::$(dof_name), $(bc_val));\n";
                    get_dirichlet *= "value[$(project_name)NodeData::$(dof_name)] = $(bc_val);\n";
                end
            end
        end
        
        # Put them in their repective bid
        if i==1
            dirichlet_bc *= "    if(" * bid_name * "){\n        " * add_dirichlet * "    }";
            dirichlet_bv *= "    if(" * bid_name * "){\n        " * get_dirichlet * "    }";
        else
            dirichlet_bc *= "    else if(" * bid_name * "){\n        " * add_dirichlet * "    }";
            dirichlet_bv *= "    else if(" * bid_name * "){\n        " * get_dirichlet * "    }";
        end
        dirichlet_bc_for_bids[i] = "                            " * add_dirichlet;
    end
    
    # if (FEQUALS(x, 0.0) and (FEQUALS(y, 0.0)) and (FEQUALS(z, 0.0))) {
    #     b.addDirichlet($(project_name)NodeData::U_1, 0.0);
    # }
    # TODO
    reference_points = "";
    
    # genfunctions to define here
    function_defs = "";
    function_decs = "";
    args ="(const TALYFEMLIB::ZEROPTV &pt, const double t)";
    finished_genfunctions = [];
    for i = 1:length(bdry_genfunctions)
        # avoid duplicates
        new_gf = true;
        for j=1:length(finished_genfunctions)
            if bdry_genfunctions[i].name == finished_genfunctions[j]
                new_gf = false;
                break;
            end
        end
        if new_gf
            push!(finished_genfunctions, bdry_genfunctions[i].name);
            str = cpp_genfunction_to_string(bdry_genfunctions[i]);
            fun = "    return "*str*";";
            
            function_defs *= "double " * bdry_genfunctions[i].name * args * "{\n";
            function_defs *= "    " * fun * "\n}\n";
            function_decs *= "        double " * bdry_genfunctions[i].name * args * ";\n";
        end
    end
    
    content = """
#pragma once

#include <TimeInfo.h>
#include <PETSc/BoundaryConditions.h>
#include <DataTypes.h>
#include "$(project_name)InputData.h"
#include "$(project_name)NodeData.h"

class $(project_name)BoundaryConditions {
public:
    $(project_name)BoundaryConditions(SubDomainBoundary *boundary, $(project_name)InputData *inputData, TimeInfo *time)
        : inputData_(inputData), time_(time), boundaries_(boundary){
        PrintInfo("Setting up boundary conditions");
    }
    
    /**
    * Method to setup Boundary conditions.
    * @param b Boundary object
    * @param position position
    */
    void setBoundaryConditions(PETSc::Boundary &b, const ZEROPTV &position){
        static const double eps = 1e-14;
        double t = time_->getCurrentTime();
$(extract_coords)
        
        // Define BIDs
$(bid_defs)
        
        bool on_bid_array[$(nbids)] = $(bid_array);
        
        /// for carved out geometry
        DENDRITE_UINT objectID = -1;
        boundaries_->generateBoundaryFlags(position, objectID);
        if (boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::SPHERE) or
            boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::CIRCLE) or
            boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::BOX) or
            boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY)){
            
            const auto &carved_geo = inputData_->carved_out_geoms_def.at(objectID);
            // Apply BC for each relevant BID and DOF
            unsigned int nbids = carved_geo.bid_V.size();
            for(int bidi = 0; bidi < nbids; bidi++){
                if(on_bid_array[carved_geo.bid_V[bidi]]){
                    // This position is part of a BID assigned to this geometry.
                    // If it is also a part of the global boundary or another geometry,
                    // Those will be ignored. Only one BC can be set for a given dof here.
                    // The first one found will be used.
$(dirichlet_bc)
                    
                    return;
                }
            }
        }
        // If it gets here, this position is not on a geometry boundary
        
        // Consider full domain boundaries if there is no geometry
$(dirichlet_bc)

$(reference_points)
        
        // This may be used at some point in the future.
        /*
        for (auto &r : inputData_->region_refine){
            if (r.forRetain){
                auto res = r.out_retain(position);
                if (res == RegionalRefine::RetainType::OUTSIDE){
                    b.addDirichlet(0, 0.0);
                    return;
                }else if (res == RegionalRefine::RetainType::ON_BOUNDARY){
                    // ??
                }
            }
        }
        */
    };
    
    /**
    * Returns a dirichlet boundary value.
    * @param value value
    * @param position position
    */
    void getBoundaryValue(const ZEROPTV &position, double* value){
        static const double eps = 1e-14;
        double t = time_->getCurrentTime();
$(extract_coords)
        
        // Define BIDs
$(bid_defs)
        
        bool on_bid_array[$(nbids)] = $(bid_array);
        
        /// for carved out geometry
        DENDRITE_UINT objectID = -1;
        boundaries_->generateBoundaryFlags(position, objectID);
        if (boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::SPHERE) or
            boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::CIRCLE) or
            boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::BOX) or
            boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY)){
            
            const auto &carved_geo = inputData_->carved_out_geoms_def.at(objectID);
            // Apply BC for each relevant BID and DOF
            unsigned int nbids = carved_geo.bid_V.size();
            for(int bidi = 0; bidi < nbids; bidi++){
                if(on_bid_array[carved_geo.bid_V[bidi]]){
                    // This position is part of a BID assigned to this geometry.
                    // If it is also a part of the global boundary or another geometry,
                    // Those will be ignored. Only one BC can be set for a given dof here.
                    // The first one found will be used.
$(dirichlet_bv)
                    
                    return;
                }
            }
        }
        // If it gets here, this position is not on a geometry boundary
        
        // Consider full domain boundaries if there is no geometry
$(dirichlet_bv)

    };
    
private:
    const $(project_name)InputData *inputData_;
    TimeInfo *time_;
    SubDomainBoundary *boundaries_;
    
    // boundary value functions /////////////////////////////////////////////////////
$(function_defs)

    /////////////////////////////////////////////////////////////////////////////////
};
"""
    println(file, content);
end

#=
The NodeData truct has enums for dofs.
=#
function dendrite_nodedata_file(var)
    config = finch_state.config;
    project_name = finch_state.project_name;
    file = add_generated_file(project_name*"NodeData.h", dir="include");
    
    # gather important numbers
    dofs_per_node = 0;
    dofs_per_loop = 0;
    dof_names = [];
    varcount = length(var);
    offset_ind = zeros(Int, varcount);
    for i=1:length(var)
        offset_ind[i] = dofs_per_node;
        dofs_per_node += var[i].total_components;
        dofs_per_loop += length(var[i].symvar);
        for j=1:length(var[i].total_components)
            push!(dof_names, string(var[i].symbol)*"_"*string(j));
        end
    end
    
    dof_enum = "";
    value_case = "";
    name_case = "";
    value_count = 0;
    for i=1:dofs_per_node
        dof_enum *= "        " * dof_names[i] * "_dofind = "*string(i-1)*",\n";
        value_case *= "            case " * dof_names[i] * "_dofind: return var_values[" * string(i-1) * "];\n";
        name_case *= "            case " * dof_names[i] * "_dofind: return \"" * dof_names[i] * "\";\n";
        value_count += 1;
    end
    # Some time steppers need extra values available to equation, so more values here.
    stepper_type = finch_state.time_stepper.type;
    if stepper_type == BDF2
        # one more for prev2
        for i=1:dofs_per_node
            j = dofs_per_node + i;
            dof_enum *= "        prev2_" * dof_names[i] * "_dofind = "*string(j-1)*",\n";
            value_case *= "            case prev2_" * dof_names[i] * "_dofind: return var_values[" * string(j-1) * "];\n";
            name_case *= "            case prev2_" * dof_names[i] * "_dofind: return \"prev2_" * dof_names[i] * "\";\n";
            value_count += 1;
        end
    end
    dof_enum *= "        " * project_name * "NODEDATA_MAX = "*string(value_count)*"\n";
    
    content = """
#pragma once

#include <exception>
#include <assert.h>
#include <DataTypes.h>
#include <talyfem/talyfem.h>

class $(project_name)NodeData {
    public:
    // number of  variables in this NodeData
    static const unsigned int NUM_VARS = $(dofs_per_node);
    // values available to equation
    double var_values[$(value_count)];
    
    $(project_name)NodeData() {
        std::memset(var_values, 0, sizeof(DENDRITE_REAL) * $(value_count));
    }
    
    enum Vars : int {
$(dof_enum)
    };

    /**
        * Returns reference to the given value in the object
        *
        * @param index the index of the desired item
        * @return reference to the desired data item
        */
    double &value(int index) {
        switch (index) {
$(value_case)
            default: throw std::runtime_error("Invalid $(project_name)NodeData index");
        }
    }

    inline double value(int index) const {
        return const_cast<$(project_name)NodeData *>(this)->value(index);
    }

    /**
        * Returns the name of the given data value in the object
        * @param index the index of the desired item
        * @return name of the specified data item
        */
    static const char *name(int index) {
        switch (index) {
$(name_case)
            default: throw std::runtime_error("Invalid $(project_name)NodeData index");
        }
    }

    /**
        * Returns the number of the data items in the object
        * @return number of the data items in the object
        */
    static int valueno() {
        return $(project_name)NODEDATA_MAX; // Should this be num_dofs instead???
    }
};
  
"""
    println(file, content);
end

#=
Refinement.
=#
function dendrite_refine_file(var)
    config = finch_state.config;
    project_name = finch_state.project_name;
    file = add_generated_file(project_name*"Refine.h", dir="include");
    
    # This will eventually contain code for driving refinement
    refinement_criteria = "//TODO";
    
    content = """
#pragma once

#include <Traversal/Refinement.h>
#include <Boundary/SubDomainBoundary.h>

#include <utility>
#include "$(project_name)InputData.h"

class $(project_name)Refine : public Refinement{
    $(project_name)InputData *inputData_;
    const DomainExtents &domainExtents_;
    SubDomainBoundary *subDomainBoundary_;
    bool doCoarsen_ = false;
    std::vector<GEOMETRY::Geometry *> ibmGeoms_;

public:
    $(project_name)Refine(DA *octDA,
                const std::vector<TREENODE> &treePart,
                const DomainExtents &domainExtents,
                $(project_name)InputData *inputData,
                SubDomainBoundary *subDomainBoundary);

    $(project_name)Refine(DA *octDA,
                const std::vector<TREENODE> &treePart,
                const DomainExtents &domainExtents,
                $(project_name)InputData *inputData,
                SubDomainBoundary *subDomainBoundary,
                std::vector<GEOMETRY::Geometry *> ibmGeoms,
                bool doCoarsen = false);
    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

    ~$(project_name)Refine() = default;
};

/**
 * @brief constructor
 * @param octDA
 * @param treePart tree partition
 * @param domainExtents domain information
 * @param inputData input data
 * @param subDomainBoundary boundary information (generateBoundaryFlags)
 */
$(project_name)Refine::$(project_name)Refine(DA *octDA,
                       const std::vector<TREENODE> &treePart,
                       const DomainExtents &domainExtents,
                       $(project_name)InputData *inputData,
                       SubDomainBoundary *subDomainBoundary)
    : Refinement(octDA, treePart, domainExtents), inputData_(inputData), domainExtents_(domainExtents), subDomainBoundary_(subDomainBoundary)
{
    this->traverse();
}
$(project_name)Refine::$(project_name)Refine(DA *octDA,
                       const std::vector<TREENODE> &treePart,
                       const DomainExtents &domainExtents,
                       $(project_name)InputData *inputData,
                       SubDomainBoundary *subDomainBoundary,
                       std::vector<GEOMETRY::Geometry *> ibmGeoms,
                       bool doCoarsen)
    : Refinement(octDA, treePart, domainExtents), inputData_(inputData), domainExtents_(domainExtents),
      subDomainBoundary_(subDomainBoundary), ibmGeoms_(std::move(ibmGeoms)), doCoarsen_(doCoarsen)
{
    this->traverse();
}

/**
 * @brief The refinement flags per element by element. If nothing is returned, its set to NO_CHANGE.
 * @param fe the element
 * @param coords vector of coords
 * @return the flags for each elements
 */
ot::OCT_FLAGS::Refine $(project_name)Refine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){
    const DomainInfo &physDomain = domainExtents_.physicalDADomain;
    const unsigned int nodeNo = m_octDA->getNumNodesPerElement();
    const DENDRITE_UINT currentLevel = this->m_level;

    /// refine walls (maximum to the refine_h level)
    const double eps = 1e-13;
    unsigned int levelForWall = inputData_->meshDef.refineLevel_base;
    const auto &refine_wall = inputData_->meshDef.refine_walls;
    if (inputData_->meshDef.refine_any_wall){
    if (refine_wall[0]){
        for (const auto &p : coords){
            if (p.x() - inputData_->meshDef.physDomain.min[0] < eps and currentLevel < inputData_->meshDef.refineLevel_channel_wall){
                levelForWall = inputData_->meshDef.refineLevel_channel_wall;
                break;
            }
        }
    }
    if (refine_wall[1]){
        for (const auto &p : coords){
            if (inputData_->meshDef.physDomain.max[0] - p.x() < eps and currentLevel < inputData_->meshDef.refineLevel_channel_wall){
                levelForWall = inputData_->meshDef.refineLevel_channel_wall;
                break;
            }
        }
    }
    if (refine_wall[2]){
        for (const auto &p : coords){
            if (p.y() - inputData_->meshDef.physDomain.min[1] < eps and currentLevel < inputData_->meshDef.refineLevel_channel_wall){
                levelForWall = inputData_->meshDef.refineLevel_channel_wall;
                break;
            }
        }
    }
    if (refine_wall[3]){
        for (const auto &p : coords){
            if (inputData_->meshDef.physDomain.max[1] - p.y() < eps and currentLevel < inputData_->meshDef.refineLevel_channel_wall){
                levelForWall = inputData_->meshDef.refineLevel_channel_wall;
                break;
            }
        }
    }
#if (DIM == 3)
    if (refine_wall[4]){
        for (const auto &p : coords){
            if (p.z() - inputData_->meshDef.physDomain.min[2] < eps and currentLevel < inputData_->meshDef.refineLevel_channel_wall){
                levelForWall = inputData_->meshDef.refineLevel_channel_wall;
                break;
            }
        }
    }
    if (refine_wall[5]){
        for (const auto &p : coords){
            if (inputData_->meshDef.physDomain.max[2] - p.z() < eps and currentLevel < inputData_->meshDef.refineLevel_channel_wall){
                levelForWall = inputData_->meshDef.refineLevel_channel_wall;
                break;
            }
        }
    }
#endif
    }

    /// region refine
    auto &rf = inputData_->region_refine;
    unsigned int maxlevelForRegion = inputData_->meshDef.refineLevel_base;
    std::vector<unsigned int> levelForRegions;
    levelForRegions.resize(rf.size());
    for (int i = 0; i < rf.size(); i++) {
        auto &r = rf[i];
        if (!r.forRetain) {
            for (const auto &p : coords){
                if (r.in_region(p)){
                    levelForRegions[i] = r.refine_region_lvl;
                    break;
                }else{
                    levelForRegions[i] = inputData_->meshDef.refineLevel_base;
                }
            }
        }
    }
    if (!levelForRegions.empty()) {
        maxlevelForRegion = *max_element(levelForRegions.begin(), levelForRegions.end());
    }
    // for region with geometry, refine the boundries.
    std::vector<unsigned int> levelForRegionsGeomBoundary;
    std::vector<int> countOfInPointsEachRegionGeom;
    levelForRegionsGeomBoundary.resize(rf.size());
    countOfInPointsEachRegionGeom.resize(rf.size());
    for (unsigned int i = 0; i < rf.size(); i++) {
        auto &r = rf[i];
        if (!r.forRetain and r.GetRefineType() == RegionalRefine::MESHOBJECT){
            for (const auto &p : coords){
                if (r.in_region(p)){
                    countOfInPointsEachRegionGeom[i]++;
                }
            }
        }
    }
    for (unsigned int i = 0; i < rf.size(); i++) {
        if (not(countOfInPointsEachRegionGeom[i] == nodeNo || countOfInPointsEachRegionGeom[i] == 0)) {
            levelForRegionsGeomBoundary[i] = (rf[i].refine_region_lvl_boundary);
        }
    }
    levelForRegionsGeomBoundary.push_back(maxlevelForRegion);
    if (!levelForRegionsGeomBoundary.empty()){
        maxlevelForRegion = *max_element(levelForRegionsGeomBoundary.begin(), levelForRegionsGeomBoundary.end());
    }
    
    /// region retain (for complete outside elements outside retain region, we don't want to refine those)
    bool outsideRetain = false;
    for (auto &r : rf){
        if (r.forRetain) {
            bool all_out = true;
            for (const auto &p : coords){
                all_out = all_out and (r.out_retain(p) == RegionalRefine::OUTSIDE);
            }
            outsideRetain = outsideRetain or all_out;
        }
    }

    DENDRITE_UINT id = -1;
    bool isObject = false;
    if (this->m_BoundaryOctant) {
        for (DENDRITE_UINT i = 0; i < m_octDA->getNumNodesPerElement(); i++) {
            subDomainBoundary_->generateBoundaryFlags(coords[i], id);
            if (subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::SPHERE) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::BOX) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::CIRCLE))
            {
                isObject = true;
                break;
            }
        }
    }
    unsigned int maxlevelForCarvedOutGeom = inputData_->meshDef.refineLevel_base;
    if (isObject) {
        maxlevelForCarvedOutGeom = std::max(inputData_->carved_out_geoms_def.at(id).refineLevel, maxlevelForCarvedOutGeom);
    }

    unsigned int maxlevelForGeom = inputData_->meshDef.refineLevel_base;
    
    // refine based on custom criteria /////////////////////////////////////////////////////
    unsigned int customCriteriaLevel = inputData_->meshDef.refineLevel_base;

$(refinement_criteria)
    
    ////////////////////////////////////////////////////////////////////////////////////////
    
    //  if (outsideRetain) {
    //    levelForWall = inputData_->meshDef.refine_l;
    //  }
    if (!doCoarsen_){
        if (currentLevel < maxlevelForGeom or
            currentLevel < maxlevelForRegion or
            currentLevel < levelForWall or
            currentLevel < maxlevelForCarvedOutGeom or
            currentLevel < customCriteriaLevel) {
            return ot::OCT_FLAGS::Refine::OCT_REFINE;
        }else{
            return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
        }
    }else{
        if (currentLevel > maxlevelForGeom and
            currentLevel > maxlevelForRegion and
            currentLevel > levelForWall and
            currentLevel > maxlevelForCarvedOutGeom and
            currentLevel > customCriteriaLevel)
        {
            return ot::OCT_FLAGS::Refine::OCT_COARSEN;
        }
    }

    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
}

"""
    println(file, content);
end

#=
Reads and maintains configuration input from config.txt or default generated values
=#
function dendrite_inputdata_file(var)
    
    # This file is now a static file. See the end of this file. 
    
    
#     project_name = finch_state.project_name;
#     file = add_generated_file(project_name*"InputData.h", dir="include");
    
#     params = finch_state.target_parameters;
#     dim = finch_state.config.dimension;
#     order = finch_state.config.basis_order_min;
#     if order > 3
#         printerr("This target only supports element order up to 3. Changing to 3.")
#         order = 3;
#     end
#     if order == 1
#         default_order = "linear";
#     elseif order == 2
#         default_order = "quadratic";
#     else # order == 3
#         default_order = "cubic";
#     end
#     default_matfree = finch_state.config.linalg_matrixfree;
#     default_refine = 3;
#     if haskey(params, :refineLevel)
#         default_refine = params[:refineLevel];
#     end
    
#     time_stepping_parts = "";
#     time_stepping_read = "";
#     if finch_state.prob.time_dependent
#         dt = finch_state.time_stepper.dt;
#         nsteps = finch_state.time_stepper.Nsteps;
#         finalT = dt*nsteps;
#         time_stepping_parts = "        // Time stepping info (can be overridden by config.txt values)\n";
#         time_stepping_parts *= "        dt = $(dt);\n";
#         time_stepping_parts *= "        unsigned int nSteps = $(nsteps);\n";
#         time_stepping_parts *= "        double finalT = $(finalT);\n";
        
#         time_stepping_read = "        // Time stepping info\n";
#         time_stepping_read *= "        if (ReadValue(\"dt\", dt)) {};\n";
#         time_stepping_read *= "        if (ReadValue(\"nSteps\", nSteps)) {};\n";
#         time_stepping_read *= "        if (ReadValue(\"finalT\", finalT)) {};\n";
#     end
    
#     content = """
# #pragma once

# #include <talyfem/input_data/input_data.h>
# #include <DataTypes.h>
# #include <point.h>

# /**
#  * Parameters for mesh
#  */
# struct MeshDef : public DomainInfo {
#     /// refinement level
#     DENDRITE_UINT refineLevel = $(default_refine);
    
#     void read_from_config(const libconfig::Setting &root) {
#         if (root.exists("refineLevel")) {
#             refineLevel= (DENDRITE_UINT) root["refineLevel"];
#         }
#         if (root.exists("min")) {
#             for (DENDRITE_UINT dim = 0; dim < $(dim); dim++) {
#                 min[dim] = (DENDRITE_REAL) root["min"][dim];
#             }
#         }else{
#             TALYFEMLIB::PrintWarning("No min defined for mesh. Using 0.0");
#             for (DENDRITE_UINT dim = 0; dim < $(dim); dim++) {
#                 min[dim] = (DENDRITE_REAL) 0.0;
#             }
#         }
#         if (root.exists("max")) {
#             for (DENDRITE_UINT dim = 0; dim < $(dim); dim++) {
#                 max[dim] = (DENDRITE_REAL) root["max"][dim];
#             }
#         }else if (root.exists("scalingFactor")) {
#             for (DENDRITE_UINT dim = 0; dim < $(dim); dim++) {
#                 max[dim] = min[dim] + (DENDRITE_REAL) root["scalingFactor"][dim];
#             }
#         }else{
#             TALYFEMLIB::PrintWarning("No max or scalingFactor defined for mesh. Using 1.0");
#             for (DENDRITE_UINT dim = 0; dim < $(dim); dim++) {
#                 max[dim] = (DENDRITE_REAL) 1.0;
#             }
#         }
#     }
# };

# class $(project_name)InputData : public TALYFEMLIB::InputData {
#     public:
#         /// Mesh definition
#         MeshDef meshDef;
#         /// Basis function
#         std::string basisFunctionStr = "$(default_order)";
#         /// Matrix  free
#         bool mfree = $(default_matfree);
#         /// Time info
#         double dt = 0.0;
#         double currentT = 0.0;
        
# $(time_stepping_parts)
        
#         SolverOptions solverOptions$(project_name);
        
#         bool dump_vec = false; // This is used for regression testing
        
#     bool ReadFromFile(const std::string &filename = std::string("config.txt")) {
#         ReadConfigFile(filename);
        
#         /// mesh size and level
#         meshDef.read_from_config(cfg.getRoot()["background_mesh"]);
        
#         /// basis function order
#         basisFunction = TALYFEMLIB::basis_string_to_enum(basisFunctionStr);
#         if (ReadValue("basisFunction", basisFunctionStr)) {
#             if (!basisFunctionStr.empty()){
#                 basisFunction = TALYFEMLIB::basis_string_to_enum(basisFunctionStr);
#             }
#         }
        
#         /// Other parameters
#         if (ReadValue("mfree", mfree)) {}
#         if (ReadValue("dump_vec", dump_vec)) {}
        
# $(time_stepping_read)
            
#         /// PETSc options
#         solverOptions$(project_name) = read_solver_options(cfg, "solver_options");
#         return true;
#     }

#     /// check if the input are valid
#     bool CheckInputData() {
#         /// Matrix free version of the code cannot have pre-conditioner
#         if (mfree) {
#             if (solverOptions$(project_name).vals.count("pc_type") == 1) {
#                 solverOptions$(project_name).vals.at("pc_type") = "none";
#                 TALYFEMLIB::PrintWarning("mfree = True, changing pc_type to 'none'");
#             }
#         }
        
#         /// Basis function order
#         if (basisFunction < 1 || basisFunction > 4){
#             basisFunction = TALYFEMLIB::basis_string_to_enum("$(default_order)");
#             TALYFEMLIB::PrintWarning("invalid basis function selection, changing order to $(default_order)");
#         }
#         return true;
#     }
# };

# """
#     println(file, content);
end

#=
The cmakelists.txt and readme.txt files
as well as the config.txt file
=#
function dendrite_build_files()
    project_name = finch_state.project_name;
    upper_name = uppercase(project_name);
    cmake_file = add_generated_file("CMakeLists.txt", make_header_text=false);
    readme_file = add_generated_file("README.txt", make_header_text=false);
    config_file = add_generated_file("config.txt", make_header_text=false);
    
    config = finch_state.config;
    prob = finch_state.prob;
    stepper = finch_state.time_stepper;
    params = finch_state.target_parameters;
    dim = config.dimension;
    
    # CMakeLists.txt
    enable2D = (dim == 2) ? "ON" : "OFF";
    enable3D = (dim == 3) ? "ON" : "OFF";
    
    content = """
## CMAKE for $(project_name) with Dendrite-KT
cmake_minimum_required(VERSION 2.8)
project($(project_name))
set(CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} -std=c++14")
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

set(CMAKE_MODULE_PATH "\${CMAKE_CURRENT_SOURCE_DIR}/dendrite-kt/cmake-modules")

# For now we just make it compulsory to have LAPACK installed.
# Later we will make it possible if LAPACK is not present to automaticall install before compiling dendro5
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "\${CMAKE_C_FLAGS} \${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} \${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "\${CMAKE_EXE_LINKER_FLAGS} \${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "\${COMPILE_FLAGS} \${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set(LINK_FLAGS "\${LINK_FLAGS} \${MPI_LINK_FLAGS}")
endif()

# # options for dendro
option(USE_64BIT_INDICES "Use 64-Bit indices. Reverts to 32-bit if turned off" ON)
option(ALLTOALLV_FIX "Use K-way all to all v" ON)
option(SPLITTER_SELECTION_FIX "Turn on Splitter Selection fix" ON)
option(DIM_2 "use the two dimentional sorting" OFF)
option(WITH_BLAS_LAPACK "build using BLAS and LAPACk" ON)
option(MANUAL_BLAS_LAPACK "configure BLAS and LAPACK Manually" OFF)
option(DENDRO_VTK_BINARY "write vtk/vtu files in binary mode " ON)
option(DENDRITE_VTU_ASCII "write vtk/vtu files in ASCII mode " OFF)
option(DENDRO_VTK_ZLIB_COMPRES "write vtk/vtu files in binary mode with zlib compression (only compatible with binary mode) " OFF)
option(BUILD_WITH_PETSC " build dendro with PETSC " ON)
option(HILBERT_ORDERING "use the Hilbert space-filling curve to order orthants" OFF)
option(BUILD_EXAMPLES "build example programs" OFF)
option(ENABLE_4D "enable 4D computation" OFF)
option(ENABLE_2D "enable 2D computation" $(enable2D))
option(ENABLE_3D "enable 3D computation" $(enable3D))
option(TENSOR "Use Tensor Operation" OFF)
option(PROFILING "Enable profiling" OFF)
option(IBM "Enable IBM Functions" ON)
option(SHBM "Enable SBM Functions" ON)

set(KWAY 128 CACHE INT 128)
set(NUM_NPES_THRESHOLD 2 CACHE INT 2)

# set the build type to release by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
        "Choose the type of build, options are: Debug Release " FORCE)
endif()

if(WITH_BLAS_LAPACK)
    add_definitions(-DWITH_BLAS_LAPACK)

    if(DEFINED ENV{MKLROOT})
        find_package(LAPACK COMPONENTS MKL REQUIRED)
        set(LAPACK_LIBRARIES \${MKL_LIBRARIES})

        if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
            set(CMAKE_EXE_LINKER_FLAGS "\${CMAKE_EXE_LINKER_FLAGS} -mkl")
        else()
            set(CMAKE_EXE_LINKER_FLAGS "\${CMAKE_EXE_LINKER_FLAGS} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm")
        endif()

        message(STATUS \${LAPACK_LIBRARIES})
    elseif(MANUAL_BLAS_LAPACK)
        if("\$ENV{BLAS}" STREQUAL "")
            message("Environment Variable BLAS is not set. Please set it to BLAS directory")
        endif()

        if("\$ENV{LAPACK}" STREQUAL "")
            message("Enviroment Variable LAPACK is note set. Please set it to LAPACK directory. ")
        endif()

        set(LAPACKE_DIR \$ENV{LAPACK}/LAPACKE)
        set(BLAS_LIBS \$ENV{BLAS}/lib)
        set(LAPACK_LIBS \$ENV{LAPACK}/lib)
        set(LAPACK_LINKER_FLAGS -llapacke -llapack -lblas -lgfortran -lquadmath)
        set(LAPACK_LIBRARIES \${LAPACK_LIBS}/liblapacke.a \${LAPACK_LIBS}/liblapack.a \${BLAS_LIBS}/libblas.a -static libgfortran.a libquadmath.a)
        set(LINK_FLAGS "\${LINK_FLAGS} \${LAPACK_LINKER_FLAGS}")
    else()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        set(LAPACK_LINKER_FLAGS -llapacke -llapack -lblas -lgfortran -lquadmath)
        set(LAPACKE_DIR \$ENV{LAPACK}/LAPACKE)
        set(LINK_FLAGS "\${LINK_FLAGS} \${LAPACK_LINKER_FLAGS}")
        find_library(LAPACKE_LIB
            NAMES lapacke lapackelib liblapacke
            HINTS "/usr/lib/"
        )
        set(LAPACK_LIBRARIES \${LAPACK_LIBRARIES} \${LAPACKE_LIB})
        message(STATUS \${LAPACK_LIBRARIES})
    endif()
endif()

if(BUILD_WITH_PETSC)
    find_package(PETSc REQUIRED)
    add_definitions(-DBUILD_WITH_PETSC)
endif()

if(TENSOR)
    add_definitions(-DTENSOR)
    message("Enabling Tensor operation")
endif()

if(DIM_2)
    add_definitions(-DDIM_2)
endif()

if(PROFILING)
    add_definitions(-DPROFILING)
endif()

if(USE_64BIT_INDICES)
    add_definitions(-DUSE_64BIT_INDICES)

    # message('Configuring 64BIT indices')
endif()

if(ALLTOALLV_FIX)
    add_definitions(-DALLTOALLV_FIX)
    add_definitions(-DKWAY=\${KWAY})
endif()

if(SPLITTER_SELECTION_FIX)
    add_definitions(-DSPLITTER_SELECTION_FIX)
    add_definitions(-DNUM_NPES_THRESHOLD=\${NUM_NPES_THRESHOLD})
endif()

if(ALLTOALL_SPARSE)
    add_definitions(-DALLTOALL_SPARSE)
endif()

if(DENDRO_VTK_BINARY)
else()
    set(DENDRO_VTK_ZLIB_COMPRES OFF)
endif()

if(DENDRITE_VTU_ASCII)
    add_definitions(-DDENDRITE_VTU_ASCII)
    message("Writing in ASCII format")
else()
    message("Writing in Binary format")
endif()

if(DENDRO_VTK_BINARY)
    add_definitions(-DDENDRO_VTU_BINARY)

    if(DENDRO_VTK_ZLIB_COMPRES)
        add_definitions(-DDENDRO_VTU_ZLIB)
    endif()
else()
    add_definitions(-DDENDRO_VTU_ASCII)
endif()

if(ENABLE_4D)
    add_definitions(-DENABLE_4D)
    message("Enabling 4D computation")
    set(ENABLE_4D,ON)
    set(ENABLE_3D,OFF)
elseif(ENABLE_2D)
    add_definitions(-DENABLE_2D)
    message("Enabling 2D computation")
    set(ENABLE_2D,ON)
    set(ENABLE_3D,OFF)
elseif(ENABLE_3D)
    add_definitions(-DENABLE_3D)
    set(ENABLE_3D,ON)
    message("Enabling 3D computation")
endif()

if(HILBERT_ORDERING)
    add_definitions(-DHILBERT_ORDERING)
endif()

if(IBM)
    add_definitions(-DIBM)
    message("Enabling IMGA Computation")
endif()

if(SHBM)
    add_definitions(-DSHBM)
    message("Enabling SBM Computation")
endif()

set(DENDRITEkT_BUILD_EXAMPLES OFF CACHE BOOL "Build Dendrite examples")
add_subdirectory(dendrite-kt "\${CMAKE_CURRENT_BINARY_DIR}/dendrite-kt")

set($(upper_name)_INC
    include/util.h
    include/SBMcalc.h
    include/$(project_name)InputData.h
    include/$(project_name)NodeData.h
    include/$(project_name)BoundaryConditions.h
    include/$(project_name)Equation.h
    include/$(project_name)Refine.h
    include/$(project_name)InputDataStructs.h
)

set($(upper_name)_SRC
    src/$(project_name).cpp)

add_executable($(project_name) \${$(upper_name)_SRC} \${$(upper_name)_INC})

target_include_directories($(project_name) PUBLIC include)
target_link_libraries($(project_name) dendriteKT dendroKT \${LAPACK_LIBRARIES} \${MPI_LIBRARIES} m)
    
""";
    println(cmake_file, content);
    
    # config.txt
    dt = string(stepper.dt);
    nSteps = string(stepper.Nsteps);
    totalT = string(stepper.dt * stepper.Nsteps);
    stype = stepper.type;
    if haskey(params, :refineLevel)
        refine_level = params[:refineLevel];
    else
        refine_level = 4;
    end
    if haskey(params, :refineLevel_channel_wall)
        refine_level_channel = params[:refineLevel_channel_wall];
    else
        refine_level_channel = refine_level;
    end
    if haskey(params, :min)
        domain_min = params[:min];
    else
        domain_min = dim == 2 ? "[0.0, 0.0]" : "[0.0, 0.0, 0.0]";
    end
    if haskey(params, :max)
        domain_max = params[:max];
    else
        domain_max = dim == 2 ? "[1.0, 1.0]" : "[1.0, 1.0, 1.0]";
    end
    if haskey(params, :refineWalls)
        refine_walls = params[:refineWalls];
    else
        refine_walls = "[false, false, false, false, false, false]";
    end
    if haskey(params, :matrixFree)
        mat_free = params[:matrixFree];
    else
        mat_free = config.linalg_matrixfree;
    end
    
    # Geometries is an array of geometry chunks, which are dicts with various info
    if haskey(params, :geometries)
        geometries = params[:geometries];
        ngeos = length(geometries);
        geometry_part = "";
        for i=1:ngeos
            # Use default values if not provided
            meshFile = get!(geometries[i], :meshFile, "UNKNOWN_MESH_FILE");
            meshName = get!(geometries[i], :meshName, "UNKNOWN_MESH_NAME");
            position = get!(geometries[i], :position, dim == 2 ? "[0.0, 0.0]" : "[0.0, 0.0, 0.0]");
            outer_bdry = get!(geometries[i], :outerBoundary, "true");
            geo_type = get!(geometries[i], :geoType, dim == 2 ? "meshobject_2d" : "meshobject");
            geo_refine_level = get!(geometries[i], :refineLevel, refine_level);
            bdry_type = get!(geometries[i], :boundaryTypes, "[]");
            bids = get!(geometries[i], :bids, "[]");
            geometry_part *= "    {\n";
            geometry_part *= "    mesh_path = \"$(meshFile)\"\n";
            geometry_part *= "    name = \"$(meshName)\"\n";
            geometry_part *= "    type = \"$(geo_type)\"\n";
            geometry_part *= "    position = $(position)\n";
            geometry_part *= "    refineLevel = $(geo_refine_level)\n";
            geometry_part *= "    bc_type_V = $(bdry_type)\n";
            geometry_part *= "    bid_V = $(bids)\n";
            geometry_part *= "    outer_boundary = $(outer_bdry)\n";
            geometry_part *= "    is_static = true\n";
            geometry_part *= "    }";
            
            if i < ngeos
                geometry_part *= ",";
            end
            geometry_part *= "\n";
        end
    else
        geometry_part = "";
    end
    
    # These are the minimum and default solver options
    solver_options = Dict([("ksp_max_it", "1000"),
                          ("ksp_type", "\"bcgs\""),
                          ("pc_type", "\"lu\""),
                          ("ksp_atol", "1e-8"),
                          ("ksp_rtol", "1e-8")]);
    # Include and solver options supplied in params
    # check these possibilities. This list needs to be expanded
    petsc_solver_ops = ["ksp_monitor", "ksp_converged_reason"];
    for op in petsc_solver_ops
        if haskey(params, Symbol(op))
            if typeof(params[Symbol(op)]) == String
                solver_options[op] = "\""*string(params[Symbol(op)])*"\"";
            else
                solver_options[op] = string(params[Symbol(op)]);
            end
        end
    end
    
    solver_options_str = "";
    for key in keys(solver_options)
        solver_options_str *= "    " * key * " = " * solver_options[key] * "\n";
    end
    
    content = """

dt = $(dt)
nSteps = $(nSteps)
totalT = $(totalT)

timeStepper="$(stype)"
solverType = "rbvms"
matrixFree = $(mat_free)

# The full background mesh
elemOrder = $(config.basis_order_min)
background_mesh = {
    refineLevel = $(refine_level)
    refineLevel_channel_wall = $(refine_level_channel)
    min = $(domain_min)
    max = $(domain_max)
    refine_walls = $(refine_walls)
}

# penalty
Cb_e = 200

# Carved out sub-mesh
geometries = (
$(geometry_part)
)

# This may be used later
region_refine = (

)

# Output
OutputStartTime = 0
OutputInterval = 1
SurfaceMonitor = [2]

#################### solver setting ####################
solver_options = {
$(solver_options_str)
}
    
""";
    println(config_file, content);
    
    # readme
    content = """
Build instructions coming soon.
""";
    println(readme_file, content);
    
end

###########################################################################################################
# Static files that will be written as is with perhaps only the project name added.

 #####   ######      ###     ######  ######   ######       #######  ######  ##        #######   #####
###   #    ##       ## ##      ##      ##    ##    ##      ##         ##    ##        ##       ###   #
  ###      ##      ##   ##     ##      ##    ##            #######    ##    ##        ######     ### 
#   ###    ##     #########    ##      ##    ##    ##      ##         ##    ##        ##       #   ###
 #####     ##    ##       ##   ##    ######   ######       ##       ######  ########  #######   #####

###########################################################################################################

function dendrite_static_files()
    project_name = finch_state.project_name;
    # These are the files that will be written here
    inputdata_file = add_generated_file(project_name*"InputData.h", dir="include");
    inputdatastructs_file = add_generated_file(project_name*"InputDataStructs.h", dir="include");
    utils_file = add_generated_file("util.h", dir="include");
    sbmcalc_file = add_generated_file("SBMcalc.h", dir="include");
    
    ###########################################################################################
    # InputData
    
    content = """
#pragma once

#include "$(project_name)InputDataStructs.h"
#include <time.h>

class $(project_name)InputData : public TALYFEMLIB::InputData{
public:
    static constexpr int nsd = DIM;

    bool Penalty = false;
    bool Special = false;
    bool Unsymmetric = false;
    bool InsideSBM = false;
    bool PenaltyNoShift = false;
    bool IfAdjointConsistency = true;
    bool IfNegativeAdjointConsistency = false;
    bool NormTheSame = false;

    bool PrintTrueSurfaceGP = false;
    bool PrintTrueSurfaceNP = false;
    bool SecondOrderTaylorQBF = true;
    bool SuperConvergencePt = false;
    bool HessianInNeumann = true;

    /// for SBM dist function calculation
    TALYFEMLIB::ZEROPTV PBoxStart;
    TALYFEMLIB::ZEROPTV PBoxEnd;
    std::vector<std::vector<ZEROPTV>> DistributePoints;
    std::vector<std::vector<int>> TriangleNumber;
    std::vector<ZEROPTV> GPPTVAll;

    /// solver (STABILIZED) and Variational multiscale based solver (RBVMS)
    enum typeSolver{
        STABILIZED = 0,
        RBVMS = 1
    };
    
    enum typeNondimension{
        FREECONV = 0,
        MIXCONV = 1
    };

    enum typeDistCalc{
        NORMAL_BASED = 0,
        NORMAL_BASED_DistributeSTL = 1,
        GP_BASED = 2
    };

    /// Declare the Solver type
    typeSolver solverType;
    
    /// Declare the dist calc type
    typeDistCalc DistCalcType;

    DENDRITE_UINT elemOrder = 1;
    bool ifMatrixFree = false;

    /// Time stepper
    std::vector<double> dt;
    std::vector<double> totalT;
    std::vector<int> nSteps;
    double OutputStartTime;
    int OutputInterval;

    /// Postprocessing
    int PostProcessingInterval = 1;
    bool printToScreen = true;
    bool writeToFile = true;

    /// weakBC parameters
    double Cb_e = 200.0;
    double Cb_f = 20.0;
    double C_DC = 0.0;
    
    /// VMS parameters may be removed if not used
    double Ci_f = 36;
    double Ci_e = 36;
    double tauM_scale = 1.0;
    double tauM_scale4Side = 1.0;
    double tauM_time = 1.0;

    /// Timestepper control
    enum typeTimestepper{
        GENERALIZED_THETA = 0,
        BDF2 = 1
    };
    /// Declare the timestepper
    typeTimestepper timestepperType;
    double thetaTimeStepping = 1;      /// Default is backward Euler.
    double thetaTimeStepping4Side = 1; /// Default is backward Euler.

    /// tolerance for the block iterations
    double blockTolerance;
    /// maximum number of iterations for block iteration
    unsigned int iterMaxBlock;

    /// postprocessing (surface indicator to monitor)
    std::vector<unsigned int> SurfaceMonitor;

    std::vector<BoundaryDef> boundary_def;
    InitialConditionDef ic_def;

    /// Curved-out geometries
    std::vector<CarvedOutGeom> carved_out_geoms_def;

    /// Solver options for PETSc are handled in these structures
    SolverOptions solverOptions;

    /// Limiter
    Limiter limiter;

    /// Debug options
    bool Debug_Integrand = false;
    bool dump_vec = false;

    /// Setup the meshDef object for subDA parameters
    MeshDef meshDef;
    std::vector<RegionalRefine> region_refine;

    ~$(project_name)InputData() = default;

    bool ReadFromFile(const std::string &filename = std::string("config.txt")){
        ReadConfigFile(filename);
        ReadValue("elemOrder", elemOrder);
        ReadValue("matrixFree", ifMatrixFree);
        
        /// SubDA (channel parameters)
        meshDef.read_from_config(cfg.getRoot()["background_mesh"]);
        if (cfg.exists("region_refine")){
            const auto &cfg_refine = cfg.getRoot()["region_refine"];
            region_refine.resize(cfg_refine.getLength());
            for (unsigned int i = 0; i < region_refine.size(); i++){
                region_refine[i].read_from_config(cfg_refine[i]);
            }
        }

        /// timestep control
        ReadVectorOrValue("dt", dt);
        ReadVectorOrValue("totalT", totalT);
        ReadVectorOrValue("nSteps", nSteps);
        
        /// Output control
        if (ReadValue("OutputStartTime", OutputStartTime)){}
        if (ReadValue("OutputInterval", OutputInterval)){}

        //if (ReadValue("PostProcessingInterval", PostProcessingInterval)){}
        //if (ReadValue("printToScreen", printToScreen)){}
        //if (ReadValue("writeToFile", writeToFile)){}

        /// Solver selection
        /// Read type of solver
        solverType = read_solver(cfg.getRoot(), "solverType");
        
        /// Timestepper option
        timestepperType = read_Timestepper(cfg.getRoot(), "timeStepper");

        /// VMS parameters
        if (solverType == typeSolver::RBVMS){
            //if (ReadValue("Ci_f", Ci_f)){};
            //if (ReadValue("Ci_e", Ci_e)){};
            //if (ReadValue("tauM_scale", tauM_scale)){};
            //if (ReadValue("tauM_scale4Side", tauM_scale4Side)){};
            //if (ReadValue("tauM_time", tauM_time)){};
        }
    
        /// WeakBC parameters
        if (ReadValue("Cb_e", Cb_e)){};
        //if (ReadValue("Cb_f", Cb_f)){};
        //if (ReadValue("C_DC", C_DC)){};

        /// read geometries
        if (cfg.exists("geometries")){
            const auto &geometries = cfg.getRoot()["geometries"];
            carved_out_geoms_def.resize(geometries.getLength());
            for (unsigned int i = 0; i < carved_out_geoms_def.size(); i++){
                carved_out_geoms_def[i].read_from_config(geometries[i]);
            }
        }

        /// Solver Options
        /// Simulation parameters and tolerances
        solverOptions = read_solver_options(cfg, "solver_options");

        /// Load surface postprocessing vector, or Nu calculation
        if (cfg.exists("SurfaceMonitor")){
            const libconfig::Setting &settings = cfg.getRoot()["SurfaceMonitor"];
            for (int i = 0; i < settings.getLength(); ++i){
                SurfaceMonitor.push_back(settings[i]);
            }
        }

        return true;
    }

    /// Function for reading a vector or a single value (stored in vector)
    template <typename T>
    void ReadVectorOrValue(const std::string &key_name, std::vector<T> &value){
        if (cfg.exists(key_name + "_V")){
            InputData::ReadVector(cfg, key_name + "_V", value);
        }else{
            double value_const;
            ReadValueRequired(key_name, value_const);
            value.push_back(value_const);
        }
    }
    
    bool CheckInputData(){
    	return true;
    }
    
    /**
    * Printout every item of inputdata for debug purpose.
    */
    void PrintInputData(){
        int rank = TALYFEMLIB::GetMPIRank();
        if (!rank){
            std::ofstream fout("InputDataOutput.txt", std::ios::app);
            time_t my_time = time(NULL);
            fout << "##############################"
                << "\\n";
            fout << ctime(&my_time);
            fout << "Total number of processor = " << TALYFEMLIB::GetMPISize() << "\\n";
            fout << "size of DendroInt " << sizeof(DendroIntL) << "\\n";
            fout << "size of PetscInt " << sizeof(PetscInt) << "\\n";

            fout << "Dimension: " << nsd << "\\n";
            fout << "basisFunctionOrder: " << elemOrder << "\\n";
            fout << "mfree: " << ifMatrixFree << "\\n";
            if (solverType == typeSolver::STABILIZED){
                fout << "solverType: STABILIZED\\n";
            }else{
                fout << "solverType: RBVMS\\n";
            }

            fout << "====== meshDef ======"
                << "\\n";
            meshDef.PrintMeshDef(fout);
            fout << "====================="
                << "\\n\\n";

            fout << "====== timestepper ======"
                << "\\n";
            PrintVector(fout, "dt", dt);
            PrintVector(fout, "totalT", totalT);
            PrintVector(fout, "nSteps", nSteps);
            if (timestepperType == typeTimestepper::GENERALIZED_THETA){
                fout << "timestepperType: GENERALIZED_THETA\\n";
            }else{
                fout << "timestepperType: BDF2\\n";
            }
            fout << "thetaTimeStepping: " << thetaTimeStepping << "\\n";
            fout << "thetaTimeStepping4Side: " << thetaTimeStepping4Side << "\\n";
            fout << "====================="
                << "\\n\\n";

            fout << "===== parameters ======="
                << "\\n";

            fout << "Cb_e: " << Cb_e << "\\n";
            fout << "====================="
                << "\\n\\n";

            fout << "\\nregion_refine: {\\n";
            for (const auto &r : region_refine)
            {
                r.PrintRegionRefineDef(fout);
                fout << "}\\n{\\n";
            }
            fout << "=========== extra stuff =========="
                << "\\n";
            // not very important
            fout << "OutputStartTime: " << OutputStartTime << "\\n";
            fout << "OutputInterval: " << OutputInterval << "\\n";
            fout << "PostProcessingInterval: " << PostProcessingInterval << "\\n";
            fout << "printToScreen: " << printToScreen << "\\n";
            fout << "writeToFile: " << writeToFile << "\\n";
            PrintVector(fout, "SurfaceMonitor", SurfaceMonitor);
            fout << "========== solver settings ==========="
                << "\\n";
            fout << "solverOptions: ["
                << "\\n";
            for (auto &val : solverOptions.vals){
                fout << val.first << " --> " << val.second << "\\n";
            }
            fout << "]"
                << "\\n";
            fout.close();
            }
        //    /// Curved-out geometries
        //    std::vector<CarvedOutGeom> carved_out_geoms_def;
    }

private:
    /// Function for reading type of Heat equation solver
    static typeSolver read_solver(libconfig::Setting &root, const char *name){
        std::string str;
        /// If nothing specified stays stabilized
        if (root.lookupValue(name, str)){
            if (str == "stabilized"){
                return STABILIZED;
            }else if (str == "rbvms"){
                return RBVMS;
            }else{
                throw TALYFEMLIB::TALYException() << "Unknown solver name: " << name << str;
            }
        }else{
            throw TALYFEMLIB::TALYException() << "Must specify solverType: stabilized or rbvms";
        }
    }

    /// Function for reading type of timestepper
    static typeTimestepper read_Timestepper(libconfig::Setting &root, const char *name){
        std::string str = "GeneralizedTheta";
        /// If nothing specified stays Generalized theta
        if (root.lookupValue(name, str)){
            if (str == "GeneralizedTheta"){
                return GENERALIZED_THETA;
            }else if (str == "BDF2"){
                return BDF2;
            }else{
                throw TALYFEMLIB::TALYException() << "Unknown timestepper: " << name << str;
            }
        }else{
            return GENERALIZED_THETA;
        }
    }
};
""";
    println(inputdata_file, content);
    ##############################################################################################
    
    ##############################################################################################
    # InputDataStructs
    
    content = """
#pragma once
// #include "util.h"
#include "IMGA/IMGA.h"
#include "$(project_name)NodeData.h"

/// read size of 3 vector as ZEROPTV
void ReadZEROPTV(const libconfig::Setting &root, const std::string &key_name, ZEROPTV &value, bool required = true){
  if (root.exists(key_name.c_str())){
    if (root[key_name.c_str()].getLength() >= DIM){
      value(0) = (double)root[key_name.c_str()][0];
      value(1) = (double)root[key_name.c_str()][1];
#if (DIM == 2)
      value(2) = 0.0;
#endif
#if (DIM == 3)
      value(2) = (double)root[key_name.c_str()][2];
#endif
    }else{
      throw TALYException() << key_name + " have size of " + std::to_string(root[key_name.c_str()].getLength());
    }
  }else if (required){
    throw TALYException() << key_name + " doesn't exist!";
  }
}

/// Function for reading a vector from root
template <typename T>
void ReadVectorRoot(const libconfig::Setting &root, const std::string &key_name, std::vector<T> &value){
  const libconfig::Setting &config_v = root[key_name.c_str()];
  if (config_v.isNumber()){
    value.push_back(config_v);
  }else{
    value.reserve(config_v.getLength());
    for (int i = 0; i < config_v.getLength(); ++i){
      value.emplace_back(config_v[i]);
    }
  }
}

/*
 * maximum value of an array
 */
double maximum(std::array<DENDRITE_REAL, DIM> &array){
  DENDRITE_REAL max = array[0];
  for (int i = 1; i < DIM; i++){
    if (array[i] > max){
      max = array[i];
    }
  }
  return max;
}

/*
 * minimum value of an array
 */
double minimum(std::array<DENDRITE_REAL, DIM> &array){
  DENDRITE_REAL min = array[0];
  for (int i = 1; i < DIM; i++){
    if (array[i] < min){
      min = array[i];
    }
  }
  return min;
}

template <typename T>
void PrintVector(std::ofstream &fstream, const std::string &name, const std::vector<T> &vec){
  int rank = TALYFEMLIB::GetMPIRank();
  if (!rank and fstream.is_open()){
    fstream << name << ": [";
    for (const auto &v : vec){
      fstream << v << ", ";
    }
    fstream << "]\\n";
  }
}

/**
 * Domain definition
 */
struct MeshDef{
  DomainInfo fullDADomain; /// The domain from which its carved out.
  DomainInfo physDomain;   /// The actual SubDAdomain

  std::array<DENDRITE_REAL, DIM> max; ///< max of domain
  std::array<DENDRITE_REAL, DIM> min; ///< min of domain

  ///< Dendro options
  DENDRITE_UINT refineLevel_base;
  DENDRITE_UINT refineLevel_channel_wall;
  /// refineLevel_interface

  bool refine_walls[6] = {false, false, false, false, false, false}; ///< whether or not fx_refine refines at min/max
  bool refine_any_wall = false;

  void read_from_config(const libconfig::Setting &root){

    refineLevel_base = (unsigned int)root["refineLevel"];

    max[0] = (DENDRITE_REAL)root["max"][0];
    max[1] = (DENDRITE_REAL)root["max"][1];
    min[0] = (DENDRITE_REAL)root["min"][0];
    min[1] = (DENDRITE_REAL)root["min"][1];
#if (DIM == 3)
    max[2] = (DENDRITE_REAL)root["max"][2];
    min[2] = (DENDRITE_REAL)root["min"][2];
#endif

    //if (!root.lookupValue("refine_walls", refine_walls)){
      // refine_walls = {false...}; // initialized to false
    //}else{
      refine_walls[0] = root["refine_walls"][0];
      refine_walls[1] = root["refine_walls"][1];
      refine_walls[2] = root["refine_walls"][2];
      refine_walls[3] = root["refine_walls"][3];
      refine_any_wall = refine_walls[0] || refine_walls[1] || refine_walls[2] || refine_walls[3];
#if (DIM == 3)
      refine_walls[4] = root["refine_walls"][4];
      refine_walls[5] = root["refine_walls"][5];
      refine_any_wall = refine_any_wall || refine_walls[4] || refine_walls[5];
#endif
    //}

    if (refine_any_wall){
      refineLevel_channel_wall = (unsigned int)root["refineLevel_channel_wall"];
      // refineLevel_pillar_wall = (unsigned int)root["refineLevel_pillar_wall"];
      if (refineLevel_base < 1 || refineLevel_base >= 31 ||
          refineLevel_channel_wall < 1 || refineLevel_channel_wall >= 31)
      {
        PrintWarning("Invalid refine_ level - should be 1 < refineLevel <= 31.");
      }
    }
    if (!refine_any_wall){
      refineLevel_channel_wall = refineLevel_base;
    }

    if (refineLevel_base > refineLevel_channel_wall){
      PrintWarning("refineLevel_base > refineLevel_channel_wall");
    }

    double minOfBox = minimum(min);
    double sizeOfBox = maximum(max);

    fullDADomain.min.fill(minOfBox);
    fullDADomain.max.fill(sizeOfBox);

    physDomain.min[0] = min[0];
    physDomain.min[1] = min[1];
    physDomain.max[0] = max[0];
    physDomain.max[1] = max[1];
#if (DIM == 3)
    physDomain.min[2] = min[2];
    physDomain.max[2] = max[2];
#endif

  }

  void PrintMeshDef(std::ofstream &fstream){
    int rank = TALYFEMLIB::GetMPIRank();
    if (!rank and fstream.is_open()){
      fstream << "max: [";
      for (int i = 0; i < DIM; i++){
        fstream << max[i] << ", ";
      }
      fstream << "]\\n";

      fstream << "fullDADomain-min: [";
      for (int i = 0; i < DIM; i++){
        fstream << fullDADomain.min[i] << ", ";
      }
      fstream << "]\\n";
      fstream << "fullDADomain-max: [";
      for (int i = 0; i < DIM; i++){
        fstream << fullDADomain.max[i] << ", ";
      }
      fstream << "]\\n";

      fstream << "physDomain-min: [";
      for (int i = 0; i < DIM; i++){
        fstream << physDomain.min[i] << ", ";
      }
      fstream << "]\\n";
      fstream << "physDomain-max: [";
      for (int i = 0; i < DIM; i++){
        fstream << physDomain.max[i] << ", ";
      }
      fstream << "]\\n";

      fstream << "refineLevel_base: " << refineLevel_base << "\\n";
      fstream << "refineLevel_channel_wall: " << refineLevel_channel_wall << "\\n";
      fstream << "refine_walls: [";
      for (int i = 0; i < DIM; i++){
        fstream << refine_walls[i*2] << ", ";
        fstream << refine_walls[i*2+1] << ", ";
      }
      fstream << "]\\n";
    }
  }
};

/**
 * Parameters for boundary condition
 */
struct BoundaryDef{
  enum Side{
    INVALID = -1,
    X_MINUS = 0,
    X_PLUS = 1,
    Y_MINUS = 2,
    Y_PLUS = 3,
    Z_MINUS = 4,
    Z_PLUS = 5,
  };
  enum Condition_Type{
    NO_GRADIENT_BC = 0,
    DIRICHLET_BC = 1,
    NEUMANN_BC = 2,
    ROBIN_BC = 3,
    WEAK_BC = 4,
    CSV_BC = 5,
    CSV_WEAK_BC = 6,
    SBM_BC = 7,
    NEUMANN_SBM_BC =8
  };

  Side side = INVALID;
  Condition_Type temperature_type = DIRICHLET_BC;
  /// for Dirichlet BC
  std::vector<DENDRITE_REAL> temperature;
  /// for Neumann BC
  std::vector<DENDRITE_REAL> flux;
  /// for robin BC
  std::vector<DENDRITE_REAL> G_vec;
  std::vector<DENDRITE_REAL> a_vec;
  
  bool ifRefine = false;

  static Side read_side_from_config(const libconfig::Setting &root){
    return str_to_side(root["side"]);
  }

  void read_from_config(const libconfig::Setting &root){
    if (root.exists("ifRefine")){
      ifRefine = (bool)root["ifRefine"];
    }
    if (root.exists("temperature_type")){
      temperature_type = read_temperature_type(root["temperature_type"]);
    }
    
    if ((temperature_type == Condition_Type::DIRICHLET_BC) or (temperature_type == Condition_Type::WEAK_BC) or (temperature_type == Condition_Type::SBM_BC)){
      ReadVectorRoot(root, "temperature", temperature);
    }

    if (temperature_type == Condition_Type::NEUMANN_BC){
      ReadVectorRoot(root, "flux", flux);
    }

    if (temperature_type == Condition_Type::ROBIN_BC){
      /// b * \\frac{\\partial u}{\\partial n} = g - a * u
      ReadVectorRoot(root, "G_constant", G_vec);
      ReadVectorRoot(root, "a_constant", a_vec);
    }
  }

  void PrintBoundaryDef(std::ofstream &fstream) const{
    int rank = TALYFEMLIB::GetMPIRank();
    if (!rank and fstream.is_open()){
      if (side == Side::INVALID){
        fstream << "side: INVALID\\n";
      }else if (side == Side::X_MINUS){
        fstream << "side: X_MINUS\\n";
      }else if (side == Side::X_PLUS){
        fstream << "side: X_PLUS\\n";
      }else if (side == Side::Y_MINUS){
        fstream << "side: Y_MINUS\\n";
      }else if (side == Side::Y_PLUS){
        fstream << "side: Y_PLUS\\n";
      }else if (side == Side::Z_MINUS){
        fstream << "side: Z_MINUS\\n";
      }else if (side == Side::Z_PLUS){
        fstream << "side: Z_PLUS\\n";
      }
      
      if (temperature_type == Condition_Type::NO_GRADIENT_BC){
        fstream << "temperature_type: NO_GRADIENT_T\\n";
      }else if (temperature_type == Condition_Type::DIRICHLET_BC){
        fstream << "temperature_type: DIRICHLET_T\\n";
      }else if (temperature_type == Condition_Type::NEUMANN_BC){
        fstream << "temperature_type: NEUMANN_T\\n";
      }else if (temperature_type == Condition_Type::ROBIN_BC){
        fstream << "temperature_type: ROBIN_T\\n";
      }else if (temperature_type == Condition_Type::WEAK_BC){
        fstream << "temperature_type: WEAK_T\\n";
      }

      PrintVector(fstream, "temperature", temperature);
      PrintVector(fstream, "flux", flux);
      PrintVector(fstream, "G_vec", G_vec);
      PrintVector(fstream, "a_vec", a_vec);
      fstream << "ifRefine: " << ifRefine << "\\n";
    }
  }

private:
  static Side str_to_side(const std::string &str){
    if (str == "x-"){
      return Side::X_MINUS;
    }else if (str == "x+"){
      return Side::X_PLUS;
    }else if (str == "y-"){
      return Side::Y_MINUS;
    }else if (str == "y+"){
      return Side::Y_PLUS;
    }else if (str == "z-"){
      return Side::Z_MINUS;
    }else if (str == "z+"){
      return Side::Z_PLUS;
    }else{
      throw TALYFEMLIB::TALYException() << "Invalid BC side";
    }
  }

  static Condition_Type read_temperature_type(const std::string &str){
    if (str == "no_gradient"){
      return Condition_Type::NO_GRADIENT_BC;
    }else if (str == "dirichlet"){
      return Condition_Type::DIRICHLET_BC;
    }else if (str == "neumann"){
      return Condition_Type::NEUMANN_BC;
    }else if (str == "robin"){
      return Condition_Type::ROBIN_BC;
    }else if (str == "weak"){
      return Condition_Type::WEAK_BC;
    }else if (str == "sbm"){
      return Condition_Type::SBM_BC;
    }else if (str == "csv"){
      return Condition_Type::CSV_BC;
    }else if (str == "csv_weak"){
      return Condition_Type::CSV_WEAK_BC;
    }else if (str == "neumann_sbm"){
        return Condition_Type::NEUMANN_SBM_BC;
    }else{
      throw TALYFEMLIB::TALYException() << "Invalid BC type for temperature";
    }
  }
};

/**
 * Parameters for initial condition
 */
struct InitialConditionDef{
  enum T_IC{
    ZERO_T = 0,
    USER_DEFINED_T = 1
  };

  T_IC t_ic_type = ZERO_T;
  double t_ic = 0.0;

  void read_from_config(const libconfig::Setting &root){
    if (root.exists("t_ic_type")){
      t_ic_type = read_temperature_ic(root["t_ic_type"]);
    }
    if (t_ic_type == T_IC::USER_DEFINED_T){
      t_ic = (double)root["t_ic"];
    }
  }

  void PrintInitialConditionDef(std::ofstream &fstream) const{
    int rank = TALYFEMLIB::GetMPIRank();
    if (!rank and fstream.is_open()){
      if (t_ic_type == T_IC::ZERO_T)
      {
        fstream << "t_ic_type: ZERO_T\\n";
      }
      else if (t_ic_type == T_IC::USER_DEFINED_T)
      {
        fstream << "t_ic_type: USER_DEFINED_T\\n";
      }
      fstream << "t_ic: " << t_ic << "\\n";
    }
  }

protected:
  static T_IC read_temperature_ic(const std::string &str){
    if (str == "zero"){
      return T_IC::ZERO_T;
    }else if (str == "user_defined"){
      return T_IC::USER_DEFINED_T;
    }else{
      throw TALYFEMLIB::TALYException() << "Invalid IC type for temperature";
    }
  }
};

/**
 * Parameters for Regional refine
 */
struct RegionalRefine{
  enum RetainType{
    OUTSIDE = 0,
    INSIDE = 1,
    ON_BOUNDARY = 2,
  };
  enum Type{
    INVALID = -1,
    CUBE = 0,
    SPHERE = 1,
    CYLINDER = 2,
    MESHOBJECT = 10,
    MESHOBJECT_2D = 11,
  };
  bool forRetain = false;
  // max refinement level
  int refine_region_lvl = 0;
  int refine_region_lvl_boundary = 0;
  std::string mesh_path;
  ZEROPTV shift;
  GEOMETRY::Geometry *geometry;
  GEOMETRY::STL *stl;
  GEOMETRY::MSH *msh;

  void read_from_config(const libconfig::Setting &root){
    if (root.exists("forRetain")){
      forRetain = (bool)root["forRetain"];
    }
    refine_type = read_type(root["type"]);
    if (forRetain){
      if (refine_type != CUBE and refine_type != MESHOBJECT){
        throw TALYException() << "Not supported retain type!";
      }
    }else{
      refine_region_lvl = (int)root["refine_region_lvl"];
      refine_region_lvl_boundary = refine_region_lvl;
      if (root.exists("refine_region_lvl_boundary"))
      {
        refine_region_lvl_boundary = (int)root["refine_region_lvl_boundary"];
      }
    }
    if (refine_type == CUBE){
      ReadZEROPTV(root, "min_c", min_corner);
      ReadZEROPTV(root, "max_c", max_corner);
    }
#if (DIM == 3)
    if (refine_type == MESHOBJECT){
      mesh_path = (const char *)root["mesh_path"];
      ReadZEROPTV(root, "shift", shift);

      stl = new GEOMETRY::STL(mesh_path, GEOMETRY::InOutTest::RAY_TRACING);
      std::array<DENDRITE_REAL, DIM> point;
      point[0] = shift[0];
      point[1] = shift[1];
      point[2] = shift[2];
      geometry = new GEOMETRY::Geometry(stl, Point<DIM>(point));
    }
#endif

#if (DIM == 2)
    if (refine_type == MESHOBJECT_2D){
      mesh_path = (const char *)root["mesh_path"];
      ReadZEROPTV(root, "shift", shift);

      msh = new GEOMETRY::MSH(mesh_path, GEOMETRY::InOutTest2D::RAY_TRACING_2D);
      std::array<DENDRITE_REAL, DIM> point;
      point[0] = shift[0];
      point[1] = shift[1];
      geometry = new GEOMETRY::Geometry(msh, Point<DIM>(point));
    }
#endif
    if (refine_type == SPHERE){
      ReadZEROPTV(root, "center", center_sphere);
      radius_sphere = (double)root["radius"];
      if (root.exists("radius_in")){
        radius_sphere_in = (double)root["radius_in"];
      }
    }
    if (refine_type == CYLINDER){
      ReadZEROPTV(root, "c1", min_corner);
      ReadZEROPTV(root, "c2", max_corner);
      radius_cylinder = (double)root["radius"];
      if (root.exists("radius_in")){
        radius_cylinder_in = (double)root["radius_in"];
      }
    }
  }

public:
  /**
   * @return true for outside retain cube (boundary nodes included (for apply BC)), false for inside retain cube
   */
  RetainType out_retain(ZEROPTV p){
    if (not(forRetain)){
      throw TALYException() << "Calling function with forRetain = false";
    }
    switch (refine_type){
    case CUBE:{
      const double eps = 1e-8;
      if (p.x() < min_corner.x() - eps ||
          p.y() < min_corner.y() - eps ||
          p.z() < min_corner.z() - eps ||
          p.x() > max_corner.x() + eps ||
          p.y() > max_corner.y() + eps ||
          p.z() > max_corner.z() + eps)
      {
        return OUTSIDE;
      }else if (p.x() > min_corner.x() + eps &&
               p.y() > min_corner.y() + eps &&
               p.z() > min_corner.z() + eps &&
               p.x() < max_corner.x() - eps &&
               p.y() < max_corner.y() - eps &&
               p.z() < max_corner.z() - eps)
      {
        return INSIDE;
      }else{
        return ON_BOUNDARY;
      }
    }
    break;
    
    case MESHOBJECT:{
      bool ifInside = geometry->ifInside(p.data());
      if (ifInside){
        return INSIDE;
      }else{
        return OUTSIDE;
      }
    }
    break;
    
    default:
      throw TALYFEMLIB::TALYException() << "Wrong type with in/out test for retain!";
    }
  }

  /**
   * @return false for outside region, true for inside region
   */
  bool in_region(ZEROPTV p){
    if (forRetain){
      throw TALYException() << "Calling function with forRetain = true";
    }
    switch (refine_type){
    case CUBE:{
      return !(p.x() < min_corner.x() ||
               p.y() < min_corner.y() ||
               p.z() < min_corner.z() ||
               p.x() > max_corner.x() ||
               p.y() > max_corner.y() ||
               p.z() > max_corner.z());
    }
    break;
    
    case SPHERE:{
      double distance = (p - center_sphere).norm();
      return (distance < radius_sphere and distance > radius_sphere_in);
    }
    break;
    
    case CYLINDER:{
      ZEROPTV AB = min_corner - max_corner;
      ZEROPTV temp = AB * (1.0 / AB.innerProduct(AB));
      ZEROPTV AP = p - max_corner;
      ZEROPTV proj_point = max_corner + temp * AP.innerProduct(AB);
      double distance_s = (proj_point - p).innerProduct(proj_point - p);
      if (distance_s < radius_cylinder * radius_cylinder and (distance_s - radius_cylinder_in * radius_cylinder_in) > -1e-6){
        double max_x = std::max(max_corner.x(), min_corner.x());
        double max_y = std::max(max_corner.y(), min_corner.y());
        double max_z = std::max(max_corner.z(), min_corner.z());
        double min_x = std::min(max_corner.x(), min_corner.x());
        double min_y = std::min(max_corner.y(), min_corner.y());
        double min_z = std::min(max_corner.z(), min_corner.z());
        return !(proj_point.x() > max_x ||
                 proj_point.y() > max_y ||
                 proj_point.z() > max_z ||
                 proj_point.x() < min_x ||
                 proj_point.y() < min_y ||
                 proj_point.z() < min_z);
      }else{
        return false;
      }
    }
    break;
    
    case MESHOBJECT:{
      bool ifInside = geometry->ifInside(p.data());
      return ifInside;
    }
    
    case MESHOBJECT_2D:{
      bool ifInside = geometry->ifInside(p.data());
      return ifInside;
    }
    break;
    
    default:
      throw TALYFEMLIB::TALYException() << "Wrong type!";
    }
  }

  Type GetRefineType() const{
    return refine_type;
  }

  void PrintRegionRefineDef(std::ofstream &fstream) const{
    int rank = TALYFEMLIB::GetMPIRank();
    if (!rank and fstream.is_open()){
      if (refine_type == Type::INVALID){
        fstream << "refine_type: INVALID\\n";
      }else if (refine_type == Type::CUBE){
        fstream << "refine_type: CUBE\\n";
        fstream << "min_corner: [";
        for (int i = 0; i < 3; i++){
          fstream << min_corner[i] << ", ";
        }
        fstream << "]\\n";
        fstream << "max_corner: [";
        for (int i = 0; i < 3; i++){
          fstream << max_corner[i] << ", ";
        }
        fstream << "]\\n";
      }else if (refine_type == Type::SPHERE){
        fstream << "refine_type: SPHERE\\n";
        fstream << "center_sphere: [";
        for (int i = 0; i < 3; i++){
          fstream << center_sphere[i] << ", ";
        }
        fstream << "]\\n";
        fstream << "radius_sphere: " << radius_sphere << "\\n";
        fstream << "radius_sphere_in: " << radius_sphere_in << "\\n";
      }else if (refine_type == Type::CYLINDER){
        fstream << "refine_type: CYLINDER\\n";
        fstream << "min_corner: [";
        for (int i = 0; i < 3; i++){
          fstream << min_corner[i] << ", ";
        }
        fstream << "]\\n";
        fstream << "max_corner: [";
        for (int i = 0; i < 3; i++){
          fstream << max_corner[i] << ", ";
        }
        fstream << "]\\n";

        fstream << "radius_cylinder: " << radius_cylinder << "\\n";
        fstream << "radius_cylinder_in: " << radius_cylinder_in << "\\n";
      }else if (refine_type == Type::MESHOBJECT){
        fstream << "refine_type: MESHOBJECT\\n";
      }else if (refine_type == Type::MESHOBJECT_2D){
        fstream << "refine_type: MESHOBJECT_2D\\n";
      }

      fstream << "shift: [";
      for (int i = 0; i < 3; i++){
        fstream << shift[i] << ", ";
      }
      fstream << "]\\n";

      fstream << "forRetain: " << forRetain << "\\n";
      fstream << "refine_region_lvl: " << refine_region_lvl << "\\n";
      fstream << "refine_region_lvl_boundary: " << refine_region_lvl_boundary << "\\n";
      fstream << "mesh_path: " << mesh_path << "\\n";
    }
  }

protected:
  Type refine_type = INVALID;
  // for sphere
  ZEROPTV center_sphere;
  double radius_sphere = 0.0;
  double radius_sphere_in = 0.0;
  // for cube
  ZEROPTV min_corner;
  ZEROPTV max_corner;
  // for cylinder
  //  ZEROPTV min_corner;
  //  ZEROPTV max_corner;
  double radius_cylinder = 0.0;
  double radius_cylinder_in = 0.0;
  static Type read_type(const std::string &str){
    if (str == "sphere"){
      return SPHERE;
    }else if (str == "cube"){
      return CUBE;
    }else if (str == "cylinder"){
      return CYLINDER;
    }else if (str == "meshobject"){
      return MESHOBJECT;
    }else if (str == "meshobject_2d"){
      return MESHOBJECT_2D;
    }else
      throw TALYFEMLIB::TALYException() << "Invalid Regional refinement type";
  }
};

/**
 * Parameters for Carved out geometry
 */
struct CarvedOutGeom{
  enum Type{
    INVALID = 0,
    SPHERE = 1,
    PILLAR = 2,
    MESHOBJECT = 3,
    MESHOBJECT_2D = 4,
    CYLINDER = 5,
    CUBE = 6,
    CIRCLE_2D_VOXEL = 7,
    SPHERE_3D_VOXEL = 8,
    BOX_2D_VOXEL = 9,
    CUBE_3D_VOXEL = 10,
  };
  enum FileFormat{
    NOTKNOWN = 0,
    GMSH = 1,
    STL = 2,
  };
  enum BCType{
    INVALID_BC = 0,
    DIRICHLET = 1,
    NEUMANN = 2,
    ROBIN = 3,
    WEAK = 4,
    SBM = 5,
    NEUMANN_SBM = 6
  };
  enum InoutTYPE{
    UNKNOWN = -1,
    PRIMITIVE = 0,
    RT = 1,
  };
  enum PrescribedMov{
    STATIONARY = 0,
    ROTATION = 1,
    TRANSALATION = 2,
  };
  ///< path to surface mesh (currently must be triangular)
  std::string mesh_path;
  ///< name of the geometry (used for printing file and etc..)
  std::string name;
  ///< level of gauss point split for triangles
  unsigned gp_level = 0;
  ///< file format depending on the file extension
  FileFormat fileformat = NOTKNOWN;
  ///< type of geometry, controls which analytical in/out function to use
  Type type = INVALID;
  ///< type for checking in_out
  InoutTYPE io_type = UNKNOWN;
  ///< refinement level for future use
  unsigned int refineLevel = 0;
  ///< initial displacement for the geometry (used only once to move geometry to desired location)
  ZEROPTV InitialDisplacement;
  ///< if postprocessing
  bool ifPostprocess = true;
  ///< retain inside or outside (default is false == retain outside)
  bool outer_boundary = false;

  /// variables depending on type...
  ZEROPTV center_of_mass;
  ///< radius for sphere or cylinder types
  DENDRITE_REAL radius = 0.0;
  ///< cylinder orientation, used for in_out test
  int cylinderDir = 2;
  double height = -1.0;
  ZEROPTV bottom_center;
  ///< first one is bottom left back corner, the second one is dimension in x-y-z
  std::vector<ZEROPTV> cube_dim{2};

  /// auto refine for IMGA loop
  GeomRefinement geomRefine;

  /// boundary conditions for geometry
  std::vector<CarvedOutGeom::BCType> bc_type_V = {};
  std::vector<unsigned int> bid_V = {};
  std::vector<DENDRITE_REAL> dirichlet_V = {};
  std::vector<DENDRITE_REAL> flux_V = {};
  std::vector<DENDRITE_REAL> G_constant_V = {};
  std::vector<DENDRITE_REAL> a_constant_V = {};

  bool is_static = true;
  std::vector<bool> translation_dof;
  bool if_buoyancy = false;
  bool if_rotation = false;
  DENDRITE_REAL acc_limiter = 1e6;
  DENDRITE_REAL rho = 1.0;

  ZEROPTV ic_vel;
  ZEROPTV ic_omega;
  DENDRITE_REAL until_time = -1.0;
  DENDRITE_REAL moving_until_time = -1.0;
  ZEROPTV bounding_min;
  ZEROPTV bounding_max;

  PrescribedMov move_type = PrescribedMov::STATIONARY;
  ZEROPTV rotation_axis = {1.0, 0.0, 0.0};
  DENDRITE_REAL rotation_speed = 0.0;
  ZEROPTV translation_speed = {0.0, 0.0, 0.0};

  std::vector<DENDRITE_REAL> getBC(int dof) const{
    if (bc_type_V[dof] == CarvedOutGeom::BCType::DIRICHLET or bc_type_V[dof] == CarvedOutGeom::BCType::WEAK or bc_type_V[dof] == CarvedOutGeom::BCType::SBM){
      return std::vector<DENDRITE_REAL>{dirichlet_V.at(dof)};
    }
    if (bc_type_V[dof] == CarvedOutGeom::BCType::NEUMANN or bc_type_V[dof] == CarvedOutGeom::BCType::NEUMANN_SBM){
      return std::vector<DENDRITE_REAL>{flux_V.at(dof)};
    }
    if (bc_type_V[dof] == CarvedOutGeom::BCType::ROBIN){
      return std::vector<DENDRITE_REAL>{G_constant_V.at(dof), a_constant_V.at(dof)};
    }
  }

  void read_from_config(const libconfig::Setting &root){
    type = str_to_type(root["type"]);
    if (type != Type::CIRCLE_2D_VOXEL and type != Type::SPHERE_3D_VOXEL and
        type != Type::BOX_2D_VOXEL and type != Type::CUBE_3D_VOXEL){
      mesh_path = (const char *)root["mesh_path"];
      if (root.exists("name")){
        name = (const char *)root["name"];
      }else{
        if (mesh_path.find_last_of('/')){
          name = mesh_path.substr(mesh_path.find_last_of('/') + 1,
                                  mesh_path.find_last_of('.') - mesh_path.find_last_of('/') - 1);
        }else{
          name = mesh_path.substr(0, mesh_path.find_last_of('.') - 1);
        }
      }
      fileformat = mesh_path_to_file_format(mesh_path);
    }else{
      name = (const char *)root["name"];
    }

    if (root.exists("gp_level")){
      gp_level = (int)root["gp_level"];
    }

    if (type != Type::CIRCLE_2D_VOXEL and type != Type::SPHERE_3D_VOXEL and
        type != Type::BOX_2D_VOXEL and type != Type::CUBE_3D_VOXEL){
      ReadZEROPTV(root, "position", InitialDisplacement);
    }
    
    io_type = type_to_iotype(type);
    if (root.exists("refineLevel")){
      refineLevel = (unsigned int)root["refineLevel"];
    }

    // for IMGA loop
    geomRefine.maxSplitIteration = 100;
    geomRefine.octantLevel = 3;
    geomRefine.ratioArea = 0.25;

    if (root.exists("ifPostprocess")){
      ifPostprocess = (bool)root["ifPostprocess"];
    }
    if (root.exists("outer_boundary")){
      outer_boundary = (bool)root["outer_boundary"];
    }
    
    assert(type == CarvedOutGeom::MESHOBJECT || type == CarvedOutGeom::MESHOBJECT_2D);

    const libconfig::Setting &bc_dof = root["bc_type_V"];
    for (int i = 0; i < bc_dof.getLength(); ++i){
      std::string a = bc_dof[i];
      CarvedOutGeom::BCType temp = str_to_bctype(bc_dof[i]);
      bc_type_V.push_back(temp);
    }
    /*
    int d_iter = 0, n_iter = 0, r_iter = 0;
    for (int i = 0; i < bc_type_V.size(); i++){
      if (bc_type_V.at(i) == CarvedOutGeom::BCType::DIRICHLET or bc_type_V.at(i) == CarvedOutGeom::BCType::WEAK or bc_type_V.at(i) == CarvedOutGeom::BCType::SBM){
        const libconfig::Setting &temp1 = root["dirichlet_V"];
        dirichlet_V.push_back(DENDRITE_REAL(temp1[d_iter++]));
        flux_V.push_back(-100.0);
        G_constant_V.push_back(-100.0);
        a_constant_V.push_back(-100.0);
      }
      if (bc_type_V.at(i) == CarvedOutGeom::BCType::NEUMANN or bc_type_V.at(i) == CarvedOutGeom::BCType::NEUMANN_SBM){
        dirichlet_V.push_back(-100.0);
        const libconfig::Setting &temp2 = root["flux_V"];
        flux_V.push_back(DENDRITE_REAL(temp2[n_iter++]));
        G_constant_V.push_back(-100.0);
        a_constant_V.push_back(-100.0);
      }
      if (bc_type_V.at(i) == CarvedOutGeom::BCType::ROBIN){
        dirichlet_V.push_back(-100.0);
        flux_V.push_back(-100.0);
        const libconfig::Setting &temp3 = root["G_constant_V"];
        const libconfig::Setting &temp4 = root["a_constant_V"];
        G_constant_V.push_back(DENDRITE_REAL(temp3[r_iter]));
        a_constant_V.push_back(DENDRITE_REAL(temp4[r_iter++]));
      }
    }
    */

  }

protected:
  void checkInput(){
    // TODO
  }

  /// read mesh type
  Type str_to_type(const std::string &str) const{
    if (str == "sphere") {
      return Type::SPHERE;
    }else if (str == "pillar"){
      return Type::PILLAR;
    }else if (str == "cylinder"){
      return Type::CYLINDER;
    }else if (str == "cube") {
      return Type::CUBE;
    }else if (str == "meshobject"){
      return Type::MESHOBJECT;
    }else if (str == "meshobject_2d"){
      return Type::MESHOBJECT_2D;
    }else if (str == "circle_2d"){
      return Type::CIRCLE_2D_VOXEL;
    }else if (str == "sphere_3d") {
      return Type::SPHERE_3D_VOXEL;
    }else if (str == "box_2d") {
      return Type::BOX_2D_VOXEL;
    }else if (str == "cube_3d"){
      return Type::CUBE_3D_VOXEL;
    }else{
      throw TALYFEMLIB::TALYException() << "Invalid geometry type '" << str
                                        << "' (expected sphere, pillar, cylinder cube or meshobject(2D) )"
                                        << "\\n for voxel, circle_2d or sphere_3d";
    }
  }

  FileFormat mesh_path_to_file_format(const std::string &str) const{
    if (str.substr(str.find_last_of('.') + 1) == "msh"){
      return FileFormat::GMSH;
    }else if (str.substr(str.find_last_of('.') + 1) == "stl"){
      return FileFormat::STL;
    } else{
      throw TALYFEMLIB::TALYException() << "Invalid file extension '" << str.substr(str.find_last_of('.') + 1)
                                        << "' (support .msh, .stl)";
    }
  }

  InoutTYPE type_to_iotype(const Type &type) const{
    if (type == Type::MESHOBJECT or type == Type::MESHOBJECT_2D){
      return RT;
    } else{
      return PRIMITIVE;
    }
  }

  /// read bc type
  BCType str_to_bctype(const std::string &str) const{
    if (str == "dirichlet"){
      return BCType::DIRICHLET;
    }else if (str == "neumann"){
      return BCType::NEUMANN;
    } else if (str == "robin"){
      return BCType::ROBIN;
    }else if (str == "weak"){
      return BCType::WEAK;
    } else if (str == "sbm") {
      return BCType::SBM;
    }else if (str == "neumann_sbm"){
        return  BCType::NEUMANN_SBM;
    } else{
      throw TALYFEMLIB::TALYException() << "Invalid bounday condition type '" << str
                                        << "' (expected dirichlet or neumann)";
    }
  }
};

/**
 * limiter
 */
struct Limiter{

  void read_from_config(const libconfig::Setting &root){
    if (root.exists("lb")){
      ReadVectorRoot(root, "lb", lb);
      if (lb.size() != $(project_name)NodeData::NUM_VARS){
        throw TALYException() << "lb has wrong size.";
      }
    }else{
      lb = std::vector<double>($(project_name)NodeData::NUM_VARS, -std::numeric_limits<PetscScalar>::infinity());
    }
    
    if (root.exists("ub")){
      ReadVectorRoot(root, "ub", ub);
      if (ub.size() != $(project_name)NodeData::NUM_VARS){
        throw TALYException() << "ub has wrong size.";
      }
    }else{
      ub = std::vector<double>($(project_name)NodeData::NUM_VARS, std::numeric_limits<PetscScalar>::infinity());
    }
    manual_limiter = true;
  }
  bool manual_limiter = false;
  std::vector<PetscScalar> lb;
  std::vector<PetscScalar> ub;
};
""";
    println(inputdatastructs_file, content);
    ##############################################################################################
    
    ##############################################################################################
    # util.h
    
    content = """
#pragma once
#include "$(project_name)NodeData.h"
#include "$(project_name)Refine.h"
#include "$(project_name)InputData.h"

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <Traversal/Analytic.h>
#include <PETSc/IO/petscVTU.h>
#include <SubDA/SubDomain.h>
#include <Boundary/SubDomainBoundary.h>
#include <Checkpoint/Checkpointer.h>
#include <talyfem/input_data/input_data.h>
#include <talyfem/talyfem.h>
#include <DataTypes.h>
#include <point.h>
#include <DendriteUtils.h>
#include <map>
#include <vector>
#include <string>
#include <PETSc/VecBounds.h>

namespace util_funcs{
  
  // This is only used by print2vtk
  void sortrows(std::vector<std::vector<double>>& matrix, int col) {
      std::sort(matrix.begin(), matrix.end(),
                [col](const std::vector<double>& lhs, const std::vector<double>& rhs) {
                    return lhs[col] > rhs[col];
                });
  }

  void print2vtk(const std::string &fPrefix, std::vector<ZEROPTV> PTV, const ZEROPTV &center){
      std::vector<ZEROPTV> PTVtemp;
      std::vector<ZEROPTV> PTVnew;

      if (TALYFEMLIB::GetMPIRank() == 0){ // if:master cpu, then:print
          ZEROPTV ptv1,ptv2;
          std::vector<std::vector<double>> AngleMap;

          PTVtemp.resize(PTV.size());
          PTVnew.resize(PTV.size());

          AngleMap.resize(PTV.size());
          for (int i =0; i<PTV.size();i+=2){

              AngleMap[i].resize(2);
              AngleMap[i+1].resize(2);

              for (int dim=0;dim<DIM;dim++){
                  if (PTV[i](dim) == PTV[i+1](dim)){
                      if (dim == 0){
                          ptv1 = {PTV[i](0),(PTV[i](1)+PTV[i+1](1))/2 + sqrt(3)*(PTV[i+1](1)-(PTV[i](1)+PTV[i+1](1))/2),0};
                          ptv2 = {PTV[i](0),(PTV[i](1)+PTV[i+1](1))/2 - sqrt(3)*(PTV[i+1](1)-(PTV[i](1)+PTV[i+1](1))/2),0};
                      } else if (dim == 1){
                          ptv1 = {(PTV[i](0)+PTV[i+1](0))/2 + sqrt(3)*(PTV[i+1](0)-(PTV[i](0)+PTV[i+1](0))/2),PTV[i](1),0};
                          ptv2 = {(PTV[i](0)+PTV[i+1](0))/2 - sqrt(3)*(PTV[i+1](0)-(PTV[i](0)+PTV[i+1](0))/2),PTV[i](1),0};
                      }
                      PTVtemp[i] = ptv1;
                      PTVtemp[i+1] = ptv2;
                      AngleMap[i][0] = i;
                      AngleMap[i][1] = atan2(PTVtemp[i].y()-center.y(),PTVtemp[i].x()-center.x());
                      AngleMap[i+1][0] = i+1;
                      AngleMap[i+1][1] = atan2(PTVtemp[i+1].y()-center.y(),PTVtemp[i+1].x()-center.x());
                  }
              }
          }

          sortrows(AngleMap,1);

          for (int i =0; i < PTVtemp.size();i++){
              PTVnew[i] = PTVtemp[AngleMap[i][0]];
          }

          char fname[256];
          snprintf(fname, sizeof(fname), "%s.vtk", fPrefix.c_str());
          std::ofstream fout(fname);

          int NGP = PTV.size();

          fout << "# vtk DataFile Version 2.0\\n";
          fout << "Created by Dendrite-KT\\n";
          fout << "ASCII\\n";
          fout << "DATASET UNSTRUCTURED_GRID\\n";
          fout << "POINTS " << NGP << " double\\n";

          for (int i = 0; i < NGP; i++) {
              fout << PTVnew[i].x() << " " << PTVnew[i].y() << " 0"
                    << "\\n";
          }
          fout << "\\n";

          fout << "CELLS " << NGP << " " << NGP * 3 << "\\n";
          for (int k = 0; k < NGP; k++){
              if (k + 1 >= NGP){
                  fout << "2 " << k << " " << k + 1 - NGP << "\\n";
              }else{
                  fout << "2 " << k << " " << k + 1 << "\\n";
              }
          }
          fout << "\\n";
          fout << "CELL_TYPES " << NGP << "\\n";
          int k = 0;
          for (int i = 0; i < NGP; i++){
              fout << "3\\n";
          }
          fout.close();
      }
  }

  DENDRITE_REAL ElementSize(const TALYFEMLIB::FEMElm &fe){
    return pow((pow(2, DIM) * fe.jacc()), (double)1 / DIM);
  }

  void TranslateGaussPoints(const std::vector<NodeAndValues<DENDRITE_REAL>> &gaussPoints, const std::vector<ZEROPTV> translations,
                            std::vector<ZEROPTV> &GPMove){
    GPMove.resize(gaussPoints.size());
    int counter = -1;
    for (const auto &gp : gaussPoints){
      counter++;
      for (int dim = 0; dim < DIM; dim++){
        GPMove[counter](dim) = gp.location[dim] + translations[counter](dim);
#if (DIM == 2)
        GPMove[counter](2) = 0;
#endif
      }
    }
  }

  void CalVec(double A, double B, double C, double x1, double y1, double &d1, double &d2){
    d1 = -A * (A * x1 + B * y1 + C) / (A * A + B * B);
    d2 = -B * (A * x1 + B * y1 + C) / (A * A + B * B);
  }

  void CalDist(double A, double B, double C, double x1, double y1, double &d){
    d = fabs(A * x1 + B * y1 + C) / sqrt(A * A + B * B);
  }

  static DENDRITE_REAL computeTriangleArea(const DENDRITE_REAL triCoords[][3]){
    /// calculate the normal and size of the surface
    TALYFEMLIB::ZEROPTV pts_0{triCoords[0][0], triCoords[0][1], triCoords[0][2]};
    TALYFEMLIB::ZEROPTV pts_1{triCoords[1][0], triCoords[1][1], triCoords[1][2]};
    TALYFEMLIB::ZEROPTV pts_2{triCoords[2][0], triCoords[2][1], triCoords[2][2]};
    TALYFEMLIB::ZEROPTV side_1, side_2;
    for (int dim = 0; dim < DIM; dim++){
      side_1(dim) = pts_1(dim) - pts_0(dim);
      side_2(dim) = pts_2(dim) - pts_0(dim);
    }

    /// this normal is either outside or inside, need to be pointed inside afterwards
    TALYFEMLIB::ZEROPTV normal;
    normal.crossProduct(side_1, side_2);
    const DENDRITE_REAL area = 0.5 * normal.SafeNormalize(); /// note: this line also normalizes normal!
    return area;
  }

  static DENDRITE_REAL computeSumTriangleArea(const ZEROPTV pt, const DENDRITE_REAL triCoords[][3]){
    /// calculate the normal and size of the surface
    std::vector<TALYFEMLIB::ZEROPTV> pts(3);
    pts[0] = {triCoords[0][0], triCoords[0][1], triCoords[0][2]};
    pts[1] = {triCoords[1][0], triCoords[1][1], triCoords[1][2]};
    pts[2] = {triCoords[2][0], triCoords[2][1], triCoords[2][2]};

    DENDRITE_REAL area = 0;

    for (int dim = 0; dim < DIM; dim++){
      pts[dim] = pt;
      TALYFEMLIB::ZEROPTV side_1, side_2;
      for (int dim = 0; dim < DIM; dim++){
        side_1(dim) = pts[1](dim) - pts[0](dim);
        side_2(dim) = pts[2](dim) - pts[0](dim);
      }

      /// this normal is either outside or inside, need to be pointed inside afterwards
      TALYFEMLIB::ZEROPTV normal;
      normal.crossProduct(side_1, side_2);
      area += 0.5 * normal.SafeNormalize(); /// note: this line also normalizes normal!
    }
    return area;
  }

  void PrintTriangleInfo(const std::vector<GEOMETRY::Triangles> m_tri){
    double minArea = 100, maxArea = -100;
    if (TALYFEMLIB::GetMPIRank() == 0){
      for (int i = 0; i < m_tri.size(); i++){
        double Area_i = computeTriangleArea(m_tri[i].triangleCoord);
        if (minArea > Area_i){
          minArea = Area_i;
        }
        if (maxArea < Area_i){
          maxArea = Area_i;
        }
      }
    }

    PrintStatus("min Triangle Area =", minArea);
    PrintStatus("max Triangle Area =", maxArea);
    PrintStatus("Triangles =", m_tri.size());
  }

  void printStl2File(const std::string &fPrefix, GEOMETRY::STL *stl){
    if (TALYFEMLIB::GetMPIRank() == 0){ // if:master cpu, then:print
      std::vector<GEOMETRY::Triangles> m_triangles = stl->getTriangles();
      char fname[256];
      sprintf(fname, "%s.csv", fPrefix.c_str());
      std::ofstream fout(fname);
      fout << "x0,y0,z0,nx0,ny0,nz0\\n";
      for (const auto &m_triangle : m_triangles){
        for (int dim = 0; dim < DIM; dim++){
          fout << m_triangle.triangleCoord[dim][0] << "," << m_triangle.triangleCoord[dim][1] << ","
               << m_triangle.triangleCoord[dim][2] << ","
               << m_triangle.normal[0] << "," << m_triangle.normal[1]
               << "," << m_triangle.normal[2] << " \\n";
        }
      }
      fout.close();
    }
  }

  void ChangingGaussPointsFormat(const std::vector<NodeAndValues<DENDRITE_REAL>> &gaussPoints,
                                 std::vector<ZEROPTV> &GPPTV){
    GPPTV.resize(gaussPoints.size());
    int counter = -1;
    for (const auto &gp : gaussPoints){
      ZEROPTV gpptv_local;
      counter++;
      for (int dim = 0; dim < DIM; dim++){
        GPPTV[counter](dim) = gp.location[dim];
#if (DIM == 2)
        GPPTV[counter](2) = 0;
#endif
      }
    }
  }

  void ChangingGaussPointsFormatWithNormal(const std::vector<NodeAndValues<DENDRITE_REAL>> &gaussPoints,
                                           std::vector<ZEROPTV> &GPPTV, std::vector<ZEROPTV> &NORMALPTV){
    GPPTV.resize(gaussPoints.size());
    NORMALPTV.resize(gaussPoints.size());
    int counter = -1;
    for (const auto &gp : gaussPoints){
      ZEROPTV gpptv_local;
      counter++;
      for (int dim = 0; dim < DIM; dim++){
        GPPTV[counter](dim) = gp.location[dim];
        NORMALPTV[counter](dim) = gp.normal[dim];
#if (DIM == 2)
        GPPTV[counter](2) = 0;
        NORMALPTV[counter](2) = 0;
#endif
      }
    }
  }

  void STL2zeroPTV(const std::vector<GEOMETRY::STL *> &stls,
                   std::vector<TALYFEMLIB::ZEROPTV> &points){
    std::vector<GEOMETRY::Triangles> m_triangles = stls[0]->getTriangles();
    int PtSize = m_triangles.size();
    points.resize(PtSize);
    for (int i = 0; i < PtSize; i++){
      for (int dim = 0; dim < DIM; dim++){
        points[i](dim) = (m_triangles[i].triangleCoord[0][dim] + m_triangles[i].triangleCoord[1][dim] + m_triangles[i].triangleCoord[2][dim]) / 3;
      }
    }
  }

  void all_together(const std::vector<TALYFEMLIB::ZEROPTV> pts_position, std::vector<TALYFEMLIB::ZEROPTV> &pts_position_all){
    int nProc = TALYFEMLIB::GetMPISize(); // be careful
    int numNodes = pts_position.size();
    std::vector<int> eachProcData(nProc);
    MPI_Allgather(&numNodes, 1, MPI_INT, eachProcData.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> disp(nProc, 0);
    for (int i = 1; i < disp.size(); i++){
      disp[i] = disp[i - 1] + eachProcData[i - 1];
    }

    int totalProcData = 0;
    for (int i = 0; i < nProc; i++){
      totalProcData += eachProcData[i];
    }

    if (TALYFEMLIB::GetMPIRank() == 0){
      pts_position_all.resize(totalProcData);
    }

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Gatherv(pts_position.data(), pts_position.size(), ZEROPTVtype, pts_position_all.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);
  }

  void ToEachProcessor(const std::vector<TALYFEMLIB::ZEROPTV> pts_position, std::vector<TALYFEMLIB::ZEROPTV> &pts_position_all){
    int nProc = TALYFEMLIB::GetMPISize(); // be careful
    int numNodes = pts_position.size();
    std::vector<int> eachProcData(nProc);
    MPI_Allgather(&numNodes, 1, MPI_INT, eachProcData.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> disp(nProc, 0);
    for (int i = 1; i < disp.size(); i++){
      disp[i] = disp[i - 1] + eachProcData[i - 1];
    }

    int totalProcData = 0;
    for (int i = 0; i < nProc; i++){
      totalProcData += eachProcData[i];
    }

    pts_position_all.resize(totalProcData);

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Allgatherv(pts_position.data(), pts_position.size(), ZEROPTVtype, pts_position_all.data(), eachProcData.data(), disp.data(), ZEROPTVtype, MPI_COMM_WORLD);
  }

  void print_Point(const std::string &fPrefix, const std::vector<TALYFEMLIB::ZEROPTV> &pts_position_all){
    if (TALYFEMLIB::GetMPIRank() == 0){ // if:master cpu, then:print
      FILE *fp = fopen(fPrefix.c_str(), "w");

#if (DIM == 2)
      fprintf(fp, "x0,y0\\n");
#endif
#if (DIM == 3)
      fprintf(fp, "x0,y0,z0\\n");
#endif

      for (int i = 0; i < pts_position_all.size(); i++){
#if (DIM == 2)
        fprintf(fp, "%.10e,%.10e\\n",
                pts_position_all[i](0), pts_position_all[i](1));
#endif
#if (DIM == 3)
        fprintf(fp, "%.10e,%.10e,%.10e\\n",
                pts_position_all[i](0), pts_position_all[i](1), pts_position_all[i](2));
#endif
      }

      fclose(fp);
    }
  }

  void print_PointandNormal(const std::string &fPrefix, const std::vector<TALYFEMLIB::ZEROPTV> &pts_position_all, const std::vector<TALYFEMLIB::ZEROPTV> &pts_normal_all){
    if (TALYFEMLIB::GetMPIRank() == 0){ // if:master cpu, then:print
      FILE *fp = fopen(fPrefix.c_str(), "w");

#if (DIM == 2)
      fprintf(fp, "x0,y0,nx0,ny0\\n");
#endif
#if (DIM == 3)
      fprintf(fp, "x0,y0,z0,nx0,ny0,nz0\\n");
#endif

      for (int i = 0; i < pts_position_all.size(); i++){
#if (DIM == 2)
        fprintf(fp, "%.10e,%.10e,%.10e,%.10e\\n",
                pts_position_all[i](0), pts_position_all[i](1), pts_normal_all[i](0), pts_normal_all[i](1));
#endif
#if (DIM == 3)
        fprintf(fp, "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\\n",
                pts_position_all[i](0), pts_position_all[i](1), pts_position_all[i](2), pts_normal_all[i](0), pts_normal_all[i](1), pts_normal_all[i](2));
#endif
      }

      fclose(fp);
    }
  }

    void Local2Global(const int &local, int &global, MPI_Comm comm){
        MPI_Reduce(&local, &global, 1, MPI_INT, MPI_SUM, 0, comm);
    }

    void Local2Global(const double &local, double &global, MPI_Comm comm){
        MPI_Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    }

    void PrintPT2File(const std::vector<ZEROPTV> &PT,const std::string &fPrefix){
        std::vector<ZEROPTV> PTGlobal;

        int nProc = TALYFEMLIB::GetMPISize(); // be careful
        int numNodes = PT.size();
        std::vector<int> eachProcData(nProc);
        MPI_Allgather(&numNodes, 1, MPI_INT, eachProcData.data(), 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<int> disp(nProc, 0);
        for (int i = 1; i < disp.size(); i++){
            disp[i] = disp[i - 1] + eachProcData[i - 1];
        }

        int totalProcData = 0;
        for (int i = 0; i < nProc; i++){
            totalProcData += eachProcData[i];
        }

        if (TALYFEMLIB::GetMPIRank() == 0){
            PTGlobal.resize(totalProcData);
        }

        MPI_Datatype ZEROPTVtype;
        MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
        MPI_Type_commit(&ZEROPTVtype);
        MPI_Gatherv(PT.data(), PT.size(), ZEROPTVtype, PTGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);

        if (TALYFEMLIB::GetMPIRank() == 0){ // if:master cpu, then:print
            FILE *fp = fopen(fPrefix.c_str(), "w");

#if (DIM == 2)
            fprintf(fp, "x,y\\n");
#endif
#if (DIM == 3)
            fprintf(fp, "x,y,z\\n");
#endif

            for (int i = 0; i < PTGlobal.size(); i++){
#if (DIM == 2)
                fprintf(fp, "%.10e,%.10e\\n",
                        PTGlobal[i](0), PTGlobal[i](1));
#endif
#if (DIM == 3)
                fprintf(fp, "%.10e,%.10e,%.10e\\n",
                        PTGlobal[i](0), PTGlobal[i](1),PTGlobal[i](2));
#endif
            }

            fclose(fp);
        }
    }

    /**
   * Calculates the L2 error of two vectors, used in block iteration.
   * @param[in] vec_u first of the two vectors for comparison
   * @param[in] vec_block second of the two vectors for comparison
   * @param[in] ndof ndof of the vector
   * @return
   */
  static double calc_error(Vec vec_u, Vec vec_block, const unsigned int ndof, MPI_Comm comm){
    PetscInt size;
    VecGetLocalSize(vec_u, &size);
    PetscInt low;
    VecGetOwnershipRange(vec_u, &low, nullptr);

    const PetscScalar *u;
    const PetscScalar *block;
    VecGetArrayRead(vec_u, &u);
    VecGetArrayRead(vec_block, &block);

    std::vector<double> totalU(ndof, 0.0);   // holds sum of u^2 at every node
    std::vector<double> totalErr(ndof, 0.0); // holds sum of (u - u_block)^2 at every node
    for (PetscInt i = 0; i < size; i++){
      const unsigned int dof = (low + i) % ndof;
      totalU[dof] += u[i] * u[i];
      double diff = u[i] - block[i];
      totalErr[dof] += diff * diff;
    }

    VecRestoreArrayRead(vec_u, &u);
    VecRestoreArrayRead(vec_block, &block);

    std::vector<double> allU(ndof);
    std::vector<double> allErrs(ndof);
    MPI_Allreduce(totalU.data(), allU.data(), ndof, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(totalErr.data(), allErrs.data(), ndof, MPI_DOUBLE, MPI_SUM, comm);
    totalU.clear();
    totalErr.clear();
    double largest = -1e16;
    for (unsigned int dof = 0; dof < ndof; dof++){
      double err = allErrs[dof] / (allU[dof] + 1e-12);
      if (std::isnan(err) || std::isinf(err)){
        throw TALYFEMLIB::TALYException() << "Infinity or NAN in block error";
      }

      if (err > largest){
        largest = err;
      }
    }
    return sqrt(largest);
  }

  /**
   * Save the octree mesh to vtk binary file.
   */
  PetscErrorCode save_timestep(DA *octDA,
                               const std::vector<TREENODE> &treePartition,
                               Vec vec,
                               unsigned int ndof,
                               const TimeInfo &ti,
                               const SubDomain &subDomain,
                               const std::string &prefix = "timestep",
                               const char **varName = nullptr /*,
                                bool saveGeometry = true,
                                bool writeData = true,
                                const std::vector<int> &dataIndex = {},
                                const std::function<double(double v_in, bool in_geom, unsigned int dof)> &f = nullptr*/
  ){
    // create directory for this timestep (if it doesnt already exist)
    char folder[PATH_MAX];
    snprintf(folder, sizeof(folder), "results_%05d", ti.getTimeStepNumber());
    int ierr = mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (ierr != 0 && errno != EEXIST){
      PrintError("Could not create folder for storing results (", strerror(errno), ").");
      return 1;
    }

    char fname[PATH_MAX];
    snprintf(fname, sizeof(fname), "%s_%05d", prefix.c_str(), ti.getTimeStepNumber());
    PETSc::petscVectopvtu(octDA, treePartition, vec, folder, fname, varName,
                          subDomain.domainExtents(), false, false, ndof);
  }

  void save_data(DA *octDA,
                 const std::vector<TREENODE> &treePartition,
                 Vec solution,
                 const $(project_name)InputData &idata,
                 const TimeInfo &ti,
                 const SubDomain &subDomain,
                 const char **varname){
    
    save_timestep(octDA, treePartition, solution, $(project_name)NodeData::NUM_VARS, ti, subDomain, "$(project_name)", varname);
  }

  void performRefinementSubDA(DA *&octDA, const std::vector<TREENODE> &treeNode, DomainExtents &domainExtents, DistTREE &dTree,
                              $(project_name)InputData &inputData, SubDomain *subdomain){
    int no_refine = 0;
    while (true){
      SubDomainBoundary subDomainBoundary(subdomain, octDA, domainExtents);

      $(project_name)Refine refine(octDA, treeNode, domainExtents, &inputData, &subDomainBoundary);
      DA *newDA = refine.getRefineSubDA(dTree);
      if (newDA == NULL){
        newDA = refine.getForceRefineSubDA(dTree);
        std::swap(newDA, octDA);
        break;
      }

      std::swap(newDA, octDA);

      delete newDA;

      subdomain->finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
      TALYFEMLIB::PrintStatus("Refinement ", no_refine, ", mesh count (node) = ", octDA->getGlobalNodeSz());
      no_refine++;
    }
  }

  void performRefinementSubDAIBM(DA *&octDA, const std::vector<TREENODE> &treeNode, DomainExtents &domainExtents, DistTREE &dTree,
                                 $(project_name)InputData &inputData, SubDomain *subdomain, std::vector<GEOMETRY::Geometry *> ibm_geoms, bool doCoarsen = false){
    int no_refine = 0;
    while (no_refine < 36){
      SubDomainBoundary subDomainBoundary(subdomain, octDA, domainExtents);
      $(project_name)Refine refine(octDA, treeNode, domainExtents, &inputData, &subDomainBoundary, ibm_geoms, doCoarsen);
      DA *newDA = refine.getRefineSubDA(dTree);
      if (newDA == NULL){
        break;
      }

      std::swap(newDA, octDA);

      delete newDA;
      subdomain->finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
      PrintStatus("Refinement ", no_refine, ", mesh count (node) = ", octDA->getGlobalNodeSz());
      no_refine++;
    }
  }

  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // These will be removed.
  // These calculations will be done in main or aren't used yet.
  
  // void refineAndIntergrid(DA *&octDA,
  //                         DistTREE &dTree,
  //                         DomainExtents &domainExtents,
  //                         $(project_name)InputData &inputData,
  //                         SubDomainBoundary *boundary,
  //                         Vec &vecHT1, Vec &vecHT2, Vec &vecHT3, Vec &vecHT4,
  //                         const char **ht_varname,
  //                         TimerGroup<MPITimer> &timers,
  //                         std::map<std::string, int> &timer_tags){
  //   int no_refine = 0;
  //   while (true){
  //     TALYFEMLIB::PrintStatus("Active comm = ",
  //                             octDA->getNpesActive(),
  //                             ", no_refine = ",
  //                             no_refine,
  //                             ", mesh count (node) = ",
  //                             octDA->getGlobalNodeSz());
                              
  //     $(project_name)Refine daRefine(octDA, dTree.getTreePartFiltered(), domainExtents, &inputData, boundary);
      
  //     DA *newDA = daRefine.getRefineSubDA(dTree);
  //     if (newDA == NULL){
  //       break;
  //     }
  //     timers.Start(timer_tags["initIntergridTransfer"]);
  //     daRefine.initPetscIntergridTransfer();
  //     timers.Stop(timer_tags["initIntergridTransfer"]);
  //     timers.Start(timer_tags["IntergridTransfer"]);
  //     daRefine.petscIntergridTransfer(newDA, dTree, vecHT1, $(project_name)NodeData::NUM_VARS);
  //     daRefine.petscIntergridTransfer(newDA, dTree, vecHT2, $(project_name)NodeData::NUM_VARS);
  //     daRefine.petscIntergridTransfer(newDA, dTree, vecHT3, $(project_name)NodeData::NUM_VARS);
  //     daRefine.petscIntergridTransfer(newDA, dTree, vecHT4, $(project_name)NodeData::NUM_VARS);
  //     timers.Stop(timer_tags["IntergridTransfer"]);
  //     daRefine.finializeIntergridTransfer();
  //     std::swap(octDA, newDA);
  //     delete newDA;
  //     no_refine++;
  //   }
  // }
  
  // /**
  //  * @brief sets vector with initial temperature. Must be allocated outside.
  //  * @param octDA
  //  * @param inputData
  //  * @param [out] Solution the final vector with temperature
  //  */
  // void setInitialConditionHT(DA *octDA, $(project_name)InputData &inputData, Vec &Solution)
  // {
  //   std::function<void(const double *, double *)> initial_condition = [&](const double *x, double *var)
  //   {
  //     /// HT only has 1 dof and should be accessed as var[0]
  //     if (inputData.dump_vec &&
  //         inputData.boundary_def[1].temperature_type == BoundaryDef::Temperature_Type::ROBIN_T)
  //     {
  //       var[0] = 100 * (1.0 - x[0] / 3);
  //     }
  //     else
  //     {
  //       var[0] = inputData.ic_def.t_ic;
  //     }
  //   };

  //   octDA->petscSetVectorByFunction(Solution, initial_condition, false, false, $(project_name)NodeData::NUM_VARS);
  // }

  // void calculateNuCarvedOut(DA *octDA, const std::vector<TREENODE> &treePart, Vec htSolution,
  //                           const DomainExtents &domainExtents, SubDomainBoundary *subDomainBoundary,
  //                           $(project_name)InputData &inputData, const TimeInfo ti)
  // {
  //   for (DENDRITE_UINT no_object = 0; no_object < inputData.carved_out_geoms_def.size(); no_object++)
  //   {
  //     DENDRITE_REAL globalNu[2];
  //     VecInfo vec(htSolution, $(project_name)NodeData::NUM_VARS, $(project_name)NodeData::TEMPERATURE);
  //     PostProcessCalc nuCalc(octDA, treePart, vec, domainExtents,
  //                            inputData.getDiffusionHT(ti.getCurrentTime()),
  //                            subDomainBoundary, no_object, PostProcessCalc::NU);
  //     nuCalc.getTotalNu(globalNu);
  //     if (TALYFEMLIB::GetMPIRank() == 0)
  //     {
  //       std::string fname = "Nu_" + std::to_string(no_object) + ".dat";
  //       //      if (surface2cal_) {
  //       //        fname = "Force_" + to_string(surface2cal_) + ".dat";
  //       //      }
  //       FILE *fp = fopen(fname.c_str(), "a");
  //       if (!fp)
  //       {
  //         throw TALYException() << "Cannot create file: " << fname;
  //       }
  //       fprintf(fp,
  //               "Timestep: %1d, Time: %.5f\\n"
  //               "Nu_normal = %.10e, Nu_penalty = %.10e\\n",
  //               ti.getTimeStepNumber(),
  //               ti.getCurrentTime(),
  //               globalNu[0], globalNu[1]);

  //       fclose(fp);
  //       //      std::ofstream fout("Nu_" + std::to_string(no_object) + ".txt", std::ios::app);
  //       //      fout << ti.getCurrentTime() << " ";
  //       //      fout << globalNu[0] << " " << globalNu[1] << " ";
  //       //      fout << "\\n";
  //       //      fout.close();
  //     }
  //   }
  // }

  // void printMinMaxValueHT(DA *octDA, Vec &sol)
  // {

  //   PetscScalar max[$(project_name)NodeData::NUM_VARS], min[$(project_name)NodeData::NUM_VARS];
  //   VecBounds<$(project_name)NodeData::NUM_VARS>::getMaxAndMinimumValues(sol, octDA->getCommActive(), max, min);
  //   for (int i = 0; i < $(project_name)NodeData::NUM_VARS; i++)
  //   {
  //     PrintStatus("dof:", i, " HT min/max: ", min[i], "/", max[i]);
  //   }
  // }

  // void setMinMaxValueHT(DA *octDA, Vec &sol, $(project_name)InputData &idata)
  // {
  //   if (idata.limiter.manual_limiter)
  //   {
  //     PetscScalar max[$(project_name)NodeData::NUM_VARS], min[$(project_name)NodeData::NUM_VARS];
  //     for (int i = 0; i < $(project_name)NodeData::NUM_VARS; i++)
  //     {
  //       min[i] = idata.limiter.lb[$(project_name)NodeData::TEMPERATURE + i];
  //       max[i] = idata.limiter.ub[$(project_name)NodeData::TEMPERATURE + i];
  //       PrintStatus("dof:", i, " Set HT min/max: ", min[i], "/", max[i]);
  //     }
  //     VecBounds<$(project_name)NodeData::NUM_VARS>::setMaxAndMinimumValues(sol, octDA->getCommActive(), max, min);
  //   }
  // }

  // void read_running_sum_step(const std::string &filename, int &running_step){
  //   if (TALYFEMLIB::GetMPIRank() == 0){
  //     std::fstream myfile(filename.c_str(), std::ios_base::in);
  //     if (myfile){
  //       myfile >> running_step;
  //       myfile.close();
  //     }
  //   }
  // }

  // void write_running_sum_step(const std::string &filename, const int &running_step){
  //   if (TALYFEMLIB::GetMPIRank() == 0){
  //     std::fstream myfile(filename.c_str(), std::ios_base::out);
  //     myfile << running_step;
  //   }
  // }

}


""";
    println(utils_file, content);
    ##############################################################################################
    
    ##############################################################################################
    # SBMcalc.h
    
    content = """
#pragma once

#include "util.h"
#include <IMGA/IMGA.h>

class SBMCalc{
private:
  double d[DIM];
  double DirichletBCValue;
  TALYFEMLIB::FEMElm fe_;
  const IMGA *imga_;
  const $(project_name)InputData *idata_;
  ZEROPTV shift_;

  ZEROPTV normal;
  const double min_domain = 0.0;
  const double max_domain = 1.0 - min_domain;
  const double min_domain_RotBox = -0.5;
  const double max_domain_RotBox = 0.5;

  /**
    * @brief function to check whether the shited point on true boundary is within the corresponding triangle
    * @param pt the carved-out based gauss point position
    * @param d the distance function between true boundary and surrogate boundary
    * @param m_triangle the corresponding triangle to check whether the point is inside
    * @param shift the initial displacement of the stl geometry
    * @return [bool] whether the point is on the plane of the triangle
    */
  bool CheckInside3DTriangle(const ZEROPTV &pt, const double (&d)[DIM], const GEOMETRY::Triangles &m_triangle, const ZEROPTV &shift);

  /**
    * @brief calculate the shortest distance from octree-based Gauss Points to triangle edges
    * @param pt the carved-out based gauss point position
    * @param m_triangle the corresponding triangle to find the distance
    * @param shift the initial displacement of the stl geometry
    * @param d [out] the distance function from octree-based Gauss Points to triangle edges
    */
  void ShortestDist2TriEdge(const ZEROPTV &pt, const GEOMETRY::Triangles &m_triangle, const ZEROPTV &shift, double (&d)[DIM]);

public:
  /**
    * @brief constructor
    * @param fe the element we use to access gauss points
    * @param idata input Data
    * @param imga imga context, and we use this to access geometry
    */
  SBMCalc(const TALYFEMLIB::FEMElm &fe, const $(project_name)InputData *idata, const IMGA *imga);
  
  /**
    * @brief calculate the distance function for different kinds of geometries
    * @param d [out] distance function
    */
  void Dist2Geo(double (&d)[DIM]);

  /**
    * @brief calculate the true normal of the SBM geometry
    * @param normal true normal of SBM geometry
    * @param d distance function
    */
  void NormalofGeo(ZEROPTV &normal, const double (&d)[DIM]);

};

SBMCalc::SBMCalc(const TALYFEMLIB::FEMElm &fe, const $(project_name)InputData *idata, const IMGA *imga)
    : fe_(fe), imga_(imga), shift_(idata->carved_out_geoms_def[0].InitialDisplacement){
  idata_ = idata;
}

bool SBMCalc::CheckInside3DTriangle(const ZEROPTV &pt, const double (&d)[DIM], const GEOMETRY::Triangles &m_triangle, const ZEROPTV &shift){
  const ZEROPTV ptMove{pt(0) + d[0] - shift(0), pt(1) + d[1] - shift(1), pt(2) + d[2] - shift(2)};

  const ZEROPTV ptA{m_triangle.triangleCoord[0][0], m_triangle.triangleCoord[0][1], m_triangle.triangleCoord[0][2]};
  const ZEROPTV ptB{m_triangle.triangleCoord[1][0], m_triangle.triangleCoord[1][1], m_triangle.triangleCoord[1][2]};
  const ZEROPTV ptC{m_triangle.triangleCoord[2][0], m_triangle.triangleCoord[2][1], m_triangle.triangleCoord[2][2]};

  ZEROPTV u;
  ZEROPTV v;
  ZEROPTV w;

  u.crossProduct(ptB - ptA, ptMove - ptA);
  v.crossProduct(ptC - ptB, ptMove - ptB);
  w.crossProduct(ptA - ptC, ptMove - ptC);

  if (u.innerProduct(v) < 0.0f){
    return false;
  }else if (u.innerProduct(w) < 0.0f){
    return false;
  }else {
    return true;
  }
}

void SBMCalc::ShortestDist2TriEdge(const ZEROPTV &pt, const GEOMETRY::Triangles &m_triangle, const ZEROPTV &shift, double (&d)[DIM]){
  std::vector<ZEROPTV> ShortestVector2Line(3);
  std::vector<ZEROPTV> ptri(3);
  for (int trinum = 0; trinum < DIM; trinum++){
    ptri[trinum] = {m_triangle.triangleCoord[trinum][0] + shift[0], m_triangle.triangleCoord[trinum][1] + shift[1], m_triangle.triangleCoord[trinum][2] + shift[2]};
  }

  ShortestVector2Line[0] = (ptri[1] - ptri[0]) * (((pt - ptri[0]).innerProduct(ptri[1] - ptri[0])) / (ptri[1] - ptri[0]).norm()) - (pt - ptri[0]);
  ShortestVector2Line[1] = (ptri[2] - ptri[1]) * (((pt - ptri[1]).innerProduct(ptri[2] - ptri[1])) / (ptri[2] - ptri[1]).norm()) - (pt - ptri[1]);
  ShortestVector2Line[2] = (ptri[0] - ptri[2]) * (((pt - ptri[2]).innerProduct(ptri[0] - ptri[2])) / (ptri[0] - ptri[2]).norm()) - (pt - ptri[2]);

  int pickNumber = 0;
  double mindist = 100;
  for (int trinum = 0; trinum < DIM; trinum++){
    if (ShortestVector2Line[trinum].norm() < mindist){
      pickNumber = trinum;
      mindist = ShortestVector2Line[trinum].norm();
    }
  }

  for (int trinum = 0; trinum < DIM; trinum++){
    d[trinum] = ShortestVector2Line[pickNumber](trinum);
  }
}

void SBMCalc::Dist2Geo(double (&d)[DIM]){
  const ZEROPTV pt = fe_.position();
  double x = pt.x();
  double y = pt.y();
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set these if specifying the distance to geometry below
  bool user_specified_distance = true;
  bool general_mesh_distance = false;
  //////////////////////////////////////////////////////////////////////////////////////////////////

#if (DIM == 2)

  if(user_specified_distance){
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Put the mesh-specific distance to geometry code here.
    // This is a circle. This code needs to be specified by the user because it is specific to the mesh
    double radius = 1.0;
    double radius_gp = sqrt((x - shift_.x()) * (x - shift_.x()) + (y - shift_.y()) * (y - shift_.y()));
    d[0] = (radius * (x - shift_.x()) / radius_gp + shift_.x()) - x;
    d[1] = (radius * (y - shift_.y()) / radius_gp + shift_.y()) - y;
    //////////////////////////////////////////////////////////////////////////////////////////////
    
  }else if(general_mesh_distance){
    
    // This will attempt to work for an unknown mesh
    auto m_lines = imga_->getGeometries()[0]->getMSH()->getLines();
    int msh_size = m_lines.size();

    double MinDist = sqrt(2) * util_funcs::ElementSize(fe_);
    for (DENDRITE_UINT i = 0; i < msh_size; i++)
    {
      if (sqrt(pow(x - m_lines[i].lineCoord[0][0] - shift_[0],2) + pow(y - m_lines[i].lineCoord[0][1] - shift_[1],2)) < MinDist)
      {
        MinDist = sqrt(pow(x - m_lines[i].lineCoord[0][0] - shift_[0],2) + pow(y - m_lines[i].lineCoord[0][1] - shift_[1],2));
        d[0] = m_lines[i].lineCoord[0][0] + shift_[0] - x;
        d[1] = m_lines[i].lineCoord[0][1] + shift_[1] - y;
      }
    }
    
  }else{
    // returning zero distance
    d[0] = 0.0;
    d[1] = 0.0;
  }

#endif
#if (DIM == 3)

  double z = pt.z();

  if(user_specified_distance){
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Put the mesh-specific distance to geometry code here.
    // This is a sphere. This code needs to be specified by the user because it is specific to the mesh
    double radius = shift_[0]; // The center is at (radius,radius,radius)
    double radius_gp = sqrt((x - radius) * (x - radius) + (y - radius) * (y - radius) + (z - radius) * (z - radius));
    d[0] = (radius * (x - radius) / radius_gp + radius) - x;
    d[1] = (radius * (y - radius) / radius_gp + radius) - y;
    d[2] = (radius * (z - radius) / radius_gp + radius) - z;
    //////////////////////////////////////////////////////////////////////////////////////////////
    
  }else if(general_mesh_distance){
    // This will attempt to work for an unknown mesh
    double MinDist = 100.0;
    ZEROPTV OnePointVector;
    ZEROPTV PickNormalVector;
    int PickGeomID = 0;
    int PickTrianleID = 0;

    switch (idata_->DistCalcType)
    {

    case $(project_name)InputData::typeDistCalc::GP_BASED:
    {
      for (const auto gp : idata_->GPPTVAll)
      {
        if (sqrt(pow(x - gp(0), 2) + pow(y - gp(1), 2) + pow(z - gp(2), 2)) < MinDist)
        {
          MinDist = sqrt(pow(x - gp(0), 2) + pow(y - gp(1), 2) + pow(z - gp(2), 2));
          d[0] = gp(0) - x;
          d[1] = gp(1) - y;
          d[2] = gp(2) - z;
        }
      }
      break;
    }

    case $(project_name)InputData::typeDistCalc::NORMAL_BASED_DistributeSTL:
    {
      for (int geoID = 0; geoID < imga_->getGeometries().size(); geoID++)
      {
        std::vector<GEOMETRY::Triangles> m_triangles = imga_->getGeometries()[geoID]->getSTL()[0].getTriangles();

        assert(idata_->PBoxEnd(0) > x and idata_->PBoxEnd(1) > y and idata_->PBoxEnd(2) > z and idata_->PBoxStart(0) < x and idata_->PBoxStart(1) < y and idata_->PBoxStart(2) < z);
        int it = 0;
        for (auto point : idata_->DistributePoints[geoID])
        {
          if (sqrt(pow(x - point.x() - shift_[0], 2) + pow(y - point.y() - shift_[1], 2) + pow(z - point.z() - shift_[2], 2)) < MinDist)
          {
            MinDist = sqrt(pow(x - point.x() - shift_[0], 2) + pow(y - point.y() - shift_[1], 2) + pow(z - point.z() - shift_[2], 2));
            for (int dim = 0; dim < DIM; dim++)
            {
              OnePointVector(dim) = m_triangles[idata_->TriangleNumber[geoID][it]].triangleCoord[0][dim] + shift_[dim] - pt(dim);
              PickNormalVector(dim) = m_triangles[idata_->TriangleNumber[geoID][it]].normal[dim];
              PickTrianleID = idata_->TriangleNumber[geoID][it];
              PickGeomID = geoID;
            }
          }
          it++;
        }
      }

      // scaling of vector
      double scale = 0.0;
      for (int dim = 0; dim < DIM; dim++)
      {
        scale += OnePointVector(dim) * PickNormalVector(dim);
      }

      for (int dim = 0; dim < DIM; dim++)
      {
        d[dim] = scale * PickNormalVector(dim);
      }

      GEOMETRY::Triangles m_triangle = imga_->getGeometries()[PickGeomID]->getSTL()[0].getTriangles()[PickTrianleID];
      if (!CheckInside3DTriangle(pt, d, m_triangle, shift_))
      {
        ShortestDist2TriEdge(pt, m_triangle, shift_, d);
      }
      break;
    }

    case $(project_name)InputData::typeDistCalc::NORMAL_BASED:
    {
      for (int geoID = 0; geoID < imga_->getGeometries().size(); geoID++)
      {
        std::vector<GEOMETRY::Triangles> m_triangles = imga_->getGeometries()[geoID]->getSTL()[0].getTriangles();
        //std::cout<<"m_triangles.size() = " << m_triangles.size() << "\n";

        for (int i = 0; i < m_triangles.size(); i++)
        {
          if (sqrt(pow(x - (m_triangles[i].triangleCoord[0][0] + m_triangles[i].triangleCoord[1][0] + m_triangles[i].triangleCoord[2][0]) / 3 - shift_[0], 2) + pow(y - (m_triangles[i].triangleCoord[0][1] + m_triangles[i].triangleCoord[1][1] + m_triangles[i].triangleCoord[2][1]) / 3 - shift_[1], 2) + pow(z - (m_triangles[i].triangleCoord[0][2] + m_triangles[i].triangleCoord[1][2] + m_triangles[i].triangleCoord[2][2]) / 3 - shift_[2], 2)) < MinDist)
          {
            MinDist = sqrt(pow(x - (m_triangles[i].triangleCoord[0][0] + m_triangles[i].triangleCoord[1][0] + m_triangles[i].triangleCoord[2][0]) / 3 - shift_[0], 2) + pow(y - (m_triangles[i].triangleCoord[0][1] + m_triangles[i].triangleCoord[1][1] + m_triangles[i].triangleCoord[2][1]) / 3 - shift_[1], 2) + pow(z - (m_triangles[i].triangleCoord[0][2] + m_triangles[i].triangleCoord[1][2] + m_triangles[i].triangleCoord[2][2]) / 3 - shift_[2], 2));

            for (int dim = 0; dim < DIM; dim++)
            {
              OnePointVector(dim) = m_triangles[i].triangleCoord[0][dim] + shift_[dim] - pt(dim);
              PickNormalVector(dim) = m_triangles[i].normal[dim];
              PickTrianleID = i;
              PickGeomID = geoID;
            }
          }
        }
      }

      // scaling of vector
      double scale = 0.0;
      for (int dim = 0; dim < DIM; dim++)
      {
        scale += OnePointVector(dim) * PickNormalVector(dim);
      }

      for (int dim = 0; dim < DIM; dim++)
      {
        d[dim] = scale * PickNormalVector(dim);
      }

      GEOMETRY::Triangles m_triangle = imga_->getGeometries()[PickGeomID]->getSTL()[0].getTriangles()[PickTrianleID];
      if (!CheckInside3DTriangle(pt, d, m_triangle, shift_))
      {
        ShortestDist2TriEdge(pt, m_triangle, shift_, d);
      }
      break;
    }
    }
    
  }else{
    // returning zero distance
    d[0] = 0.0;
    d[1] = 0.0;
    d[2] = 0.0;
  }
#endif
}

void SBMCalc::NormalofGeo(ZEROPTV &normal, const double (&d)[DIM]){
  double R2 = 0;

  for (int dim = 0; dim < DIM; dim++){
    R2 += pow(d[dim], 2);
  }
  for (int dim = 0; dim < DIM; dim++){
    normal(dim) = -d[dim] / sqrt(R2);
  }
}

""";
    println(sbmcalc_file, content);
    ##############################################################################################
end