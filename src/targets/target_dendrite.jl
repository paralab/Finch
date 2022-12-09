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
        dendrite_inputdata_file(var);
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
        
        code = gen_func_name * "(p, t)";
        
    elseif op === :KNOWN_VAR
        # code = "p_data_->valueFEM(fe, dofind)";
        if IR.args[4] == 1
            code = "p_data_->valueFEM(fe, " * string(IR.args[2]-1) * ")";
        else
            code = "p_data_->valueFEM(fe, " * string(IR.args[2]-1) * " + eq_idata->NUM_VARS)";
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
        #     var[0] = sin(M_PI * x[0]) * sin(M_PI * x[1]) * sin(M_PI * x[2]);
        #     var[1] = ...
        # };
        solution_step *= "    std::function<void(const double *, double *)> initial_condition = [](const double *pt, double *var) {\n";
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
        solution_step *= "    TALYFEMLIB::PrintStatus(\"TIme step progress (%) 0\");\n";
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
        solution_step *= "    // Beginning time steps\n";
        solution_step *= "    unsigned int currentStep = 0;\n";
        solution_step *= "    double currentT = 0.0;\n";
        solution_step *= "    double dt = idata.dt;\n";
        solution_step *= "    unsigned int nSteps = idata.nSteps;\n";
        solution_step *= step_loop;
        
    else # not time dependent
        # A single solve step
        solution_step = "    // Solve it\n    " * project_name * "Solver->solve();\n"
    end # solution_step
    
    # Output
    output_part = "    // Output to vtk\n";
    output_part *= "    static const char *varname[]{\"" * dof_names[1] * "\"";
    for i=2:length(dof_names)
        output_part *= ", \"" * dof_names[i] * "\"";
    end
    output_part *= "};\n";
    if prob.time_dependent
        output_part *= "    petscVectopvtu(octDA, prev_solution, \"$(project_name)\", varname, physDomain, false, false, $(project_name)NodeData::NUM_VARS);\n";
    else
        output_part *= "    petscVectopvtu(octDA, $(project_name)Solver->getCurrentSolution(), \"$(project_name)\", varname, physDomain, false, false, $(project_name)NodeData::NUM_VARS);\n";
    end
    
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
#include <PostProcessing/postProcessing.h>
#include <Traversal/Analytic.h>
#include <IO/VTU.h>
// Petsc
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <PETSc/IO/petscVTU.h>

// generated
#include <$(project_name)Equation.h>
#include <$(project_name)BoundaryConditions.h>
#include <$(project_name)NodeData.h>
#include <$(project_name)InputData.h>

using namespace PETSc;
int main(int argc, char* argv[]) {
    /// initialize
    dendrite_init(argc, argv);
    int rank = TALYFEMLIB::GetMPIRank();
    
    /// read parameters from config.txt
    heat3dInputData idata;
    std::ifstream configFile("config.txt");
    if (configFile.good()) {
        if (!idata.ReadFromFile()) {
            throw std::runtime_error("[ERR] Error reading input data, check the config file!");
        }
        if (!idata.CheckInputData()) {
            throw std::runtime_error("[ERR] Problem with input data, check the config file!");
        }
    } else{
        TALYFEMLIB::PrintStatus("No config.txt file found. Copy the generated config.txt into the build directory");
        return -1;
    }
    
    // Configuration variables
    unsigned int eleOrder = idata.basisFunction;
    unsigned int level = idata.meshDef.refineLevel;
    bool mfree = idata.mfree;
    
    // Print status
    TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
    TALYFEMLIB::PrintStatus("DIM =  ", DIM);
    TALYFEMLIB::PrintStatus("Element order ", eleOrder);
    TALYFEMLIB::PrintStatus("Refinement level ", level);
    TALYFEMLIB::PrintStatus("Mfree ", mfree);
    
    // Make the tree
    std::vector<TREENODE> treePart;
    DA *octDA = createRegularDA(treePart, level, eleOrder);
    
    // Problem specification
    static const unsigned int NUM_VARS = $(dofs_per_node);
    DomainInfo physDomain;
    physDomain.min = idata.meshDef.min;
    physDomain.max = idata.meshDef.max;
    
    // Time steps (define even if no time stepping)
    unsigned int currentStep = 0;
    double currentT = 0.0;
    double dt = idata.dt;
    unsigned int nSteps = idata.nSteps;
    
    auto $(project_name)Eq = new TalyEquation<$(project_name)Equation, $(project_name)NodeData>(octDA, physDomain, NUM_VARS, false, &idata);
    
    // Setup PETSc solver
    LinearSolver *$(project_name)Solver = setLinearSolver($(project_name)Eq, octDA, NUM_VARS, mfree);
    /// apply solver parameters from config.txt
    idata.solverOptions$(project_name).apply_to_petsc_options("-$(project_name)_");
    KSP m_ksp = $(project_name)Solver->ksp();
    KSPSetOptionsPrefix(m_ksp, "$(project_name)_");
    KSPSetFromOptions(m_ksp);
    
    // Boundary conditions
    $(project_name)BoundaryConditions $(project_name)BC(&idata);
    $(project_name)Solver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary {
        Boundary b;
        $(project_name)BC.get$(project_name)BoundaryCondition(b, pos, currentT);
        return b;
    });
    
$(solution_step)

$(output_part)
    
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
#include <DataTypes.h>
#include <Basis/MatVec.h>
#include <Basis/Vec.h>
#include <Basis/Mat.h>
class $(project_name)Equation : public TALYFEMLIB::CEquation<$(project_name)NodeData> {
    
    public:
    explicit $(project_name)Equation(const $(project_name)InputData * idata)
        : TALYFEMLIB::CEquation<$(project_name)NodeData>(false, TALYFEMLIB::kAssembleGaussPoints) {
        eq_idata = idata;
        dt = idata->dt;
        currentT = 0.0;
    }
    
    // Do nothing. Needed by talyfem? //////////////////////
    void Solve(double dt, double t) override {
        assert(false);
    }
    void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                    TALYFEMLIB::ZEROARRAY<double> &be) override {
        assert(false);
    }
    ///////////////////////////////////////////////////////
    
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
    
    void advanceTime(const double dti){
        currentT += dti;
    }
    
    protected:
    const $(project_name)InputData *eq_idata;
    double dt;
    double currentT;
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
    
    if config.dimension == 1
        extract_coords = """
        double x = pos.x();
        Point<1> domainMin(inputData_->meshDef.min);
        Point<1> domainMax(inputData_->meshDef.max);
"""
    elseif config.dimension == 2
        extract_coords = """
        double x = pos.x();
        double y = pos.y();
        Point<2> domainMin(inputData_->meshDef.min);
        Point<2> domainMax(inputData_->meshDef.max);
"""
    elseif config.dimension == 3
        extract_coords = """
        double x = pos.x();
        double y = pos.y();
        double z = pos.z();
        Point<3> domainMin(inputData_->meshDef.min);
        Point<3> domainMax(inputData_->meshDef.max);
"""
    end
    
    # bool on_BID_1 = (fabs(y - domainMax.x(1)) < eps);
    bid_defs = "";
    nbids = length(prob.bid);
    for i=1:nbids
        bid_name = "on_bid_"*string(prob.bid[i]);
        bid_location = cpp_bid_def(prob.bid_def[i]);
        bid_defs *= "    bool " * bid_name * " = " * bid_location * ";\n";
    end
    
    # genfunctions
    bdry_genfunctions = [];
    
    # if on_BID_1 {
    #     b.addDirichlet($(project_name)NodeData::U_1, 0.0);
    # } else if on_BID_2 {
    #     b.addDirichlet($(project_name)NodeData::U_1, 0.0);
    # }
    dirichlet_bc = "";
    for i=1:nbids
        bid_name = "on_bid_"*string(prob.bid[i]);
        # One line for each var with a dirichlet bdry on this bid
        add_dirichlet = "";
        for vi = 1:length(var)
            var_ind = var[vi].index;
            # One line for each component
            for ci=1:var[vi].total_components
                if prob.bc_type[var_ind] == DIRICHLET
                    # Constant or genfunction values
                    if typeof(prob.bc_func[var_ind][ci]) <: Number
                        bc_val = string(prob.bc_func[var_ind][ci]);
                    elseif typeof(prob.bc_func[var_ind][ci]) == GenFunction
                        bc_val = prob.bc_func[var_ind][ci].name * "(pos, t)";
                        push!(bdry_genfunctions, prob.bc_func[var_ind][ci]);
                    else # uh oh
                        printerr("Boundary values that are not constant or [x,y,z,t] gen functions must be entered manually.");
                        bc_val = "0.0";
                    end
                    
                    dof_name = string(var[vi].symbol) * "_" * string(ci) * "_dofind"; # u_1_dofind
                    add_dirichlet *= "b.addDirichlet($(project_name)NodeData::$(dof_name), $(bc_val));\n";
                end
            end
        end
        
        # Put them in their repective bid
        if i==1
            dirichlet_bc *= "    if(" * bid_name * "){\n        " * add_dirichlet * "    }";
        else
            dirichlet_bc *= " else if(" * bid_name * "){\n        " * add_dirichlet * "    }";
        end
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
            str = cpp_genfunction_to_string(bdry_genfunctions[i]);
            fun = "    return "*str*";";
            
            function_defs *= "double $(project_name)BoundaryConditions::" * bdry_genfunctions[i].name * args * "{\n";
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
    private:
        const $(project_name)InputData *inputData_;
    public:
        $(project_name)BoundaryConditions(const $(project_name)InputData *idata);
        void get$(project_name)BoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos, const double t);
$(function_decs)
};
$(project_name)BoundaryConditions::$(project_name)BoundaryConditions(const $(project_name)InputData *idata) {
    inputData_ = idata;
}

void $(project_name)BoundaryConditions::get$(project_name)BoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos, const double t) {
    static constexpr double eps = 1e-14;
$(extract_coords)
    
$(bid_defs)
    
$(dirichlet_bc)
    
$(reference_points)
}

$(function_defs)
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

class $(project_name)NodeData {
    public:
    // number of  variables in this NodeData
    static const unsigned int NUM_VARS = $(dofs_per_node);
    // values available to equation
    double var_values[$(value_count)];

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
        return $(project_name)NODEDATA_MAX;
    }
};
  
"""
    println(file, content);
end

#=
Stores configuration input from config.txt or default generated values
The basics are:
basisFunction = "linear" # "quadratic" "cubic" "hermite"
mfree = false # "true"

Time stepping:
dt = 0.1
nSteps = 10
finalT = 1,0

Mesh options: (depends on MeshDef struct)
background_mesh = {
  refineLevel = 4
  min = [0.0, 0.0]   
  max = [1.0, 1.0] 
  # max or scalingFactor = [1.0, 1.0] 
  # You can pass either of these two but not both. InputData needs to be updated accordingly 
}

PETSc options: (any other PETSc options can be added as needed)
solver_options = {
  ksp_max_it = 2000
  ksp_type = "bcgs"
  pc_type = "asm"
  ksp_atol = 1e-7
  ksp_rtol = 1e-10
  ksp_monitor = ""
  ksp_converged_reason = ""
  
  # nonlinear snes
  snes_atol = 1e-8
  snes_rtol = 1e-8
  snes_stol = 1e-10
  snes_max_it = 40
  snes_max_funcs = 80000
  snes_monitor = ""
  snes_converged_reason = ""
  
  # multigrid
  ksp_type = "bcgs"
  pc_type = "gamg"
  pc_gamg_asm_use_agg = True
  mg_levels_ksp_type = "bcgs"
  mg_levels_pc_type = "asm"
  mg_levels_sub_pc_type = "lu"
  mg_levels_ksp_max_it = 50
  
  # direct solver
  ksp_type = "bcgs"
  pc_type = "lu"
  pc_factor_mat_solver_package = "mumps"
  pc_factor_levels = 3
}
=#
function dendrite_inputdata_file(var)
    project_name = finch_state.project_name;
    file = add_generated_file(project_name*"InputData.h", dir="include");
    
    params = finch_state.target_parameters;
    DIM = finch_state.config.dimension;
    order = finch_state.config.basis_order_min;
    if order > 3
        printerr("This target only supports element order up to 3. Changing to 3.")
        order = 3;
    end
    if order == 1
        default_order = "linear";
    elseif order == 2
        default_order = "quadratic";
    else # order == 3
        default_order = "cubic";
    end
    default_matfree = finch_state.config.linalg_matrixfree;
    default_refine = 3;
    if haskey(params, :refineLevel)
        default_refine = params[:refineLevel];
    end
    
    time_stepping_parts = "";
    time_stepping_read = "";
    if finch_state.prob.time_dependent
        dt = finch_state.time_stepper.dt;
        nsteps = finch_state.time_stepper.Nsteps;
        finalT = dt*nsteps;
        time_stepping_parts = "        // Time stepping info (can be overridden by config.txt values)\n";
        time_stepping_parts *= "        dt = $(dt);\n";
        time_stepping_parts *= "        unsigned int nSteps = $(nsteps);\n";
        time_stepping_parts *= "        double finalT = $(finalT);\n";
        
        time_stepping_read = "        // Time stepping info\n";
        time_stepping_read *= "        if (ReadValue(\"dt\", dt)) {};\n";
        time_stepping_read *= "        if (ReadValue(\"nSteps\", nSteps)) {};\n";
        time_stepping_read *= "        if (ReadValue(\"finalT\", finalT)) {};\n";
    end
    
    content = """
#pragma once

#include <talyfem/input_data/input_data.h>
#include <DataTypes.h>
#include <point.h>

/**
 * Parameters for mesh
 */
struct MeshDef : public DomainInfo {
    /// refinement level
    DENDRITE_UINT refineLevel = $(default_refine);
    
    void read_from_config(const libconfig::Setting &root) {
        if (root.exists("refineLevel")) {
            refineLevel= (DENDRITE_UINT) root["refineLevel"];
        }
        if (root.exists("min")) {
            for (DENDRITE_UINT dim = 0; dim < $(DIM); dim++) {
                min[dim] = (DENDRITE_REAL) root["min"][dim];
            }
        }else{
            TALYFEMLIB::PrintWarning("No min defined for mesh. Using 0.0");
            for (DENDRITE_UINT dim = 0; dim < $(DIM); dim++) {
                min[dim] = (DENDRITE_REAL) 0.0;
            }
        }
        if (root.exists("max")) {
            for (DENDRITE_UINT dim = 0; dim < $(DIM); dim++) {
                max[dim] = (DENDRITE_REAL) root["max"][dim];
            }
        }else if (root.exists("scalingFactor")) {
            for (DENDRITE_UINT dim = 0; dim < $(DIM); dim++) {
                max[dim] = min[dim] + (DENDRITE_REAL) root["scalingFactor"][dim];
            }
        }else{
            TALYFEMLIB::PrintWarning("No max or scalingFactor defined for mesh. Using 1.0");
            for (DENDRITE_UINT dim = 0; dim < $(DIM); dim++) {
                max[dim] = (DENDRITE_REAL) 1.0;
            }
        }
    }
};

class $(project_name)InputData : public TALYFEMLIB::InputData {
    public:
        /// Mesh definition
        MeshDef meshDef;
        /// Basis function
        std::string basisFunctionStr = "$(default_order)";
        /// Matrix  free
        bool mfree = $(default_matfree);
        /// Time info
        double dt = 0.0;
        double currentT = 0.0;
        
$(time_stepping_parts)
        
        SolverOptions solverOptions$(project_name);
        
        bool dump_vec = false; // This is used for regression testing
        
    bool ReadFromFile(const std::string &filename = std::string("config.txt")) {
        ReadConfigFile(filename);
        
        /// mesh size and level
        meshDef.read_from_config(cfg.getRoot()["background_mesh"]);
        
        /// basis function order
        basisFunction = TALYFEMLIB::basis_string_to_enum(basisFunctionStr);
        if (ReadValue("basisFunction", basisFunctionStr)) {
            if (!basisFunctionStr.empty()){
                basisFunction = TALYFEMLIB::basis_string_to_enum(basisFunctionStr);
            }
        }
        
        /// Other parameters
        if (ReadValue("mfree", mfree)) {}
        if (ReadValue("dump_vec", dump_vec)) {}
        
$(time_stepping_read)
            
        /// PETSc options
        solverOptions$(project_name) = read_solver_options(cfg, "solver_options");
        return true;
    }

    /// check if the input are valid
    bool CheckInputData() {
        /// Matrix free version of the code cannot have pre-conditioner
        if (mfree) {
            if (solverOptions$(project_name).vals.count("pc_type") == 1) {
                solverOptions$(project_name).vals.at("pc_type") = "none";
                TALYFEMLIB::PrintWarning("mfree = True, changing pc_type to 'none'");
            }
        }
        
        /// Basis function order
        if (basisFunction < 1 || basisFunction > 4){
            basisFunction = TALYFEMLIB::basis_string_to_enum("$(default_order)");
            TALYFEMLIB::PrintWarning("invalid basis function selection, changing order to $(default_order)");
        }
        return true;
    }
};

"""
    println(file, content);
end

#=
The cmakelists.txt and readme.txt files
as well as the config.txt file
=#
function dendrite_build_files()
    project_name = finch_state.project_name;
    cmake_file = add_generated_file("CMakeLists.txt", make_header_text=false);
    readme_file = add_generated_file("README.txt", make_header_text=false);
    config_file = add_generated_file("config.txt", make_header_text=false);
    
end

###########################################################################################################
# Static files that will be written as is

 #####   ######      ###     ######  ######   ######       #######  ######  ##        #######   #####
###   #    ##       ## ##      ##      ##    ##    ##      ##         ##    ##        ##       ###   #
  ###      ##      ##   ##     ##      ##    ##            #######    ##    ##        ######     ### 
#   ###    ##     #########    ##      ##    ##    ##      ##         ##    ##        ##       #   ###
 #####     ##    ##       ##   ##    ######   ######       ##       ######  ########  #######   #####

###########################################################################################################

function dendrite_static_files()
    # TODO
end