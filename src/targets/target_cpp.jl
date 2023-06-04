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
        cpp_write_static_files();
        # Build and info files
        cpp_build_files(finch_state.target_parameters);
        
        # The generated code
        cpp_main_file(var, IR);
        
        cpp_genfunction_file();
        cpp_boundary_file(var);
        cpp_solve_file(var);
        cpp_utils_file(var, config);
        cpp_output_file(var);
    end
    
    # All procs will write their mesh file
    cpp_mesh_file();
    
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
                    code *= generate_from_IR_external(IR.index[i], IRtypes) * "-1";
                elseif typeof(IR.index[i]) <: Int
                    code *= string(IR.index[i]-1);
                else
                    code *= string(IR.index[i]) * "-1";
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
                
                return indent * "// GLOBAL_ARRAY_ASSIGNMENT: "*string(IR.args[1].label);
            end
            
            # Note that if the symbol (IR.args[1]) is not yet on the declared list, 
            # it will be here.
            if typeof(IR.args[1]) <: IR_data_node
                typestring = cpp_type_name(IRtypes.name[IR.args[1].type]);
                # arrays need a *
                if length(IR.args[1].size) > 0
                    typestring *= "*";
                end
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
            
            # Note that if the symbol (IR.args[1]) is not yet on the declared list, 
            # it will be here.
            if typeof(IR.args[2]) <: IR_data_node
                typestring = cpp_type_name(IRtypes.name[IR.args[2].type]);
                name = IR.args[2].label;
                # check if it has been declared
                foundit = false;
                for d in declared
                    if name === d
                        foundit = true;
                        break;
                    end
                end
                if !foundit
                    code = indent * typestring * " " * generate_from_IR_external(IR.args[2], IRtypes);
                    push!(declared, name);
                else
                    code = indent * generate_from_IR_external(IR.args[2], IRtypes);
                end
                
            elseif typeof(IR.args[2]) <: IR_part
                code = indent * generate_from_IR_external(IR.args[2], IRtypes);
            else
                code = indent * string(IR.args[2]);
            end
            
            code *= " "*math_op*"= ";
            if typeof(IR.args[2]) <: IR_part
                code *= generate_from_IR_external(IR.args[2], IRtypes);
            else
                code *= string(IR.args[2]);
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
        # "finch::genfunction_"*string(cval)*"(x, y, z, t, nid)"
        gen_func_name = finch_state.coefficients[IR.args[2]].value[IR.args[3]].name;
        
        code = "finch::" * gen_func_name * "(x, y, z, t, nodeID)";
        
    elseif op === :KNOWN_VAR
        # code = "variables[" * string(IR.args[2]) * "].values[nodeID]";
        code = "SORRY_THIS_IS_NOT_GOING_TO_WORK_YET_SORRY_SORRY"
        
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
        # args[2] is the I matrix, args[3] is the J matrix
        # code = indent * "set_matrix_indices!(" * generate_from_IR_external(IR.args[2], IRtypes) * ", " * 
        #             generate_from_IR_external(IR.args[3], IRtypes) * ", dofs_per_node, mesh)";
        code = indent * "// INIT_MATRIX_IJ_FV named op not ready";
        
    elseif op === :GLOBAL_FORM_MATRIX
        # // Petsc assembles the system
        code =  indent * "// Petsc assembles the system\n" *
                indent * "PetscCall( VecAssemblyBegin(global_vector) );\n" * 
                indent * "PetscCall( VecAssemblyEnd(global_vector) );\n" * 
                indent * "PetscCall( MatAssemblyBegin(global_matrix, MAT_FINAL_ASSEMBLY) );\n" * 
                indent * "PetscCall( MatAssemblyEnd(global_matrix, MAT_FINAL_ASSEMBLY) );\n";
        
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
        # int solve_code = finch::solve(lhs_matrix, rhs_vector, solution);
        code = indent * "int solve_code = finch::solve(global_matrix, global_vector, solution);";
        
    elseif op === :GLOBAL_DISTRIBUTE_VECTOR
        # code = indent * generate_from_IR_external(IR.args[2], IRtypes) * " = distribute_solution("*generate_from_IR_external(IR.args[2], IRtypes) * ", nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config);"
        code = indent * "// GHOST_DISTRIBUTE_VECTOR named op not ready";
        
    elseif op === :GATHER_VARS
        # place variable arrays in global vector
        code = indent * "// TODO: When initial conditions are needed, copy variable values into global vector";
        
    elseif op === :BDRY_TO_VECTOR
        # FV_copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node);
        # if length(IR.args) < 3
        #     code = indent * "copy_bdry_vals_to_vector(var, "* generate_from_IR_external(IR.args[2], IRtypes) *", mesh, dofs_per_node);";
        # else
        #     code = indent * "copy_bdry_vals_to_vector("* generate_from_IR_external(IR.args[3], IRtypes) *", "* 
        #                     generate_from_IR_external(IR.args[2], IRtypes) *", mesh, dofs_per_node);";
        # end
        code = indent * "// BDRY_TO_VECTOR named op not ready";
        
    elseif op === :BDRY_TO_VAR
        # copy_bdry_vals_to_variables(var, solution, mesh, dofs_per_node, true)
        # if length(IR.args) < 4
        #     code = indent * "copy_bdry_vals_to_variables(var, "* generate_from_IR_external(IR.args[2], IRtypes) *
        #                     ", mesh, dofs_per_node, "* generate_from_IR_external(IR.args[3], IRtypes) *");";
        # else
        #     code = indent * "copy_bdry_vals_to_variables("* generate_from_IR_external(IR.args[4], IRtypes) *", "* 
        #                     generate_from_IR_external(IR.args[2], IRtypes) *
        #                     ", mesh, dofs_per_node, "* generate_from_IR_external(IR.args[3], IRtypes) *");";
        # end
        code = indent * "// BDRY_TO_VAR named op not ready";
        
    elseif op === :SCATTER_VARS
        # place solution in variable arrays
        code = indent * "// TODO: When more than one variable, distribute global solution to variable arrays";
        
    elseif op === :LOCAL2GLOBAL
        # put elemental matrix and vector in global system
        code = """
        // Place elemental matrix and vector in global versions
        for (int nodei = 0; nodei < mesh.nodes_per_element[eid]; nodei++) { // Loop over rows
            unsigned long this_nodeid = mesh.loc2glb[eid][nodei];
            int8_t this_bid = mesh.nodebid[this_nodeid]; // The BID for the node corresponding to these rows.
            for (int dofi = 0; dofi < dofs_per_node; dofi++) {               //
                row_id = mesh.partition2global_n[this_nodeid] * dofs_per_node + dofi;
                // For the vector
                local_row = nodei * dofs_per_node + dofi;
                vec_indices[local_row] = row_id;
                
                // If it is a boundary node with a boundary condition, this row is special.
                if(this_bid >= 0 && finch::get_bc_type(this_bid) > 0){
                    // Only do this once to avoid the add_values/insert_values issue
                    // Only set if this node belongs to this partition.
                    if(mesh.node_owner[this_nodeid] == partition_index && this_bid < 100){
                        for(int dimi=0; dimi<mesh.dimension; dimi++){
                            node_x[dimi] = mesh.allnodes[this_nodeid*mesh.dimension + dimi];
                        }
                        // For now, only Dirichlet -> identity row
                        PetscCall( MatSetValue(lhs_matrix, row_id, row_id, 1.0, ADD_VALUES) );
                        vec_values[local_row] = finch::evaluate_bc(node_x, this_bid, mesh.partition2global_n[this_nodeid]);
                        
                        mesh.nodebid[this_nodeid] = 127;
                    }else{
                        // It has already been set or will be by another proc, do nothing
                    }
                    
                }else{ // No boundary condition
                    tmp_index = 0;
                    for (int nodej = 0; nodej < mesh.nodes_per_element[eid]; nodej++) { // Loop over columns
                        for (int dofj = 0; dofj < dofs_per_node; dofj++) {               //
                            col_indices[tmp_index] = mesh.partition2global_n[mesh.loc2glb[eid][nodej]] * dofs_per_node + dofj;
                            mat_values[tmp_index] = element_matrix[local_row * (nodes_per_element*dofs_per_node) + nodej*dofs_per_node + dofj]; // change from block dofs to interleaved dofs
                            tmp_index++;
                        }
                    }
                    // Is it faster to set this one row at a time, or the whole thing at once?
                    // This is one row at a time.
                    PetscCall( MatSetValues(lhs_matrix, 1, &row_id, col_indices.size(), col_indices.data(), mat_values.data(), ADD_VALUES) );
                    
                    vec_values[local_row] = element_vector[dofi * nodes_per_element + nodei]; // change from block dofs to interleaved dofs
                }
            }// dof loop
        }// node loop
        
        // Set vector values
        PetscCall( VecSetValues(rhs_vector, vec_indices.size(), vec_indices.data(), vec_values.data(), ADD_VALUES) );
        """
        
    elseif op === :LOCAL2GLOBAL_VEC
        # put elemental vector in global system
        code = """
        // Place elemental vector in global version
        for (int nodei = 0; nodei < mesh.nodes_per_element[eid]; nodei++) { // Loop over rows
            unsigned long this_nodeid = mesh.loc2glb[eid][nodei];
            int8_t this_bid = mesh.nodebid[this_nodeid]; // The BID for the node corresponding to these rows.
            for (int dofi = 0; dofi < dofs_per_node; dofi++) {               //
                row_id = mesh.partition2global_n[this_nodeid] * dofs_per_node + dofi;
                // For the vector
                local_row = nodei * dofs_per_node + dofi;
                vec_indices[local_row] = row_id;
                
                // If it is a boundary node with a boundary condition, this row is special.
                if(this_bid >= 0 && finch::get_bc_type(this_bid) > 0){
                    // Only do this once to avoid the add_values/insert_values issue
                    // Only set if this node belongs to this partition.
                    if(mesh.node_owner[this_nodeid] == partition_index && this_bid < 100){
                        for(int dimi=0; dimi<mesh.dimension; dimi++){
                            node_x[dimi] = mesh.allnodes[this_nodeid*mesh.dimension + dimi];
                        }
                        // For now, only Dirichlet -> identity row
                        vec_values[local_row] = finch::evaluate_bc(node_x, this_bid, mesh.partition2global_n[this_nodeid]);
                        
                        mesh.nodebid[this_nodeid] = 127;
                    }else{
                        // It has already been set or will be by another proc, do nothing
                    }
                    
                }else{ // No boundary condition
                    vec_values[local_row] = element_vector[dofi * nodes_per_element + nodei]; // change from block dofs to interleaved dofs
                }
            }// dof loop
        }// node loop
        
        // Set vector values
        PetscCall( VecSetValues(rhs_vector, vec_indices.size(), vec_indices.data(), vec_values.data(), ADD_VALUES) );
        """
        
    elseif op === :ADD_GLOBAL_VECTOR_AND_NORM
        # u = u + du 
        # and find absolute and relative norms of du
        # args are: u, du, abs_residual, rel_residual
        # code = generate_from_IR_external(IntermediateRepresentation.generate_residual_norms_and_update(
        #     IR.args[2], IR.args[3], IR.args[4], IR.args[5]), IRtypes, indent);
        code = indent * "// ADD_GLOBAL_VECTOR_AND_NORM named op not ready";
        
    elseif op === :UPDATE_GLOBAL_VECTOR_AND_NORM
        # b = a
        # and find absolute and relative norms of (a-b)
        # args are: a, b, abs_residual, rel_residual
        # code = generate_from_IR_external(IntermediateRepresentation.generate_difference_norms_and_update(
        #     IR.args[2], IR.args[3], IR.args[4], IR.args[5]), IRtypes, indent);
        code = indent * "// UPDATE_GLOBAL_VECTOR_AND_NORM named op not ready";
        
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
        for(unsigned int row=1; row<=nodes_per_element; row++){
$init_lines
            for(unsigned int col=1; col<=nodes_per_element; col++){
$compute_lines
            }
        }";
    
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
        
        code = generate_from_IR_external(generate_linalg_TDM_product(IR.args[2], IR.args[3], IR.args[4], i_index, j_index, k_index), IRtypes);
            
    elseif op === :LINALG_Tv
        # Tcode = generate_from_IR_external(IR.args[2], IRtypes) * "[col + (row-1)*qnodes_per_element]";
        # vcode = generate_from_IR_external(IntermediateRepresentation.apply_indexed_access(IR.args[3], [:col], IRtypes), IRtypes);
        # code = "*(" * Tcode *", "* vcode * ")";
        
        i_index = :row;
        j_index = :col;
        
        code = generate_from_IR_external(generate_linalg_Tv_product(IR.args[2], IR.args[3], i_index, j_index), IRtypes);
    
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
    elseif T == "Int64"
        return "int64_t";
    elseif T == "Float32"
        return "float";
    elseif T == "Float64"
        return "double";
    elseif T == "CustomFloat"
        return "double";
    elseif T == "CustomInt"
        return "int";
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

function cpp_genfunction_to_string(genfun)
    newex = cpp_swap_symbol(:pi, :M_PI, genfun.expr); # swap pi for M_PI
    newex = cpp_change_math_ops(newex); # change operators to match C++
    s = string(newex);
    ns = replace(s, r"([\d)])([(A-Za-z])" => s"\1*\2"); # explicitly multiply with "*" (2*x not 2x)
    return ns;
end

#######################################################
# Write code files

#######  ######  ##        #######   #####
##         ##    ##        ##       ###   #
#######    ##    ##        ######     ### 
##         ##    ##        ##       #   ###
##       ######  ########  #######   #####

#######################################################

function cpp_main_file(var, IR)
    project_name = finch_state.project_name;
    file = add_generated_file(project_name*".cpp", dir="src");
    hfile = add_generated_file(project_name*".hpp", dir="include");
    
    if !(typeof(var) <:Array)
        var = [var];
    end
    
    # gather important numbers
    varcount = 1;
    dofs_per_node = var[1].total_components;
    dofs_per_loop = length(var[1].symvar);
    offset_ind = [0];
    dof_names = [];
    if length(var) > 1
        varcount = length(var);
        offset_ind = zeros(Int, varcount);
        for i=2:length(var)
            offset_ind[i] = dofs_per_node;
            dofs_per_node += var[i].total_components;
            dofs_per_loop += length(var[i].symvar);
            for i=1:length(var[i].symvar)
                push!(dof_names, string(var[1].symvar[i]));
            end
        end
    end
    offset_ind_str = cpp_gen_string(offset_ind);
    
    # Make a list of all needed timers
    timer_names = [];
    # Output file for timing info
    timer_file = project_name*"_timing.txt";
    
    content = """
/**
* Main file for $(project_name).
*
* See readme.txt for instructions
*/ 

// Option parameters are defined here
#define finch_KSP_TYPE KSPCG
#define finch_KSP_RTOL 1e-6
#define finch_KSP_ATOL 1e-12
#define finch_KSP_DTOL PETSC_DEFAULT
#define finch_KSP_MAXITERS 1000
#define finch_PC_TYPE PCJACOBI

#include \"$(project_name).hpp\"
finch::Mesh mesh;
finch::Refel refel;
finch::GeometricFactors geometric_factors;

#ifndef PetscCall
void PetscCall(PetscErrorCode e){
    // handle Petsc errors TODO
}
#endif

// This file contains the solve function.
// int solve(Mat lhs, Vec rhs, Vec solution)
#include "finch_solve.cpp"

#include "finch_utils.cpp"

// The main function 
int main(int argc, char* argv[]) {
    
    // Init PETSc and MPI
    PetscCall( PetscInitialize(&argc, &argv, NULL, NULL) );
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status Stat;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    // A tiny number that is approximately zero
    // This should be scaled when used for geometry
    const double eps = 1E-14;
    
    // Number of dofs per node
    const unsigned int dofs_per_node = $(dofs_per_node);
    
    // timing variables
    double total_time = 0.0;
    double setup_time = 0.0;
    double assemble_time = 0.0;
    double solve_time = 0.0;
    const char* timer_filename = \"$(timer_file)\"; 
    
    total_time = MPI_Wtime();
    setup_time = total_time;
    
    // For now partition index is same as rank.
    int num_partitions = size;
    int partition_index = rank;
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Import mesh, refel and geometric factors precomputed in Finch
    ////////////////////////////////////////////////////////////////////////////////////////////
    mesh.import_mesh("MeshData", num_partitions, partition_index);
    refel.import_refel("RefelData");
    geometric_factors.import_geometric_factors("GeoData", mesh.dimension, num_partitions, partition_index);
    
    if (!rank) {
        std::cout << "============ Imported mesh, geo factors and refel ================\\n";
        std::cout << "Read "<<refel.dimension<<"D refel with "<<refel.Nfaces<<" faces and "<<refel.Np<<" nodes.\\n";
        std::cout << "Read "<<mesh.dimension<<"D mesh with "<<mesh.nel_global<<" elements and "<<mesh.nnodes_global<<" nodes.\\n";
        std::cout << "Read geo factors with "<<geometric_factors.vals_per_element<<" values per element.\\n";
    }
    
    // Define needed symbols
    // Useful symbols for FEM
    double* Q = refel.Q;
    double* wg = refel.wg;
    
    // Prepare some useful numbers
    // dofs_per_node = $(dofs_per_node); // set above
    const unsigned int dofs_per_loop = $(dofs_per_loop);
    unsigned int* dof_offsets = $(offset_ind_str);
    unsigned int nnodes_partition = mesh.nnodes_local;
    unsigned int nnodes_global = mesh.nnodes_global;
    unsigned int num_elements = mesh.nel_owned;
    unsigned int num_elements_global = mesh.nel_global;
    unsigned int num_elements_ghost = mesh.nel_ghost;
    unsigned int num_faces = mesh.nface_owned + mesh.nface_ghost;
    
    unsigned int dofs_global = dofs_per_node * nnodes_global;
    unsigned int dofs_partition = dofs_per_node * nnodes_partition;
    unsigned int fv_dofs_global = dofs_per_node * num_elements_global;
    unsigned int fv_dofs_partition = dofs_per_node * (num_elements + num_elements_ghost);
    
    unsigned int nodes_per_element = refel.Np;
    unsigned int qnodes_per_element = refel.Nqp;
    unsigned int faces_per_element = refel.Nfaces;
    unsigned int nodes_per_face = refel.Nfp[1];
    unsigned int dofs_per_element = dofs_per_node * nodes_per_element;
    unsigned int local_system_size = dofs_per_loop * nodes_per_element;
    
    unsigned int nnodes_owned = mesh.nnodes_local - mesh.nnodes_borrowed;
    unsigned int dofs_owned = nnodes_owned * dofs_per_node;
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Create PETSc objects
    ////////////////////////////////////////////////////////////////////////////////////////////
    //unsigned long dofs_global = dofs_per_node * mesh.nnodes_global;
    //unsigned long dofs_partition = dofs_per_node * (mesh.nnodes_local - mesh.nnodes_borrowed);
    //unsigned int* dofs_per_element = new unsigned int[mesh.nel_local];
    //for (unsigned int e = 0; e < mesh.nel_local; e++) {
    //    dofs_per_element[e] = mesh.nodes_per_element[e] * dofs_per_node;
    //}
    
    // create rhs and solution vectors
    Vec rhs_vector, solution, exact;
    PetscCall( VecCreate(comm, &rhs_vector) );
    PetscCall( VecCreate(comm, &solution) );
    if (size > 1) {
        PetscCall( VecSetType(rhs_vector, VECMPI) );
        PetscCall( VecSetType(solution, VECMPI) );
    } else {
        PetscCall( VecSetType(rhs_vector, VECSEQ) );
        PetscCall( VecSetType(solution, VECSEQ) );
    }
    PetscCall( VecSetSizes(rhs_vector, dofs_owned, dofs_global) );
    PetscCall( VecSetSizes(solution, dofs_owned, dofs_global) );
    PetscCall( VecSet(rhs_vector, 0.0) );
    PetscCall( VecSet(solution, 0.0) );
    // exact solution
    PetscCall( VecDuplicate(solution, &exact) );
    
    // Create lhs matrix
    Mat lhs_matrix;
    int matrix_bands = dofs_per_element * (mesh.faces_per_element[0] + 1);
    PetscCall( MatCreate(comm, &lhs_matrix) );
    PetscCall( MatSetSizes(lhs_matrix, dofs_owned , dofs_owned , dofs_global, dofs_global) );
    if (size > 1) {
        PetscCall( MatSetType(lhs_matrix, MATMPIAIJ) );
        PetscCall( MatMPIAIJSetPreallocation(lhs_matrix, matrix_bands, PETSC_NULL, matrix_bands, PETSC_NULL) );
    } else {
        PetscCall( MatSetType(lhs_matrix, MATSEQAIJ) );
        PetscCall( MatSeqAIJSetPreallocation(lhs_matrix, matrix_bands, PETSC_NULL) );
    }
    // eals with an allocation error
    PetscCall( MatSetOption(lhs_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE) );
    // PetscCall( MatSetOption(lhs_matrix,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE) );
    PetscCall( MatZeroEntries(lhs_matrix) );
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Elemental assembly
    ////////////////////////////////////////////////////////////////////////////////////////////
    assemble_time = MPI_Wtime();
    setup_time = assemble_time - setup_time;
    
    std::vector<PetscScalar> mat_values(dofs_per_element);
    std::vector<PetscInt> col_indices(dofs_per_element);
    std::vector<PetscScalar> vec_values(dofs_per_element);
    std::vector<PetscInt> vec_indices(dofs_per_element);
    PetscInt row_id;
    int tmp_index, local_row;
    
    // double* element_matrix = new double [dofs_per_element * dofs_per_element]; // elemental matrix
    // double* element_vector = new double[dofs_per_element]; // elemental vector
    double* node_x = new double[mesh.dimension]; // nodal coordinates
    double t = 0.0;
    
    // hackity hack
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    unsigned long eid = 0;
    int index_offset = 0;
    unsigned long nodeID = 0;
    
    // Begin code generated from IR //////////////////////////////////////////////////////////////////////
    
    """*
    ######################################################################################################
    # This is the part that is generated from IR
    ######################################################################################################
    generate_from_IR_external(IR, indent="    ") * 
    """
    
    // End code generated from IR //////////////////////////////////////////////////////////////////////
    
    //delete [] element_matrix;
    //delete [] element_vector;
    //delete [] node_x;
    
    // Petsc assembles the system
    PetscCall( VecAssemblyBegin(rhs_vector) );
    PetscCall( VecAssemblyEnd(rhs_vector) );
    PetscCall( MatAssemblyBegin(lhs_matrix, MAT_FINAL_ASSEMBLY) );
    PetscCall( MatAssemblyEnd(lhs_matrix, MAT_FINAL_ASSEMBLY) );
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Solve
    ////////////////////////////////////////////////////////////////////////////////////////////
    solve_time = MPI_Wtime();
    assemble_time = solve_time - assemble_time;
    
    int solve_code = finch::solve(lhs_matrix, rhs_vector, solution);
    
    double end_time = MPI_Wtime();
    solve_time = end_time - solve_time;
    total_time = end_time - total_time;
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Output
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Output timing
    if (rank == 0) {
        std::cout << "-              setup time = " << setup_time << "\\n";
        std::cout << "- elemental assembly time = " << assemble_time << "\\n";
        std::cout << "-              solve time = " << solve_time << "\\n";
        std::cout << "-              total time = " << total_time << "\\n";
        
        // timing output file, open in append mode
        std::ofstream timing_file;
        timing_file.open(timer_filename, std::fstream::app);
        timing_file << "procs = " << size << ", elements = " << num_elements_global << ", dofs = " << dofs_global << "\\n";
        timing_file << "setup    assemble    solve    total\\n";
        timing_file << setup_time << assemble_time << solve_time << total_time << "\\n";
        timing_file.close();
    }
    
    // Output data 
    #include "finch_output.cpp"
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Cleanup
    ////////////////////////////////////////////////////////////////////////////////////////////
    // delete[] dofs_per_element;
    
    PetscCall( VecDestroy(&solution) );
    PetscCall( VecDestroy(&rhs_vector) );
    PetscCall( MatDestroy(&lhs_matrix) );
    
    PetscCall( PetscFinalize() );
    
    return 0;
}   
"""
    println(file, content);
    
    hppcontent = """
    #pragma once
    
    // Basics
    #include <iostream>
    #include <fstream>
    #include <math.h>
    #include <stdio.h>
    #include <time.h>
    #include <functional>
    
    // Parallel
    #include <mpi.h>
    
    // PETSc
    #include <petsc.h>
    #include <petscvec.h>
    #include <petscmat.h>
    #include <petscksp.h>
    #include <petscerror.h>
    
    // Finch
    #include "finch_mesh.hpp"
    #include "finch_geometry.hpp"
    #include "finch_functions.hpp"
    #include "finch_bdry.hpp"
    
    extern finch::Mesh mesh;
    extern finch::Refel refel;
    extern finch::GeometricFactors geometric_factors;    
"""
    println(hfile, hppcontent);
end

#=
The mesh file contains any code related to setting up the mesh.
The refeldata and meshdata file contains all of the data from the Refel and Grid structs
These are to be read into cpp by a function in the mesh file.
=#
function cpp_mesh_file()
    project_name = finch_state.project_name;
    config = finch_state.config;
    grid_data = finch_state.grid_data;
    geometric_factors = finch_state.geo_factors;
    refel = finch_state.refel;
    
    # First write the binary files for grid_data and refel
    part_suffix = "";
    if config.num_partitions > 1
        part_suffix = "_p" * string(config.partition_index);
    end
    meshfile = add_generated_file("MeshData"*part_suffix, make_header_text=false);
    Finch.CodeGenerator.write_grid_to_file(meshfile, grid_data);
    
    geofile = add_generated_file("GeoData"*part_suffix, make_header_text=false);
    Finch.CodeGenerator.write_geometric_factors_to_file(geofile, geometric_factors);
    
    if config.proc_rank == 0
        refelfile = add_generated_file("RefelData", make_header_text=false);
        Finch.CodeGenerator.write_refel_to_file(refelfile, refel);
    end
    
    # The files to read these in are in the static files section at the end of this file.
end

# Functions for boundary conditions, coeficients, etc.
function cpp_genfunction_file()
    project_name = finch_state.project_name;
    file = add_generated_file("finch_functions.cpp", dir="src");
    hfile = add_generated_file("finch_functions.hpp", dir="include");
    
    genfunctions = finch_state.genfunctions;
    
    function_defs = "";
    indent = "";
    args = ["x", "y", "z", "t", "nid"];
    argtypes = ["const double", "const double", "const double", "const double", "const unsigned long"];
    rettype = "double";
    function_declarations = "";
    for i = 1:length(genfunctions)
        str = cpp_genfunction_to_string(genfunctions[i]);
        fun = ["return "*str*";"];
        
        lines = cpp_function_def(indent, "finch::"*genfunctions[i].name, args, argtypes, rettype, fun);
        for j=1:length(lines)
            function_defs *= lines[j] * "\n";
        end
        
        function_declarations *= "    double "*genfunctions[i].name*"(const double x, const double y, const double z, const double t, const unsigned long nid);\n";
    end
    
    content = """
/*
Genfunctions for things like coefficients, boundary conditions etc.
*/
#include "finch_functions.hpp"

$function_defs

"""
    println(file, content);
    
    content = """
/*
Genfunctions for things like coefficients, boundary conditions etc.
*/
#pragma once
#include <math.h>
#include <vector>
#include <cstdint>

namespace finch{
$function_declarations
}
"""
    println(hfile, content);
    
end

# This file has boundary condition parts.
function cpp_boundary_file(var)
    project_name = finch_state.project_name;
    file = add_generated_file("finch_bdry.cpp", dir="src");
    hfile = add_generated_file("finch_bdry.hpp", dir="include");
    
    config = finch_state.config;
    prob = finch_state.prob;
    
    bids = prob.bid[1,:];
    default_bc = "0.0";
    bctype_const = Dict(NO_BC=>0, DIRICHLET=>1, NEUMANN=>3, FLUX=>4);
    
    bid_select = "";
    bctype_select = "
    if(bid == 127){ // This is a special number that will be skipped, but acts like dirichlet
        return 1;
    }else ";
    # TODO: This assumes one scalar variable only
    if !(typeof(var) <: Array)
        var = [var];
    end
    for b in bids
        sb = string(b);
        bid_select *= "    if(bid == $sb){\n"
        bctype_select *= "    if(bid == $sb){\n"
        
        bid_select *= "        // BC type: "*prob.bc_type[var[1].index, b]*"\n"
        bctype_select *= "        // BC type: "*prob.bc_type[var[1].index, b]*"\n"
        bctype_select *= "        return "*string(bctype_const[prob.bc_type[var[1].index, b]])*";\n"
        
        bfunc = prob.bc_func[var[1].index, b];
        if typeof(bfunc[1]) <: Number
            bid_select *= "        return "*string(bfunc[1])*";\n"
        elseif typeof(bfunc[1]) == GenFunction
            if config.dimension == 1
                bid_select *= "        return finch::"*string(bfunc[1].name)*"(x[0], 0.0, 0.0, 0.0, nid);\n"
            elseif config.dimension == 2
                bid_select *= "        return finch::"*string(bfunc[1].name)*"(x[0], x[1], 0.0, 0.0, nid);\n"
            else
                bid_select *= "        return finch::"*string(bfunc[1].name)*"(x[0], x[1], x[2], 0.0, nid);\n"
            end
            
        else
            println("unexpected BC value: "*string(bfunc[1]));
        end
        
        bid_select *= "    }else ";
        bctype_select *= "    }else ";
    end
    bid_select *= 
"    {
        return $default_bc ; // default value for missing BC
    }";
    bctype_select *= 
"    {
        return 0; // default value for missing BC
    }";
    
    content = """
/*
Functions for evaluating boundary conditions.
*/
#pragma once
#include <string>
#include "finch_functions.hpp"

namespace finch{
    double evaluate_bc(const double* x, const int bid, const unsigned long nid);
    int get_bc_type(const int bid);
}    
""";
    
    println(hfile, content)
    
    content = """
/*
Functions for evaluating boundary conditions.
*/
#include <iostream>
#include <sstream>
#include <stdio.h>
#include "finch_bdry.hpp"

double finch::evaluate_bc(const double* x, const int bid, const unsigned long nid){
$bid_select
}

// BC types: 0=NO_BC, 1=DIRICHLET, 2=NEUMANN, 3=ROBIN, 4=FLUX
int finch::get_bc_type(const int bid){
$bctype_select
}
""";
    
    println(file, content)
end

# Output code to be directly included in the main file.
function cpp_output_file(var)
    project_name = finch_state.project_name;
    file = add_generated_file("finch_output.cpp", dir="src");
    
    point_data_part = """
    
        std::vector<PetscInt> point_indices (dofs_per_node);
        std::vector<PetscScalar> point_values (dofs_per_node);
    """;
    cell_data_part = "";
    pvtu_includes = "";
    
    if typeof(var) <: Array
        vararray = var;
    else
        vararray = [var];
    end
    
    dof_offset = 0;
    for vi=1:length(vararray)
        comps = length(vararray[vi].symvar);
        compsstr = string(comps);
        vname = string(vararray[vi].symbol);
        offsetstr = string(dof_offset);
        if vararray[vi].location == CELL
            #TODO
        else
            point_data_part *= """
        
        dof_offset = $offsetstr;
        comps = $compsstr;
        vtufile << "        <DataArray type=\\"Float64\\" Name=\\"$vname\\" NumberOfComponents=\\"$compsstr\\" format=\\"ascii\\">\\n";
        for(int ni=0; ni<num_points; ni++){
            for(int d=0; d<comps; d++){
                vtufile << solution_values[ni * comps + dof_offset] << " ";
            }
            vtufile << "\\n";
        }
        vtufile << "        </DataArray>\\n";
        """
        pvtu_includes *= """    vtufile << "        <PDataArray type=\\"Float64\\" Name=\\"$vname\\" NumberOfComponents=\\"$compsstr\\"/>\\n\";"""
        pvtu_includes *= "\n";
        end
        
        dof_offset += comps;
    end
    
    content = 
"""
// Write the vtu file
std::ofstream vtufile;
std::string vtufilename = "$project_name";
vtufilename = vtufilename + "_" + std::to_string(size) + "_" + std::to_string(rank) + ".vtu";
vtufile.open(vtufilename, std::fstream::out);

unsigned long num_points = mesh.nnodes_local;
unsigned long num_cells = mesh.nel_local;

// If the mesh is partitioned some borrowed node data will not be available.
// Make a new data vector and scatter all of the local node data here.
Vec expanded_solution;
VecScatter scatter_ctx; 
IS from_IS, to_IS;
PetscScalar *solution_values; // This array will hold the values of all local nodes in local order
int *from_indices = new int[num_points];
int *to_indices = new int[num_points];
for(int ni=0; ni<num_points; ni++){
    from_indices[ni] = mesh.partition2global_n[ni];
    to_indices[ni] = ni;
}
VecCreateSeq(PETSC_COMM_SELF, num_points, &expanded_solution);
ISCreateGeneral(PETSC_COMM_SELF, num_points, from_indices, PETSC_COPY_VALUES, &from_IS);
ISCreateGeneral(PETSC_COMM_SELF, num_points, to_indices, PETSC_COPY_VALUES, &to_IS);

VecScatterCreate(solution, from_IS, expanded_solution, to_IS, &scatter_ctx);
VecScatterBegin(scatter_ctx, solution, expanded_solution, INSERT_VALUES, SCATTER_FORWARD);
VecScatterEnd(scatter_ctx, solution, expanded_solution, INSERT_VALUES, SCATTER_FORWARD);
VecGetArray(expanded_solution, &solution_values);

ISDestroy(&from_IS);
ISDestroy(&to_IS);
VecScatterDestroy(&scatter_ctx);

int tmp = 0;
unsigned long offset = 0;
int dof_offset, comps;

// cell types: ??
// int nodes_per_element = mesh.vertices_per_element[0];
int cell_type;
if(mesh.dimension == 1){
    cell_type = 3;
}else if(mesh.dimension == 2){
    if(nodes_per_element == 3){
        cell_type = 5;
    }else if(nodes_per_element == 4){
        cell_type = 9;
    }
}else{
    if(nodes_per_element == 4){
        cell_type = 10;
    }else if(nodes_per_element == 8){
        cell_type = 12;
    }
}

// header
vtufile << "<?xml version=\\"1.0\\" encoding=\\"utf-8\\"?>\\n";
vtufile << "<VTKFile type=\\"UnstructuredGrid\\" version=\\"1.0\\" byte_order=\\"LittleEndian\\">\\n";
vtufile << "  <UnstructuredGrid>\\n";
vtufile << "    <Piece NumberOfPoints=\\""+std::to_string(num_points)+"\\" NumberOfCells=\\""+std::to_string(num_cells)+"\\">\\n";

// Points
vtufile << "      <Points>\\n";
vtufile << "        <DataArray type=\\"Float64\\" Name=\\"Points\\" NumberOfComponents=\\"3\\" format=\\"ascii\\">\\n";

for(int ni=0; ni<num_points; ni++){
    vtufile << "          ";
    vtufile << mesh.allnodes[ni*mesh.dimension] << " ";
    vtufile << ((mesh.dimension > 1) ? mesh.allnodes[ni*mesh.dimension+1] : 0.0) << " ";
    vtufile << ((mesh.dimension > 2) ? mesh.allnodes[ni*mesh.dimension+2] : 0.0);
    vtufile << "\\n";
}

vtufile << "        </DataArray>\\n";
vtufile << "      </Points>\\n";

// Cells
vtufile << "      <Cells>\\n";
vtufile << "        <DataArray type=\\"Int32\\" Name=\\"connectivity\\" format=\\"ascii\\">\\n";

for(int ci=0; ci<num_cells; ci++){
    vtufile << "          ";
    for(int ni=0; ni<nodes_per_element; ni++){
        vtufile << mesh.glbvertex[ci][ni] << " ";
    }
    vtufile << "\\n";
}
vtufile << "        </DataArray>\\n";

vtufile << "        <DataArray type=\\"Int32\\" Name=\\"offsets\\" format=\\"ascii\\">\\n";
tmp = 0;
offset = 0;
for(int ci=0; ci<num_cells; ci++){
    if(tmp == 0){ vtufile << "          "; }
    offset += nodes_per_element;
    vtufile << offset << " ";
    tmp += 1;
    if(tmp == 20 || ci == num_cells-1){ 
        vtufile << "\\n"; 
        tmp = 0; 
    }
}
vtufile << "        </DataArray>\\n";

vtufile << "        <DataArray type=\\"UInt8\\" Name=\\"types\\" format=\\"ascii\\">\\n";
tmp = 0;
for(int ci=0; ci<num_cells; ci++){
    if(tmp == 0){ vtufile << "          "; }
    vtufile << cell_type << " ";
    tmp += 1;
    if(tmp == 20 || ci == num_cells-1){ 
        vtufile << "\\n"; 
        tmp = 0; 
    }
}
vtufile << "        </DataArray>\\n";
vtufile << "      </Cells>\\n";

// Point data
vtufile << "      <PointData>\\n";

"""*
point_data_part *
"""

vtufile << "      </PointData>\\n";

// Cell data
vtufile << "      <CellData>\\n";

"""*
cell_data_part *
"""

vtufile << "      </CellData>\\n";

vtufile << "    </Piece>\\n";
vtufile << "  </UnstructuredGrid>\\n";
vtufile << "</VTKFile>\\n";

vtufile.close();

// If the mesh is partitioned also write a .pvtu file to tie them together
if(size > 1 && rank == 0){
    vtufilename = \"$project_name\";
    vtufilename = vtufilename + \"_\" + std::to_string(size) + \".pvtu\";
    vtufile.open(vtufilename, std::fstream::out);
    
    vtufile << \"<?xml version=\\"1.0\\" encoding=\\"utf-8\\"?>\\n\";
    vtufile << \"<VTKFile type=\\"PUnstructuredGrid\\" version=\\"1.0\\" byte_order=\\"LittleEndian\\">\\n\";
    vtufile << \"<PUnstructuredGrid GhostLevel=\\"0\\">\\n\";
    for(int i=0; i<size; i++){
        vtufile << \"    <Piece Source=\\"$project_name\" << \"_\" << size << \"_\" << i << \".vtu\\"/>\\n\";
    }
    vtufile << \"    <PPoints>\\n\";
    vtufile << \"        <PDataArray type=\\"Float64\\" Name=\\"Points\\" NumberOfComponents=\\"3\\"/>\\n\";
    vtufile << \"    </PPoints>\\n\";
    vtufile << \"    <PPointData>\\n\";
""" *
pvtu_includes *
"""
    vtufile << \"    </PPointData>\\n\";
    vtufile << \"</PUnstructuredGrid>\\n\";
    vtufile << \"</VTKFile>\\n\";
    
    vtufile.close();
}
""";
    
    println(file, content);
end

# Solve code to be directly included in the main file.
function cpp_solve_file(var)
    project_name = finch_state.project_name;
    file = add_generated_file("finch_solve.cpp", dir="src");
    
    content = 
"""
namespace finch {
    PetscErrorCode solve(Mat lhs, Vec rhs, Vec out) {
        // Krylov solver (PETSc)
        KSP ksp;
        // preconditioner (PETSc)
        PC  pc;
        
        // Set up the KSP.
        PetscCall( KSPCreate(MPI_COMM_WORLD, &ksp) );
        PetscCall( KSPSetType(ksp, finch_KSP_TYPE) );
        PetscCall( KSPSetFromOptions(ksp) );
        PetscCall( KSPSetTolerances(ksp, finch_KSP_RTOL, finch_KSP_ATOL, finch_KSP_DTOL, finch_KSP_MAXITERS) );
        
        // Set the matrix in the KSP
        PetscCall( KSPSetOperators(ksp, lhs, lhs) );
        
        // Set up the preconditioner
        PetscCall( KSPGetPC(ksp,&pc) );
        PetscCall( PCSetType(pc, finch_PC_TYPE) );
        PetscCall( PCSetFromOptions(pc) );
        
        // Solve the system
        PetscCall( KSPSolve(ksp, rhs, out) );
        
        // clean up KSP
        PetscCall( KSPDestroy(&ksp) );
        return 0;
    }
}
""";
    
    println(file, content);
end

# Utilities such as build_derivative_matrix
function cpp_utils_file(var, config)
    project_name = finch_state.project_name;
    file = add_generated_file("finch_utils.cpp", dir="src");
    
    # build_derivative_matrix is called like:
    # build_derivative_matrix(refel, geometric_factors, 1, eid, 0, RQ1);
    # build_derivative_matrix(refel, geometric_factors, 2, eid, 0, RQ2);
    # build_derivative_matrix(refel, geometric_factors, 3, eid, 0, RQ3);
    #
    # build_derivative_matrix(Refel refel, GeometricFactors geometric_factors, int direction, unsigned int eid, int type, double* RQn)
    
    # If constant jacobian, the row index for r_x etc. will just be 0
    rowind = "row";
    if length(finch_state.geo_factors.J[1].rx) == 1 # constant jacobian
        rowind = "0";
    end
    
    if config.dimension == 1
        derivmat = 
"""
        // Build derivative matrices:
        // RQn are quadrature matrices for the derivatives of the basis functions
        // with Jacobian factors. They are made like this.
        // |RQ1|   | rx || Qx |
        
        double* d_r;
        int idx;
        if(type == 0){
            d_r = refel.Qr;
        }else{
            d_r = refel.Ddr;
        }
        int nqnodes = refel.Nqp;
        int nnodes = refel.Np;
        
        for(int row=0; row<nqnodes; row++){
            for(int col=0; col<nnodes; col++){
                idx = row*nnodes + col;
                RQn[idx] = geometric_factors.rx[eid][$rowind] * d_r[idx];
            }
        }
"""
    elseif config.dimension == 2
        derivmat = 
"""
        // Build derivative matrices:
        // RQn are quadrature matrices for the derivatives of the basis functions
        // with Jacobian factors. They are made like this.
        // |RQ1|   | rx sx|| Qx |
        // |RQ2| = | ry sy|| Qy |
        
        double* r_xyz;
        double* s_xyz;
        double* d_r;
        double* d_s;
        int idx;
        if(direction == 1){
            r_xy = geometric_factors.rx[eid];
            s_xy = geometric_factors.sx[eid];
        }else{
            r_xy = geometric_factors.ry[eid];
            s_xy = geometric_factors.sy[eid];
        }
        if(type == 0){
            d_r = refel.Qr;
            d_s = refel.Qs;
        }else{
            d_r = refel.Ddr;
            d_s = refel.Dds;
        }
        int nqnodes = refel.Nqp;
        int nnodes = refel.Np;
        
        for(int row=0; row<nqnodes; row++){
            for(int col=0; col<nnodes; col++){
                idx = row*nnodes + col;
                RQn[idx] = r_xy[$rowind] * d_r[idx] + s_xy[$rowind] * d_s[idx];
            }
        }
"""
    else
        derivmat = 
"""
        // Build derivative matrices:
        // RQn are quadrature matrices for the derivatives of the basis functions
        // with Jacobian factors. They are made like this.
        // |RQ1|   | rx sx tx || Qx |
        // |RQ2| = | ry sy ty || Qy |
        // |RQ3|   | rz sz tz || Qz |
        
        double* r_xyz;
        double* s_xyz;
        double* t_xyz;
        double* d_r;
        double* d_s;
        double* d_t;
        int idx;
        if(direction == 1){
            r_xyz = geometric_factors.rx[eid];
            s_xyz = geometric_factors.sx[eid];
            t_xyz = geometric_factors.tx[eid];
        }else if(direction == 2){
            r_xyz = geometric_factors.ry[eid];
            s_xyz = geometric_factors.sy[eid];
            t_xyz = geometric_factors.ty[eid];
        }else{
            r_xyz = geometric_factors.rz[eid];
            s_xyz = geometric_factors.sz[eid];
            t_xyz = geometric_factors.tz[eid];
        }
        if(type == 0){
            d_r = refel.Qr;
            d_s = refel.Qs;
            d_t = refel.Qt;
        }else{
            d_r = refel.Ddr;
            d_s = refel.Dds;
            d_t = refel.Ddt;
        }
        int nqnodes = refel.Nqp;
        int nnodes = refel.Np;
        
        for(int row=0; row<nqnodes; row++){
            for(int col=0; col<nnodes; col++){
                idx = row*nnodes + col;
                RQn[idx] = r_xyz[$rowind] * d_r[idx] + s_xyz[$rowind] * d_s[idx] + t_xyz[$rowind] * d_t[idx];
            }
        }
"""
    end
    
    content = 
"""
// Build the derivative matrix RQn or RDn depending on type value.
// RQn computes the derivative from nodal values and interpolates to quadrature nodes.
// RDn computes the derivative from nodal values without interpolating.
// The n represents (x,y,z) depending on direction value.
void build_derivative_matrix(finch::Refel refel, finch::GeometricFactors geometric_factors, int direction, unsigned int eid, int type, double* RQn){
    $derivmat
}
""";
    
    println(file, content);
end

# Writes files to build the project and a readme.
function cpp_build_files(params)
    project_name = finch_state.project_name;
    cmakefile = add_generated_file("CMakeLists.txt", make_header_text=false);
    readmefile = add_generated_file("readme.txt", make_header_text=false);
    
    config = finch_state.config;
    param_val = Dict(true=>"ON", false=>"OFF");
    
    content = 
"""    
cmake_minimum_required(VERSION 3.7)

# set name of project, which can be used by \${PROJECT_NAME}
project($(project_name))

# SET THESE for your PETSc installation if cmake can't find them automatically
# set (PETSC_DIR /path/to/petsc)
# set (PETSC_ARCH arch-linux-c-debug)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)

# find and load settings from external project, if not found then stop with error message
find_package(MPI REQUIRED)
# find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


set(LAPACK_LINKER_FLAGS -llapacke -llapack -lblas -lgfortran -lquadmath)
set(LAPACKE_DIR \$ENV{LAPACK}/LAPACKE)
set(LINK_FLAGS "\${LINK_FLAGS} \${LAPACK_LINKER_FLAGS}")
set(LAPACK_LIBRARIES \${LAPACK_LIBRARIES} \${LAPACKE_LIB})
message(STATUS \${LAPACK_LIBRARIES})
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message("\${CMAKE_CXX_COMPILER_ID} compiler detected adding -mkl flag for BLAS LAPACK")
    set (CMAKE_C_FLAGS "\${CMAKE_C_FLAGS} -mkl")
    set (CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} -mkl")
endif()

# when OpenMP is found
# if(OpenMP_FOUND)
#     set (CMAKE_C_FLAGS "\${CMAKE_C_FLAGS} \${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} \${OpenMP_CXX_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "\${CMAKE_EXE_LINKER_FLAGS} \${OpenMP_EXE_LINKER_FLAGS}")
# endif()

if(MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "\${COMPILE_FLAGS} \${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set(LINK_FLAGS "\${LINK_FLAGS} \${MPI_LINK_FLAGS}")
endif()

set(INCLUDE_FILES )

set(SOURCE_FILES
    src/finch_mesh.cpp
    src/finch_geometry.cpp
    src/finch_functions.cpp
    src/finch_bdry.cpp)

# cmake options, which will be visible at ccmake ../
option(BUILD_WITH_PETSC "Build code with the petsc" ON)

# if BUILD_WITH_PETSC ON , #define BUILD_WITH_PETSC
if(BUILD_WITH_PETSC)
    list(APPEND CMAKE_MODULE_PATH "\${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")
    find_package(PETSc REQUIRED)
    add_definitions(-DBUILD_WITH_PETSC)
endif(BUILD_WITH_PETSC)

#set (CMAKE_C_FLAGS "\${CMAKE_C_FLAGS} -qopt-report=5 -qopt-report-phase=vec -qopt-report-file=stdout")
#set (CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} -qopt-report=5 -qopt-report-phase=vec -qopt-report-file=stdout")

set(EIGEN_HEADER_DIR .)

add_executable($(project_name) src/$(project_name).cpp include/$(project_name).hpp \${INCLUDE_FILES} \${SOURCE_FILES})
target_include_directories($(project_name) PUBLIC include)
target_include_directories($(project_name) PRIVATE \${MPI_INCLUDE_PATH})
target_include_directories($(project_name) PRIVATE \${EIGEN_HEADER_DIR})
target_link_libraries($(project_name) \${MPI_LIBRARIES} m)

if(BUILD_WITH_PETSC)
    target_include_directories($(project_name) PUBLIC \${PETSC_INCLUDES})
    target_link_libraries($(project_name) \${PETSC_LIBRARIES})
endif()
""";
    
    println(cmakefile, content);
    
    content = """
Readme for $(project_name) for c++.
Follow these steps to build this project.

1. Go to this directory.
2. Copy the CMake modules folder from the github repo: src/targets/cmake-modules
3. Change the PETSC_DIR and PETSC_ARCH in CMakeLists.txt for your PETSc installation.
4. Create a build directory here and compile with CMake.
   For example:
    \$ mkdir build
    \$ cd build
    \$ ccmake ..   (then follow the CMake procedures: "c c g")
    \$ make
5. Move all of the MeshData, GeoData and RefelData files from this directory 
   into the build directory.
6. Run with something like: mpirun -n """*string(config.num_partitions)*""" ./$(project_name)\n
7. Timing results will be in $(project_name)_timing.txt
   Data output will be in the .vtu/.pvtu files. Open the .pvtu file in Paraview.

NOTES:
- The number of processes must equal the number used when generating the MeshData files in Finch.
  The separate mesh partitions and geodata are suffixed with _pn where n is the partition index.
""";
    println(readmefile, content);
    
end


###########################################################################################################
# Symbolic to code layer generation

  ######      #####     ######    #######       ##           ###     ##    ##  #######  ######
##     ##   ###   ###   ##   ##   ##            ##          ## ##     ##  ##   ##       ##   ##
##          ##     ##   ##    ##  ######        ##         ##   ##     ####    ######   ######
##     ##   ###   ###   ##   ##   ##            ##        #########     ##     ##       ##  ##
  ######      #####     ######    #######       #######  ##       ##    ##     #######  ##   ##

###########################################################################################################

# If needed, build derivative matrices
function cpptarget_build_derivative_matrices()
    config = finch_state.config;
    
    row_ind = "row";
    if length(finch_state.geo_factors.J[1].rx) == 1 # constant jacobian
        row_ind = "0";
    end
    if config.dimension == 1
        todelete = ["RQ1"];
        code = 
"
        // Build derivative matrices:
        // RQn are quadrature matrices for the derivatives of the basis functions
        // with Jacobian factors. They are made like this.
        // |RQ1|   | rx || Qx |

        double *RQ1 = new double[refel.Np * refel.Nqp];
        for(int row=0; row<nqnodes; row++){
            for(int col=0; col<nnodes; col++){
                int idx = row*nnodes + col;
                RQ1[idx] = geometric_factors.rx[eid]["*row_ind*"] * refel.Qr[idx];
            }
        }
";
    elseif config.dimension == 2
        todelete = ["RQ1", "RQ2"];
        code = 
"
        // Build derivative matrices:
        // RQn are quadrature matrices for the derivatives of the basis functions
        // with Jacobian factors. They are made like this.
        // |RQ1|   | rx sx|| Qx |
        // |RQ2| = | ry sy|| Qy |

        double *RQ1 = new double[refel.Np * refel.Nqp];
        double *RQ2 = new double[refel.Np * refel.Nqp];
        for(int row=0; row<nqnodes; row++){
            for(int col=0; col<nnodes; col++){
                int idx = row*nnodes + col;
                RQ1[idx] = geometric_factors.rx[eid]["*row_ind*"] * refel.Qr[idx] + geometric_factors.sx[eid]["*row_ind*"] * refel.Qs[idx];
                RQ2[idx] = geometric_factors.ry[eid]["*row_ind*"] * refel.Qr[idx] + geometric_factors.sy[eid]["*row_ind*"] * refel.Qs[idx];
            }
        }
";
    elseif config.dimension == 3
        todelete = ["RQ1", "RQ2", "RQ3"];
        code = 
"
        // Build derivative matrices:
        // RQn are quadrature matrices for the derivatives of the basis functions
        // with Jacobian factors. They are made like this.
        // |RQ1|   | rx sx tx || Qx |
        // |RQ2| = | ry sy ty || Qy |
        // |RQ3|   | rz sz tz || Qz |
        
        double *RQ1 = new double[refel.Np * refel.Nqp];
        double *RQ2 = new double[refel.Np * refel.Nqp];
        double *RQ3 = new double[refel.Np * refel.Nqp];
        for(int row=0; row<nqnodes; row++){
            for(int col=0; col<nnodes; col++){
                int idx = row*nnodes + col;
                RQ1[idx] = geometric_factors.rx[eid]["*row_ind*"] * refel.Qr[idx] + geometric_factors.sx[eid]["*row_ind*"] * refel.Qs[idx] + geometric_factors.tx[eid]["*row_ind*"] * refel.Qt[idx];
                RQ2[idx] = geometric_factors.ry[eid]["*row_ind*"] * refel.Qr[idx] + geometric_factors.sy[eid]["*row_ind*"] * refel.Qs[idx] + geometric_factors.ty[eid]["*row_ind*"] * refel.Qt[idx];
                RQ3[idx] = geometric_factors.rz[eid]["*row_ind*"] * refel.Qr[idx] + geometric_factors.sz[eid]["*row_ind*"] * refel.Qs[idx] + geometric_factors.tz[eid]["*row_ind*"] * refel.Qt[idx];
            }
        }
";
    end
    return (code, todelete);
end

# Allocate, compute, or fetch all needed values
function cpptarget_prepare_needed_values(entities, var, lorr, vors)
    to_delete = []; # arrays allocated with new that need deletion
    used_names = []; # Don't want to duplicate variables
    code = "";
    coef_loop = "";
    coef_interp = "";
    coef_interp_names = [];
    
    # Determine if derivative matrices will be required
    need_derivs = false;
    for i=1:length(entities)
        if length(entities[i].derivs) > 0
            need_derivs = true;
            break;
        end
    end
    if need_derivs
        (dcode, todel) = cpptarget_build_derivative_matrices();
        code *= dcode;
        append!(to_delete, todel);
    end
    
    for i=1:length(entities)
        cname = CodeGenerator.make_entity_name(entities[i]);
        # Make sure it hasn't been prepared yet.
        # Since entities with the same name should have the same value, just skip.
        for used in used_names
            if used == cname
                continue;
            end
        end
        push!(used_names, cname);
        if CodeGenerator.is_test_function(entities[i])
            # Assign it a transpose quadrature matrix
            if length(entities[i].derivs) > 0
                xyzchar = ["x","y","z"];
                for di=1:length(entities[i].derivs)
                    code *= "        double * " * cname * " = RQ"*string(entities[i].derivs[di])*"; // d/d"*xyzchar[entities[i].derivs[di]]*" of test function\n";
                end
            else
                code *= "        double * " * cname * " = refel.Q; // test function.\\n";
            end
        elseif CodeGenerator.is_unknown_var(entities[i], var) && lorr == LHS
            if length(entities[i].derivs) > 0
                xyzchar = ["x","y","z"];
                for di=1:length(entities[i].derivs)
                    code *= "        double * " * cname * " = RQ"*string(entities[i].derivs[di])*"; // d/d"*xyzchar[entities[i].derivs[di]]*" of trial function\n";
                end
            else
                code *= "        double * " * cname * " = refel.Q; // trial function.\\n";
            end
        else
            # Is coefficient(number or function) or variable(array)?
            (ctype, cval) = CodeGenerator.get_coef_val(entities[i]);
            if ctype == -1
                # It was a special symbol like dt. These should already be available.
            elseif ctype == 0
                # It was a number. Do nothing.
            elseif ctype == 1 # a constant wrapped in a coefficient
                # This generates something like: double coef_k_i = 4;
                if length(entities[i].derivs) > 0
                    code *= "        double " * cname * " = 0.0; // NOTE: derivative applied to constant coefficient = 0\n";
                else
                    code *= "        double " * cname * " = " * string(cval) * ";\n";
                end
                
            elseif ctype == 2 # a coefficient function
                if vors == "volume"
                    push!(to_delete, "NODAL"*cname);
                    push!(to_delete, cname);
                    code *= "        double* NODAL"*cname*" = new double[nnodes]; // value at nodes\n"; # at nodes
                    code *= "        double* "*cname*" = new double[nqnodes];     // value at quadrature points\n"; # at quadrature points
                    push!(coef_interp_names, cname);
                    coef_loop *= "            NODAL"*cname*"[ni] = finch::genfunction_"*string(cval)*"(x, y, z, t, nid);\n";
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                        # for(int row=0; row<nqnodes; row++){
                        #     Qcoef_f_1[row] = 0.0;
                        #     for(int col=0; col<nnodes; col++){
                        #         Qcoef_f_1[row] += refel.Q[row*nnodes + col] * coef_f_1[col];
                        #     }
                        # }
                    if length(entities[i].derivs) > 0
                        coef_interp *= "                "*cname * "[row] += RQ"*string(entities[i].derivs[di])*"[row*nnodes + col] * NODAL"*cname*"[col];\n";
                    else
                        coef_interp *= "                "*cname * "[row] += refel.Q[row*nnodes + col] * NODAL"*cname*"[col];\n";
                    end
                else
                    #TODO surface
                end
                
            elseif ctype == 3 # a known variable value
                # This should only occur for time dependent problems
                # Use the values from the previous time step
                if vors == "volume"
                    push!(to_delete, cname);
                    push!(to_delete, "NODAL"*cname);
                    code *= "        double* NODAL"*cname*" = new double[nnodes]; // value at nodes\n"; # at nodes
                    code *= "        double* "*cname*" = new double[nqnodes];     // value at quadrature points\n"; # at quadrature points
                    push!(coef_interp_names, cname);
                    coef_loop *= "            NODAL"*cname*"[ni] = 0.0; // TODO extract variable values\n";
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        coef_interp *= "                "*cname * "[row] += RQ"*string(entities[i].derivs[di])*"[row*nnodes + col] * NODAL"*cname*"[col];\n";
                    else
                        coef_interp *= "                "*cname * "[row] += refel.Q[row*nnodes + col] * NODAL"*cname*"[col];\n";
                    end
                else
                    #TODO surface
                end
                
            end
        end # if coefficient
    end # entity loop
    
    # Loop to compute coefficients at nodes
    if length(coef_loop) > 2
        code *= 
"
        // Loop to compute coefficients at nodes.
        for(int ni=0; ni<nnodes; ni++){
            int nid = (mesh.loc2glb[eid][ni]-1);
            double x = mesh.allnodes[nid*dim];
            double y = mesh.allnodes[nid*dim+1];
            double z = mesh.allnodes[nid*dim+2];
            double t = 0.0;
            
"*coef_loop*"
        }
";
    
        # Interpolate
        zerocode = "";
        for cn in coef_interp_names
            zerocode *= "            "*cn*"[row] = 0.0;\n";
        end
        code *= 
"
        // Interpolate at quadrature points and apply derivatives if needed.
        for(int row=0; row<nqnodes; row++){
"*zerocode*"
            for(int col=0; col<nnodes; col++){
"*coef_interp*"
            }
        }
";
    end
    
    return (code, to_delete);
end

function cpptarget_make_elemental_computation(terms, var, lorr, vors)
    # Here is where I make some assumption about the form of the expression.
    # Since it was expanded by the parser it should look like a series of terms: t1 + t2 + t3...
    # Where each term is multiplied by one test function component, and if LHS, involves one unknown component.
    # The submatrix modified by a term is determined by these, so go through the terms and divide them
    # into their submatrix expressions. 
    # Each term will look something like 
    # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
    # RHS: test_part * (weight_part .* coef_part)
    code = "";
    
    dofs_per_node = 0;
    dof_names = [];
    if typeof(var) <:Array
        dofs_per_node = length(var[1].symvar);
        for i=1:length(var.symvar)
            push!(dof_names, string(var[1].symvar[i]));
        end
        for vi=2:length(var)
            dofs_per_node = dofs_per_node + length(var[vi].symvar);
            for i=1:length(var.symvar)
                push!(dof_names, string(var[vi].symvar[i]));
            end
        end
    else
        dofs_per_node = length(var.symvar);
        for i=1:length(var.symvar)
            push!(dof_names, string(var.symvar[i]));
        end
    end
    
    # For constant jacobian, this step can be skipped
    detj_update = "";
    if length(finch_state.geo_factors.J[1].rx) > 1 # not constant jacobian
        if lorr==LHS
            if vors == "volume"
                detj_update = "detj = geometric_factors.detJ[eid][i];";
            else
               #TODO 
            end
        else
            if vors == "volume"
                detj_update = "detj = geometric_factors.detJ[eid][col];";
            else
               #TODO 
            end
        end
    end
    
    # Separate the factors of each term into test, trial, coef and form the calculation
    if dofs_per_node > 1
        # # Submatrices or subvectors for each component
        # if lorr == LHS
        #     submatrices = Array{String, 2}(undef, dofs_per_node, dofs_per_node);
        # else # RHS
        #     submatrices = Array{String, 1}(undef, dofs_per_node);
        # end
        # for smi=1:length(submatrices)
        #     submatrices[smi] = "";
        # end
        
        # if typeof(var) <: Array
        #     for vi=1:length(var) # variables
        #         # Process the terms for this variable
        #         for ci=1:length(terms[vi]) # components
        #             for i=1:length(terms[vi][ci])
        #                 (term_result, test_ind, trial_ind) = generate_term_calculation_cg_julia(terms[vi][ci][i], var, lorr, vors);
                        
        #                 # println(terms)
        #                 # println(terms[vi])
        #                 # println(terms[vi][ci])
        #                 # println(terms[vi][ci][i])
        #                 # println(term_result * " : "*string(test_ind)*", "*string(trial_ind))
                        
        #                 # Find the appropriate submatrix for this term
        #                 submati = offset_ind[vi] + test_ind;
        #                 submatj = trial_ind;
        #                 if lorr == LHS
        #                     submat_ind = submati + dofs_per_node * (submatj-1);
        #                 else
        #                     submat_ind = submati;
        #                 end
                        
                        
        #                 if length(submatrices[submat_ind]) > 1
        #                     submatrices[submat_ind] *= " .+ " * term_result;
        #                 else
        #                     submatrices[submat_ind] = term_result;
        #                 end
        #             end
        #         end
                
        #     end # vi
            
        # else # only one variable
        #     # Process the terms for this variable
        #     for ci=1:length(terms) # components
        #         for i=1:length(terms[ci])
        #             (term_result, test_ind, trial_ind) = generate_term_calculation_cg_julia(terms[ci][i], var, lorr, vors);
                    
        #             # Find the appropriate submatrix for this term
        #             if lorr == LHS
        #                 submat_ind = test_ind + dofs_per_node * (trial_ind-1);
        #             else
        #                 submat_ind = test_ind;
        #             end
                    
        #             if length(submatrices[submat_ind]) > 1
        #                 submatrices[submat_ind] *= " + " * term_result;
        #             else
        #                 submatrices[submat_ind] = term_result;
        #             end
        #         end
        #     end
            
        # end
        
        # # Put the submatrices together into element_matrix or element_vector
        # if lorr == LHS
        #     for emi=1:dofs_per_node
        #         for emj=1:dofs_per_node
        #             if length(submatrices[emi, emj]) > 1
        #                 rangei = "("*string(emi-1)*"*refel.Np + 1):("*string(emi)*"*refel.Np)";
        #                 rangej = "("*string(emj-1)*"*refel.Np + 1):("*string(emj)*"*refel.Np)";
        #                 code *= "element_matrix["*rangei*", "*rangej*"] = " * submatrices[emi,emj] * "\n";
        #             end
        #         end
        #     end
        #     code *= "return element_matrix;\n"
            
        # else # RHS
        #     for emi=1:dofs_per_node
        #         if length(submatrices[emi]) > 1
        #             rangei = "("*string(emi-1)*"*refel.Np + 1):("*string(emi)*"*refel.Np)";
        #             code *= "element_vector["*rangei*"] = " * submatrices[emi] * "\n";
        #         end
        #     end
        #     code *= "return element_vector;\n"
        # end
        
    else # one dof
        terms = terms[1];
        
        # Zero the output
        if lorr == LHS
            code *= 
        "
        // zero the output
        for(int row=0; row<nnodes; row++){
            for(int col=0; col<nnodes; col++){
                ke[row*nnodes + col] = 0.0;
            }
        }
        ";
        else
            code *= 
        "
        // zero the output
        for(int row=0; row<nnodes; row++){
            be[row] = 0.0;
        }
        ";
        end
        
        #process each term
        result = "";
        for i=1:length(terms)
            (term_result, test_ind, trial_ind) = cpptarget_generate_term_calculation(terms[i], var, lorr, vors);
            
            if i > 1
                result *= " + " * term_result;
            else
                result = term_result;
            end
        end
        
        if lorr == LHS
            code *= 
"
        // compute the elemental matrix
        for(int i=0; i<nqnodes; i++){
            for(int row=0; row<nnodes; row++){
                for(int col=0; col<nnodes; col++){
                    "*detj_update*"
                    ke[row*nnodes + col] += "* result *";
                }
            }
        }
";
        else
            code *= 
"
        // compute the elemental vector
        for(int col=0; col<nqnodes; col++){
            for(int row=0; row<nnodes; row++){
                "*detj_update*"
                be[row] += "* result *";
            }
        }
";
        end
    end
    
    return code;
end

function cpptarget_generate_term_calculation(term, var, lorr, vors)
    
    result = "";
    test_ind = 1;
    trial_ind = 1;
    
    if lorr == LHS
        (test_part, trial_part, coef_part, test_ind, trial_ind) = CodeGenerator.separate_factors(term, var);
        # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
        # RQ1[i*nqnodes + row] * refel.wg[i] * detj * RQ1[i*nqnodes + col]
        if !(coef_part === nothing)
            result = string(cpp_change_math_ops(CodeGenerator.replace_entities_with_symbols(test_part))) * "[i*nqnodes + row] * refel.wg[i] * detj * (" * 
                    string(cpp_change_math_ops(CodeGenerator.replace_entities_with_symbols(coef_part, index="i"))) * ") * " * 
                    string(cpp_change_math_ops(CodeGenerator.replace_entities_with_symbols(trial_part))) * "[i*nqnodes + col]";
        else # no coef_part
            result = string(cpp_change_math_ops(CodeGenerator.replace_entities_with_symbols(test_part))) * "[i*nqnodes + row] * refel.wg[i] * detj * " * 
                    string(cpp_change_math_ops(CodeGenerator.replace_entities_with_symbols(trial_part))) * "[i*nqnodes + col]";
        end
    else
        (test_part, trial_part, coef_part, test_ind, trial_ind) = CodeGenerator.separate_factors(term);
        # RHS: test_part * (weight_part .* coef_part)
        if !(coef_part === nothing)
            result = string(cpp_change_math_ops(CodeGenerator.replace_entities_with_symbols(test_part))) * "[col*nnodes + row] * refel.wg[col] * detj * " * 
                    string(cpp_change_math_ops(CodeGenerator.replace_entities_with_symbols(coef_part, index="col"))) * "";
        else
            result = string(cpp_change_math_ops(CodeGenerator.replace_entities_with_symbols(test_part))) * " * refel.wg[col] * detj";
        end
    end
    
    return (result, test_ind, trial_ind);
end

###########################################################################################################
# Static files that will be written as is

 #####   ######      ###     ######  ######   ######       #######  ######  ##        #######   #####
###   #    ##       ## ##      ##      ##    ##    ##      ##         ##    ##        ##       ###   #
  ###      ##      ##   ##     ##      ##    ##            #######    ##    ##        ######     ### 
#   ###    ##     #########    ##      ##    ##    ##      ##         ##    ##        ##       #   ###
 #####     ##    ##       ##   ##    ######   ######       ##       ######  ########  #######   #####

###########################################################################################################

function cpp_write_static_files()
    # Files include:
    # - finch_mesh.cpp / finch_mesh.hpp
    # - finch_geometry.cpp / finch_geometry.hpp
    
meshcpp = """
/*
Mesh and reference element data structures that are imported from files.
Finch exports binary files containing the mesh(grid_data in Finch), and refel.

For partitioned meshes, each partition is placed in a separate file with "_pn" appended for partition n.
*/
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include "finch_mesh.hpp"

// File IO check macro
#define read_check(s) if(!s){ std::cout << "Error: Failed to read some data from the file. Exiting.\\n"; exit(0); }

finch::Mesh::Mesh(){
    dimension = 0;
    nel_global = 0UL;
    nel_local = 0UL;
    nnodes_global = 0UL;
    nnodes_local = 0UL;
    allnodes = nullptr;
    num_bids = 0;
    nodes_per_bid = nullptr;
    faces_per_bid = nullptr;
    bdry = nullptr;
    bdry_face = nullptr;
    bdry_normal = nullptr;
    bids = nullptr;
    nodebid = nullptr;
    nodes_per_element = nullptr;
    vertices_per_element = nullptr;
    faces_per_element = nullptr;
    loc2glb = nullptr;
    glbvertex = nullptr;
    num_faces = 0;
    face2glb = nullptr;
    element2face = nullptr;
    face2element = nullptr;
    face_normal = nullptr;
    face_refel_index = nullptr;
    face_bid = nullptr;
    is_subgrid = false;
    elemental_order = nullptr;
    nel_owned = 0UL;
    nel_ghost = 0UL;
    nface_owned = 0UL;
    nface_ghost = 0UL;
    nnodes_borrowed = 0UL;
    element_owner = nullptr;
    node_owner = nullptr;
    partition2global_e = nullptr;
    partition2global_n = nullptr;
    global_bdry_index = nullptr;
    num_neighbor_partitions = 0;
    neighboring_partitions = nullptr;
    ghost_counts = nullptr;
    ghost_index = nullptr;
    
    is_bdry_node = nullptr;
}

finch::Mesh::~Mesh(){
    // There is an issue with deletion. Skip the deconstructor for now.
    return;
    if (nnodes_local < 1){
        return;
    }
    delete [] allnodes;
    delete [] nodes_per_bid;
    delete [] faces_per_bid;
    for (unsigned long i=0; i<num_bids; i++){
        delete [] bdry[i];
    }
    delete [] bdry;
    for (unsigned long i=0; i<num_bids; i++){
        delete [] bdry_face[i];
    }
    delete [] bdry_face;
    for (unsigned long i=0; i<num_bids; i++){
        delete [] bdry_normal[i];
    }
    delete [] bdry_normal;
    delete [] bids;
    delete [] nodebid;
    delete [] nodes_per_element;
    delete [] vertices_per_element;
    delete [] faces_per_element;
    for (unsigned long i=0; i<nel_local; i++){
        delete [] loc2glb[i];
    }
    delete [] loc2glb;
    for (unsigned long i=0; i<nel_local; i++){
        delete [] glbvertex[i];
    }
    delete [] glbvertex;
    for (unsigned long i=0; i<num_faces; i++){
        delete [] face2glb[i];
    }
    delete [] face2glb;
    for (unsigned long i=0; i<nel_local; i++){
        delete [] element2face[i];
    }
    delete [] element2face;
    delete [] face2element;
    delete [] face_normal;
    delete [] face_refel_index;
    delete [] face_bid;
    delete [] elemental_order;
    
    if(num_neighbor_partitions > 0){
        
        delete [] partition2global_e;
        if(nnodes_borrowed > 0){// FE only
            delete [] partition2global_n;
            delete [] node_owner;
            delete [] global_bdry_index;
        }
        if(nel_ghost > 0){// FV only
            delete [] element_owner;
            delete [] neighboring_partitions;
            delete [] ghost_counts;
            for (unsigned long i=0; i<num_neighbor_partitions; i++){
                delete [] ghost_index[i];
            }
            delete [] ghost_index;
        }
        
    }
    
    delete [] is_bdry_node;
}

void finch::Mesh::import_mesh(std::string filename, int num_partitions, int my_partition){
    std::ifstream file;
    if(num_partitions > 1){
        // Partitioned meshes are stored in separate files for each partition.
        // Their names have a "_pn" on the end of the filename. (n=partition number)
        std::stringstream newname;
        newname << filename << "_p" << my_partition;
        file.open(newname.str(), std::ios::binary);
        if(!file){
            std::cout << "Error: Partition "<<my_partition<<" couldn't open Finch mesh file: " << newname.str() << ". Exiting.\\n";
            exit(0);
        }
        
    }else{
        file.open(filename, std::ios::binary);
        if(!file){
            std::cout << "Error: couldn't open Finch mesh file: " << filename << ". Exiting.\\n";
            exit(0);
        }
    }
    
    // These will temporarily hold the read values
    char in8[8];
    char in1[1];
    unsigned long count, size;
    
    read_check(file.read(in8, 8));
    dimension = ((int64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    nel_local = ((uint64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    nnodes_local = ((uint64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    nel_global = ((uint64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    nnodes_global = ((uint64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    int max_nodes_per_element = ((int64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    int max_vertices_per_element = ((int64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    int max_faces_per_element = ((int64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    num_faces = ((uint64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    int max_nodes_per_face = ((int64_t*)in8)[0];
    
    // Nodes
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    allnodes = new double[count];
    read_check(file.read((char*)allnodes, count*size));
    
    // Boundary
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    num_bids = count;
    if(num_bids > 0){
        nodes_per_bid = new unsigned long[num_bids];
        faces_per_bid = new unsigned long[num_bids];
        bdry = new unsigned long*[num_bids];
        bdry_face = new unsigned long*[num_bids];
        bdry_normal = new double*[num_bids];
        bids = new int[num_bids];
        
        for(int i=0; i<num_bids; i++){
            read_check(file.read((char*)&count, 8));
            read_check(file.read((char*)&size, 8));
            nodes_per_bid[i] = count;
            bdry[i] = new unsigned long[count];
            read_check(file.read((char*)bdry[i], count*size));
        }
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        for(int i=0; i<num_bids; i++){
            read_check(file.read((char*)&count, 8));
            read_check(file.read((char*)&size, 8));
            faces_per_bid[i] = count;
            bdry_face[i] = new unsigned long[count];
            read_check(file.read((char*)bdry_face[i], count*size));
        }
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        for(int i=0; i<num_bids; i++){
            read_check(file.read((char*)&count, 8));
            read_check(file.read((char*)&size, 8));
            bdry_normal[i] = new double[count];
            read_check(file.read((char*)bdry_normal[i], count*size));
        }
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        read_check(file.read((char*)bids, count*size));
        
    }else{
        // bdry info is all empty, but these zeros need to be skipped
        read_check(file.read((char*)&count, 8)); // bdry_face
        read_check(file.read((char*)&size, 8));
        read_check(file.read((char*)&count, 8)); // bdry_normal
        read_check(file.read((char*)&size, 8));
        read_check(file.read((char*)&count, 8)); // bids
        read_check(file.read((char*)&size, 8));
    }
    // node bids
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    nodebid = new int[count];
    long tmp_nodebid = 0;
    for(unsigned long i=0; i<count; i++){
        read_check(file.read((char*)&tmp_nodebid, 8));
        nodebid[i] = tmp_nodebid;
    }
    
    // Elements
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    // The data has the same number for all elements, but extras may be zeros. TODO resize each loc2glb
    if(count > 0){
        loc2glb = new unsigned long*[nel_local];
        nodes_per_element = new int[nel_local];
        if(nel_local*max_nodes_per_element != count){
            std::cout << "Error: element node count mismatch.\\n";
            std::cout << "loc2glb size " << count << "\\n";
            std::cout << "nel_local " << nel_local << "\\n";
            std::cout << "nodes per element " << max_nodes_per_element << "\\n";
            exit(0);
        }
        for(unsigned long i=0; i<nel_local; i++){
            nodes_per_element[i] = max_nodes_per_element;
            loc2glb[i] = new unsigned long[nodes_per_element[i]];
            read_check(file.read((char*)loc2glb[i], nodes_per_element[i]*size));
        }
    }else{
        // There are no elements?
        std::cout << "Error: elemental node map is empty." << count << "\\n";
        exit(0);
    }
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    // The data has the same number for all elements, but extras may be zeros. TODO resize each glbvertex
    if(count > 0){
        glbvertex = new unsigned long*[nel_local];
        vertices_per_element = new int[nel_local];
        if(nel_local*max_vertices_per_element != count){
            std::cout << "Error: element vertex count mismatch.\\n";
            exit(0);
        }
        for(unsigned long i=0; i<nel_local; i++){
            vertices_per_element[i] = max_vertices_per_element;
            glbvertex[i] = new unsigned long[vertices_per_element[i]];
            read_check(file.read((char*)glbvertex[i], vertices_per_element[i]*size));
        }
    }else{
        // There are no elements?
        std::cout << "Error: elemental vertex map is empty.\\n";
        exit(0);
    }
    
    // Faces
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    if(count > 0){
        face2glb = new unsigned long*[num_faces];
        nodes_per_face = new int[num_faces];
        if(num_faces*max_nodes_per_face != count){
            std::cout << "Error: face node count mismatch.\\n";
            exit(0);
        }
        for(unsigned long i=0; i<num_faces; i++){
            nodes_per_face[i] = max_nodes_per_face ;
            face2glb[i] = new unsigned long[nodes_per_face[i]];
            read_check(file.read((char*)face2glb[i], nodes_per_face[i]*size));
        }
    }else{
        // There are no faces?
        std::cout << "Error: face map is empty.\\n";
        exit(0);
    }
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    if(count > 0){
        element2face = new unsigned long*[nel_local];
        faces_per_element = new int[nel_local];
        for(unsigned long i=0; i<nel_local; i++){
            faces_per_element[i] = max_faces_per_element;
            element2face[i] = new unsigned long[faces_per_element[i]];
            read_check(file.read((char*)element2face[i], faces_per_element[i]*size));
        }
    }else{
        std::cout << "Error: element2face map is empty.\\n";
        exit(0);
    }
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    face2element = new unsigned long[num_faces * 2];
    read_check(file.read((char*)face2element, num_faces*2*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    face_normal = new double[num_faces * dimension];
    read_check(file.read((char*)face_normal, num_faces*dimension*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    face_refel_index = new unsigned long[num_faces*2];
    read_check(file.read((char*)face_refel_index, num_faces*2*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    face_bid = new long[num_faces];
    read_check(file.read((char*)face_bid, num_faces*size));
    
    // For partitioned meshes
    read_check(file.read(in8, 8));
    is_subgrid = ((int64_t*)in8)[0] > 0;
    
    // elemental loop order
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    elemental_order = new unsigned long[count];
    read_check(file.read((char*)elemental_order, count*size));
    
    if(is_subgrid){
        read_check(file.read(in8, 8));
        nel_owned = ((int64_t*)in8)[0];
        read_check(file.read(in8, 8));
        nel_ghost = ((int64_t*)in8)[0];
        read_check(file.read(in8, 8));
        nface_owned = ((int64_t*)in8)[0];
        read_check(file.read(in8, 8));
        nface_ghost = ((int64_t*)in8)[0];
        
        read_check(file.read(in8, 8));
        nnodes_borrowed = ((int64_t*)in8)[0];
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        partition2global_e = new unsigned long[count];
        read_check(file.read((char*)partition2global_e, count*size));
        
        if(nel_ghost == 0){// FE only
            read_check(file.read((char*)&count, 8));
            read_check(file.read((char*)&size, 8));
            partition2global_n = new unsigned long[nnodes_local];
            read_check(file.read((char*)partition2global_n, count*size));
            
            read_check(file.read((char*)&count, 8));
            read_check(file.read((char*)&size, 8));
            node_owner = new unsigned long[nnodes_local];
            read_check(file.read((char*)node_owner, count*size));
            
            read_check(file.read((char*)&count, 8));
            read_check(file.read((char*)&size, 8));
            global_bdry_index = new int8_t[nnodes_global];
            read_check(file.read((char*)global_bdry_index, count*size));
        }
        
        if(nel_ghost > 0){// FV only
            read_check(file.read((char*)&count, 8));
            read_check(file.read((char*)&size, 8));
            element_owner = new unsigned long[count];
            read_check(file.read((char*)element_owner, count*size));
            
            read_check(file.read(in8, 8));
            num_neighbor_partitions = ((int64_t*)in8)[0];
            
            if(num_neighbor_partitions > 0){
                neighboring_partitions = new unsigned long[num_neighbor_partitions];
                ghost_counts = new unsigned long[num_neighbor_partitions];
                ghost_index = new unsigned long*[num_neighbor_partitions];
                for(unsigned long i=0; i<num_neighbor_partitions; i++){
                    read_check(file.read((char*)&count, 8));
                    read_check(file.read((char*)&size, 8));
                    ghost_index[i] = new unsigned long[count];
                    read_check(file.read((char*)ghost_index[i], count*size));
                }
            }
        }
        
    }else{ // Not partitioned
        nel_owned = nel_local;
        nel_ghost = 0;
        nface_owned = num_faces;
        nface_ghost = 0;
        nnodes_global = nnodes_local;
        nnodes_borrowed = 0;
        
        // partition2global_e = new unsigned long[nel_local];
        partition2global_n = new unsigned long[nnodes_local];
        node_owner = new unsigned long[nnodes_local];
        for(unsigned long i=0; i<nnodes_local; i++){
            partition2global_n[i] = i;
            node_owner[i] = my_partition;
        }
    }
    
    is_bdry_node = new bool[nnodes_local];
    for(unsigned long i=0; i<nnodes_local; i++){
        is_bdry_node[i] = false;
    }
    for(int i=0; i<num_bids; i++){
        for(int j=0; j<nodes_per_bid[i]; j++){
            is_bdry_node[bdry[i][j]] = true;
        }
    }
    
    file.close();
}

void display_array(double *a, int cols, int rows){
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            std::cout << "	" << a[i*cols + j];
        }
        std::cout << "\\n";
    }
    std::cout << "\\n";
}
void display_array(unsigned long *a, int cols, int rows){
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            std::cout << "	" << a[i*cols + j];
        }
        std::cout << "\\n";
    }
    std::cout << "\\n";
}

void finch::Mesh::display(int until){
    int maxel = until;
    int maxnode = until;
    if(maxel>nel_local){ maxel = nel_local;}
    if(maxnode>nnodes_local){ maxnode = nnodes_local;}
    
    std::cout << "Mesh info (only showing first " << until << " parts)\\n";
    std::cout << "dimension: " << dimension << "\\n";
    std::cout << "nel_local: " << nel_local << "\\n";
    std::cout << "nnodes_local: " << nnodes_local << "\\n";
    std::cout << "nodes: \\n";
    display_array(allnodes, dimension, maxnode);
    std::cout << "element partition2global: \\n";
    display_array(partition2global_e, maxel, 1);
    std::cout << "nnodes_borrowed: " << nnodes_borrowed << "\\n";
    std::cout << "node partition2global: \\n";
    display_array(partition2global_n, maxnode, 1);
    
}


finch::Refel::Refel(){
    dimension = 0;
    element_type = finch::ElementType::UNSUPPORTED;
    N = 0;
    Np = 0;
    Nqp = 0;
    Nfaces = 0;
    Nfp = nullptr;
    r = nullptr;
    wr = nullptr;
    g = nullptr;
    wg = nullptr;
    V = nullptr;
    invV = nullptr;
    gradV = nullptr;
    Vg = nullptr;
    invVg = nullptr;
    gradVg = nullptr;
    Q = nullptr;
    Qr = nullptr;
    Qs = nullptr;
    Qt = nullptr;
    Ddr = nullptr;
    Dds = nullptr;
    Ddt = nullptr;
    face2local = nullptr;
    surf_r = nullptr;
    surf_wr = nullptr;
    surf_g = nullptr;
    surf_wg = nullptr;
    surf_V = nullptr;
    surf_gradV = nullptr;
    surf_Vg = nullptr;
    surf_gradVg = nullptr;
    surf_Q = nullptr;
    surf_Qr = nullptr;
    surf_Qs = nullptr;
    surf_Qt = nullptr;
    surf_Ddr = nullptr;
    surf_Dds = nullptr;
    surf_Ddt = nullptr;
}

finch::Refel::~Refel(){
    // There is an issue with deletion. Skip the deconstructor for now.
    return;
    if (Np < 1){
        return;
    }
    Np = 0;
    delete [] Nfp;
    delete [] r;
    delete [] wr;
    delete [] g;
    delete [] wg;
    delete [] V;
    delete [] invV;
    delete [] gradV;
    delete [] Vg;
    delete [] invVg;
    delete [] gradVg;
    delete [] Q;
    delete [] Qr;
    delete [] Ddr;
    if(dimension > 1){
        delete [] Qs;
        delete [] Dds;
    }
    if(dimension > 2){
        delete [] Qt;
        delete [] Ddt;
    }
    for (int i=0; i<Nfaces; i++)
        delete [] face2local[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_r[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_wr[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_g[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_wg[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_V[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_gradV[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_Vg[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_gradVg[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_Q[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_Qr[i];
    for (int i=0; i<Nfaces; i++)
        delete [] surf_Ddr[i];
    if(dimension > 1){
        for (int i=0; i<Nfaces; i++)
            delete [] surf_Qs[i];
        for (int i=0; i<Nfaces; i++)
            delete [] surf_Dds[i];
            
        delete [] surf_Qs;
        delete [] surf_Dds;
    }
    if(dimension > 2){
        for (int i=0; i<Nfaces; i++)
            delete [] surf_Qt[i];
        for (int i=0; i<Nfaces; i++)
            delete [] surf_Ddt[i];
        
        delete [] surf_Qt;
        delete [] surf_Ddt;
    }
    delete [] face2local;
    delete [] surf_r;
    delete [] surf_wr;
    delete [] surf_g;
    delete [] surf_wg;
    delete [] surf_V;
    delete [] surf_gradV;
    delete [] surf_Vg;
    delete [] surf_gradVg;
    delete [] surf_Q;
    delete [] surf_Qr;
    delete [] surf_Ddr;
}

void finch::Refel::import_refel(std::string filename){
    std::ifstream file;
    file.open(filename, std::ios::binary);
    if(!file){
        std::cout << "Error: couldn't open Finch refel file: " << filename << ". Exiting.\\n";
        exit(0);
    }
    
    // These will temporarily hold the read values
    char in8[8];
    char in1[1];
    unsigned long count, count2, size;
    
    read_check(file.read(in8, 8));
    dimension = ((int64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    N = ((int64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    Np = ((int64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    Nqp = ((int64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    Nfaces = ((int64_t*)in8)[0];
    
    element_type = finch::ElementType::UNSUPPORTED;
    if(dimension == 1){
        element_type = finch::ElementType::LINE;
    }else if(dimension == 2 && Nfaces == 3){
        element_type = finch::ElementType::TRI;
    }else if(dimension == 2 && Nfaces == 4){
        element_type = finch::ElementType::QUAD;
    }else if(dimension == 3 && Nfaces == 4){
        element_type = finch::ElementType::TET;
    }else if(dimension == 3 && Nfaces == 6){
        element_type = finch::ElementType::HEX;
    }
    if(element_type == finch::ElementType::UNSUPPORTED){
        std::cout << "Error: Unsupported element type: dimension = " << dimension << ", num_faces = " << Nfaces << ". Exiting.\\n";
        exit(0);
    }
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    Nfp = new unsigned long[count];
    read_check(file.read((char*)Nfp, count*size));
    
    // Nodes and vandermonde
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    r = new double[count];
    read_check(file.read((char*)r, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    wr = new double[count];
    read_check(file.read((char*)wr, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    g = new double[count];
    read_check(file.read((char*)g, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    wg = new double[count];
    read_check(file.read((char*)wg, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    V = new double[count];
    read_check(file.read((char*)V, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    gradV = new double[count];
    read_check(file.read((char*)gradV, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    invV = new double[count];
    read_check(file.read((char*)invV, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    Vg = new double[count];
    read_check(file.read((char*)Vg, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    gradVg = new double[count];
    read_check(file.read((char*)gradVg, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    invVg = new double[count];
    read_check(file.read((char*)invVg, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    Q = new double[count];
    read_check(file.read((char*)Q, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    Qr = new double[count];
    read_check(file.read((char*)Qr, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    if(dimension > 1){
        Qs = new double[count];
        read_check(file.read((char*)Qs, count*size));
    }
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    if(dimension > 2){
        Qt = new double[count];
        read_check(file.read((char*)Qt, count*size));
    }
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    Ddr = new double[count];
    read_check(file.read((char*)Ddr, count*size));
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    if(dimension > 1){
        Dds = new double[count];
        read_check(file.read((char*)Dds, count*size));
    }
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    if(dimension > 2){
        Ddt = new double[count];
        read_check(file.read((char*)Ddt, count*size));
    }
    
    // Surface
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    if(count != Nfaces){
        // ignore surface data
        return;
    }
    face2local = new unsigned long*[Nfaces];
    surf_r = new double*[Nfaces];
    surf_wr = new double*[Nfaces];
    surf_g = new double*[Nfaces];
    surf_wg = new double*[Nfaces];
    surf_V = new double*[Nfaces];
    surf_gradV = new double*[Nfaces];
    surf_Vg = new double*[Nfaces];
    surf_gradVg = new double*[Nfaces];
    surf_Q = new double*[Nfaces];
    surf_Qr = new double*[Nfaces];
    surf_Ddr = new double*[Nfaces];
    if(dimension > 1){
        surf_Qs = new double*[Nfaces];
        surf_Dds = new double*[Nfaces];
    }
    if(dimension > 2){
        surf_Qt = new double*[Nfaces];
        surf_Ddt = new double*[Nfaces];
    }
    
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        face2local[i] = new unsigned long[count2];
        read_check(file.read((char*)face2local[i], count2*size));
    }
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_r[i] = new double[count2];
        read_check(file.read((char*)surf_r[i], count2*size));
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_wr[i] = new double[count2];
        read_check(file.read((char*)surf_wr[i], count2*size));
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_g[i] = new double[count2];
        read_check(file.read((char*)surf_g[i], count2*size));
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_wg[i] = new double[count2];
        read_check(file.read((char*)surf_wg[i], count2*size));
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_V[i] = new double[count2];
        read_check(file.read((char*)surf_V[i], count2*size));
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_gradV[i] = new double[count2];
        read_check(file.read((char*)surf_gradV[i], count2*size));
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_Vg[i] = new double[count2];
        read_check(file.read((char*)surf_Vg[i], count2*size));
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_gradVg[i] = new double[count2];
        read_check(file.read((char*)surf_gradVg[i], count2*size));
    }
    
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_Q[i] = new double[count2];
        read_check(file.read((char*)surf_Q[i], count2*size));
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_Qr[i] = new double[count2];
        read_check(file.read((char*)surf_Qr[i], count2*size));
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        if(count2 > 0){
            surf_Qs[i] = new double[count2];
            read_check(file.read((char*)surf_Qs[i], count2*size));
        }
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        if(count2 > 0){
            surf_Qt[i] = new double[count2];
            read_check(file.read((char*)surf_Qt[i], count2*size));
        }
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        surf_Ddr[i] = new double[count2];
        read_check(file.read((char*)surf_Ddr[i], count2*size));
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        if(count2 > 0){
            surf_Dds[i] = new double[count2];
            read_check(file.read((char*)surf_Dds[i], count2*size));
        }
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    for(int i=0; i<count; i++){
        read_check(file.read((char*)&count2, 8));
        read_check(file.read((char*)&size, 8));
        if(count2 > 0){
            surf_Ddt[i] = new double[count2];
            read_check(file.read((char*)surf_Ddt[i], count2*size));
        }
    }
    
    file.close();
}
"""
meshhpp = """
/*
Utilities for interfacing aMat with mesh and reference element data from Finch.
Finch exports binary files containing the mesh(grid_data in Finch), and refel.

For partitioned meshes, each partition is placed in a separate file with "_pn" appended for partition n.
*/
#pragma once
#include <string>

namespace finch{
    
    enum class ElementType{UNSUPPORTED, LINE, TRI, QUAD, TET, HEX};       // element types
    
    // The local partition of the mesh
    struct Mesh{
        int dimension;                      // dimension
        unsigned long nel_global, nel_local;// number of elements (global/local)
        unsigned long nnodes_global, nnodes_local;// number of nodes
        unsigned long *elemental_order; // order for the elemental loop
        // nodes
        double        *allnodes;        // coordinates of all nodes ([x1, y1, z1, x2, y2, z2,...])
        // boundary
        int           num_bids;         // Number of boundary IDs
        unsigned long *nodes_per_bid;   // Number of nodes for each bid
        unsigned long *faces_per_bid;   // Number of faces for each bid
        unsigned long **bdry;           // Indices of boundary nodes for each BID (bdry[bid][nodes])*note:array of arrays
        unsigned long **bdry_face;      // Indices of faces touching each BID (bdryface[bid][faces])*note:array of arrays
        double        **bdry_normal;    // Normal vector for boundary nodes for each BID (bdrynorm[bid][dim, nodes])*note:array of arrays
        int           *bids;            // BID corresponding to arrays of bdrynodes
        int           *nodebid;         // BID of each node
        // elements
        int           *nodes_per_element;// Number of nodes per element
        int           *vertices_per_element;// Number of vertex nodes per element
        int           *faces_per_element;// Number of faces per element
        unsigned long **loc2glb;        // local to global map for each element's nodes (size is (Np, nel))
        unsigned long **glbvertex;      // global indices of each elements' vertices (size is (Nvertex, nel))
        // faces
        unsigned long num_faces;        // Number of faces
        int           *nodes_per_face;  // Number of nodes per face
        unsigned long **face2glb;       // local to global map for faces (size is (Nfp, G, Nfaces))
        unsigned long **element2face;   // face indices for each element (size is (Nfaces, nel))
        unsigned long *face2element;    // elements on both sides of a face, 0=boundary (size is (2, Nfaces))
        double        *face_normal;     // normal vector for each face
        unsigned long *face_refel_index;// Index for face within the refel for each side
        long          *face_bid;        // BID of each face (0=interior face)
        // partition info
        bool          is_subgrid;       // Is this a partition of a greater grid? should be true
        unsigned long nel_owned;        // Number of elements owned by this partition
        unsigned long nel_ghost;        // Number of ghost elements
        unsigned long nface_owned;      // Number of faces owned by this partition
        unsigned long nface_ghost;      // Number of ghost faces that are not owned
        unsigned long nnodes_borrowed;  // Number of nodes shared with other partitions that are not owned by this one.
        unsigned long  *element_owner;  // The partition of each ghost element's owner or -1 if locally owned
        unsigned long  *node_owner;     // The partition of each ghost node's owner
        unsigned long *partition2global_e;// Map from partition elements to global mesh element index
        unsigned long *partition2global_n;// Map from partition nodes to global node index
        int8_t        *global_bdry_index;// boundary ID for all Global nodes 
        unsigned long  num_neighbor_partitions;// number of partitions that share ghosts with this.
        unsigned long  *neighboring_partitions;// IDs of neighboring partitions
        unsigned long  *ghost_counts;          // How many ghosts for each neighbor
        unsigned long **ghost_index;          // Lists of ghost elements to send/recv for each neighbor
        
        // Things not imported from Finch, but set herein.
        bool          *is_bdry_node;    // True if this node is a bdry node
        
        Mesh();
        ~Mesh();
        
        void import_mesh(std::string filename, int num_partitions, int my_partition);
        
        void display(int until);
    };
    
    // Reference element
    struct Refel{
        int dimension;      // dimension
        ElementType element_type;   // type of element from the enum above
        int N;          // order of polynomials
        int Np;         // number of nodes
        int Nqp;        // number of quadrature nodes
        int Nfaces;         // number of faces
        unsigned long *Nfp;   // number of nodes per face
        
        double *r, *wr;     // GLL nodes and quadrature weights
        double *g, *wg;     // Gauss nodes and quadrature weights
        
        double *V, *invV;   // Vandermonde matrix of basis at r and inverse
        double *gradV;      // grad of V
        
        double *Vg, *invVg; // Vandermonde matrix of basis at g and inverse
        double *gradVg;     // grad of Vg
        
        // Useful quadrature matrices for the volume integrals
        // Use these like transpose(Q) * diag(wg) * Q  to compute integral(phi*phi, dx)
        // Derivative versions will also require geometric factors.
        double *Q;          // quadrature matrix: like Vg*invV
        double *Qr;         // quadrature of derivative matrix: like gradVg*invV
        double *Qs;         //
        double *Qt;         //
        
        double *Ddr;        // Derivatives at the elemental nodes, not quadrature nodes
        double *Dds;        //
        double *Ddt;        //
        
        // Surface versions for surface quadrature
        unsigned long **face2local;// local indices for face nodes for each face
        
        double **surf_r, **surf_wr;// nodes and weights on surface
        double **surf_g, **surf_wg;// Gauss nodes and quadrature weights on surface
        
        double **surf_V;          // Vandermonde matrix of basis at r
        double **surf_gradV;      // grad of V
        
        double **surf_Vg;         // Vandermonde matrix of basis at g 
        double **surf_gradVg;     // grad of Vg
        
        double **surf_Q;          // quadrature matrix: like Vg*invV
        double **surf_Qr;         // quadrature of derivative matrix: like gradVg*invV
        double **surf_Qs;         //
        double **surf_Qt;         //
        
        double **surf_Ddr;        // Derivatives at the elemental nodes, not quadrature nodes
        double **surf_Dds;        //
        double **surf_Ddt;        //
        
        Refel();
        ~Refel();
        
        void import_refel(std::string filename);
    };
}
"""

geocpp = """
/*
Utilities related to geometric factors.
*/
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include "finch_geometry.hpp"

// File IO check macro
#define read_check(s) if(!s){ std::cout << "Error: Failed to read some data from the geometry file. Exiting.\\n"; exit(0); }

finch::GeometricFactors::GeometricFactors(){
    dimension = 0;
    nel = 0;
    vals_per_element = 0;
    constant_jacobian = false;
    detJ = nullptr;
    rx = nullptr;
    ry = nullptr;
    rz = nullptr;
    sx = nullptr;
    sy = nullptr;
    sz = nullptr;
    tx = nullptr;
    ty = nullptr;
    tz = nullptr;
}

finch::GeometricFactors::~GeometricFactors(){
    // There is an issue with deletion. Skip the deconstructor for now.
    return;
    if(dimension > 0){
        for(unsigned long i=0; i<nel; i++){
            delete [] detJ[i];
            delete [] rx[i];
            if(dimension > 1){
                delete [] ry[i];
                delete [] sx[i];
                delete [] sy[i];
            }
            if(dimension > 2){
                delete [] rz[i];
                delete [] sz[i];
                delete [] tx[i];
                delete [] ty[i];
                delete [] tz[i];
            }
        }
        
        delete [] detJ;
        delete [] rx;
        
        if(dimension > 1){
            delete [] ry;
            delete [] sx;
            delete [] sy;
        }
        
        if(dimension > 2){
            delete [] rz;
            delete [] sz;
            delete [] tx;
            delete [] ty;
            delete [] tz;
        }
    }
}

void finch::GeometricFactors::import_geometric_factors(std::string filename, int dim, int num_partitions, int my_partition){
    dimension = dim;
    
    std::ifstream file;
    if(num_partitions > 1){
        // Partitioned meshes are stored in separate files for each partition.
        // Their names have a "_pn" on the end of the filename. (n=partition number)
        std::stringstream newname;
        newname << filename << "_p" << my_partition;
        file.open(newname.str(), std::ios::binary);
        if(!file){
            std::cout << "Error: Partition "<<my_partition<<" couldn't open Finch geofacs file: " << newname.str() << ". Exiting.\\n";
            exit(0);
        }
        
    }else{
        file.open(filename, std::ios::binary);
        if(!file){
            std::cout << "Error: couldn't open Finch geofacs file: " << filename << ". Exiting.\\n";
            exit(0);
        }
    }
    
    // These will temporarily hold the read values
    char in8[8];
    char in1[1];
    unsigned long count, size;
    
    read_check(file.read(in1, 1));
    constant_jacobian = (bool)(in8[0]);
    
    read_check(file.read(in8, 8));
    nel = ((uint64_t*)in8)[0];
    
    read_check(file.read(in8, 8));
    vals_per_element = ((uint64_t*)in8)[0]; // if constant jacobian, this should be 1
    if(vals_per_element > 1 && constant_jacobian){
        // This could be an error, but just switch constant_jacobian
        constant_jacobian = false;
    }
    
    // detJ
    detJ = new double*[nel];
    for(unsigned long i=0; i<nel; i++){
        detJ[i] = new double[vals_per_element];
    }
    read_check(file.read((char*)&count, 8));
    read_check(file.read((char*)&size, 8));
    if(count != nel*vals_per_element){
        std::cout << "Error: Reading geometric factors, number of values for jacobian incorrect: " << nel*vals_per_element << " vs. " << count << ". Exiting.\\n";
        exit(0);
    }
    for(unsigned long i=0; i<nel; i++){
        read_check(file.read((char*)detJ[i], vals_per_element*size));
    }
    
    // J
    rx= new double*[nel];
    if(dimension > 1){
        ry= new double*[nel];
        sx= new double*[nel];
        sy= new double*[nel];
    }
    if(dimension > 2){
        rz= new double*[nel];
        sz= new double*[nel];
        tx= new double*[nel];
        ty= new double*[nel];
        tz= new double*[nel];
    }
    for(unsigned long i=0; i<nel; i++){
        rx[i] = new double[vals_per_element];
        if(dimension > 1){
            ry[i] = new double[vals_per_element];
            sx[i] = new double[vals_per_element];
            sy[i] = new double[vals_per_element];
        }
        if(dimension > 2){
            rz[i] = new double[vals_per_element];
            sz[i] = new double[vals_per_element];
            tx[i] = new double[vals_per_element];
            ty[i] = new double[vals_per_element];
            tz[i] = new double[vals_per_element];
        }
    }
    for(unsigned long i=0; i<nel; i++){
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        read_check(file.read((char*)rx[i], count*size));
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        if(dimension>1){ read_check(file.read((char*)ry[i], count*size)); }
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        if(dimension>2){ read_check(file.read((char*)rz[i], count*size)); }
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        if(dimension>1){ read_check(file.read((char*)sx[i], count*size)); }
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        if(dimension>1){ read_check(file.read((char*)sy[i], count*size)); }
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        if(dimension>2){ read_check(file.read((char*)sz[i], count*size)); }
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        if(dimension>2){ read_check(file.read((char*)tx[i], count*size)); }
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        if(dimension>2){ read_check(file.read((char*)ry[i], count*size)); }
        
        read_check(file.read((char*)&count, 8));
        read_check(file.read((char*)&size, 8));
        if(dimension>2){ read_check(file.read((char*)tz[i], count*size)); }
    }
    
}
"""

geohpp = """
/*
Functions for evaluating boundary conditions.
*/
#pragma once
#include <string>

namespace finch{
    
    struct GeometricFactors{
        int dimension;
        unsigned long nel;
        unsigned long vals_per_element;
        bool constant_jacobian;
        double **detJ;
        double **rx, **ry, **rz;
        double **sx, **sy, **sz;
        double **tx, **ty, **tz;
        
        GeometricFactors();
        ~GeometricFactors();
        
        void import_geometric_factors(std::string filename, int dim, int num_partitions, int my_partition); // When importing values computed in Finch.
    };
}
"""

    meshcppfile = add_generated_file("finch_mesh.cpp", dir="src");
    println(meshcppfile, meshcpp);

    meshhppfile = add_generated_file("finch_mesh.hpp", dir="include");
    println(meshhppfile, meshhpp);

    geocppfile = add_generated_file("finch_geometry.cpp", dir="src");
    println(geocppfile, geocpp);

    geohppfile = add_generated_file("finch_geometry.hpp", dir="include");
    println(geohppfile, geohpp);

end