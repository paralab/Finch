#=
Functions for building an IR from the DSL for FEM type problems.
These consist of assembling a global matrix/vector from local ones.
Quadrature is done with matrices from the refel.

allocate
first loop
    init all
    matrix
    vector
build mat
(solve+place) -or -
(time stepper
    loops
        eval parts needed by rhs only
        lhs matrix if t-dependent ??TODO
        vector
    solve+place
)
=#
function build_IR_fem(lhs_vol, lhs_surf, rhs_vol, rhs_surf, var, indices, config, prob, time_stepper)
    dimension = config.dimension;
    # Count variables, dofs, and store offsets
    varcount = 1;
    dofsper = 0;
    dofsper_loop = 0;
    offset_ind = [0];
    if typeof(var) <:Array
        varcount = length(var);
        offset_ind = zeros(Int, varcount);
        dofsper = var[1].total_components;
        dofsper_loop = length(var[1].symvar);
        for i=2:length(var)
            offset_ind[i] = dofsper;
            dofsper = dofsper + var[i].total_components;
            dofsper_loop = dofsper_loop + length(var[i].symvar);
        end
    else
        var = [var]; # put it in an array for consistency
        offset_ind = zeros(Int, 1);
        dofsper = var[1].total_components;
        dofsper_loop = length(var[1].symvar);
    end
    
    IRtypes = IR_entry_types();
    
    # These will hold the IR
    allocate_block = IR_block_node([],"allocation");
    prepare_block = IR_block_node([],"prepare");
    derivmat_block = IR_block_node([],"deriv mat");
    time_prepare_block = IR_block_node([], "prepare");
    time_derivmat_block = IR_block_node([], "deriv mat");
    matrix_block = IR_block_node([],"elemental matrix");
    vector_block = IR_block_node([],"elemental vector");
    bdry_block = IR_block_node([],"boundary");
    toglobal_block = IR_block_node([],"local to global");
    time_bdry_block = IR_block_node([],"boundary");
    time_toglobal_block = IR_block_node([],"local to global");
    
    # Allocate the global matrix and vector
    push!(allocate_block.parts, IR_comment_node("Allocate global matrix(IJV form) and vector."));
    allocatedNZ = IR_data_node(IRtypes.float_64_data, :allocated_nonzeros);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        allocatedNZ,
        IR_operation_node(IRtypes.math_op, [:*, :num_elements, :dofs_per_element, :dofs_per_element])
        ]));
    globalmat_I = IR_data_node(IRtypes.int_64_data, :global_matrix_I, [:allocated_nonzeros], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        globalmat_I,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.int_64_data, :allocated_nonzeros])
        ]));
    globalmat_J = IR_data_node(IRtypes.int_64_data, :global_matrix_J, [:allocated_nonzeros], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        globalmat_J,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.int_64_data, :allocated_nonzeros])
        ]));
    globalmat_V = IR_data_node(IRtypes.float_64_data, :global_matrix_V, [:allocated_nonzeros], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        globalmat_V,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :allocated_nonzeros])
        ]));
    
    # Global vectors
    globalvec = IR_data_node(IRtypes.float_64_data, :global_vector, [:dofs_global], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        globalvec,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])
        ]));
    
    solvec = IR_data_node(IRtypes.float_64_data, :solution, [:dofs_global], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        solvec,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])
        ]));
        
    # I and J should be initialized to 1 for julia
    # How about c++ ?? do we need a named op for init of IJV?
    push!(allocate_block.parts, IR_comment_node("I and J vectors should init as ones"));
    push!(allocate_block.parts, IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY, globalmat_I, 1, :allocated_nonzeros]));
    push!(allocate_block.parts, IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY, globalmat_J, 1, :allocated_nonzeros]));
    
    # Allocate the elemental matrix and vector
    push!(allocate_block.parts, IR_comment_node("Allocate elemental matrix and vector."));
    elementmat = IR_data_node(IRtypes.float_64_data, :element_matrix, [:local_system_size, :local_system_size], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        elementmat,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :local_system_size, :local_system_size])
        ]));
    
    elementvec = IR_data_node(IRtypes.float_64_data, :element_vector, [:local_system_size], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        elementvec,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :local_system_size])
        ]));
    
    # bdry done flag for each node
    push!(allocate_block.parts, IR_comment_node("Boundary done flag for each node."));
    bdry_done = IR_data_node(IRtypes.int_64_data, :bdry_done, [:nnodes_global], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        bdry_done,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.int_64_data, :nnodes_global])
        ]));
    
    # coefficient prep
    # a list of all entities and rhs only ones
    all_entities = [];
    rhs_entities = [];
    counts = zeros(Int,4); # how many entities for each piece 
    
    # LHS volume
    if !(lhs_vol === nothing)
        entities = extract_entities(lhs_vol);
        append!(all_entities, entities);
        counts[1] = length(entities);
    end
    
    # RHS volume
    if !(rhs_vol === nothing)
        entities = extract_entities(rhs_vol);
        append!(all_entities, entities);
        append!(rhs_entities, entities);
        counts[2] = length(entities);
    end
    
    # LHS surface
    if !(lhs_surf === nothing)
        entities = extract_entities(lhs_surf);
        append!(all_entities, entities);
        counts[3] = length(entities);
    end
    
    # RHS surface
    if !(rhs_surf === nothing)
        entities = extract_entities(rhs_surf);
        append!(all_entities, entities);
        append!(rhs_entities, entities);
        counts[4] = length(entities);
    end
    
    # full loop
    (allocate, coef) = prepare_coefficient_values(all_entities, var, dimension, counts);
    push!(allocate_block.parts, IR_comment_node("Allocate coefficient vectors."));
    append!(allocate_block.parts, allocate);
    push!(prepare_block.parts, IR_comment_node("Evaluate coefficients."));
    append!(prepare_block.parts, coef);
    push!(prepare_block.parts, IR_operation_node(IRtypes.assign_op, [:detj, 
        IR_operation_node(IRtypes.member_op, [:geometric_factors, 
            IR_data_node(IRtypes.float_64_data, :detJ, [:?,:?],  [1, :eid])])]))
    
    # derivative matrix prep
    (allocate, deriv) = prepare_derivative_matrices(all_entities, counts, var);
    if length(allocate)>0 
        push!(allocate_block.parts, IR_comment_node("Allocate for derivative matrices."));
        append!(allocate_block.parts, allocate);
        push!(deriv, IR_comment_node("Prepare derivative matrices."));
        derivmat_block.parts = deriv;
    end
    
    # rhs only loop
    counts[1] = 0; # no lhs entities
    counts[3] = 0;
    (allocate, rhscoef) = prepare_coefficient_values(rhs_entities, var, dimension, counts);
    push!(time_prepare_block.parts, IR_comment_node("Evaluate coefficients."));
    append!(time_prepare_block.parts, rhscoef);
    push!(time_prepare_block.parts, IR_operation_node(IRtypes.assign_op, [:detj, 
        IR_operation_node(IRtypes.member_op, [:geometric_factors, 
            IR_data_node(IRtypes.float_64_data, :detJ, [:?,:?],  [1, :eid])])]))
    
    (allocate, rhsderiv) = prepare_derivative_matrices(rhs_entities, counts, var);
    if length(rhsderiv)>0 
        push!(rhsderiv, IR_comment_node("Prepare derivative matrices."));
        time_derivmat_block.parts = rhsderiv;
    end
    
    # computation
    if !(lhs_vol === nothing)
        lhsvol_terms = process_terms(lhs_vol);
        push!(matrix_block.parts, make_elemental_computation_fem(lhsvol_terms, var, dofsper, offset_ind, LHS, "volume"));
    end
    if !(rhs_vol === nothing)
        rhsvol_terms = process_terms(rhs_vol);
        push!(vector_block.parts, make_elemental_computation_fem(rhsvol_terms, var, dofsper, offset_ind, RHS, "volume"));
    end
    if !(lhs_surf === nothing) 
        lhssurf_terms = process_terms(lhs_surf);
        push!(matrix_block.parts, make_elemental_computation_fem(lhssurf_terms, var, dofsper, offset_ind, LHS, "surface"));
    end
    if !(rhs_surf === nothing)
        rhssurf_terms = process_terms(rhs_surf);
        push!(vector_block.parts, make_elemental_computation_fem(rhssurf_terms, var, dofsper, offset_ind, RHS, "surface"));
    end
    
    # bdry block 
    push!(bdry_block.parts, IR_comment_node("Apply boundary conditions."));
    push!(bdry_block.parts, IR_operation_node(IRtypes.function_op, [:apply_boundary_conditions_elemental, 
        :var, :eid, :mesh, :refel, :geometric_factors, :prob, :t, :element_matrix, :element_vector, :bdry_done, :index_offset]));
    push!(time_bdry_block.parts, IR_comment_node("Apply boundary conditions."));
    push!(time_bdry_block.parts, IR_operation_node(IRtypes.function_op, [:apply_boundary_conditions_elemental_rhs, 
        :var, :eid, :mesh, :refel, :geometric_factors, :prob, :t, :element_vector, :bdry_done, :index_offset]));
    
    # add to global sysem
    push!(toglobal_block.parts, IR_comment_node("Place elemental parts in global system."));
    # push!(toglobal_block.parts, generate_local_to_global_fem(dofsper));
    push!(toglobal_block.parts, IR_operation_node(IRtypes.named_op, [:LOCAL2GLOBAL, dofsper]));
    push!(time_toglobal_block.parts, IR_comment_node("Place elemental parts in global system."));
    # push!(time_toglobal_block.parts, generate_local_to_global_fem(dofsper, vec_only=true));
    push!(time_toglobal_block.parts, IR_operation_node(IRtypes.named_op, [:LOCAL2GLOBAL_VEC, dofsper]));
    
    # assembly loop
    assembly_loop = generate_assembly_loop_fem(var, indices);
    rhs_assembly_loop = generate_assembly_loop_fem(var, indices);
    
    # find the innermost assembly loop
    inner_loop = assembly_loop;
    rhs_inner_loop = rhs_assembly_loop;
    while length(inner_loop.body.parts) > 0 && typeof(inner_loop.body.parts[1]) == IR_loop_node
        inner_loop = inner_loop.body.parts[1];
        rhs_inner_loop = rhs_inner_loop.body.parts[1];
    end
    
    # fill the loop
    append!(inner_loop.body.parts, [
        derivmat_block,
        prepare_block,
        matrix_block,
        vector_block,
        bdry_block,
        toglobal_block
    ])
    append!(rhs_inner_loop.body.parts, [
        time_derivmat_block,
        time_prepare_block,
        vector_block,
        time_bdry_block,
        time_toglobal_block
    ])
    
    # Package it all in one piece to add to master
    if prob.time_dependent
        if prob.nonlinear
            # The time stepping loop to be altered
            (t_allocate, step_loop) = generate_time_stepping_loop_fem(time_stepper, rhs_assembly_loop, false);
            push!(allocate_block.parts, t_allocate);
            
            # Pieces needed for nonlinear iteration
            nl_change_abs = IR_data_node(IRtypes.float_64_data, :nl_change_abs);
            nl_change_rel = IR_data_node(IRtypes.float_64_data, :nl_change_rel);
            oldsolvec = IR_data_node(IRtypes.float_64_data, :solution_old, [:dofs_global], []);
            zero_el_vec = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY,
                IR_data_node(IRtypes.float_64_data, :global_vector, [:dofs_global], []), 0, :dofs_global]);
            zero_bdry_done = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY,
                IR_data_node(IRtypes.int_64_data, :bdry_done, [:nnodes_global], []), 0, :nnodes_global]);
            push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
                oldsolvec,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])
            ]));
            
            # The nonlinear iteration loop will contain the body of the time step loop
            nonlinear_loop = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [nl_change_abs, 1]),
                IR_operation_node(IRtypes.assign_op, [nl_change_rel, 1]),
                # The nonlinear iteration
                wrap_in_timer(:nonlinear_iteration, IR_loop_node(IRtypes.while_loop, :while_loop, :nl_iter, 0, 
                    IR_operation_node(IRtypes.math_op, [:&&, 
                        IR_operation_node(IRtypes.math_op, [:<, :nl_iter, prob.max_iters]),
                        IR_operation_node(IRtypes.math_op, [:>, nl_change_rel, prob.relative_tol]),
                        IR_operation_node(IRtypes.math_op, [:>, nl_change_abs, prob.absolute_tol])
                    ]), 
                    IR_block_node([
                        # The body of the time step loop
                        step_loop.body,
                        
                        wrap_in_timer(:update_vars, IR_block_node([
                            # oldsolvec = oldsolvec + solvec (u=u+du)
                            # and find absolute and relative norms of solvec
                            IR_operation_node(IRtypes.named_op, [:ADD_GLOBAL_VECTOR_AND_NORM, oldsolvec, solvec, nl_change_abs, nl_change_rel]),
                            # put the updated values into vars
                            IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, oldsolvec, :nl_var]),
                        ]))
                    ])
                ))
                
            ],"nonlinear iteration")
            
            # Put the nonlinear iteration back inside the time step loop
            step_loop.body = nonlinear_loop;
            
            compute_block = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [:t, 0]),
                IR_operation_node(IRtypes.assign_op, [:dt, time_stepper.dt]),
                IR_comment_node("Initial loop to build matrix"),
                wrap_in_timer(:first_assembly, assembly_loop),
                IR_operation_node(IRtypes.named_op, [:GLOBAL_FINALIZE]),
                IR_operation_node(IRtypes.named_op, [:GATHER_SOLUTION, oldsolvec, :nl_var]),
                
                IR_comment_node("###############################################"),
                IR_comment_node("Time stepping loop"),
                step_loop
            ],"computation")
            
        else # linear
            (t_allocate, step_loop) = generate_time_stepping_loop_fem(time_stepper, rhs_assembly_loop);
            # step_loop.body = IR_block_node([rhs_assembly_loop]); # TODO add solve and scatter
            push!(allocate_block.parts, t_allocate);
            compute_block = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [:t, 0]),
                IR_operation_node(IRtypes.assign_op, [:dt, time_stepper.dt]),
                IR_comment_node("Initial loop to build matrix"),
                wrap_in_timer(:first_assembly, assembly_loop),
                IR_operation_node(IRtypes.named_op, [:GLOBAL_FINALIZE]),
                IR_operation_node(IRtypes.named_op, [:GATHER_SOLUTION, solvec]),
                
                IR_comment_node("###############################################"),
                IR_comment_node("Time stepping loop"),
                step_loop
            ],"computation");
        end
        
    else # stationary
        if prob.nonlinear
            # # Newton's method will look like:
            # # solve for DELTAu ( = -F/F' = A\b)
            # # u = u + DELTAu
            # # check norm(DELTAu) for convergence 
            # #
            # nl_change_abs = IR_data_node(IRtypes.float_64_data, :nl_change_abs);
            # nl_change_rel = IR_data_node(IRtypes.float_64_data, :nl_change_rel);
            # oldsolvec = IR_data_node(IRtypes.float_64_data, :solution_old, [:dofs_global], []);
            # zero_el_vec = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY,
            #     IR_data_node(IRtypes.float_64_data, :global_vector, [:dofs_global], []), 0, :dofs_global]);
            # zero_bdry_done = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY,
            #     IR_data_node(IRtypes.int_64_data, :bdry_done, [:nnodes_global], []), 0, :nnodes_global]);
            # push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            #     oldsolvec,
            #     IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])
            # ]));
            # compute_block = IR_block_node([
            #     IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_64_data, :t), 0]),
            #     IR_operation_node(IRtypes.assign_op, [nl_change_abs, 1]),
            #     IR_operation_node(IRtypes.assign_op, [nl_change_rel, 1]),
            #     # build vector from existing sol vars for initial guess
            #     IR_operation_node(IRtypes.named_op, [:GATHER_SOLUTION, oldsolvec, :nl_var]),
                
            #     # The nonlinear iteration
            #     wrap_in_timer(:nonlinear_iteration, IR_loop_node(IRtypes.while_loop, :while_loop, :nl_iter, 0, 
            #         IR_operation_node(IRtypes.math_op, [:&&, 
            #             IR_operation_node(IRtypes.math_op, [:<, :nl_iter, prob.max_iters]),
            #             IR_operation_node(IRtypes.math_op, [:>, nl_change_rel, prob.relative_tol]),
            #             IR_operation_node(IRtypes.math_op, [:>, nl_change_abs, prob.absolute_tol])
            #         ]), 
            #         IR_block_node([
            #             # zero vec and bdry (May need mat as well depending on target?)
            #             zero_el_vec,
            #             zero_bdry_done,
            #             # compute
            #             wrap_in_timer(:assembly, assembly_loop),
            #             IR_operation_node(IRtypes.named_op, [:GLOBAL_FINALIZE]),
            #             wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, solvec, :global_matrix, :global_vector])),
            #             wrap_in_timer(:update_vars, IR_block_node([
            #                 # oldsolvec = oldsolvec + solvec (u=u+du)
            #                 # and find absolute and relative norms of solvec
            #                 IR_operation_node(IRtypes.named_op, [:ADD_GLOBAL_VECTOR_AND_NORM, oldsolvec, solvec, nl_change_abs, nl_change_rel]),
            #                 # put the updated values into vars
            #                 IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, oldsolvec, :nl_var]),
            #             ]))
            #         ])
            #     ))
                
            # ],"computation")
            
            
            ###### A different approach using Taylor approx for nonlinear terms
            nl_change_abs = IR_data_node(IRtypes.float_64_data, :nl_change_abs);
            nl_change_rel = IR_data_node(IRtypes.float_64_data, :nl_change_rel);
            oldsolvec = IR_data_node(IRtypes.float_64_data, :solution_old, [:dofs_global], []);
            zero_el_vec = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY,
                IR_data_node(IRtypes.float_64_data, :global_vector, [:dofs_global], []), 0, :dofs_global]);
            zero_bdry_done = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY,
                IR_data_node(IRtypes.int_64_data, :bdry_done, [:nnodes_global], []), 0, :nnodes_global]);
            push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
                oldsolvec,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])
                ]));
            compute_block = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_64_data, :t), 0]),
                IR_operation_node(IRtypes.assign_op, [nl_change_abs, 1]),
                IR_operation_node(IRtypes.assign_op, [nl_change_rel, 1]),
                # The nonlinear iteration
                wrap_in_timer(:nonlinear_iteration, IR_loop_node(IRtypes.while_loop, :while_loop, :nl_iter, 0, 
                    IR_operation_node(IRtypes.math_op, [:&&, 
                        IR_operation_node(IRtypes.math_op, [:<, :nl_iter, prob.max_iters]),
                        IR_operation_node(IRtypes.math_op, [:>, nl_change_rel, prob.relative_tol]),
                        IR_operation_node(IRtypes.math_op, [:>, nl_change_abs, prob.absolute_tol])
                    ]), 
                    IR_block_node([
                        # zero vec and bdry (May need mat as well depending on target?)
                        zero_el_vec,
                        zero_bdry_done,
                        # compute
                        wrap_in_timer(:assembly, assembly_loop),
                        IR_operation_node(IRtypes.named_op, [:GLOBAL_FINALIZE]),
                        wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, solvec, :global_matrix, :global_vector])),
                        wrap_in_timer(:scatter, IR_block_node([
                            IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, solvec, :var]),
                            IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, solvec, :nl_var]),
                        ])),
                        # norm of change and update old
                        # generate_difference_norms_and_update(solvec, oldsolvec, nl_change_abs, nl_change_rel))
                        wrap_in_timer(:norm_and_update, IR_operation_node(IRtypes.named_op, [
                            :UPDATE_GLOBAL_VECTOR_AND_NORM, solvec, oldsolvec, nl_change_abs, nl_change_rel]))
                    ])
                ))
                
            ],"computation")
            
        else # linear
            compute_block = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [:t, 0]),
                wrap_in_timer(:assembly, assembly_loop),
                IR_operation_node(IRtypes.named_op, [:GLOBAL_FINALIZE]),
                wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, solvec, :global_matrix, :global_vector])),
                wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, solvec])),
            ],"computation")
        end
    end
    
    # Put them all together in a master block
    master_block = IR_block_node([
        allocate_block,
        compute_block
    ],"master");
    
    return master_block;
end

function prepare_derivative_matrices(entities, counts, var)
    IRtypes = IR_entry_types();
    deriv_part = Vector{IR_part}(undef,0); # Derivative matrix building inside elemental loop
    allocate_part = Vector{IR_part}(undef,0); # Allocations that are done outside of the loops
    
    needed_derivative_matrices = fill(false, 8); # 1,2,3 = x,y,z quadrature points, 5,6,7 = nodes
    
    # Loop over entities to check for derivatives
    for i=1:(counts[1]+counts[2])
        # Only check for derivatives
        if length(entities[i].derivs) > 0
            for di=1:length(entities[i].derivs)
                needed_derivative_matrices[entities[i].derivs[di]] = true;
            end
        end
    end
    
    for i=(counts[1]+counts[2]+1):length(entities)
        if is_test_function(entities[i]) || is_unknown_var(entities[i], var)
            # Only check for derivatives
            if length(entities[i].derivs) > 0
                for di=1:length(entities[i].derivs)
                    needed_derivative_matrices[entities[i].derivs[di]] = true;
                end
            end
        else  # Is a coefficient(number or function) or variable(array)?
            # If derivatives are needed, must evaluate at surface nodes, not quadrature points
            if length(entities[i].derivs) > 0
                for di=1:length(entities[i].derivs)
                    needed_derivative_matrices[entities[i].derivs[di] + 4] = true; # This is the RDn matrix for nodal derivatives
                end
            end
        end # if coefficient
    end # entity loop
    
    # Derivative matrices
    for i=1:3
        # allocate
        if needed_derivative_matrices[i]
            push!(allocate_part, IR_operation_node(
                IRtypes.assign_op,[
                IR_data_node(IRtypes.float_64_data, Symbol("RQ"*string(i)), [:qnodes_per_element, :nodes_per_element], []),
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :qnodes_per_element, :nodes_per_element])
                ]));
        end
        if needed_derivative_matrices[i+4]
            push!(allocate_part, IR_operation_node(
                IRtypes.assign_op,[
                IR_data_node(IRtypes.float_64_data, Symbol("RD"*string(i)), [:nodes_per_element, :nodes_per_element], []),
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :nodes_per_element, :nodes_per_element])
                ]));
        end
        
        # Build derivative matrices
        if needed_derivative_matrices[i]
            push!(deriv_part, 
                IR_operation_node(IRtypes.function_op, [:build_derivative_matrix, :refel, :geometric_factors, i, :eid, 0, Symbol("RQ"*string(i))]));
        end
        if needed_derivative_matrices[i+4]
            push!(deriv_part, 
                IR_operation_node(IRtypes.function_op, [:build_derivative_matrix, :refel, :geometric_factors, i, :eid, 1, Symbol("RD"*string(i))]));
        end
    end
    
    return (allocate_part, deriv_part);
end

# Allocate, compute, or fetch all needed values
function prepare_coefficient_values(entities, var, dimension, counts; rhs_only=false)
    IRtypes = IR_entry_types();
    row_col_matrix_index = IR_operation_node(IRtypes.named_op, [:ROWCOL_TO_INDEX, :row, :col, :nodes_per_element]);
    col_row_matrix_index = IR_operation_node(IRtypes.named_op, [:ROWCOL_TO_INDEX, :col, :row, :nodes_per_element]);
    
    # These parts will be returned
    allocate_part = Vector{IR_part}(undef,0); # Allocations that are done outside of the loops
    coef_part = Vector{IR_part}(undef,0); # Coefficient evaluation/preparation inside elemental loop
    deriv_part = Vector{IR_part}(undef,0); # Derivative matrix building inside elemental loop
    
    # coef_loop_body = Vector{IR_part}(undef,0); # Loop that evaluates coefficient values
    coef_init_loop_body = [];
    coef_interp_loop_body = [];
    
    # Check to see if derivative matrices are needed
    needed_derivative_matrices = fill(false, 8); # 1,2,3 = x,y,z quadrature points, 5,6,7 = nodes
    need_normals = false;
    need_vol_coef_loop = false;
    need_vol_interp_loop = false;
    need_surf_coef_loop = false;
    need_surf_interp_loop = false;
    
    unique_entity_names = []; # avoid duplicate names
    
    # This loop evaluates coefficient functions INSIDE THE ELEMENTAL LOOP
    # If these are to be precomputed, this must change.
    # First get x,y,z,t,nodeID
    coef_loop_body = [
        # nodeID = mesh.loc2glb[eid, ni]
        IR_operation_node(IRtypes.assign_op,[
            :nodeID, 
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.int_64_data, :loc2glb, 
                                                        [:nodes_per_element, :num_elements], [:ni, :eid])])]),
        # node_offset = (nodeID-1)*dimension+1
        IR_operation_node(IRtypes.assign_op,[
            :node_offset, 
            IR_operation_node(IRtypes.math_op, [:+, 
                IR_operation_node(IRtypes.math_op, [:*,
                    IR_operation_node(IRtypes.math_op, [:-, :nodeID, 1]),
                    dimension]),
                1
            ])]),
        # t = the current step time
        # IR_statement_node(IRtypes.assign_statement, :t, 0.0),
        # x = mesh.allnodes[nodeID*dimension]
        IR_operation_node(IRtypes.assign_op, [
            :x, 
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.float_64_data, :allnodes, 
                                                        [:dimension, :nnodes_partition], [1,:nodeID])])])
    ];
    if dimension > 1
        # y = mesh.allnodes[nodeID*dimension+1]
        push!(coef_loop_body, IR_operation_node(IRtypes.assign_op, [
            :y, 
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.float_64_data, :allnodes, 
                                                        [:dimension, :nnodes_partition], [2,:nodeID])])]))
    else
        # y = 0
        push!(coef_loop_body, IR_operation_node(IRtypes.assign_op, [:y, 0]));
    end
    if dimension > 2
        # z = mesh.allnodes[nodeID*dimension+2]
        push!(coef_loop_body, IR_operation_node(IRtypes.assign_op, [
            :z, 
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.float_64_data, :allnodes, 
                                                        [:dimension, :nnodes_partition], [3,:nodeID])])]))
    else
        # z = 0
        push!(coef_loop_body, IR_operation_node(IRtypes.assign_op, [:z, 0]));
    end
    
    # Loop over entities to perpare for each one
    for i=1:length(entities)
        if i <= counts[1] 
            lorr = LHS; vors = "volume";
        elseif i <= (counts[1]+counts[2]) 
            lorr = RHS; vors = "volume";
        elseif i <= (counts[1]+counts[2]+counts[3]) 
            lorr = LHS; vors = "surface";
        else 
            lorr = RHS; vors = "surface";
        end
         
        if is_test_function(entities[i])
            # Do nothing
        elseif is_unknown_var(entities[i], var) && lorr == LHS
            # Do nothing
        else  # It is a coefficient(number or function) or known variable(array)
            cname = make_entity_name(entities[i]);
            is_unique = true;
            for n in unique_entity_names
                if cname == n && is_unique
                    is_unique = false;
                    break;
                end
            end
            if is_unique
                push!(unique_entity_names, cname); # 
            else
                continue;
            end
            
            (ctype, cval) = get_coef_val(entities[i]);
            if ctype == -1
                # It was a special symbol like dt or normal
                if entities[i].name == "FACENORMAL1" || entities[i].name == "FACENORMAL2"
                    need_normals = true;
                end
                
            elseif ctype == 0
                # It was a number, do nothing?
                
            elseif ctype == 1 # a constant wrapped in a coefficient will be replaced by a number
                push!(coef_part, IR_operation_node(IRtypes.assign_op,[
                    IR_data_node(IRtypes.float_64_data, Symbol(cname), [], []), cval]))
                
            elseif ctype == 2 || ctype == 4 # a coefficient function or indexed coefficient function
                # Build the index IR
                if typeof(entities[i].index) <: Array
                    # It is an indexed variable
                    if length(entities[i].index) == 1
                        indstr = "INDEX_VAL_"*entities[i].index[1];
                        index_IR = Symbol(indstr);
                    else
                        # There is more than one index. Need to form an expression for it.
                        indstr = "(INDEX_VAL_"*entities[i].index[1];
                        index_IR = Symbol(indstr);
                        indices = variables[cval].indexer;
                        for indi=2:length(entities[i].index)
                            indstr *= " + ("*string(length(indices[indi-1].range))*"*(INDEX_VAL_"*entities[i].index[indi]*"-1)";
                            this_ind = "INDEX_VAL_"*entities[i].index[indi];
                            index_IR = IR_operation_node(IRtypes.math_op, [:+, index_IR,
                                IR_operation_node(IRtypes.math_op, [:*, length(indices[indi-1].range),
                                    IR_operation_node(IRtypes.math_op, [:-, Symbol(this_ind), 1])])]);
                        end
                        for indi=1:length(entities[i].index)
                            indstr *= ")";
                        end
                    end
                    
                else
                    indstr = string(entities[i].index);
                    index_IR = entities[i].index;
                end
                if vors == "volume"
                    need_vol_interp_loop = true;
                    need_vol_coef_loop = true;
                    # Need to allocate NP 
                    np_allocate = IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :nodes_per_element]);
                    nqp_allocate = IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :qnodes_per_element]);
                    nodal_coef_name = "NODAL"*cname;
                    # For the nodal values
                    push!(allocate_part, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_64_data, Symbol(nodal_coef_name), [np_allocate], []),
                        np_allocate]));
                    # For the interpolated/differentiated quadrature values
                    push!(allocate_part, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_64_data, Symbol(cname), [nqp_allocate], []),
                        nqp_allocate]));
                    
                    coef_index = get_coef_index(entities[i]);
                    
                    push!(coef_loop_body, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_64_data, Symbol(nodal_coef_name), [np_allocate], [:ni]),
                        IR_operation_node(IRtypes.named_op, [:COEF_EVAL, coef_index, index_IR, :x, :y, :z, :t, :nodeID])]));
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        quad_coef_node = IR_data_node(IRtypes.float_64_data, Symbol(cname), [nqp_allocate], [:col]);
                        nodal_coef_node = IR_data_node(IRtypes.float_64_data, Symbol(nodal_coef_name), [np_allocate], [:row]);
                        push!(coef_init_loop_body, IR_operation_node(IRtypes.assign_op, [quad_coef_node, 0.0]));
                        for di=1:length(entities[i].derivs)
                            deriv_quad_mat = IR_data_node(IRtypes.float_64_data, Symbol("RQ"*string(entities[i].derivs[di])), 
                                                            [:qnodes_per_element, :nodes_per_element], [col_row_matrix_index]);
                            push!(coef_interp_loop_body, IR_operation_node(IRtypes.assign_op,[
                                quad_coef_node,
                                IR_operation_node(IRtypes.math_op, [:(+), quad_coef_node, IR_operation_node(IRtypes.math_op, [:(*), 
                                                                    deriv_quad_mat, nodal_coef_node])])
                                ]));
                        end
                    else
                        quad_coef_node = IR_data_node(IRtypes.float_64_data, Symbol(cname), [nqp_allocate], [:col]);
                        nodal_coef_node = IR_data_node(IRtypes.float_64_data, Symbol(nodal_coef_name), [np_allocate], [:row]);
                        # refelQ = IR_operation_node(IRtypes.member_op, [:refel, IR_data_node(IRtypes.array_data, :Q, [row_col_matrix_index])]);
                        refelQ = IR_data_node(IRtypes.float_64_data, :Q, [:nodes_per_element, :nodes_per_element], [col_row_matrix_index]);
                        push!(coef_init_loop_body, IR_operation_node(IRtypes.assign_op, [quad_coef_node, 0.0]));
                        push!(coef_interp_loop_body, IR_operation_node(IRtypes.assign_op,[
                            quad_coef_node,
                            IR_operation_node(IRtypes.math_op, [:(+), quad_coef_node, IR_operation_node(IRtypes.math_op, [:(*), refelQ, nodal_coef_node])])
                            ]));
                    end
                    
                else # surface
                    # If derivatives are needed, must evaluate at all volume nodes.
                    if length(entities[i].derivs) > 0
                        for di=1:length(entities[i].derivs)
                            # TODO
                        end
                        
                    else # no derivatives
                        # TODO
                    end
                    # Interpolate at surface quadrature nodes
                    # TODO
                end
                
            elseif ctype == 3 # a known variable value
                # Build the index IR
                if typeof(entities[i].index) <: Array
                    # It is an indexed variable
                    if length(entities[i].index) == 1
                        indstr = "INDEX_VAL_"*entities[i].index[1];
                        index_IR = Symbol(indstr);
                    else
                        # There is more than one index. Need to form an expression for it.
                        indstr = "(INDEX_VAL_"*entities[i].index[1];
                        index_IR = Symbol(indstr);
                        indices = variables[cval].indexer;
                        for indi=2:length(entities[i].index)
                            indstr *= " + ("*string(length(indices[indi-1].range))*"*(INDEX_VAL_"*entities[i].index[indi]*"-1)";
                            this_ind = "INDEX_VAL_"*entities[i].index[indi];
                            index_IR = IR_operation_node(IRtypes.math_op, [:+, index_IR,
                                IR_operation_node(IRtypes.math_op, [:*, length(indices[indi-1].range),
                                    IR_operation_node(IRtypes.math_op, [:-, Symbol(this_ind), 1])])]);
                        end
                        for indi=1:length(entities[i].index)
                            indstr *= ")";
                        end
                    end
                    
                else
                    indstr = string(entities[i].index);
                    index_IR = entities[i].index;
                end
                if vors == "volume"
                    need_vol_interp_loop = true;
                    need_vol_coef_loop = true;
                    # Need to allocate NP 
                    np_allocate = IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :nodes_per_element]);
                    nqp_allocate = IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :qnodes_per_element]);
                    nodal_coef_name = "NODAL"*cname;
                    # For the nodal values
                    push!(allocate_part, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_64_data, Symbol(nodal_coef_name), [np_allocate], []),
                        np_allocate]));
                    # For the interpolated/differentiated quadrature values
                    push!(allocate_part, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_64_data, Symbol(cname), [nqp_allocate], []),
                        nqp_allocate]));
                    
                    push!(coef_loop_body, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_64_data, Symbol(nodal_coef_name), [np_allocate], [:ni]),
                        IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, index_IR, :nodeID])]));
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        quad_coef_node = IR_data_node(IRtypes.float_64_data, Symbol(cname), [nqp_allocate], [:col]);
                        nodal_coef_node = IR_data_node(IRtypes.float_64_data, Symbol(nodal_coef_name), [np_allocate], [:row]);
                        push!(coef_init_loop_body, IR_operation_node(IRtypes.assign_op, [quad_coef_node, 0.0]));
                        for di=1:length(entities[i].derivs)
                            deriv_quad_mat = IR_data_node(IRtypes.float_64_data, Symbol("RQ"*string(entities[i].derivs[di])), 
                                                            [:qnodes_per_element, :nodes_per_element], [col_row_matrix_index]);
                            push!(coef_interp_loop_body, IR_operation_node(IRtypes.assign_op,[
                                quad_coef_node,
                                IR_operation_node(IRtypes.math_op, [:(+), quad_coef_node, IR_operation_node(IRtypes.math_op, [:(*), deriv_quad_mat, nodal_coef_node])])
                                ]));
                        end
                    else
                        quad_coef_node = IR_data_node(IRtypes.float_64_data, Symbol(cname), [nqp_allocate], [:col]);
                        nodal_coef_node = IR_data_node(IRtypes.float_64_data, Symbol(nodal_coef_name), [np_allocate], [:row]);
                        # refelQ = IR_operation_node(IRtypes.member_op, [:refel, IR_data_node(IRtypes.array_data, :Q, [row_col_matrix_index])]);
                        refelQ = IR_data_node(IRtypes.float_64_data, :Q, [:nodes_per_element, :nodes_per_element], [col_row_matrix_index]);
                        push!(coef_init_loop_body, IR_operation_node(IRtypes.assign_op, [quad_coef_node, 0.0]));
                        push!(coef_interp_loop_body, IR_operation_node(IRtypes.assign_op,[
                            quad_coef_node,
                            IR_operation_node(IRtypes.math_op, [:(+), quad_coef_node, IR_operation_node(IRtypes.math_op, [:(*), refelQ, nodal_coef_node])])
                            ]));
                    end
                else # surface
                    # TODO
                end
            # elseif ctype == 4 # an indexed coefficient
                # TODO
            end
        end # if coefficient
    end # entity loop
    
    # Build the coef eval loop
    coef_eval_loop = IR_loop_node(IRtypes.space_loop, :nodes, :ni, 1, :nodes_per_element, coef_loop_body);
    
    push!(coef_init_loop_body, IR_loop_node(IRtypes.space_loop, :nodes, :row, 1, :nodes_per_element, coef_interp_loop_body))
    coef_interp_loop = IR_loop_node(IRtypes.space_loop, :qnodes, :col, 1, :qnodes_per_element, coef_init_loop_body);
    
    # If there are no needed coefficients, return empty lists
    if need_vol_coef_loop
        push!(coef_part, coef_eval_loop);
    end
    if need_vol_interp_loop
        push!(coef_part, coef_interp_loop);
    end
    
    return (allocate_part, coef_part);
end

function make_elemental_computation_fem(terms, var, dofsper, offset_ind, lorr, vors)
    # Here is where I make some assumption about the form of the expression.
    # Since it was expanded by the parser it should look like a series of terms: t1 + t2 + t3...
    # Where each term is multiplied by one test function component, and if LHS, involves one unknown component.
    # The submatrix modified by a term is determined by these, so go through the terms and divide them
    # into their submatrix expressions. 
    # Each term will look something like 
    # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
    # RHS: test_part * (weight_part .* coef_part)
    
    # LINALG_MATRIX_BLOCK and LINALG_VECTOR_BLOCK are special named ops that include these args
    # matrix: n_blocks, blockwidth, matrixname, i,j,term_IR, k,l,term_IR, ... up to n_blocks
    # vector: similar but one index per block
    
    IRtypes = IR_entry_types();
    
    # This will be returned
    compute_block = IR_block_node([],"elemental compute");
    
    # Separate the factors of each term into test, trial, coef and form the calculation
    if dofsper > 1
        # Submatrices or subvectors for each component
        if lorr == LHS
            submatrices = Array{Vector{IR_part}, 2}(undef, dofsper, dofsper);
            for i=1:dofsper
                for j=1:dofsper
                    submatrices[i,j] = [];
                end
            end
        else # RHS
            submatrices = Array{Vector{IR_part}, 1}(undef, dofsper);
            for j=1:dofsper
                submatrices[j] = [];
            end
        end
        
        if typeof(var) <: Array
            for vi=1:length(var) # variables
                # Process the terms for this variable
                for ci=1:length(terms[vi]) # components
                    for i=1:length(terms[vi][ci])
                        (test_part, trial_part, coef_part, test_ind, trial_ind) = generate_term_calculation_fem(terms[vi][ci][i], var, lorr, vors);
                        
                        # Turn these three parts into an expression like A'DB or A'Dv = A'd
                        # Where D is a diagonal matrix specified by a vector in the IR
                        if lorr == LHS
                            term_IR = IR_operation_node(IRtypes.named_op, [:LINALG_TDM, test_part, coef_part, trial_part]);
                        else
                            term_IR = IR_operation_node(IRtypes.named_op, [:LINALG_Tv, test_part, coef_part]);
                        end
                        
                        # Find the appropriate submatrix for this term
                        submati = offset_ind[vi] + test_ind;
                        submatj = trial_ind;
                        if lorr == LHS
                            submat_ind = submati + dofsper * (submatj-1);
                        else
                            submat_ind = submati;
                        end
                        
                        push!(submatrices[submat_ind], term_IR);
                    end
                end
                
            end # vi
            
        else # only one variable
            # Process the terms for this variable
            for ci=1:length(terms) # components
                for i=1:length(terms[ci])
                    (test_part, trial_part, coef_part, test_ind, trial_ind) = generate_term_calculation_fem(terms[ci][i], var, lorr, vors);
                    
                    # Turn these three parts into an expression like A'DB or A'Dv = A'd
                    # Where D is a diagonal matrix specified by a vector in the IR
                    if lorr == LHS
                        term_IR = IR_operation_node(IRtypes.named_op, [:LINALG_TDM, test_part, coef_part, trial_part]);
                    else
                        term_IR = IR_operation_node(IRtypes.named_op, [:LINALG_Tv, test_part, coef_part]);
                    end
                    
                    # Find the appropriate submatrix for this term
                    if lorr == LHS
                        submat_ind = test_ind + dofsper * (trial_ind-1);
                    else
                        submat_ind = test_ind;
                    end
                    
                    push!(submatrices[submat_ind], term_IR);
                end
            end
            
        end
        
        # Put the submatrices together into element_matrix or element_vector
        num_nonzero_blocks = 0;
        
        if lorr == LHS
            linalg_matrix_block_args = [];
            push!(linalg_matrix_block_args, :LINALG_MATMAT_BLOCKS);
            push!(linalg_matrix_block_args, 0);
            push!(linalg_matrix_block_args, :nodes_per_element);
            push!(linalg_matrix_block_args, :element_matrix);
            for smi=1:dofsper
                for smj=1:dofsper
                    submat_ind = smj + (smi-1)*dofsper;
                    if length(submatrices[smi,smj]) > 0
                        if length(submatrices[smi,smj]) > 1
                            new_term_vec = [];
                            push!(new_term_vec, :+);
                            append!(new_term_vec, submatrices[smi,smj]);
                            submat_rhs = IR_operation_node(IRtypes.math_op, new_term_vec);
                        else
                            submat_rhs = submatrices[smi,smj][1];
                        end
                        # push!(compute_block.parts, IR_operation_node(IRtypes.named_op, [
                        #                             :LINALG_MATMAT_BLOCKS, 
                        #                             smi, smj, :nodes_per_element, :element_matrix, submat_rhs]));
                        push!(linalg_matrix_block_args, smi);
                        push!(linalg_matrix_block_args, smj);
                        push!(linalg_matrix_block_args, submat_rhs);
                        
                        num_nonzero_blocks += 1;
                    end
                end
            end
            linalg_matrix_block_args[2] = num_nonzero_blocks;
            push!(compute_block.parts, IR_operation_node(IRtypes.named_op, linalg_matrix_block_args));
            
        else # RHS
            linalg_vector_block_args = [];
            push!(linalg_vector_block_args, :LINALG_MATVEC_BLOCKS);
            push!(linalg_vector_block_args, 0);
            push!(linalg_vector_block_args, :nodes_per_element);
            push!(linalg_vector_block_args, :element_vector);
            for smj=1:dofsper
                if length(submatrices[smj]) > 0
                    if length(submatrices[smj]) > 1
                        new_term_vec = [];
                        push!(new_term_vec, :+);
                        append!(new_term_vec, submatrices[smj]);
                        submat_rhs = IR_operation_node(IRtypes.math_op, new_term_vec);
                    else
                        submat_rhs = submatrices[smj][1];
                    end
                    
                    # push!(compute_block.parts, IR_operation_node(IRtypes.named_op, [
                    #                             :LINALG_MATVEC_BLOCKS, 
                    #                             smj, :nodes_per_element, :element_vector, submat_rhs]));
                    push!(linalg_vector_block_args, smj);
                    push!(linalg_vector_block_args, submat_rhs);
                    
                    num_nonzero_blocks += 1;
                end
            end
            linalg_vector_block_args[2] = num_nonzero_blocks;
            push!(compute_block.parts, IR_operation_node(IRtypes.named_op, linalg_vector_block_args));
        end
        
        
    else # one dof
        terms = terms[1][1];
        term_vec = Vector{IR_part}(undef,0);
        
        #process each term
        for i=1:length(terms)
            (test_part, trial_part, coef_part, test_ind, trial_ind) = generate_term_calculation_fem(terms[i], var, lorr, vors);
            
            # Turn these three parts into an expression like A'DB or A'Dv = A'd
            # Where D is a diagonal matrix specified by a vector in the IR
            # Note T=transpose matrix, D=diagonal matrix, M=matrix, v=vector, t=transpose vector
            if lorr == LHS
                term_IR = IR_operation_node(IRtypes.named_op, [:LINALG_TDM, test_part, coef_part, trial_part]);
            else
                term_IR = IR_operation_node(IRtypes.named_op, [:LINALG_Tv, test_part, coef_part]);
            end
            
            push!(term_vec, term_IR);
            
        end
        if length(term_vec) > 1
            new_term_vec = [];
            push!(new_term_vec, :+);
            append!(new_term_vec, term_vec);
            combined_rhs = IR_operation_node(IRtypes.math_op, new_term_vec);
            if lorr == LHS
                push!(compute_block.parts, IR_operation_node(IRtypes.named_op, [
                    :LINALG_MATMAT_BLOCKS, 1, :nodes_per_element, :element_matrix, 1, 1, combined_rhs]));
            else
                push!(compute_block.parts, IR_operation_node(IRtypes.named_op, [ 
                    :LINALG_MATVEC_BLOCKS, 1, :nodes_per_element, :element_vector, 1, combined_rhs]));
            end
        else
            if lorr == LHS
                push!(compute_block.parts, IR_operation_node(IRtypes.named_op, [
                    :LINALG_MATMAT_BLOCKS, 1, :nodes_per_element, :element_matrix, 1, 1, term_vec[1]]));
            else
                push!(compute_block.parts, IR_operation_node(IRtypes.named_op, [ 
                    :LINALG_MATVEC_BLOCKS, 1, :nodes_per_element, :element_vector, 1, term_vec[1]]));
            end
        end
        
    end
    
    return compute_block;
end

# This takes a term expression that should have a form like test_part * (coef_parts) * trial_part
# The parts are separated, test and trial parts are translated into quadrature matrices (refel.Q or RQn)
# and they are returned as IR_parts
function generate_term_calculation_fem(term, var, lorr, vors)
    IRtypes = IR_entry_types();
    
    if lorr == LHS
        (test_ex, trial_ex, coef_ex, test_ind, trial_ind) = separate_factors(term, var);
    else
        (test_ex, trial_ex, coef_ex, test_ind, trial_ind) = separate_factors(term);
    end
    
    test_negative = false;
    if (typeof(test_ex) == Expr) && (test_ex.head == :call && (test_ex.args[1] == :- || test_ex.args[1] == :.-) && length(test_ex.args) == 2)
        test_ex = test_ex.args[2];
        test_negative = true;
    end
    trial_negative = false;
    if (typeof(trial_ex) == Expr) && (trial_ex.head == :call && (trial_ex.args[1] == :- || trial_ex.args[1] == :.-) && length(trial_ex.args) == 2)
        trial_ex = trial_ex.args[2];
        trial_negative = true;
        if test_negative
            test_negative = false;
            trial_negative = false;
        end
    end
    
    # Determine the matrix corresponding to test and trial
    if typeof(test_ex) == SymEntity
        if length(test_ex.derivs) > 0
            deriv_index = test_ex.derivs[1];
            test_part = IR_data_node(IRtypes.float_64_data, Symbol("RQ"*string(deriv_index)), [:qnodes_per_element, :nodes_per_element], []);
        else
            # test_part = IR_operation_node(IRtypes.member_op, [:refel, IR_data_node(IRtypes.matrix_data, :Q, [:qnodes_per_element, :nodes_per_element])]);
            test_part = IR_data_node(IRtypes.float_64_data, :Q, [:qnodes_per_element, :nodes_per_element], []);
        end
    else
        test_part = nothing;
    end
    if typeof(trial_ex) == SymEntity
        if length(trial_ex.derivs) > 0
            deriv_index = trial_ex.derivs[1];
            trial_part = IR_data_node(IRtypes.float_64_data, Symbol("RQ"*string(deriv_index)), [:qnodes_per_element, :nodes_per_element], []);
        else
            # trial_part = IR_operation_node(IRtypes.member_op, [:refel, IR_data_node(IRtypes.matrix_data, :Q, [:qnodes_per_element, :nodes_per_element])]);
            trial_part = IR_data_node(IRtypes.float_64_data, :Q, [:qnodes_per_element, :nodes_per_element], []);
        end
    else
        trial_part = nothing;
    end
    
    # Turn the coefficient part into IR
    wg_part = IR_data_node(IRtypes.float_64_data, :wg, [:qnodes_per_element], []);
    if false # non-constant J
        detj_part = IR_data_node(IRtypes.float_64_data, :detj, [:qnodes_per_element], []);
    else
        detj_part = IR_data_node(IRtypes.float_64_data, :detj);
    end
    if !(coef_ex === nothing)
        coef_part = arithmetic_expr_to_IR(coef_ex);
        if trial_negative || test_negative
            # coef_part = IR_operation_node(IRtypes.math_op, [:*, IR_operation_node(IRtypes.member_op, [:refel, wg_part]), :detj, coef_part, -1]);
            coef_part = IR_operation_node(IRtypes.math_op, [:*, wg_part, detj_part, coef_part, -1]);
        else
            # coef_part = IR_operation_node(IRtypes.math_op, [:*, IR_operation_node(IRtypes.member_op, [:refel, wg_part]), :detj, coef_part]);
            coef_part = IR_operation_node(IRtypes.math_op, [:*, wg_part, detj_part, coef_part]);
        end
        
    else
        if trial_negative || test_negative
            # coef_part = IR_operation_node(IRtypes.math_op, [:*, IR_operation_node(IRtypes.member_op, [:refel, wg_part]), :detj, -1]);
            coef_part = IR_operation_node(IRtypes.math_op, [:*, wg_part, detj_part, -1]);
        else
            # coef_part = IR_operation_node(IRtypes.math_op, [:*, IR_operation_node(IRtypes.member_op, [:refel, wg_part]), :detj]);
            coef_part = IR_operation_node(IRtypes.math_op, [:*, wg_part, detj_part]);
        end
        
    end
    
    return (test_part, trial_part, coef_part, test_ind, trial_ind);
end

# This part inserts the elemental matrix and vector into the global one.
# Since this is target dependent, it should not be generated automatically,
# but rather through a named operation.
function generate_local_to_global_fem(dofs_per_node; vec_only=false)
    IRtypes = IR_entry_types();
    result_block = IR_block_node([]);
    vector_loop_body = IR_block_node([]);
    matrix_loop_body = IR_block_node([]);
    if dofs_per_node == 1
        # next_ind = +(1, *(*(-(eid, 1), nodes_per_element), nodes_per_element));
        # for ni = 1:nodes_per_element
        #     glb_i = mesh.loc2glb[ni, eid];
        #     global_vector[glb_i] = +(global_vector[glb_i], element_vector[ni]);
        #     for nj = 1:nodes_per_element
        #         glb_j = mesh.loc2glb[nj, eid];
        #         global_matrix_I[next_ind] = glb_i;
        #         global_matrix_J[next_ind] = glb_j;
        #         global_matrix_V[next_ind] = element_matrix[ni, nj];
        #         next_ind = +(next_ind, 1);
        #     end
        # end
        
        # only needed if not vec_only
        if !vec_only
            push!(result_block.parts, IR_operation_node(IRtypes.assign_op, [
                :next_ind,
                IR_operation_node(IRtypes.math_op, [:+, 1,
                    IR_operation_node(IRtypes.math_op, [:*, 
                        IR_operation_node(IRtypes.math_op, [:-, :eid, 1]), :nodes_per_element, :nodes_per_element])])
            ]));
        end
        
        push!(vector_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
            :glb_i,
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.int_64_data, :loc2glb, [:nodes_per_element, :num_elements], [:ni, :eid])])
        ]));
        push!(vector_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.float_64_data, :global_vector, [:dofs_global], [:glb_i]),
            IR_operation_node(IRtypes.math_op, [:+, 
                IR_data_node(IRtypes.float_64_data, :global_vector, [:dofs_global], [:glb_i]),
                IR_data_node(IRtypes.float_64_data, :element_vector, [:dofs_per_element], [:ni])])
        ]));
        
        # only to the matrix part if not vec_only
        if !vec_only
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                :glb_j,
                IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.int_64_data, :loc2glb, [:nodes_per_element, :num_elements], [:nj, :eid])])
            ]));
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.int_64_data, :global_matrix_I, [:allocated_nonzeros], [:next_ind]),
                :glb_i
            ]));
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.int_64_data, :global_matrix_J, [:allocated_nonzeros], [:next_ind]),
                :glb_j
            ]));
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :global_matrix_V, [:allocated_nonzeros], [:next_ind]),
                IR_data_node(IRtypes.float_64_data, :element_matrix, [:dofs_per_element, :dofs_per_element], [:ni, :nj])
            ]));
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                :next_ind,
                IR_operation_node(IRtypes.math_op, [:+, :next_ind, 1])
            ]));
            
            push!(vector_loop_body.parts, IR_loop_node(IRtypes.space_loop, :nodes, :nj, 1, :nodes_per_element, matrix_loop_body));
        end
        
        toglobal_loop = IR_loop_node(IRtypes.space_loop, :nodes, :ni, 1, :nodes_per_element, vector_loop_body);
        push!(result_block.parts, toglobal_loop)
        
    else
        # next_ind = +(1, *(*(-(eid, 1), dofs_per_element), dofs_per_element)) + nodes_per_element*nodes_per_element*index_offset;
        # for ni = 1:nodes_per_element
        #     glb_i = mesh.loc2glb[ni, eid];
        #     for di = 1:dofs_per_node
        #         glb_dofi = (glb_i - 1)*dofs_per_node + di + index_offset;
        #         global_vector[glb_dofi] = +(global_vector[glb_dofi], element_vector[(di-1)*nodes_per_element + ni]);
        #         for nj = 1:nodes_per_element
        #             glb_j = mesh.loc2glb[nj, eid];
        #             for dj = 1:dofs_per_node
        #                 glb_dofj = (glb_j - 1)*dofs_per_node + dj + index_offset;
        #                 global_matrix_I[next_ind] = glb_dofi;
        #                 global_matrix_J[next_ind] = glb_dofj;
        #                 global_matrix_V[next_ind] = element_matrix[(di-1)*nodes_per_element + ni, (dj-1)*nodes_per_element + nj];
        #                 next_ind = +(next_ind, 1);
        #             end
        #         end
        #     end
        # end
        vector_dof_loop_body = IR_block_node([]);
        matrix_dof_loop_body = IR_block_node([]);
        
        # only needed if not vec_only
        if !vec_only
            push!(result_block.parts, IR_operation_node(IRtypes.assign_op, [
                :next_ind,
                IR_operation_node(IRtypes.math_op, [:+, 1,
                    IR_operation_node(IRtypes.math_op, [:*, 
                        IR_operation_node(IRtypes.math_op, [:-, :eid, 1]), :dofs_per_element, :dofs_per_element]),
                    IR_operation_node(IRtypes.math_op, [:*, :nodes_per_element, :nodes_per_element, :index_offset])
                ])
            ]));
        end
        
        push!(vector_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
            :glb_i,
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.int_64_data, :loc2glb, [:nodes_per_element, :num_elements], [:ni, :eid])])
        ]));
        push!(vector_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
            :glb_dofi,
            generate_the_one_pattern(:glb_i, 1, :dofs_per_node, IR_operation_node(IRtypes.math_op, [:+, :di, :index_offset]))
        ]));
        push!(vector_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.float_64_data, :global_vector, [:dofs_global], [:glb_dofi]),
            IR_operation_node(IRtypes.math_op, [:+, 
                IR_data_node(IRtypes.float_64_data, :global_vector, [:dofs_global], [:glb_dofi]),
                IR_data_node(IRtypes.float_64_data, :element_vector, [:dofs_per_element], [generate_the_one_pattern(:di, 1, :nodes_per_element, :ni)])])
        ]));
        
        # only needed if not vec_only
        if !vec_only
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                :glb_j,
                IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.int_64_data, :loc2glb, [:nodes_per_element, :num_elements], [:nj, :eid])])
            ]));
            push!(matrix_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                :glb_dofj,
                generate_the_one_pattern(:glb_j, 1, :dofs_per_node, IR_operation_node(IRtypes.math_op, [:+, :dj, :index_offset]))
            ]));
            push!(matrix_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.int_64_data, :global_matrix_I, [:allocated_nonzeros], [:next_ind]),
                :glb_dofi
            ]));
            push!(matrix_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.int_64_data, :global_matrix_J, [:allocated_nonzeros], [:next_ind]),
                :glb_dofj
            ]));
            push!(matrix_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :global_matrix_V, [:allocated_nonzeros], [:next_ind]),
                IR_data_node(IRtypes.float_64_data, :element_matrix, [:dofs_per_element, :dofs_per_element], [
                    generate_the_one_pattern(:di, 1, :nodes_per_element, :ni),
                    generate_the_one_pattern(:dj, 1, :nodes_per_element, :nj)])
            ]));
            push!(matrix_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                :next_ind,
                IR_operation_node(IRtypes.math_op, [:+, :next_ind, 1])
            ]));
            
            push!(matrix_loop_body.parts, IR_loop_node(IRtypes.space_loop, :nodes, :dj, 1, :dofs_per_loop, matrix_dof_loop_body));
            push!(vector_dof_loop_body.parts, IR_loop_node(IRtypes.space_loop, :nodes, :nj, 1, :nodes_per_element, matrix_loop_body));
        end
        
        push!(vector_loop_body.parts, IR_loop_node(IRtypes.space_loop, :nodes, :di, 1, :dofs_per_loop, vector_dof_loop_body));
        toglobal_loop = IR_loop_node(IRtypes.space_loop, :nodes, :ni, 1, :nodes_per_element, vector_loop_body);
        push!(result_block.parts, toglobal_loop)
        
        toglobal_loop = IR_loop_node(IRtypes.space_loop, :nodes, :ni, 1, :nodes_per_element, vector_loop_body);
    end
    
    return result_block;
end

# Generate the assembly loop structures and insert the content
function generate_assembly_loop_fem(var, indices)
    IRtypes = IR_entry_types();
    
    # Make names for all of the index variables like INDEX_VAL_something
    index_names = [];
    elements_included = false;
    ind_shift = 0;
    for i=1:length(indices)
        if indices[i] == "elements" || indices[i] == "cells"
            elements_included = true;
            push!(index_names, "elements");
        else
            push!(index_names, IR_data_node(IRtypes.int_64_data, Symbol("INDEX_VAL_"*string(indices[i].symbol))));
        end
    end
    
    # If elements were not included, make them the outermost loop
    if !elements_included && length(index_names) > 0
        index_names = ["elements"; index_names];
        ind_shift = 1;
    end
    
    # Make an offset that is set at the beginning of the loop body
    # We have to assume that unknown variables have the same set of indices.
    # If they don't, we need much more complicated logic.
    index_offset = 0; # offset in variable values for this index
    if var[1].indexer === nothing
        # no indexer. 
    else
        # coul be in an array, or just a single indexer
        if typeof(var[1].indexer) <: Array
            index_offset = Symbol("INDEX_VAL_"*string(var[1].indexer[1].symbol));
            prev_ind = length(var[1].indexer[1].range);
            for i=2:length(var[1].indexer)
                index_offset = IR_operation_node(IRtypes.math_op, [:+, index_offset,
                    IR_operation_node(IRtypes.math_op, [:*, prev_ind, 
                        IR_operation_node(IRtypes.math_op, [:-, Symbol("INDEX_VAL_"*string(var[1].indexer[i].symbol)), 1])])]);
                prev_ind *= length(var[1].indexer[i].range);
            end
            index_offset = IR_operation_node(IRtypes.math_op, [:-, index_offset, 1]);
        else
            index_offset = IR_operation_node(IRtypes.math_op, [:-, Symbol("INDEX_VAL_"*string(var[1].indexer.symbol)), 1]);
        end
    end
    
    # Placeholder computation that will be replaced by the actual one
    placeholder = IR_block_node([
        IR_operation_node(IRtypes.assign_op, [:index_offset, index_offset]);
    ]);
    
    # generate the loop structures
    if length(index_names) > 0
        # The innermost loop holds placeholder
        if index_names[end] == "elements"
            assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :eid, 1, :num_elements, placeholder);
        else
            assembly_loop = IR_loop_node(IRtypes.index_loop, indices[end].symbol, 
                            index_names[end].label, indices[end].range[1], indices[end].range[end], placeholder);
        end
        
        # work outwards nesting assembly_loop
        for i=(length(index_names)-1):-1:1
            if index_names[i] == "elements"
                assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :eid, 1, :num_elements, assembly_loop);
            elseif i > ind_shift
                assembly_loop = IR_loop_node(IRtypes.index_loop, indices[i-ind_shift].symbol, 
                                index_names[i].label, indices[i-ind_shift].range[1], 
                                indices[i-ind_shift].range[end], assembly_loop);
            end
        end
    else # only an element loop
        assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :eid, 1, :num_elements, placeholder);
    end
    
    return assembly_loop
end

# Generates the time stepping loop using the supplied assembly block and stepper
function generate_time_stepping_loop_fem(stepper, assembly, include_var_update=true)
    IRtypes = IR_entry_types();
    tloop_body = IR_block_node([]);
    zero_el_vec = IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY, IR_data_node(IRtypes.float_64_data, :global_vector, [:dofs_global], []), 0, :dofs_global]);
    zero_bdry_done = IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY, IR_data_node(IRtypes.int_64_data, :bdry_done, [:nnodes_global], []), 0, :nnodes_global]);
    # If pre or post-step functions have been set, include those
    if !(prob.pre_step_function === nothing)
        pre_step_call = IR_operation_node(IRtypes.function_op, [:pre_step_function]);
    else
        pre_step_call = IR_comment_node("No pre-step function specified");
    end
    if !(prob.post_step_function === nothing)
        post_step_call = IR_operation_node(IRtypes.function_op, [:post_step_function]);
    else
        post_step_call = IR_comment_node("No post-step function specified");
    end
    var_update = IR_comment_node("");
    allocate_block = IR_block_node([]);
    if stepper.stages < 2
        if include_var_update
            var_update = wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :solution]));
        end
        
        # # This part would be used if the solution is updated like sol = sol + dt*()
        # sol_i = IR_data_node(IRtypes.array_data, :solution, [:update_i]);
        # dsol_i = IR_data_node(IRtypes.array_data, :d_solution, [:update_i]);
        # update_loop = IR_loop_node(IRtypes.dof_loop, :dof, :update_i, 1, :dofs_global, IR_block_node([
        #     IR_operation_node(IRtypes.assign_op, [
        #         sol_i,
        #         IR_operation_node(IRtypes.math_op, [:+, sol_i, IR_operation_node(IRtypes.math_op, [:*, time_stepper.dt, dsol_i])])
        #     ])
        # ]))
        # tloop_body.parts = [
        #     wrap_in_timer(:step_assembly, assembly),
        #     wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :d_solution, :global_matrix, :global_vector])),
        #     wrap_in_timer(:update_sol, update_loop),
        #     wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :solution])),
        #     IR_operation_node(IRtypes.assign_op, [
        #         :t,
        #         IR_operation_node(IRtypes.math_op, [:+, :t, stepper.dt])
        #     ])
        # ];
        
        # This part is used if the update is coded into the system sol = A\b (not sol = sol + dt*(A\b))
        tloop_body.parts = [
            zero_el_vec,
            zero_bdry_done,
            pre_step_call,
            wrap_in_timer(:step_assembly, assembly),
            wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector])),
            var_update,
            post_step_call,
            IR_operation_node(IRtypes.assign_op, [
                :t,
                IR_operation_node(IRtypes.math_op, [:+, :t, stepper.dt])
            ])
        ];
        
        time_loop = IR_loop_node(IRtypes.time_loop, :time, :ti, 1, stepper.Nsteps, tloop_body);
        
    else # multistage explicit steppers
        # LSRK4 is a special case, low storage
        if stepper.type == LSRK4
            # Low storage RK4: 
            # p0 = u
            #   ki = ai*k(i-1) + dt*f(p(i-1), t+ci*dt)
            #   pi = p(i-1) + bi*ki
            # u = p5
            if include_var_update
                var_update = wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :tmppi]));
            end
            na = length(stepper.a);
            nb = length(stepper.b);
            nc = length(stepper.c);
            tmp_pi = IR_data_node(IRtypes.float_64_data, :tmppi, [:dofs_global], []);
            tmp_ki = IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], []);
            piki_loop_one = IR_loop_node(IRtypes.space_loop, :dofs, :piki_i, 1, :dofs_global, IR_block_node([
                # tmpki .= stepper.dt .* sol;
                # tmppi .= tmppi + stepper.b[rki].*tmpki;
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], [:piki_i]),
                    IR_operation_node(IRtypes.math_op, [:*, :dt, IR_data_node(IRtypes.float_64_data, :solution, [:dofs_global], [:piki_i])])
                ]),
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_64_data, :tmppi, [:dofs_global], [:piki_i]),
                    IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.float_64_data, :tmppi, [:dofs_global], [:piki_i]), 
                        IR_operation_node(IRtypes.math_op, [:*, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :b, [nb], [:rki])]), 
                            IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], [:piki_i])])
                    ])
                ])
            ]))
            piki_loop_two = IR_loop_node(IRtypes.space_loop, :dofs, :piki_i, 1, :dofs_global, IR_block_node([
                # tmpki .= stepper.a[rki].*tmpki + stepper.dt .* sol;
                # tmppi .= tmppi + stepper.b[rki].*tmpki;
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], [:piki_i]),
                    IR_operation_node(IRtypes.math_op, [:+, 
                        IR_operation_node(IRtypes.math_op, [:*, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :a, [na], [:rki])]), 
                            IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], [:piki_i])]),
                        IR_operation_node(IRtypes.math_op, [:*, :dt, IR_data_node(IRtypes.float_64_data, :solution, [:dofs_global], [:piki_i])])
                    ])
                ]),
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_64_data, :tmppi, [:dofs_global], [:piki_i]),
                    IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.float_64_data, :tmppi, [:dofs_global], [:piki_i]), 
                        IR_operation_node(IRtypes.math_op, [:*, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :b, [nb], [:rki])]), 
                            IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], [:piki_i])])
                    ])
                ])
            ]))
            piki_condition = IR_conditional_node(IR_operation_node(IRtypes.math_op, [:<, :rki, 2]),
                IR_block_node([piki_loop_one]),
                IR_block_node([piki_loop_two]));
            
            stage_loop_body = IR_block_node([
                # last_t = t;
                # t = last_t + stepper.c[rki]*stepper.dt;
                IR_operation_node(IRtypes.assign_op, [
                    :t,
                    IR_operation_node(IRtypes.math_op, [:+, :t, 
                        IR_operation_node(IRtypes.math_op, [:*, :dt, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :c, [nc], [:rki])])])])
                ]),
                # assemble
                zero_el_vec,
                zero_bdry_done,
                pre_step_call,
                wrap_in_timer(:step_assembly, assembly),
                # solve
                wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector])),
                # copy_bdry_vals_to_variables(var, sol, grid_data, dofs_per_node, true);
                IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :solution, :var, :true]),
                # update tmppi and tmpki
                piki_condition,
                # copy_bdry_vals_to_vector(var, tmppi, grid_data, dofs_per_node);
                IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :tmppi]),
                var_update,
                post_step_call
            ]);
            stage_loop = IR_loop_node(IRtypes.time_loop, :stages, :rki, 1, stepper.stages, stage_loop_body);
            
            push!(tloop_body.parts, IR_operation_node(IRtypes.named_op, [:GATHER_SOLUTION, :tmppi]));
            push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:last_t, :t]));
            push!(tloop_body.parts, stage_loop);
            push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:t, IR_operation_node(IRtypes.math_op, [:+, :last_t, :dt])]));
            
            allocate_block = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [
                    tmp_pi,
                    IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])]),
                IR_operation_node(IRtypes.assign_op, [
                    tmp_ki,
                    IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])]),
            ]);
            time_loop = IR_loop_node(IRtypes.time_loop, :time, :ti, 1, stepper.Nsteps, tloop_body);
            
        else
            # Explicit multi-stage methods: 
            # x = x + dt*sum(bi*ki)
            # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
            if include_var_update
                var_update = wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :tmpresult]));
                var_update2 = wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :last_result]));
            else
                var_update2 = IR_comment_node("");
            end
            na = length(stepper.a);
            nb = length(stepper.b);
            nc = length(stepper.c);
            tmp_last = IR_data_node(IRtypes.float_64_data, :last_result, [:dofs_global], []);
            tmp_result = IR_data_node(IRtypes.float_64_data, :tmpresult, [:dofs_global], []);
            tmp_ki = IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], []);
            
            update_ki_loop_one = IR_loop_node(IRtypes.space_loop, :dofs, :k, 1, :dofs_global, IR_block_node([
                # # Update the values in vars to be used in the next stage
                # for k=1:dofs_global
                #     tmpki[k,stage] = solution[k];
                #     tmpresult[k] = last_result[k];
                #     for j=1:stage
                #         tmpresult[k] += stepper.dt * stepper.a[stage+1, j] * tmpki[k,j];
                #     end
                # end
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], [:k, :rki]),
                    IR_data_node(IRtypes.float_64_data, :solution, [:dofs_global], [:k])
                ]),
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_64_data, :tmpresult, [:dofs_global], [:k]),
                    IR_data_node(IRtypes.float_64_data, :last_result, [:dofs_global], [:k])
                ]),
                IR_loop_node(IRtypes.time_loop, :stage, :j, 1, :rki, IR_block_node([
                    IR_operation_node(IRtypes.assign_op, [
                        IR_data_node(IRtypes.float_64_data, :tmpresult, [:dofs_global], [:k]),
                        IR_operation_node(IRtypes.math_op, [:+,
                            IR_data_node(IRtypes.float_64_data, :tmpresult, [:dofs_global], [:k]),
                            IR_operation_node(IRtypes.math_op, [:*, :dt, 
                                IR_operation_node(IRtypes.member_op, [:time_stepper, 
                                    IR_data_node(IRtypes.float_64_data, :a, [na], [IR_operation_node(IRtypes.math_op, [:+, :rki, 1]), :j])]),
                                IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], [:k, :j])
                            ])
                        ])
                    ])
                ]))
            ]))
            update_ki_loop_two = IR_loop_node(IRtypes.space_loop, :dofs, :k, 1, :dofs_global, IR_block_node([
                # for k=1:dofs_global
                #     tmpki[k,stage] = solution[k];
                # end
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], [:k, :rki]),
                    IR_data_node(IRtypes.float_64_data, :solution, [:dofs_global], [:k])
                ])
            ]))
            update_ki_condition = IR_conditional_node(IR_operation_node(IRtypes.math_op, [:<, :rki, stepper.stages]),
                IR_block_node([
                    update_ki_loop_one,
                    # copy_bdry_vals_to_vector(var, tmpresult, grid_data, dofs_per_node);
                    IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :tmpresult]),
                    # place_sol_in_vars(var, tmpresult, stepper);
                    var_update
                ]),
                IR_block_node([update_ki_loop_two]));
            
            stage_loop_body = IR_block_node([
                # last_t = t;
                # t = last_t + stepper.c[rki]*stepper.dt;
                IR_operation_node(IRtypes.assign_op, [
                    :t,
                    IR_operation_node(IRtypes.math_op, [:+, :t, 
                        IR_operation_node(IRtypes.math_op, [:*, :dt, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :c, [nc], [:rki])])])])
                ]),
                # assemble
                zero_el_vec,
                zero_bdry_done,
                pre_step_call,
                wrap_in_timer(:step_assembly, assembly),
                # solve
                wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector])),
                # copy_bdry_vals_to_variables(var, sol, grid_data, dofs_per_node, true);
                IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :solution, :var, :true]),
                # update tmpki
                update_ki_condition,
                post_step_call
            ]);
            stage_loop = IR_loop_node(IRtypes.time_loop, :stages, :rki, 1, stepper.stages, stage_loop_body);
            
            # last_result[k] += stepper.dt * stepper.b[rki] * tmpki[k, rki];
            combine_loop_body = IR_block_node([
                IR_loop_node(IRtypes.time_loop, :stages, :rki, 1, stepper.stages, IR_block_node([
                    IR_operation_node(IRtypes.assign_op, [
                        IR_data_node(IRtypes.float_64_data, :last_result, [:dofs_global], [:k]),
                        IR_operation_node(IRtypes.math_op, [:+,
                            IR_data_node(IRtypes.float_64_data, :last_result, [:dofs_global], [:k]), 
                            IR_operation_node(IRtypes.math_op, [:*, :dt, 
                                IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :b, [nb], [:rki])]),
                                IR_data_node(IRtypes.float_64_data, :tmpki, [:dofs_global], [:k, :rki])
                            ])
                        ])
                    ])
                ]))
            ])
            combine_loop = IR_loop_node(IRtypes.space_loop, :dofs, :k, 1, :dofs_global, combine_loop_body);
            
            push!(tloop_body.parts, IR_operation_node(IRtypes.named_op, [:GATHER_SOLUTION, :last_result]));
            push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:last_t, :t]));
            push!(tloop_body.parts, stage_loop);
            push!(tloop_body.parts, combine_loop);
            # copy_bdry_vals_to_vector(var, last_result, grid_data, dofs_per_node);
            push!(tloop_body.parts, IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :last_result]))
            # place_sol_in_vars(var, last_result, stepper);
            push!(tloop_body.parts, var_update2);
            push!(tloop_body.parts, post_step_call);
            # update time
            push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:t, IR_operation_node(IRtypes.math_op, [:+, :last_t, :dt])]));
            
            allocate_block = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [
                    tmp_last,
                    IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])]),
                IR_operation_node(IRtypes.assign_op, [
                    tmp_result,
                    IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])]),
                IR_operation_node(IRtypes.assign_op, [
                    tmp_ki,
                    IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global, stepper.stages])])
            ]);
            time_loop = IR_loop_node(IRtypes.time_loop, :time, :ti, 1, stepper.Nsteps, tloop_body);
            
        end
    end
    
    return (allocate_block, time_loop);
end

# TDM(A,b,C) = transpose(A) * diagm(b) * C
# P_ij =  A'_ik * b_k * C_kj = A_ki * b_k * C_kj
# return the IR for A[k,i] * b[k] * C[k,j]
# A,b,C are symbols or IR_data_node
# i,j,k are symbols, numbers, IR_part
function generate_linalg_TDM_product(A, b, C, i, j, k)
    IRtypes = IR_entry_types();
    if typeof(A) == IR_data_node
        A_part = IR_data_node(IRtypes.float_64_data, A.label, [:?,:?], [k,i]);
    elseif typeof(A) <: IR_part
        A_part = apply_indexed_access(A, [k,i], IRtypes);
    else
        A_part = IR_data_node(IRtypes.float_64_data, A, [:?,:?], [k,i]);
    end
    if typeof(b) == IR_data_node
        b_part = IR_data_node(IRtypes.float_64_data, b.label, [:?], [k]);
    elseif typeof(b) <: IR_part
        b_part = apply_indexed_access(b, [k], IRtypes);
    else
        b_part = IR_data_node(IRtypes.float_64_data, b, [:?], [k]);
    end
    if typeof(C) == IR_data_node
        C_part = IR_data_node(IRtypes.float_64_data, C.label, [:?,:?], [k,j]);
    elseif typeof(C) <: IR_part
        C_part = apply_indexed_access(C, [k,j], IRtypes);
    else
        C_part = IR_data_node(IRtypes.float_64_data, C, [:?,:?], [k,j]);
    end
    
    return IR_operation_node(IRtypes.math_op, [:*, A_part, b_part, C_part]);
end

# Tv(A,b) = transpose(A) * b
# P_i =  A'_ij * b_j = A_ji * b_j
# return the IR for A[j,i] * b[k]
# A,b are symbols or IR_data_node
# i,j are symbols, numbers, IR_parts
function generate_linalg_Tv_product(A, b, i, j)
    IRtypes = IR_entry_types();
    if typeof(A) == IR_data_node
        A_part = IR_data_node(IRtypes.float_64_data, A.label, [:?,:?], [j,i]);
    elseif typeof(A) <: IR_part
        A_part = apply_indexed_access(A, [j,i], IRtypes);
    else
        A_part = IR_data_node(IRtypes.float_64_data, A, [:?,:?], [j,i]);
    end
    if typeof(b) == IR_data_node
        b_part = IR_data_node(IRtypes.float_64_data, b.label, [:?], [j]);
    elseif typeof(b) <: IR_part
        b_part = apply_indexed_access(b, [j], IRtypes);
    else
        b_part = IR_data_node(IRtypes.float_64_data, b, [:?], [j]);
    end
    
    return IR_operation_node(IRtypes.math_op, [:*, A_part, b_part]);
end

# For array a, b and scalar c,d do
# a = a + b
# c = norm(b)
# d = c / norm(a)
function generate_residual_norms_and_update(vec_a, vec_b, val_c, val_d)
    IRtypes = IR_entry_types();
    # tmp = (vec_a[norm_i] - vec_b[norm_i])
    vec_a_i = apply_indexed_access(vec_a, [:norm_i], IRtypes);
    vec_b_i = apply_indexed_access(vec_b, [:norm_i], IRtypes);
    
    block = IR_block_node([
        # zero c
        IR_operation_node(IRtypes.assign_op, [val_c, 0]),
        IR_operation_node(IRtypes.assign_op, [val_d, 1e-20]),
        # loop over vector
        IR_loop_node(IRtypes.index_loop, :dofs, :norm_i, 1, :dofs_global, IR_block_node([
            # update a = a + b
            IR_operation_node(IRtypes.assign_op, [vec_a_i, 
                IR_operation_node(IRtypes.math_op, [:+, vec_a_i, vec_b_i])
            ]),
            # add to norm val_c = val_c + (vec_a[norm_i])^2
            IR_operation_node(IRtypes.assign_op, [val_c, 
                IR_operation_node(IRtypes.math_op, [:+, val_c, 
                    IR_operation_node(IRtypes.math_op, [:*, vec_b_i, vec_b_i])
                ])
            ]),
            IR_operation_node(IRtypes.assign_op, [val_d, 
                IR_operation_node(IRtypes.math_op, [:+, val_d, 
                    IR_operation_node(IRtypes.math_op, [:*, vec_a_i, vec_a_i])
                ])
            ])
        ])),
        # sqrt(val_c)
        # c / sqrt(val_d)
        IR_operation_node(IRtypes.assign_op, [val_c, IR_operation_node(IRtypes.math_op, [:sqrt, val_c])]),
        IR_operation_node(IRtypes.assign_op, [val_d, IR_operation_node(IRtypes.math_op, [:sqrt, val_d])]),
        IR_operation_node(IRtypes.assign_op, [val_d, IR_operation_node(IRtypes.math_op, [:/, val_c, val_d])])
    ])
    
    return block;
end

# For array a, b and scalar c,d do
# c = norm(a)
# d = c / norm(b)
function generate_residual_norms(vec_a, vec_b, val_c, val_d)
    IRtypes = IR_entry_types();
    # tmp = (vec_a[norm_i] - vec_b[norm_i])
    vec_a_i = apply_indexed_access(vec_a, [:norm_i], IRtypes);
    vec_b_i = apply_indexed_access(vec_b, [:norm_i], IRtypes);
    
    block = IR_block_node([
        # zero c
        IR_operation_node(IRtypes.assign_op, [val_c, 0]),
        IR_operation_node(IRtypes.assign_op, [val_d, 1e-20]),
        # loop over vector
        IR_loop_node(IRtypes.index_loop, :dofs, :norm_i, 1, :dofs_global, IR_block_node([
            # add to norm val_c = val_c + (vec_a[norm_i])^2
            IR_operation_node(IRtypes.assign_op, [val_c, 
                IR_operation_node(IRtypes.math_op, [:+, val_c, 
                    IR_operation_node(IRtypes.math_op, [:*, vec_a_i, vec_a_i])
                ])
            ]),
            IR_operation_node(IRtypes.assign_op, [val_d, 
                IR_operation_node(IRtypes.math_op, [:+, val_d, 
                    IR_operation_node(IRtypes.math_op, [:*, vec_b_i, vec_b_i])
                ])
            ])
        ])),
        # sqrt(val_c)
        # c / sqrt(val_d)
        IR_operation_node(IRtypes.assign_op, [val_c, IR_operation_node(IRtypes.math_op, [:sqrt, val_c])]),
        IR_operation_node(IRtypes.assign_op, [val_d, IR_operation_node(IRtypes.math_op, [:sqrt, val_d])]),
        IR_operation_node(IRtypes.assign_op, [val_d, IR_operation_node(IRtypes.math_op, [:/, val_c, val_d])])
    ])
    
    return block;
end

# For array a, b and scalar c, d do
# c = norm(a-b) and copy a into b
# d = c / norm(a)
function generate_difference_norms_and_update(vec_a, vec_b, val_c, val_d)
    IRtypes = IR_entry_types();
    # tmp = (vec_a[norm_i] - vec_b[norm_i])
    vec_a_i = apply_indexed_access(vec_a, [:norm_i], IRtypes);
    vec_b_i = apply_indexed_access(vec_b, [:norm_i], IRtypes);
    dif = IR_operation_node(IRtypes.math_op, [:-, vec_a_i, vec_b_i]);
    
    block = IR_block_node([
        # zero c
        IR_operation_node(IRtypes.assign_op, [val_c, 0]),
        IR_operation_node(IRtypes.assign_op, [val_d, 0]),
        # loop over vector
        IR_loop_node(IRtypes.index_loop, :dofs, :norm_i, 1, :dofs_global, IR_block_node([
            # add to norm val_c = val_c + (vec_a[norm_i] - vec_b[norm_i])
            IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_64_data, :tmp), dif]),
            IR_operation_node(IRtypes.assign_op, [val_c, 
                IR_operation_node(IRtypes.math_op, [:+, val_c, 
                    IR_operation_node(IRtypes.math_op, [:*, :tmp, :tmp])
                ])
            ]),
            IR_operation_node(IRtypes.assign_op, [val_d, 
                IR_operation_node(IRtypes.math_op, [:+, val_d, 
                    IR_operation_node(IRtypes.math_op, [:*, vec_a_i, vec_a_i])
                ])
            ]),
            # update b
            IR_operation_node(IRtypes.assign_op, [vec_b_i, vec_a_i])
        ])),
        # sqrt(val_c)
        IR_operation_node(IRtypes.assign_op, [val_c, IR_operation_node(IRtypes.math_op, [:sqrt, val_c])]),
        IR_operation_node(IRtypes.assign_op, [val_d, IR_operation_node(IRtypes.math_op, [:sqrt, val_d])]),
        IR_operation_node(IRtypes.assign_op, [val_d, IR_operation_node(IRtypes.math_op, [:/, val_c, val_d])])
    ])
    
    return block;
end

# creates something like (a-b)*c + d
function generate_the_one_pattern(a, b, c, d)
    IRtypes = IR_entry_types();
    result = IR_operation_node(IRtypes.math_op, [:+, 
                d,
                IR_operation_node(IRtypes.math_op, [:*,
                    c, 
                    IR_operation_node(IRtypes.math_op, [:-, a, b])
                    ])
                ]);
    
    return result;
end