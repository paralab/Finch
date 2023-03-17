#=
Functions for building an IR from the symbolic expressions for FVM type problems.

allocate
first loop
    matrix (if implicit)
build mat (if implicit)

time stepper
    loops
        eval parts needed by rhs only
        lhs matrix if t-dependent ??TODO
        vector
    solve+place (if implicit)

=#
function build_IR_fvm(input_exprs, var, indices, config, prob, time_stepper, fv_info)
    lhs_vol = input_exprs[1];
    rhs_vol = input_exprs[2];
    lhs_surf = input_exprs[3];
    rhs_surf = input_exprs[4];
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
    
    # Options
    average_coefficients = false; # true to integrate and average coefficients over cells and faces.
    
    # A matrix only needs to be created if the stepper is implicit
    need_matrix = time_stepper.implicit;
    
    IRtypes = IR_entry_types();
    
    # These will hold the IR
    allocate_block = IR_block_node([], "allocation");
    vol_coef_block = IR_block_node([], "volume coefficients");
    surf_coef_block = IR_block_node([], "surface coefficients");
    
    source_block = IR_block_node([], "surface integral");
    flux_block = IR_block_node([], "volume integral");
    toglobal_block = IR_block_node([], "local to global");
    
    # Allocate the global matrix and vector
    if need_matrix
        push!(allocate_block.parts, IR_comment_node("Allocate global matrix(IJV form)"));
        # allocated size of nonzeros is dofs*dofs*(nel + nfaces*4)
        allocatedNZ = IR_data_node(IRtypes.int_data, :allocated_nonzeros, [], []);
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            allocatedNZ,
            IR_operation_node(IRtypes.math_op, [:*, :dofs_per_node, :dofs_per_node,
                IR_operation_node(IRtypes.math_op, [:+, :num_elements, IR_operation_node(IRtypes.math_op, [:*, :num_faces,4])])
            ])
        ]));
        globalmat_I = IR_data_node(IRtypes.int_data, :global_matrix_I, [:allocated_nonzeros], []);
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            globalmat_I,
            IR_operation_node(IRtypes.allocate_op, [IRtypes.int_data, :allocated_nonzeros])
            ]));
        globalmat_J = IR_data_node(IRtypes.int_data, :global_matrix_J, [:allocated_nonzeros], []);
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            globalmat_J,
            IR_operation_node(IRtypes.allocate_op, [IRtypes.int_data, :allocated_nonzeros])
            ]));
        globalmat_V = IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocated_nonzeros], []);
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            globalmat_V,
            IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :allocated_nonzeros])
            ]));
        
        # # I and J indices only need to be set once, so init them here
        # push!(allocate_block.parts, IR_comment_node("Global matrix I,J values are set once."));
        # push!(allocate_block.parts, IR_operation_node(IRtypes.named_op, [:INIT_MATRIX_IJ_FV, :global_matrix_I, :global_matrix_J]));
        
        # next_nonzero_index = 1;
        push!(allocate_block.parts, IR_comment_node("I,J,V vectors for the sparse matrix are filled as needed."));
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.int_data, :next_nonzero_index), 1]));
        
        # Elemental matrices
        push!(allocate_block.parts, IR_comment_node("Allocate elemental matrix"));
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.float_data, :source_mat, [:dofs_per_loop, :dofs_per_node], []),
            IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :dofs_per_loop, :dofs_per_node])
        ]));
        flux_mat_cols = IR_operation_node(IRtypes.math_op, [:*, :dofs_per_node, :faces_per_element, 2]);
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, flux_mat_cols], []),
            IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :dofs_per_loop, flux_mat_cols])
        ]));
        flux_mat_tmp_cols = IR_operation_node(IRtypes.math_op, [:*, :dofs_per_node, 2]);
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.float_data, :flux_mat_tmp, [:dofs_per_loop, flux_mat_tmp_cols], []),
            IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :dofs_per_loop, flux_mat_tmp_cols])
        ]));
    end
    
    # Global vectors
    push!(allocate_block.parts, IR_comment_node("Allocate global vectors."));
    globalvec = IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        globalvec,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :fv_dofs_partition])
        ]));
    gsolvec = IR_data_node(IRtypes.float_data, :global_solution, [:fv_dofs_global], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        gsolvec,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :fv_dofs_global])
    ]));
    solvec = IR_data_node(IRtypes.float_data, :solution, [:fv_dofs_partition], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        solvec,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :fv_dofs_partition])
        ]));
    
    # elemental flux and source
    push!(allocate_block.parts, IR_comment_node("Allocate elemental source and flux."));
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        IR_data_node(IRtypes.float_data, :source, [:dofs_per_loop], []),
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :dofs_per_loop])
        ]));
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        IR_data_node(IRtypes.float_data, :flux, [:dofs_per_loop], []),
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :dofs_per_loop])
        ]));
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        IR_data_node(IRtypes.float_data, :flux_tmp, [:dofs_per_loop], []),
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :dofs_per_loop])
        ]));
    
    
    # bdry done flag for each face
    push!(allocate_block.parts, IR_comment_node("Boundary done flag for each face."));
    bdry_done = IR_data_node(IRtypes.int_data, :bdry_done, [:num_faces], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        bdry_done,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.int_data, :num_faces])
    ]));
    
    # Face flux done flag
    push!(allocate_block.parts, IR_comment_node("Flux done flag for each face so that it is not done twice."));
    face_flux_done = IR_data_node(IRtypes.boolean_data, :face_flux_done, [:num_faces], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        face_flux_done,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.boolean_data, :num_faces])
    ]));
    
    # index_values
    if length(indices) > 0
        max_index_tag = 0;
        for i=1:length(indices)
            if typeof(indices[i]) == Indexer
                max_index_tag = max(max_index_tag, indices[i].tag);
            end
        end
        push!(allocate_block.parts, IR_comment_node("index values to be passed to BCs if needed"));
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.int_data, :index_values, [max_index_tag], []),
            IR_operation_node(IRtypes.allocate_op, [IRtypes.int_data, max_index_tag])
        ]));
        
    else
        push!(allocate_block.parts, IR_comment_node("No indexed variables"));
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.int_data, :index_values, [0], []),
            IR_operation_node(IRtypes.allocate_op, [IRtypes.int_data, 0])
        ]));
    end
    
    #########################################################
    # If order > 1, parent-child maps are used
    #########################################################
    if fv_info.fluxOrder > 1 && !(fv_info.parentMaps === nothing)
        # # patch_size = fv_info.parentMaps.patch_size;
        # # patch_cells = zeros(Int, patch_size)
        # push!(allocate_block.parts, 
        #         IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.int_data, :patch_size),
        #             IR_operation_node(IRtypes.member_op, [:fv_info, 
        #                 IR_operation_node(IRtypes.member_op, [:parentMaps, 
        #                     IR_data_node(IRtypes.int_data, :patch_size)
        #         ])])]));
        # push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        #     IR_data_node(IRtypes.int_data, :patch_cells, [:patch_size], []),
        #     IR_operation_node(IRtypes.allocate_op, [IRtypes.int_data, :patch_size])
        # ]));
        
        # num_parents = fv_info.parentMaps.num_parents;
        # num_children = fv_info.parentMaps.children_per_parent;
        # num_child_faces = fv_info.parentMaps.faces_per_parent;
        parent_part = IR_block_node([
            IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.int_data, :num_parents),
                IR_operation_node(IRtypes.member_op, [:fv_info, 
                    IR_operation_node(IRtypes.member_op, [:parentMaps, 
                        IR_data_node(IRtypes.int_data, :num_parents)
            ])])]),
            IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.int_data, :num_children),
                IR_operation_node(IRtypes.member_op, [:fv_info, 
                    IR_operation_node(IRtypes.member_op, [:parentMaps, 
                        IR_data_node(IRtypes.int_data, :children_per_parent)
            ])])]),
            IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.int_data, :num_child_faces),
                IR_operation_node(IRtypes.member_op, [:fv_info, 
                    IR_operation_node(IRtypes.member_op, [:parentMaps, 
                        IR_data_node(IRtypes.int_data, :faces_per_parent)
            ])])])
        ])
        need_parent = true;
        
    else # no parent-child
        parent_part = IR_comment_node("No parent-child mesh needed");
        need_parent = false;
    end
    
    #########################################################
    # coefficient prep
    #########################################################
    # a list of all entities
    all_entities = [];
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
        counts[4] = length(entities);
    end
    
    #########################################################
    # Prepare coefficient and other elemental values
    #########################################################
    (vol_coef, surf_coef) = prepare_coefficient_values_fvm(all_entities, var, dimension, counts, fv_info, need_parent);
    
    # volume coefficients (for source)
    push!(vol_coef_block.parts, IR_comment_node("Evaluate volume coefficients."));
    append!(vol_coef_block.parts, vol_coef.parts);
    push!(vol_coef_block.parts, IR_operation_node(IRtypes.assign_op, [:volume, 
        IR_operation_node(IRtypes.member_op, [:geometric_factors, 
            IR_data_node(IRtypes.float_data, :volume, [:num_elements], [:eid])])]))
    
    # surface coefficients (for flux)
    push!(surf_coef_block.parts, IR_comment_node("Evaluate surface coefficients."));
    append!(surf_coef_block.parts, surf_coef.parts);
    push!(surf_coef_block.parts, IR_operation_node(IRtypes.assign_op, [:area, 
        IR_operation_node(IRtypes.member_op, [:geometric_factors, 
            IR_data_node(IRtypes.float_data, :area, [:num_faces], [:fid])])]))
    push!(surf_coef_block.parts, IR_operation_node(IRtypes.assign_op, [:area_over_volume, 
        IR_operation_node(IRtypes.math_op, [:(/), :area, :volume])]))
    
    #########################################################
    # source block
    #########################################################
    push!(source_block.parts, vol_coef_block);
    push!(source_block.parts, IR_comment_node("Compute source terms (volume integral)"));
    # first the vector(rhs or explicit)
    if !(rhs_vol === nothing) && !is_empty_expression(rhs_vol)
        rhsvol_terms = process_terms(rhs_vol);
        push!(source_block.parts, make_elemental_computation_fvm(rhsvol_terms, var, dofsper, offset_ind, RHS, "volume"));
    else
        
    end
    # The matrix is only for implicit
    if need_matrix && !(lhs_vol === nothing) && !is_empty_expression(lhs_vol)
        lhsvol_terms = process_terms(lhs_vol);
        push!(source_block.parts, make_elemental_computation_fvm(lhsvol_terms, var, dofsper, offset_ind, LHS, "volume"));
    end
    
    #########################################################
    # flux block (includes face loop)
    #########################################################
    push!(flux_block.parts, IR_comment_node("Compute flux terms (surface integral) in a face loop"));
    fid = :fid # fid = IR_data_node(IRtypes.int_data, :fid, [], []); # Face ID
    fbid = IR_data_node(IRtypes.int_data, :fbid, [], []); # Face BID
    leftel = IR_data_node(IRtypes.int_data, :left_el, [], []); # left element
    rightel = IR_data_node(IRtypes.int_data, :right_el, [], []); # right element
    in_side = IR_data_node(IRtypes.int_data, :in_side, [], []);
    out_side = IR_data_node(IRtypes.int_data, :out_side, [], []);
    # zero flux_tmp inside the face loop
    if dofsper == 1
        zero_flux_t = IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_data, :flux_tmp, [:dofs_per_loop], [1]), 0.0]);
    else
        zero_flux_t = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY, IR_data_node(IRtypes.float_data, :flux_tmp, [:dofs_per_loop], []), 0.0, :dofs_per_loop]);
    end
    if need_matrix
        zero_flux_mat_t = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY, :flux_mat_tmp, 0.0, :dofs_per_loop, flux_mat_tmp_cols]);
    else
        zero_flux_mat_t = IR_comment_node("");
    end
    # The face loop
    face_loop_body = IR_block_node([
        # zero flux_tmp arrays
        zero_flux_t,
        zero_flux_mat_t,
        # fid = fv_grid.element2face[i, eid]; # face ID 
        IR_operation_node(IRtypes.assign_op, [:fid,
            IR_operation_node(IRtypes.member_op, [:mesh, 
                IR_data_node(IRtypes.int_data, :element2face, [:faces_per_element, :num_elements], [:fi, :eid])])
        ]),
        # # Only compute flux once for each face
        # IR_conditional_node(IR_data_node(IRtypes.boolean_data, :face_flux_done, [:num_faces], [:fid]),
        #     IR_block_node([IR_operation_node(IRtypes.named_op, [:CONTINUE])])
        # ),
        # fbid = fv_grid.facebid[fid]; # BID of this face
        IR_operation_node(IRtypes.assign_op, [fbid,
            IR_operation_node(IRtypes.member_op, [:mesh, 
                IR_data_node(IRtypes.int_data, :facebid, [:num_faces], [:fid])])
        ]),
        # (leftel, rightel) = fv_grid.face2element[:,fid];
        IR_operation_node(IRtypes.assign_op, [leftel,
            IR_operation_node(IRtypes.member_op, [:mesh, 
                IR_data_node(IRtypes.int_data, :face2element, [2, :num_faces], [1, :fid])])
        ]),
        IR_operation_node(IRtypes.assign_op, [rightel,
            IR_operation_node(IRtypes.member_op, [:mesh, 
                IR_data_node(IRtypes.int_data, :face2element, [2, :num_faces], [2, :fid])])
        ]),
        # if (eid == rightel || rightel == 0) (neighbor = leftel) else (neighbor = rightel)
        IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :eid, rightel]),
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [:neighbor, leftel]),
                IR_operation_node(IRtypes.assign_op, [in_side, 2]),
                IR_operation_node(IRtypes.assign_op, [out_side, 1])
            ]),
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [in_side, 1]),
                IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), rightel, 0]),
                    IR_block_node([
                        IR_operation_node(IRtypes.assign_op, [:neighbor, leftel]),
                        IR_operation_node(IRtypes.assign_op, [out_side, 1])
                    ]),
                    IR_block_node([
                        IR_operation_node(IRtypes.assign_op, [:neighbor, rightel]),
                        IR_operation_node(IRtypes.assign_op, [out_side, 2])
                    ])
                )
            ])
        ),
        
    ]);
    
    push!(face_loop_body.parts, surf_coef_block);
    push!(face_loop_body.parts, IR_comment_node("Compute flux terms (surface integral)"));
    # first the vector(rhs or explicit)
    if !(rhs_surf === nothing) && !is_empty_expression(rhs_surf)
        rhssurf_terms = process_terms(rhs_surf);
        push!(face_loop_body.parts, make_elemental_computation_fvm(rhssurf_terms, var, dofsper, offset_ind, RHS, "surface"));
    end
    
    # The matrix is only for implicit
    if need_matrix && !(lhs_surf === nothing)  && !is_empty_expression(lhs_surf)
        lhssurf_terms = process_terms(lhs_surf);
        push!(face_loop_body.parts, make_elemental_computation_fvm(lhssurf_terms, var, dofsper, offset_ind, LHS, "surface"));
    end
    
    boundary_condition_block = IR_block_node([], "boundary conditions");
    push!(boundary_condition_block.parts, IR_comment_node("Apply boundary conditions"));
    if need_matrix && !(lhs_surf === nothing)  && !is_empty_expression(lhs_surf)
        push!(boundary_condition_block.parts, IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>), :fbid, 0]),
            IR_block_node([IR_operation_node(IRtypes.function_op, [
                :apply_boundary_conditions_face, 
                :var, :eid, :fid, :fbid, :mesh, :refel, :geometric_factors, :fv_info, :prob, :t, :dt, :flux_mat_tmp, :flux_tmp, :bdry_done, :index_offset, :index_values])
            ])));
    else
        push!(boundary_condition_block.parts, IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>), :fbid, 0]),
            IR_block_node([IR_operation_node(IRtypes.function_op, [
                :apply_boundary_conditions_face_rhs, 
                :var, :eid, :fid, :fbid, :mesh, :refel, :geometric_factors, :fv_info, :prob, :t, :dt, :flux_tmp, :bdry_done, :index_offset, :index_values])
            ])));
    end
    push!(face_loop_body.parts, boundary_condition_block);
    
    # Add flux_tmp to flux for this element
    if dofsper_loop == 1
        # add RHS flux_tmp to flux after handling BCs
        push!(face_loop_body.parts, IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_data, :flux, [:dofs_per_loop], [1]), 
            IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.float_data, :flux, [:dofs_per_loop], [1]), 
                IR_operation_node(IRtypes.math_op, [:*, IR_data_node(IRtypes.float_data, :flux_tmp, [:dofs_per_loop], [1]), :area_over_volume])])]));
                
        # add LHS flux_tmp to flux if needed
        if need_matrix && !(lhs_surf === nothing)  && !is_empty_expression(lhs_surf)
            push!(face_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, flux_mat_cols], [1, :fi]), 
                IR_operation_node(IRtypes.math_op, [:+, 
                    IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, flux_mat_cols], [1, :fi]), 
                    IR_operation_node(IRtypes.math_op, [:*, 
                        IR_data_node(IRtypes.float_data, :flux_mat_tmp, [:dofs_per_loop, flux_mat_tmp_cols], [1, 1]), 
                        :area_over_volume])])]));
            flux_mat_index = IR_operation_node(IRtypes.math_op, [:+, :faces_per_element, :fi]); 
            push!(face_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, flux_mat_cols], [1, flux_mat_index]), 
                IR_operation_node(IRtypes.math_op, [:+, 
                    IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, flux_mat_cols], [1, flux_mat_index]), 
                    IR_operation_node(IRtypes.math_op, [:*, 
                        IR_data_node(IRtypes.float_data, :flux_mat_tmp, [:dofs_per_loop, flux_mat_tmp_cols], [1, 2]), 
                        :area_over_volume])])]));
        end
        
    else # dofs > 1
        push!(face_loop_body.parts, IR_loop_node(IRtypes.dof_loop, :dofs, :dofi, 1, :dofs_per_loop, IR_block_node([
            IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_data, :flux, [:dofs_per_loop], [:dofi]), 
            IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.float_data, :flux, [:dofs_per_loop], [:dofi]), 
                IR_operation_node(IRtypes.math_op, [:*, IR_data_node(IRtypes.float_data, :flux_tmp, [:dofs_per_loop], [:dofi]), :area_over_volume])])])
        ])));
        # add LHS flux_tmp to flux if needed
        if need_matrix && !(lhs_surf === nothing)  && !is_empty_expression(lhs_surf)
            flux_mat_index1 = IR_operation_node(IRtypes.math_op, [:+, :dofj,
                IR_operation_node(IRtypes.math_op, [:*, :dofs_per_node, 
                    IR_operation_node(IRtypes.math_op, [:-, :fi, 1])])
            ]);
            flux_mat_index2 = IR_operation_node(IRtypes.math_op, [:+, :dofj,
                IR_operation_node(IRtypes.math_op, [:*, :dofs_per_node, 
                    IR_operation_node(IRtypes.math_op, [:-, :fi, 1])]),
                IR_operation_node(IRtypes.math_op, [:*, :faces_per_element, :dofs_per_node])
            ]);
            dofjplus = IR_operation_node(IRtypes.math_op, [:+, :dofj, :dofs_per_node]);
            push!(face_loop_body.parts, IR_loop_node(IRtypes.dof_loop, :dofs, :dofi, 1, :dofs_per_loop, IR_block_node([
                IR_loop_node(IRtypes.dof_loop, :dofs, :dofj, 1, :dofs_per_node, IR_block_node([
                    IR_operation_node(IRtypes.assign_op, [
                        IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, flux_mat_cols], [:dofi, flux_mat_index1]), 
                        IR_operation_node(IRtypes.math_op, [:+, 
                            IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, flux_mat_cols], [:dofi, flux_mat_index1]), 
                            IR_operation_node(IRtypes.math_op, [:*, 
                                IR_data_node(IRtypes.float_data, :flux_mat_tmp, [:dofs_per_loop, flux_mat_tmp_cols], [:dofi, :dofj]), 
                                :area_over_volume])])]),
                        IR_operation_node(IRtypes.assign_op, [
                            IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, flux_mat_cols], [:dofi, flux_mat_index2]), 
                            IR_operation_node(IRtypes.math_op, [:+, 
                                IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, flux_mat_cols], [:dofi, flux_mat_index2]), 
                                IR_operation_node(IRtypes.math_op, [:*, 
                                    IR_data_node(IRtypes.float_data, :flux_mat_tmp, [:dofs_per_loop, flux_mat_tmp_cols], [:dofi, dofjplus]), 
                                    :area_over_volume])])])
                ]))
            ])));
        end
    end
    
    # Place the body in the face loop
    push!(flux_block.parts, IR_loop_node(IRtypes.space_loop, :faces, :fi, 1, :faces_per_element, face_loop_body));
    
    #########################################################
    # add to global sysem
    #########################################################
    push!(toglobal_block.parts, IR_comment_node("Place elemental parts in global system."));
    push!(toglobal_block.parts, generate_local_to_global_fvm(dofsper, offset_ind, vec_only=!need_matrix));
    
    #########################################################
    # assembly loop
    #########################################################
    # The elemental loop will be different for parent-child maps
    assembly_loop = generate_assembly_loop_fvm(var, indices, need_parent);
    
    # find the innermost assembly loop
    inner_loop = assembly_loop;
    while length(inner_loop.body.parts) > 0 && typeof(inner_loop.body.parts[1]) == IR_loop_node
        inner_loop = inner_loop.body.parts[1];
    end
    
    # zero elemental source and flux
    if dofsper == 1
        zero_source = IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_data, :source, [:dofs_per_loop], [1]), 0.0]);
        zero_flux = IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_data, :flux, [:dofs_per_loop], [1]), 0.0]);
    else
        zero_source = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY, IR_data_node(IRtypes.float_data, :source, [:dofs_per_loop], []), 0.0, :dofs_per_loop]);
        zero_flux = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY, IR_data_node(IRtypes.float_data, :flux, [:dofs_per_loop], []), 0.0, :dofs_per_loop]);
    end
    
    # fill the elemental loop
    if need_matrix
        zero_source_mat = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY, IR_data_node(IRtypes.float_data, :source_mat, [:dofs_per_loop, :dofs_per_node], []), 0.0, :dofs_per_loop, :dofs_per_node]);
        zero_flux_mat = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY, IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :flux_mat_cols], []), 0.0, :dofs_per_loop, flux_mat_cols]);
        
        append!(inner_loop.body.parts, [
            zero_source,
            zero_flux,
            zero_source_mat,
            zero_flux_mat,
            source_block,
            flux_block,
            toglobal_block
        ])
        
    else
        append!(inner_loop.body.parts, [
            zero_source,
            zero_flux,
            source_block,
            flux_block,
            toglobal_block
        ])
    end
    # Wrap the assembly loop in a block to label it
    assembly_loop = IR_block_node([assembly_loop], "assembly loop");
    
    #########################################################
    # time loop
    #########################################################
    need_time_loop = prob.time_dependent; # This should be true for FV
    # Exchange ghosts if needed
    ghost_exchange = IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>), :num_partitions, 1]),
                        IR_block_node([IR_operation_node(IRtypes.named_op, [:GHOST_EXCHANGE_FV, 0])]));
    if need_time_loop
        stepper_dt = IR_operation_node(IRtypes.member_op, [:time_stepper, :dt]);
        if need_matrix # implicit steppers
            step_loop = generate_time_stepping_loop_fvm(time_stepper, assembly_loop, prob);
            compute_block = IR_block_node([
                parent_part,
                ghost_exchange,
                IR_operation_node(IRtypes.named_op, [:GATHER_VARS, solvec]),
                IR_operation_node(IRtypes.assign_op, [:t, 0.0]),
                IR_operation_node(IRtypes.assign_op, [:dt, stepper_dt]),
                IR_comment_node("###############################################"),
                IR_comment_node("Time stepping loop"),
                wrap_in_timer(:time_steps, step_loop)
            ]);
            
        else # explicit steppers (no matrix)
            step_loop = generate_time_stepping_loop_fvm(time_stepper, assembly_loop, prob);
            compute_block = IR_block_node([
                parent_part,
                ghost_exchange,
                IR_operation_node(IRtypes.named_op, [:GATHER_VARS, solvec]),
                IR_operation_node(IRtypes.assign_op, [:t, 0.0]),
                IR_operation_node(IRtypes.assign_op, [:dt, stepper_dt]),
                IR_comment_node("###############################################"),
                IR_comment_node("Time stepping loop"),
                wrap_in_timer(:time_steps, step_loop)
            ]);
        end
        
    else # no time stepping
        # Do we need this?
        compute_block = IR_block_node([IR_comment_node("FVM only supported for time dependent problems.")])
    end
    
    #########################################################
    # Put them all together in a master block
    #########################################################
    master_block = IR_block_node([
        wrap_in_timer(:allocate, allocate_block),
        compute_block
    ]);
    
    return master_block;
end

# Compute or fetch all needed values
function prepare_coefficient_values_fvm(entities, var, dimension, counts, fv_info, need_parent)
    IRtypes = IR_entry_types();
    row_col_matrix_index = IR_operation_node(IRtypes.named_op, [:ROWCOL_TO_INDEX, :row, :col, :nodes_per_element]);
    col_row_matrix_index = IR_operation_node(IRtypes.named_op, [:ROWCOL_TO_INDEX, :col, :row, :nodes_per_element]);
    
    FACENORMAL1_1 = IR_data_node(IRtypes.float_data, :FACENORMAL1_1, [], []);
    FACENORMAL1_2 = IR_data_node(IRtypes.float_data, :FACENORMAL1_2, [], []);
    FACENORMAL1_3 = IR_data_node(IRtypes.float_data, :FACENORMAL1_3, [], []);
    FACENORMAL2_1 = IR_data_node(IRtypes.float_data, :FACENORMAL2_1, [], []);
    FACENORMAL2_2 = IR_data_node(IRtypes.float_data, :FACENORMAL2_2, [], []);
    FACENORMAL2_3 = IR_data_node(IRtypes.float_data, :FACENORMAL2_3, [], []);
    fid = IR_data_node(IRtypes.int_data, :fid, [], []); # Face ID
    leftel = IR_data_node(IRtypes.int_data, :left_el, [], []); # left element
    rightel = IR_data_node(IRtypes.int_data, :right_el, [], []); # right element
    neighbor = IR_data_node(IRtypes.int_data, :neighbor, [], []);
    
    # These parts will be returned
    vol_block = IR_block_node([]);
    surf_block = IR_block_node([]);
    
    unique_entity_names = []; # avoid duplicate names
    
    # First get x,y,z for cell centers
    cell_coords = IR_block_node([]);
    need_cell_coords = false;
    push!(cell_coords.parts, IR_operation_node(IRtypes.assign_op, [:x, 
        IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :cellCenters, [dimension, :num_elements], [1,:eid])])]));
    if dimension > 1
        push!(cell_coords.parts, IR_operation_node(IRtypes.assign_op, [:y, 
            IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :cellCenters, [dimension, :num_elements], [2,:eid])])]));
    else
        push!(cell_coords.parts, IR_operation_node(IRtypes.assign_op, [:y, 0.0]));
    end
    if dimension > 2
        push!(cell_coords.parts, IR_operation_node(IRtypes.assign_op, [:z, 
            IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :cellCenters, [dimension, :num_elements], [3,:eid])])]));
    else
        push!(cell_coords.parts, IR_operation_node(IRtypes.assign_op, [:z, 0.0]));
    end
    
    # And the same for face centers
    face_coords = IR_block_node([]);
    need_face_coords = false;
    push!(face_coords.parts, IR_operation_node(IRtypes.assign_op, [:x, 
        IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :faceCenters, [dimension, :num_faces], [1,:fid])])]));
    if dimension > 1
        push!(face_coords.parts, IR_operation_node(IRtypes.assign_op, [:y, 
            IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :faceCenters, [dimension, :num_faces], [2,:fid])])]));
    else
        push!(face_coords.parts, IR_operation_node(IRtypes.assign_op, [:y, 0.0]));
    end
    if dimension > 2
        push!(face_coords.parts, IR_operation_node(IRtypes.assign_op, [:z, 
            IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :faceCenters, [dimension, :num_faces], [3,:fid])])]));
    else
        push!(face_coords.parts, IR_operation_node(IRtypes.assign_op, [:z, 0.0]));
    end
    
    # Are face normals needed?
    need_normals = false;
    normal_part = IR_block_node([]);
    # FACENORMAL1[1] -> FACENORMAL1_1 = mesh.facenormals[1,fid];
    push!(normal_part.parts, IR_operation_node(IRtypes.assign_op, [FACENORMAL1_1, 
        IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.float_data, :facenormals, [dimension, :num_faces], [1,:fid])])]));
    if dimension > 1
        push!(normal_part.parts, IR_operation_node(IRtypes.assign_op, [FACENORMAL1_2, 
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.float_data, :facenormals, [dimension, :num_faces], [2,:fid])])]));
    else
        # push!(normal_part.parts, IR_operation_node(IRtypes.assign_op, [FACENORMAL1_2, 0.0]));
    end
    if dimension > 2
        push!(normal_part.parts, IR_operation_node(IRtypes.assign_op, [FACENORMAL1_3, 
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.float_data, :facenormals, [dimension, :num_faces], [3,:fid])])]));
    else
        # push!(normal_part.parts, IR_operation_node(IRtypes.assign_op, [FACENORMAL1_3, 0.0]));
    end
    
    if dimension == 1
        push!(normal_part.parts, IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :eid, leftel]),
            # normal is correct, make FACENORMAL2 = -FACENORMAL1
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_1, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_1])])]),
            # normal is reversed
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_1, FACENORMAL1_1]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL1_1, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_1])])])
        ))
    elseif dimension == 2
        push!(normal_part.parts, IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :eid, leftel]),
            # normal is correct, make FACENORMAL2 = -FACENORMAL1
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_1, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_1])]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_2, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_2])])]),
            # normal is reversed
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_1, FACENORMAL1_1]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_2, FACENORMAL1_2]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL1_1, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_1])]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL1_2, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_2])])])
        ))
    elseif dimension == 3
        push!(normal_part.parts, IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :eid, leftel]),
            # normal is correct, make FACENORMAL2 = -FACENORMAL1
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_1, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_1])]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_2, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_2])]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_3, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_3])])]),
            # normal is reversed
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_1, FACENORMAL1_1]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_2, FACENORMAL1_2]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL2_3, FACENORMAL1_3]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL1_1, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_1])]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL1_2, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_2])]),
                IR_operation_node(IRtypes.assign_op, [FACENORMAL1_3, IR_operation_node(IRtypes.math_op, [:-, FACENORMAL1_3])])])
        ))
    end
    
    # Is the distance between cell centers needed for derivatives?
    # normal scaled by distance between cell centers"
    # cell_dist = sqrt( (xn-xe)^2 + ()^2 + ()^2 )
    # cell_dx = cell_dist * FACENORMAL1_1
    need_deriv_dist = false;
    # get cell_xe, cell_xn, ...
    cell_coords_2 = IR_block_node([
            IR_operation_node(IRtypes.assign_op, [:cell_xe, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :cellCenters, [dimension, :num_elements], [1,:eid])])
            ]),
            IR_operation_node(IRtypes.assign_op, [:cell_xn, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :cellCenters, [dimension, :num_elements], [1,:neighbor])])
            ])
        ]);
    if dimension > 1
        append!(cell_coords_2.parts, [
            IR_operation_node(IRtypes.assign_op, [:cell_ye, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :cellCenters, [dimension, :num_elements], [2,:eid])])
            ]),
            IR_operation_node(IRtypes.assign_op, [:cell_yn, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :cellCenters, [dimension, :num_elements], [2,:neighbor])])
            ])
        ])
    end
    if dimension > 2
        append!(cell_coords_2.parts, [
            IR_operation_node(IRtypes.assign_op, [:cell_ze, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :cellCenters, [dimension, :num_elements], [3,:eid])])
            ]),
            IR_operation_node(IRtypes.assign_op, [:cell_zn, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_data, :cellCenters, [dimension, :num_elements], [3,:neighbor])])
            ])
        ])
    end
    # find the square distance (xn-xe)^2 + ...
    square_dist = IR_operation_node(IRtypes.math_op, [:*,
        IR_operation_node(IRtypes.math_op, [:-,
            IR_data_node(IRtypes.float_data, :cell_xn, [], []),
            IR_data_node(IRtypes.float_data, :cell_xe, [], [])
        ]),
        IR_operation_node(IRtypes.math_op, [:-,
            IR_data_node(IRtypes.float_data, :cell_xn, [], []),
            IR_data_node(IRtypes.float_data, :cell_xe, [], [])
        ])
    ])
    if dimension > 1
        square_dist = IR_operation_node(IRtypes.math_op, [:+, square_dist, 
            IR_operation_node(IRtypes.math_op, [:*,
                IR_operation_node(IRtypes.math_op, [:-,
                    IR_data_node(IRtypes.float_data, :cell_yn, [], []),
                    IR_data_node(IRtypes.float_data, :cell_ye, [], [])
                ]),
                IR_operation_node(IRtypes.math_op, [:-,
                    IR_data_node(IRtypes.float_data, :cell_yn, [], []),
                    IR_data_node(IRtypes.float_data, :cell_ye, [], [])
                ])
            ])
        ])
    end
    if dimension > 2
        push!(square_dist.args, 
            IR_operation_node(IRtypes.math_op, [:*,
                IR_operation_node(IRtypes.math_op, [:-,
                    IR_data_node(IRtypes.float_data, :cell_zn, [], []),
                    IR_data_node(IRtypes.float_data, :cell_ze, [], [])
                ]),
                IR_operation_node(IRtypes.math_op, [:-,
                    IR_data_node(IRtypes.float_data, :cell_zn, [], []),
                    IR_data_node(IRtypes.float_data, :cell_ze, [], [])
                ])
            ])
        )
    end
    # Put it all together
    deriv_dist = IR_block_node([IR_comment_node("Normal scaled by distance between cell centers"),
        cell_coords_2,
        IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.float_data, :cell_dist, [], []),
            IR_operation_node(IRtypes.math_op, [:sqrt, square_dist])
        ]),
        IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.float_data, :dxyz_1, [], []),
            IR_operation_node(IRtypes.math_op, [:*, :cell_dist, FACENORMAL1_1])
        ])
    ]);
    if dimension > 1
        push!(deriv_dist.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :dxyz_2, [], []),
                IR_operation_node(IRtypes.math_op, [:*, :cell_dist, FACENORMAL1_2])
            ]))
    end
    if dimension > 2
        push!(deriv_dist.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :dxyz_3, [], []),
                IR_operation_node(IRtypes.math_op, [:*, :cell_dist, FACENORMAL1_3])
            ]))
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
         
        if is_unknown_var(entities[i], var) && lorr == LHS
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
                if entities[i].name == "FACENORMAL1"
                    need_normals = true;
                    # FACENORMAL1_n = normal[n]
                elseif entities[i].name == "FACENORMAL2"
                    need_normals = true;
                    # FACENORMAL2_n = -normal[n]
                end
                
            elseif ctype == 0
                # It was a number, do nothing?
                
            elseif ctype == 1 # a constant wrapped in a coefficient will be replaced by a number
                # cval holds the number
                if vors == "volume"
                    push!(vol_block.parts, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_data, Symbol(cname), [], []), cval]))
                else
                    push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_data, Symbol(cname), [], []), cval]))
                end
                
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
                        indices = finch_state.variables[cval].indexer;
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
                    need_cell_coords = true;
                    coef_index = get_coef_index(entities[i]);
                    
                    push!(vol_block.parts, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_data, Symbol(cname), [], []),
                        IR_operation_node(IRtypes.named_op, [:COEF_EVAL, coef_index, index_IR, :x, :y, :z, :t, :eid, 0, :index_values])]));
                    # If derivatives are needed, I'm not sure what to do here
                    
                else # surface
                    need_face_coords = true;
                    coef_index = get_coef_index(entities[i]);
                    
                    push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_data, Symbol(cname), [], []),
                        IR_operation_node(IRtypes.named_op, [:COEF_EVAL, coef_index, index_IR, :x, :y, :z, :t, :eid, :fid, :index_values])]));
                    # If derivatives are needed, I'm not sure what to do here
                end
                
            elseif ctype == 3 # a known variable value
                cellside = get_face_side_info(entities[i]); # in code_generator_utils.jl
                if cellside == 1
                    cellsidename = :in_side;
                elseif cellside == 2
                    cellsidename = :out_side;
                else
                    cellsidename = 3;
                end
                
                # Make an index string for indexed variables
                if typeof(entities[i].index) <: Array
                    # It is an indexed variable
                    if length(entities[i].index) == 1
                        indsymbol = Symbol("INDEX_VAL_"*entities[i].index[1]);
                    else
                        # There is more than one index. Need to form an expression for it.
                        indsymbol = "(INDEX_VAL_"*entities[i].index[1];
                        indices = finch_state.variables[cval].indexer;
                        for indi=2:length(entities[i].index)
                            indsymbol *= " + ("*string(length(indices[indi-1].range))*"*(INDEX_VAL_"*entities[i].index[indi]*"-1)";
                        end
                        for indi=1:length(entities[i].index)
                            indsymbol *= ")";
                        end
                    end
                    
                else # not indexed variable
                    indsymbol = entities[i].index;
                end
                
                if vors == "volume"
                    # It is possible that the variable is defined at nodes, not cells
                    if finch_state.variables[cval].discretization == FV
                        # finch_state.variables[cval].values[index,eid]
                        push!(vol_block.parts, IR_operation_node(IRtypes.assign_op,[
                            IR_data_node(IRtypes.float_data, Symbol(cname), [], []),
                            IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :eid])
                        ]));
                    else
                        # sum(finch_state.variables[cval].values[index,mesh.loc2glb[:,eid]]) / nodes_per_element
                        # TODO
                    end
                    
                else # surface
                    if finch_state.variables[cval].discretization == FV
                        # Need to reconstruct it using neighboring cells
                        if need_parent
                            # cname = FV_reconstruct_value(   )
                            reconstruct_args = [:FV_reconstruct_face_value,
                            IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval]),
                            indsymbol, :fid, cellsidename, :mesh, :fv_info];
                            append!(reconstruct_args, entities[i].derivs);
                            
                            push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                                IR_data_node(IRtypes.float_data, Symbol(cname)),
                                IR_operation_node(IRtypes.function_op, reconstruct_args)
                            ]));
                            
                        else # no parent
                            if length(entities[i].derivs) > 0
                                need_normals = true;
                                need_deriv_dist = true;
                                # cname = finch_state.variables[cval].values[index, els2] - finch_state.variables[cval].values[index, els1]
                                # if (els1 != els2 && abs(normal[entities[i].derivs[1]]) > 1e-10)
                                #     cname = cname / dxyz[entities[i].derivs[1]]
                                # else
                                #     cname = 0
                                # end
                                if entities[i].derivs[1] == 1
                                    normal_n = FACENORMAL1_1;
                                elseif entities[i].derivs[1] == 2
                                    normal_n = FACENORMAL1_2;
                                else
                                    normal_n = FACENORMAL1_3;
                                end
                                push!(surf_block.parts, IR_conditional_node(
                                    IR_operation_node(IRtypes.math_op, [:(&&), 
                                        IR_operation_node(IRtypes.math_op, [:(!=), :eid, :neighbor]),
                                        IR_operation_node(IRtypes.math_op, [:(>), 
                                            IR_operation_node(IRtypes.function_op, [:abs, normal_n]),
                                            1e-10])
                                    ]), 
                                    IR_block_node([
                                        IR_operation_node(IRtypes.assign_op,[
                                            IR_data_node(IRtypes.float_data, Symbol(cname), [], []),
                                            IR_operation_node(IRtypes.math_op, [:/,
                                                IR_operation_node(IRtypes.math_op, [:-,
                                                    IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :neighbor]),
                                                    IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :eid])
                                                ]),
                                                IR_data_node(IRtypes.float_data, Symbol("dxyz_"*string(entities[i].derivs[1])), [], [])
                                            ])
                                        ])
                                    ]), 
                                    IR_block_node([
                                        IR_operation_node(IRtypes.assign_op,[
                                            IR_data_node(IRtypes.float_data, Symbol(cname), [], []), 0.0])
                                    ])
                                ));
                                
                            else
                                if cellside == 0 || cellside == 3
                                    # No side was specified or central approx, so use the average
                                    # cname = 0.5 * (finch_state.variables[cval].values[index,eid] + finch_state.variables[cval].values[index,neighbor])
                                    push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                                        IR_data_node(IRtypes.float_data, Symbol(cname), [], []),
                                        IR_operation_node(IRtypes.math_op, [:*, 0.5,
                                            IR_operation_node(IRtypes.math_op, [:+,
                                                IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :neighbor]),
                                                IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :eid])
                                            ])
                                        ])
                                    ]))
                                elseif cellside == 1
                                    # cname = finch_state.variables[cval].values[index, whichone]
                                    push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                                        IR_data_node(IRtypes.float_data, Symbol(cname)),
                                        IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :eid])
                                    ]))
                                elseif cellside == 2
                                    # cname = finch_state.variables[cval].values[index, whichone]
                                    push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                                        IR_data_node(IRtypes.float_data, Symbol(cname)),
                                        IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :neighbor])
                                    ]))
                                # elseif cellside == 3 # central
                                    # same as 0 for this case
                                    # cname = 0.5 * (finch_state.variables[cval].values[index,eid] + finch_state.variables[cval].values[index,neighbor])
                                elseif cellside == 4 # neighborhood
                                    # This is a special case only for callback functions.
                                    # Rather than representing a value, it constructs a Neighborhood object to be passed.
                                    # code *= cname * " = Finch.Neighborhood(els, cellx, [Finch.finch_state.variables["*string(cval)*"].values["*indstr*", els[1]], Finch.finch_state.variables["*string(cval)*"].values["*indstr*", els[2]]]);\n";
                                    # TODO
                                end
                            end
                        end
                        
                    else # defined at nodes
                        # TODO
                    end
                end # surf/vol
            end # ctype
        end # if coefficient
    end # entity loop
    
    # include these pieces if needed
    if need_cell_coords
        append!(cell_coords.parts, vol_block.parts);
        vol_block.parts = cell_coords.parts;
    end
    if need_deriv_dist
        append!(deriv_dist.parts, surf_block.parts);
        surf_block.parts = deriv_dist.parts;
    end
    if need_face_coords
        append!(face_coords.parts, surf_block.parts);
        surf_block.parts = face_coords.parts;
    end
    if need_normals
        append!(normal_part.parts, surf_block.parts);
        surf_block.parts = normal_part.parts;
    end
    
    return (vol_block, surf_block);
end

function make_elemental_computation_fvm(terms, var, dofsper, offset_ind, lorr, vors)
    # Here is where I make some assumption about the form of the expression.
    # Since it was expanded by the parser it should look like a series of terms: t1 + t2 + t3...
    # Where each term, if LHS, involves one unknown component and possibly some coefficients.
    
    IRtypes = IR_entry_types();
    
    dofsper_loop = length(var[1].symvar);
    for i=2:length(var)
        dofsper_loop = dofsper_loop + length(var[i].symvar);
    end
    
    # This will be returned
    compute_block = IR_block_node([]);
    
    # source or flux being computed?
    if vors == "surface"
        if lorr == LHS
            if dofsper > 1
                result_size = [:dofs_per_loop, :?];
                
            else # one dof
                result_size = [:dofs_per_loop, 2];
            end
            result_part = IR_data_node(IRtypes.float_data, :flux_mat_tmp, result_size, [])
            
        else # RHS
            result_part = IR_data_node(IRtypes.float_data, :flux_tmp, [:dofs_per_loop], [])
        end
        
    else # volume
        if lorr == LHS
            result_size = [:dofs_per_loop, :dofs_per_node];
            result_part = IR_data_node(IRtypes.float_data, :source_mat, result_size, [])
            
        else # RHS
            result_size = [:dofs_per_loop];
            result_part = IR_data_node(IRtypes.float_data, :source, result_size, [])
        end
    end
    
    # Separate the factors of each term into test, trial, coef and form the calculation
    if dofsper > 1
        # Submatrices or subvectors for each component
        if lorr == LHS
            submatrices_L = Array{Vector, 2}(undef, dofsper, dofsper);
            submatrices_R = Array{Vector, 2}(undef, dofsper, dofsper);
            for i=1:dofsper
                for j=1:dofsper
                    submatrices_L[i,j] = [];
                    submatrices_R[i,j] = [];
                end
            end
        else # RHS
            submatrices = Array{Vector, 1}(undef, dofsper);
            for j=1:dofsper
                submatrices[j] = [];
            end
        end
        
        for vi=1:length(var) # variables
            # Process the terms for this variable
            for ci=1:length(terms[vi]) # components
                for i=1:length(terms[vi][ci])
                    (term_IR, term_IR2, var_ind) = generate_term_calculation_fvm(terms[vi][ci][i], var, lorr, vors);
                    
                    # Find the appropriate subvector for this term
                    subveci = offset_ind[vi] + ci;
                    subvecj = var_ind;
                    if lorr == LHS
                        subvec_ind = subveci + dofsper * (subvecj-1);
                    else
                        subvec_ind = subveci;
                    end
                    
                    if lorr == LHS
                        if !(term_IR == 0)
                            push!(submatrices_L[subveci, subvecj], term_IR);
                        end
                        if !(term_IR2 == 0)
                            push!(submatrices_R[subveci, subvecj], term_IR2);
                        end
                    else
                        push!(submatrices[subveci], term_IR);
                    end
                end
            end
            
        end # vi
        
        # Put the submatrices together
        num_nonzero_blocks = 0;
        
        if lorr == LHS
            num_nonzero_blocks_L = 0;
            num_nonzero_blocks_R = 0;
            linalg_matrix_block_args_L = [:LINALG_MAT_BLOCKS, 0, 1, result_part];
            linalg_matrix_block_args_R = [:LINALG_MAT_BLOCKS, 0, 1, result_part];
            for smi=1:dofsper
                for smj=1:dofsper
                    if length(submatrices_L[smi,smj]) > 0
                        if length(submatrices_L[smi,smj]) > 1
                            new_term_vec_L = [];
                            push!(new_term_vec_L, :+);
                            append!(new_term_vec_L, submatrices_L[smi,smj]);
                            submat_rhs_L = IR_operation_node(IRtypes.math_op, new_term_vec_L);
                        else
                            submat_rhs_L = submatrices_L[smi,smj][1];
                        end
                        smjplus = IR_operation_node(IRtypes.math_op, [:+, smj, :index_offset]);
                        push!(linalg_matrix_block_args_L, smi);
                        push!(linalg_matrix_block_args_L, smjplus);
                        push!(linalg_matrix_block_args_L, submat_rhs_L);
                        
                        num_nonzero_blocks_L += 1;
                    end
                    if length(submatrices_R[smi,smj]) > 0
                        if length(submatrices_R[smi,smj]) > 1
                            new_term_vec_R = [];
                            push!(new_term_vec_R, :+);
                            append!(new_term_vec_R, submatrices_R[smi,smj]);
                            submat_rhs_R = IR_operation_node(IRtypes.math_op, new_term_vec_R);
                        else
                            submat_rhs_R = submatrices_R[smi,smj][1];
                        end
                        smjplus = IR_operation_node(IRtypes.math_op, [:+, smj, :index_offset, :dofs_per_node])
                        push!(linalg_matrix_block_args_R, smi);
                        push!(linalg_matrix_block_args_R, smjplus);
                        push!(linalg_matrix_block_args_R, submat_rhs_R);
                        
                        num_nonzero_blocks_R += 1;
                    end
                end
            end
            linalg_matrix_block_args_L[2] = num_nonzero_blocks_L;
            linalg_matrix_block_args_R[2] = num_nonzero_blocks_R;
            push!(compute_block.parts, IR_operation_node(IRtypes.named_op, linalg_matrix_block_args_L));
            push!(compute_block.parts, IR_operation_node(IRtypes.named_op, linalg_matrix_block_args_R));
            
        else # RHS
            linalg_vector_block_args = [];
            push!(linalg_vector_block_args, :LINALG_VEC_BLOCKS);
            push!(linalg_vector_block_args, 0);
            push!(linalg_vector_block_args, 1); # number of entries per dof is 1 for FVM
            push!(linalg_vector_block_args, result_part);
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
        if lorr == LHS
            terms = terms[1][1];
            
            if vors == "surface"
                term_vec = [];
                term_vec2 = [];
                
                #process each term
                for i=1:length(terms)
                    (term_IR, term_IR2, var_ind) = generate_term_calculation_fvm(terms[i], var, lorr, vors);
                    
                    push!(term_vec, term_IR);
                    push!(term_vec2, term_IR2);
                end
                if length(term_vec) > 1
                    new_term_vec_L = [];
                    new_term_vec_R = [];
                    push!(new_term_vec_L, :+);
                    push!(new_term_vec_R, :+);
                    append!(new_term_vec_L, term_vec);
                    append!(new_term_vec_R, term_vec2);
                    result_rhs_L = IR_operation_node(IRtypes.math_op, new_term_vec_L);
                    result_rhs_R = IR_operation_node(IRtypes.math_op, new_term_vec_R);
                    
                else
                    result_rhs_L = term_vec[1];
                    result_rhs_R = term_vec2[1];
                end
                
                # NO linalg_matrix_block stuff, just set the values directly
                push!(compute_block.parts, IR_operation_node(IRtypes.assign_op, [apply_indexed_access(result_part, [1,1], IRtypes), result_rhs_L]));
                push!(compute_block.parts, IR_operation_node(IRtypes.assign_op, [apply_indexed_access(result_part, [1,2], IRtypes), result_rhs_R]));
                
            else # volume
                term_vec = [];
                
                #process each term
                for i=1:length(terms)
                    (term_IR, term_IR2, var_ind) = generate_term_calculation_fvm(terms[i], var, lorr, vors);
                    push!(term_vec, term_IR);
                end
                if length(term_vec) > 1
                    new_term_vec = [];
                    push!(new_term_vec, :+);
                    append!(new_term_vec, term_vec);
                    result_rhs = IR_operation_node(IRtypes.math_op, new_term_vec);
                    
                else
                    result_rhs = term_vec[1];
                end
                
                # NO linalg_matrix_block stuff, just set the values directly
                push!(compute_block.parts, IR_operation_node(IRtypes.assign_op, [apply_indexed_access(result_part, [1,1], IRtypes), result_rhs]));
                
            end
            
        else # RHS
            terms = terms[1][1];
            term_vec = Vector{Union{IR_part,Number}}(undef,0);
            
            #process each term
            for i=1:length(terms)
                (term_IR, term_IR2, var_ind) = generate_term_calculation_fvm(terms[i], var, lorr, vors);
                push!(term_vec, term_IR);
            end
            if length(term_vec) > 1
                new_term_vec = [];
                push!(new_term_vec, :+);
                append!(new_term_vec, term_vec);
                result_rhs = IR_operation_node(IRtypes.math_op, new_term_vec);
                
            else
                result_rhs = term_vec[1];
            end
            
            # NO linalg_vector_block stuff, just set the value directly
            push!(compute_block.parts, IR_operation_node(IRtypes.assign_op, [apply_indexed_access(result_part, [1], IRtypes), result_rhs]));
        end
    end
    
    return compute_block;
end

# This takes a term expression that should have a form like var_part * coef_parts
function generate_term_calculation_fvm(term, var, lorr, vors)
    IRtypes = IR_entry_types();
    result1 = 0;
    result2 = 0;
    
    # Note: separate_factors return test and trial info, but FV will not have any test functions.
    # var_part refers to the unknown variable part.
    # example:
    # 0.1*D1__u_1*_FACENORMAL1_1  ->  0.1 * normal[1]                    on LHS
    #                             ->  0.1 * (coef_D1xu_1 .* normal[1])   on RHS
    if lorr == LHS
        (test_part, var_part, coef_part, test_ind, var_ind) = separate_factors(term, var);
        which_side = get_face_side_info(var_part); # 0=none, 1/2=side 1/2, 3=average, 4=neighborhood
        
        if coef_part === nothing # Only an unknown variable part
            if vors == "surface"
                # This is tricky. First make two copies of the coef_part.
                # Since the coef_part could actually contain a variable hidden in an expression like a conditional,
                # Set the correct side to 1 and the other side to 0, or if averaging, set both to 0.5
                var_side1 = copy(var_part);
                var_side2 = copy(var_part);
                var_side1 = replace_lhs_surface_var_entities(var_side1, var, 1, false);
                var_side2 = replace_lhs_surface_var_entities(var_side2, var, 2, false);
                if which_side == 0 # nothing specified. 
                    result1 = arithmetic_expr_to_IR(var_side1);
                    result2 = arithmetic_expr_to_IR(var_side2);
                elseif which_side == 1
                    result1 = arithmetic_expr_to_IR(var_side1);
                elseif which_side == 2
                    result2 = arithmetic_expr_to_IR(var_side2);
                elseif which_side == 3 # average them
                    exp1 = :(0.5*a);
                    exp2 = :(0.5*a);
                    exp1.args[3] = var_side1;
                    exp2.args[3] = var_side2;
                    result1 = arithmetic_expr_to_IR(exp1);
                    result2 = arithmetic_expr_to_IR(exp2);
                else
                    # neighborhoods are not supported for implicit steppers yet
                end
            else
                result1 = arithmetic_expr_to_IR(replace_lhs_surface_var_entities(var_part, var, 0, false));
            end
            
        else # a coefficient part exists
            combined_exp = :(a*b);
            combined_exp.args[2] = coef_part;
            combined_exp.args[3] = var_part;
            if vors == "surface"
                # This is tricky. First make two copies of the coef_part.
                # Since the coef_part could actually contain a variable hidden in an expression like a conditional,
                # Set the correct side to 1 and the other side to 0, or if no side, set both to 0.5
                var_side1 = copy(combined_exp);
                var_side2 = copy(combined_exp);
                var_side1 = replace_lhs_surface_var_entities(var_side1, var, 1, false);
                var_side2 = replace_lhs_surface_var_entities(var_side2, var, 2, false);
                if which_side == 0 # Nothing specified
                    result1 = arithmetic_expr_to_IR(var_side1);
                    result2 = arithmetic_expr_to_IR(var_side2);
                elseif which_side == 1
                    result1 = arithmetic_expr_to_IR(var_side1);
                elseif which_side == 2
                    result2 = arithmetic_expr_to_IR(var_side2);
                elseif which_side == 3 # average them
                    exp1 = :(0.5*a);
                    exp2 = :(0.5*a);
                    exp1.args[3] = var_side1;
                    exp2.args[3] = var_side2;
                    result1 = arithmetic_expr_to_IR(exp1);
                    result2 = arithmetic_expr_to_IR(exp2);
                else
                    # neighborhoods are not supported for implicit steppers yet
                end
            else
                result1 = arithmetic_expr_to_IR(replace_lhs_surface_var_entities(combined_exp, var, 0, false));
            end
        end
        
    else # RHS
        (test_part, var_part, coef_part, test_ind, var_ind) = separate_factors(term);
        # RHS only has coef_part
        if !(coef_part === nothing)
            result1 = arithmetic_expr_to_IR(coef_part);
        else
            result1 = 0;
        end
    end
    
    return (result1, result2, var_ind);
end

# This part inserts the elemental matrix and vector into the global one
function generate_local_to_global_fvm(dofs_per_node, offset_ind; vec_only=false)
    IRtypes = IR_entry_types();
    result_block = IR_block_node([]);
    next_nonzero = IR_data_node(IRtypes.int_data, :next_nonzero_index);
    fid = IR_data_node(IRtypes.int_data, :fid, [], []); # Face ID
    leftel = IR_data_node(IRtypes.int_data, :left_el, [], []); # left element
    rightel = IR_data_node(IRtypes.int_data, :right_el, [], []); # right element
    neighbor = IR_data_node(IRtypes.int_data, :neighbor, [], []);
    if dofs_per_node == 1
        if !vec_only # implicit stepper
            # global_vector[eid] = source + flux + solution[eid];
            result_block = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], [:eid]),
                    IR_operation_node(IRtypes.math_op, [:+,
                        IR_data_node(IRtypes.float_data, :source, [:dofs_per_loop], [1]),
                        IR_data_node(IRtypes.float_data, :flux, [:dofs_per_loop], [1]),
                        IR_data_node(IRtypes.float_data, :solution, [:fv_dofs_global], [:eid])
                    ])
                ])
            ])
        else
            # global_vector[eid] = source + flux;
            result_block = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], [:eid]),
                    IR_operation_node(IRtypes.math_op, [:+,
                        IR_data_node(IRtypes.float_data, :source, [:dofs_per_loop], [1]),
                        IR_data_node(IRtypes.float_data, :flux, [:dofs_per_loop], [1])
                    ])
                ])
            ])
        end
        
        if !vec_only # matrix also needed
            # global_matrix_I[next_nonzero_index] = eid;
            # global_matrix_J[next_nonzero_index] = eid;
            # global_matrix_V[next_nonzero_index] = source_mat[1, 1] + 1;
            push!(result_block.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :global_matrix_I, [:allocatedNZ], [next_nonzero]), :eid
            ]));
            push!(result_block.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :global_matrix_J, [:allocatedNZ], [next_nonzero]), :eid
            ]));
            push!(result_block.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [next_nonzero]),
                IR_operation_node(IRtypes.math_op, [:+,
                    IR_data_node(IRtypes.float_data, :source_mat, [:dofs_per_loop, :dofs_per_node], [1,1]),
                    1
                ])
            ]));
            push!(result_block.parts, IR_operation_node(IRtypes.assign_op, [next_nonzero,
                IR_operation_node(IRtypes.math_op, [:+, next_nonzero, 1])
            ]));
            
            neighbor_col = IR_operation_node(IRtypes.math_op, [:+, :fi, :faces_per_element]);
            face_loop_body = IR_block_node([
                #     fid = mesh.element2face[fi, eid]
                IR_operation_node(IRtypes.assign_op, [fid,
                    IR_operation_node(IRtypes.member_op, [:mesh, 
                        IR_data_node(IRtypes.int_data, :element2face, [:faces_per_element, :num_elements], [:fi, :eid])])
                ]),
                #     left_el = mesh.face2element[1, fid]
                IR_operation_node(IRtypes.assign_op, [leftel,
                    IR_operation_node(IRtypes.member_op, [:mesh, 
                        IR_data_node(IRtypes.int_data, :face2element, [2, :num_faces], [1, :fid])])
                ]),
                #     right_el = mesh.face2element[2, fid]
                IR_operation_node(IRtypes.assign_op, [rightel,
                    IR_operation_node(IRtypes.member_op, [:mesh, 
                        IR_data_node(IRtypes.int_data, :face2element, [2, :num_faces], [2, :fid])])
                ]),
                #     if (eid == right_el)
                #         neighbor = left_el
                #     else
                #         neighbor = right_el
                #     end
                IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :eid, rightel]),
                    IR_block_node([ # true block
                        IR_operation_node(IRtypes.assign_op, [neighbor, leftel])
                    ]),
                    IR_block_node([ # else block
                        IR_operation_node(IRtypes.assign_op, [:neighbor, rightel])
                    ])
                ),
                #     # diagonal left
                #     global_matrix_I[next_nonzero_index] = eid
                #     global_matrix_J[next_nonzero_index] = eid
                #     global_matrix_V[next_nonzero_index] = flux_mat[1,fi];
                IR_comment_node("Diagonal, left element"),
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_data, :global_matrix_I, [:allocatedNZ], [next_nonzero]), :eid
                ]),
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_data, :global_matrix_J, [:allocatedNZ], [next_nonzero]), :eid
                ]),
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [next_nonzero]),
                    IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :?], [1, :fi])
                ]),
                IR_operation_node(IRtypes.assign_op, [next_nonzero,
                    IR_operation_node(IRtypes.math_op, [:+, next_nonzero, 1])
                ]),
                #     if neighbor > 0
                IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>), neighbor, 0]),
                    IR_block_node([ # true block
                        #         # off-diagonal left
                        #         global_index_start = global_index_start + 1; # + dofs_squared
                        #         global_matrix_V[global_index_start] = flux_mat[1,right_index];
                        IR_operation_node(IRtypes.assign_op, [
                            IR_data_node(IRtypes.float_data, :global_matrix_I, [:allocatedNZ], [next_nonzero]), :eid
                        ]),
                        IR_operation_node(IRtypes.assign_op, [
                            IR_data_node(IRtypes.float_data, :global_matrix_J, [:allocatedNZ], [next_nonzero]), :neighbor
                        ]),
                        IR_operation_node(IRtypes.assign_op, [
                            IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [next_nonzero]),
                            IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :?], [1,neighbor_col])
                        ]),
                        IR_operation_node(IRtypes.assign_op, [next_nonzero,
                            IR_operation_node(IRtypes.math_op, [:+, next_nonzero, 1])
                        ]),
                        # #         # diagonal right
                        # #         global_index_start = global_index_start + 1; # + dofs_squared*2
                        # #         global_matrix_V[global_index_start] = -flux_mat[1,right_index];
                        # IR_comment_node("Diagonal, right element"),
                        # IR_operation_node(IRtypes.assign_op, [global_i_st,
                        #     IR_operation_node(IRtypes.math_op, [:+, global_i_st, 1])
                        # ]),
                        # IR_operation_node(IRtypes.assign_op, [
                        #     IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [global_i_st]),
                        #     IR_operation_node(IRtypes.math_op, [:-,
                        #         IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :?], [1,right_index])
                        #     ])
                        # ]),
                        # #         # off-diagonal right
                        # #         global_index_start = global_index_start + 1; # + dofs_squared*3
                        # #         global_matrix_V[global_index_start] = -flux_mat[1,left_index];
                        # IR_comment_node("Off-diagonal, right element"),
                        # IR_operation_node(IRtypes.assign_op, [global_i_st,
                        #     IR_operation_node(IRtypes.math_op, [:+, global_i_st, 1])
                        # ]),
                        # IR_operation_node(IRtypes.assign_op, [
                        #     IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [global_i_st]),
                        #     IR_operation_node(IRtypes.math_op, [:-,
                        #         IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :?], [1,left_index])
                        #     ])
                        # ])
                    ])
                )
            ])
            # Add to global matrix
            # for fi = 1:faces_per_element
            push!(result_block.parts, IR_comment_node("Add to global matrix."));
            push!(result_block.parts, IR_loop_node(IRtypes.space_loop, :faces, :fi, 1, :faces_per_element, face_loop_body));
        end
        
    else # more than one dof
        
        if !vec_only # implicit stepper
            # for i=1:dofs_per_loop
            #     row_index = index_offset + i + dofs_per_node * (eid-1)
            #     global_vector[dofid] = sourcevec[i] + fluxvec[i] + solution[dofid];
            #     
            #     for j=1:dofs_per_node
            #         column_index = (j + (dofs_per_node * (eid - 1)));
            #         global_matrix_V[matind] = source_mat[i, j];
            #         if i==j
            #             global_matrix_V[matind] = global_matrix_V[matind] + 1
            # end
            result_block = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.int_data, :dofs_squared, [], []),
                    IR_operation_node(IRtypes.math_op, [:*, :dofs_per_node, :dofs_per_node])
                ])
                IR_loop_node(IRtypes.space_loop, :dofs, :i, 1, :dofs_per_loop, IR_block_node([
                    IR_operation_node(IRtypes.assign_op, [
                        IR_data_node(IRtypes.int_data, :row_index, [], []),
                        IR_operation_node(IRtypes.math_op, [:+,
                            :index_offset,
                            :i,
                            IR_operation_node(IRtypes.math_op, [:*,
                                :dofs_per_node,
                                IR_operation_node(IRtypes.math_op, [:-,
                                    IR_data_node(IRtypes.int_data, :eid, [], []),
                                    1
                                ])
                            ])
                        ])
                    ]),
                    IR_operation_node(IRtypes.assign_op, [
                        IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], [:row_index]),
                        IR_operation_node(IRtypes.math_op, [:+,
                            IR_data_node(IRtypes.float_data, :source, [:dofs_per_loop], [:i]),
                            IR_data_node(IRtypes.float_data, :flux, [:dofs_per_loop], [:i]),
                            IR_data_node(IRtypes.float_data, :solution, [:fv_dofs_global], [:row_index])
                        ])
                    ]),
                    
                    IR_loop_node(IRtypes.space_loop, :dofs, :j, 1, :dofs_per_node, IR_block_node([
                        # if abs(source_mat[i, j]) > 0.0 || (i==(j-index_offset)) 
                        IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(||), 
                                            IR_operation_node(IRtypes.math_op, [:(==), :i, 
                                                IR_operation_node(IRtypes.math_op, [:-, :j, :index_offset])]),
                                            IR_operation_node(IRtypes.math_op, [:(>),
                                                IR_operation_node(IRtypes.math_op, [:(abs),
                                                    IR_data_node(IRtypes.float_data, :source_mat, [:dofs_per_loop, :dofs_per_loop], [:i, :j])]),
                                                    0.0
                                                ])
                                            ]),
                            IR_block_node([
                                # column_index = (j + (dofs_per_node * (eid - 1)));
                                IR_operation_node(IRtypes.assign_op, [
                                    IR_data_node(IRtypes.int_data, :column_index, [], []),
                                    IR_operation_node(IRtypes.math_op, [:+,
                                        :j,
                                        IR_operation_node(IRtypes.math_op, [:*,
                                            :dofs_per_node,
                                            IR_operation_node(IRtypes.math_op, [:-,
                                                IR_data_node(IRtypes.int_data, :eid, [], []),
                                                1
                                            ])
                                        ])
                                    ])
                                ]),
                                # global_matrix_I[next_nonzero_index] = row_index;
                                # global_matrix_J[next_nonzero_index] = column_index;
                                # global_matrix_V[next_nonzero_index] = source_mat[i, j];
                                IR_operation_node(IRtypes.assign_op, [
                                    IR_data_node(IRtypes.float_data, :global_matrix_I, [:allocatedNZ], [next_nonzero]), :row_index
                                ]),
                                IR_operation_node(IRtypes.assign_op, [
                                    IR_data_node(IRtypes.float_data, :global_matrix_J, [:allocatedNZ], [next_nonzero]), :column_index
                                ]),
                                IR_operation_node(IRtypes.assign_op, [
                                    IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [next_nonzero]),
                                    IR_data_node(IRtypes.float_data, :source_mat, [:dofs_per_loop, :dofs_per_loop], [:i, :j])
                                ]),
                                
                                # if i==(j-index_offset) add one
                                IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :i, 
                                                    IR_operation_node(IRtypes.math_op, [:-, :j, :index_offset])]),
                                    IR_block_node([
                                        IR_operation_node(IRtypes.assign_op, [
                                            IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [next_nonzero]),
                                            IR_operation_node(IRtypes.math_op, [:+,
                                            IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [next_nonzero]),
                                                1
                                            ])
                                        ])
                                    ])
                                ),
                                IR_operation_node(IRtypes.assign_op, [next_nonzero,
                                    IR_operation_node(IRtypes.math_op, [:+, next_nonzero, 1])
                                ]),
                            ])
                        ),
                        
                    ]))
                ]))
            ])
        else
            # for i=1:dofs_per_node
            #     row_index = i + dofs_per_node * (eid-1)
            #     global_vector[row_index] = sourcevec[i] + fluxvec[i];
            # end
            result_block = IR_block_node([
                IR_loop_node(IRtypes.space_loop, :dofs, :i, 1, :dofs_per_loop, IR_block_node([
                    IR_operation_node(IRtypes.assign_op, [
                        IR_data_node(IRtypes.int_data, :row_index, [], []),
                        IR_operation_node(IRtypes.math_op, [:+,
                            :index_offset,
                            :i,
                            IR_operation_node(IRtypes.math_op, [:*,
                                :dofs_per_node,
                                IR_operation_node(IRtypes.math_op, [:-,
                                    IR_data_node(IRtypes.int_data, :eid, [], []),
                                    1
                                ])
                            ])
                        ])
                    ]),
                    IR_operation_node(IRtypes.assign_op, [
                        IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], [:row_index]),
                        IR_operation_node(IRtypes.math_op, [:+,
                            IR_data_node(IRtypes.float_data, :source, [:dofs_per_loop], [:i]),
                            IR_data_node(IRtypes.float_data, :flux, [:dpfs_per_loop], [:i])
                        ])
                    ])
                ]))
            ])
        end
        
            
        if !vec_only # matrix also
            # Note: source part was done above. Just do flux now.
            LR_ind_offset = IR_operation_node(IRtypes.math_op, [:+,
                IR_operation_node(IRtypes.math_op, [:*,
                    IR_operation_node(IRtypes.math_op, [:-, :fi, 1]),
                    :dofs_per_node
                ]),
                1
            ])
            face_loop_body = IR_block_node([
                #     fid = mesh.element2face[fi, eid]
                IR_operation_node(IRtypes.assign_op, [fid,
                    IR_operation_node(IRtypes.member_op, [:mesh, 
                        IR_data_node(IRtypes.int_data, :element2face, [:faces_per_element, :num_elements], [:fi, :eid])])
                ]),
                #     left_el = mesh.face2element[1, fid]
                IR_operation_node(IRtypes.assign_op, [leftel,
                    IR_operation_node(IRtypes.member_op, [:mesh, 
                        IR_data_node(IRtypes.int_data, :face2element, [2, :num_faces], [1, :fid])])
                ]),
                #     right_el = mesh.face2element[2, fid]
                IR_operation_node(IRtypes.assign_op, [rightel,
                    IR_operation_node(IRtypes.member_op, [:mesh, 
                        IR_data_node(IRtypes.int_data, :face2element, [2, :num_faces], [2, :fid])])
                ]),
                #     if (eid == right_el)
                #         neighbor = left_el
                #     else
                #         neighbor = right_el
                #     end
                IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :eid, rightel]),
                    IR_block_node([ # true block
                        IR_operation_node(IRtypes.assign_op, [neighbor, leftel])
                    ]),
                    IR_block_node([ # else block
                        IR_operation_node(IRtypes.assign_op, [:neighbor, rightel])
                    ])
                ),
                #     # diagonal left
                IR_comment_node("Diagonal, left element"),
                # for i=1:dofs_per_loop
                #   row_index = (index_offset + i + (dofs_per_node * (eid - 1)));
                #   for j=1:dofs_per_node
                #         column_index = (j + (dofs_per_node * (eid - 1)));
                #         global_matrix_I[next_nonzero_index] = row_index;
                #         global_matrix_J[next_nonzero_index] = column_index;
                #         global_matrix_V[next_nonzero_index] = flux_mat[i, j + dofs_per_node*(fi-1)];
                #         next_nonzero_index += 1;
                IR_loop_node(IRtypes.space_loop, :dofs, :i, 1, :dofs_per_loop, IR_block_node([
                    #   row_index = (index_offset + i + (dofs_per_node * (eid - 1)));
                    IR_operation_node(IRtypes.assign_op, [
                        IR_data_node(IRtypes.int_data, :row_index, [], []),
                        IR_operation_node(IRtypes.math_op, [:+,
                            :index_offset,
                            :i,
                            IR_operation_node(IRtypes.math_op, [:*,
                                :dofs_per_node,
                                IR_operation_node(IRtypes.math_op, [:-,
                                    IR_data_node(IRtypes.int_data, :eid, [], []),
                                    1
                                ])
                            ])
                        ])
                    ]),
                    IR_loop_node(IRtypes.space_loop, :dofs, :j, 1, :dofs_per_node, IR_block_node([
                        # if abs(flux_mat[i, j + dofs_per_node*(fi-1)]) > 0.0
                        IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>),
                                                IR_operation_node(IRtypes.math_op, [:(abs),
                                                    IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :dofs_per_loop], [:i,
                                                        IR_operation_node(IRtypes.math_op, [:+,
                                                            :j,
                                                            IR_operation_node(IRtypes.math_op, [:*,
                                                                :dofs_per_node,
                                                                IR_operation_node(IRtypes.math_op, [:-,
                                                                    IR_data_node(IRtypes.int_data, :fi, [], []),
                                                                    1
                                                                ])
                                                            ])
                                                        ])
                                                    ])
                                                ]),
                                                0.0
                                            ]),
                            IR_block_node([
                                # column_index = (j + (dofs_per_node * (eid - 1)));
                                IR_operation_node(IRtypes.assign_op, [
                                    IR_data_node(IRtypes.int_data, :column_index, [], []),
                                    IR_operation_node(IRtypes.math_op, [:+,
                                        :j,
                                        IR_operation_node(IRtypes.math_op, [:*,
                                            :dofs_per_node,
                                            IR_operation_node(IRtypes.math_op, [:-,
                                                IR_data_node(IRtypes.int_data, :eid, [], []),
                                                1
                                            ])
                                        ])
                                    ])
                                ]),
                                # global_matrix_I[next_nonzero_index] = row_index;
                                # global_matrix_J[next_nonzero_index] = column_index;
                                # global_matrix_V[next_nonzero_index] = flux_mat[i, j + dofs_per_node*(fi-1)];
                                # next_nonzero_index += 1;
                                IR_operation_node(IRtypes.assign_op, [
                                    IR_data_node(IRtypes.float_data, :global_matrix_I, [:allocatedNZ], [next_nonzero]), :row_index
                                ]),
                                IR_operation_node(IRtypes.assign_op, [
                                    IR_data_node(IRtypes.float_data, :global_matrix_J, [:allocatedNZ], [next_nonzero]), :column_index
                                ]),
                                IR_operation_node(IRtypes.assign_op, [
                                    IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [next_nonzero]),
                                    IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :dofs_per_loop], [:i,
                                        IR_operation_node(IRtypes.math_op, [:+,
                                            :j,
                                            IR_operation_node(IRtypes.math_op, [:*,
                                                :dofs_per_node,
                                                IR_operation_node(IRtypes.math_op, [:-,
                                                    IR_data_node(IRtypes.int_data, :fi, [], []),
                                                    1
                                                ])
                                            ])
                                        ])
                                    ])
                                ]),
                                IR_operation_node(IRtypes.assign_op, [next_nonzero,
                                    IR_operation_node(IRtypes.math_op, [:+, next_nonzero, 1])
                                ])
                            ])
                        )
                    ]))
                ])),
                #     if neighbor > 0
                IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>), neighbor, 0]),
                    IR_block_node([ # true block
                        #         # off-diagonal left
                        IR_comment_node("Off-diagonal, left element"),
                        # for i=1:dofs_per_loop
                        #   row_index = (index_offset + i + (dofs_per_node * (eid - 1)));
                        #   for j=1:dofs_per_node
                        #         column_index = (j + (dofs_per_node * (neighbor - 1)));
                        #         global_matrix_I[next_nonzero_index] = row_index;
                        #         global_matrix_J[next_nonzero_index] = column_index;
                        #         global_matrix_V[next_nonzero_index] = flux_mat[i, j + dofs_per_node*(fi-1) + dofs_per_node*faces_per_element];
                        #         next_nonzero_index += 1;
                        IR_loop_node(IRtypes.space_loop, :dofs, :i, 1, :dofs_per_loop, IR_block_node([
                            #   row_index = (index_offset + i + (dofs_per_node * (eid - 1)));
                            IR_operation_node(IRtypes.assign_op, [
                                IR_data_node(IRtypes.int_data, :row_index, [], []),
                                IR_operation_node(IRtypes.math_op, [:+,
                                    :index_offset,
                                    :i,
                                    IR_operation_node(IRtypes.math_op, [:*,
                                        :dofs_per_node,
                                        IR_operation_node(IRtypes.math_op, [:-,
                                            IR_data_node(IRtypes.int_data, :eid, [], []),
                                            1
                                        ])
                                    ])
                                ])
                            ]),
                            IR_loop_node(IRtypes.space_loop, :dofs, :j, 1, :dofs_per_node, IR_block_node([
                                # if abs(flux_mat[i, j + dofs_per_node*(fi-1) + dofs_per_node*faces_per_element]) > 0.0
                                IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>),
                                                        IR_operation_node(IRtypes.math_op, [:(abs),
                                                            IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :dofs_per_loop], [:i,
                                                                IR_operation_node(IRtypes.math_op, [:+,
                                                                    :j,
                                                                    IR_operation_node(IRtypes.math_op, [:*,
                                                                        :dofs_per_node,
                                                                        IR_operation_node(IRtypes.math_op, [:+,
                                                                            IR_data_node(IRtypes.int_data, :fi, [], []),
                                                                            :faces_per_element,
                                                                            -1
                                                                        ])
                                                                    ])
                                                                ])
                                                            ])
                                                        ]),
                                                        0.0
                                                    ]),
                                    IR_block_node([
                                        # column_index = (j + (dofs_per_node * (neighbor - 1)));
                                        IR_operation_node(IRtypes.assign_op, [
                                            IR_data_node(IRtypes.int_data, :column_index, [], []),
                                            IR_operation_node(IRtypes.math_op, [:+,
                                                :j,
                                                IR_operation_node(IRtypes.math_op, [:*,
                                                    :dofs_per_node,
                                                    IR_operation_node(IRtypes.math_op, [:-,
                                                        IR_data_node(IRtypes.int_data, :neighbor, [], []),
                                                        1
                                                    ])
                                                ])
                                            ])
                                        ]),
                                        # global_matrix_I[next_nonzero_index] = row_index;
                                        # global_matrix_J[next_nonzero_index] = column_index;
                                        # global_matrix_V[next_nonzero_index] = flux_mat[i, j + dofs_per_node*(fi-1) + dofs_per_node*faces_per_element];
                                        # next_nonzero_index += 1;
                                        IR_operation_node(IRtypes.assign_op, [
                                            IR_data_node(IRtypes.float_data, :global_matrix_I, [:allocatedNZ], [next_nonzero]), :row_index
                                        ]),
                                        IR_operation_node(IRtypes.assign_op, [
                                            IR_data_node(IRtypes.float_data, :global_matrix_J, [:allocatedNZ], [next_nonzero]), :column_index
                                        ]),
                                        IR_operation_node(IRtypes.assign_op, [
                                            IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [next_nonzero]),
                                            IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :dofs_per_loop], [:i,
                                                IR_operation_node(IRtypes.math_op, [:+,
                                                    :j,
                                                    IR_operation_node(IRtypes.math_op, [:*,
                                                        :dofs_per_node,
                                                        IR_operation_node(IRtypes.math_op, [:+,
                                                            IR_data_node(IRtypes.int_data, :fi, [], []),
                                                            :faces_per_element,
                                                            -1
                                                        ])
                                                    ])
                                                ])
                                            ])
                                        ]),
                                        IR_operation_node(IRtypes.assign_op, [next_nonzero,
                                            IR_operation_node(IRtypes.math_op, [:+, next_nonzero, 1])
                                        ])
                                    ])
                                )
                            ]))
                        ])),
                        # #         # diagonal right
                        # IR_comment_node("Diagonal, right element"),
                        # # for i=1:dofs_per_loop
                        # #   for j=1:dofs_per_node
                        # #     global_matrix_V[mat_ind] = -flux_mat[i,right_index + j - 1];
                        # IR_loop_node(IRtypes.space_loop, :dofs, :i, 1, :dofs_per_loop, IR_block_node([
                        #     IR_loop_node(IRtypes.space_loop, :dofs, :j, 1, :dofs_per_node, IR_block_node([
                        #         IR_operation_node(IRtypes.assign_op, [
                        #             IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [:mat_ind]),
                        #             IR_operation_node(IRtypes.math_op, [:-,
                        #                 IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :?], [:i,IR_operation_node(IRtypes.math_op, [:+,right_index,:j,-1])])
                        #             ])
                        #         ]),
                        #         IR_operation_node(IRtypes.assign_op, [
                        #             IR_data_node(IRtypes.int_data, :mat_ind, [], []),
                        #             IR_operation_node(IRtypes.math_op, [:+,
                        #                 IR_data_node(IRtypes.int_data, :mat_ind, [], []),
                        #                 1
                        #             ])
                        #         ])
                        #     ]))
                        # ])),
                        # #         # off-diagonal right
                        # IR_comment_node("Off-diagonal, right element"),
                        # # for i=1:dofs_per_loop
                        # #   for j=1:dofs_per_node
                        # #     global_matrix_V[mat_ind] = -flux_mat[i,left_index + j - 1];
                        # IR_loop_node(IRtypes.space_loop, :dofs, :i, 1, :dofs_per_loop, IR_block_node([
                        #     IR_loop_node(IRtypes.space_loop, :dofs, :j, 1, :dofs_per_node, IR_block_node([
                        #         IR_operation_node(IRtypes.assign_op, [
                        #             IR_data_node(IRtypes.float_data, :global_matrix_V, [:allocatedNZ], [:mat_ind]),
                        #             IR_operation_node(IRtypes.math_op, [:-,
                        #                 IR_data_node(IRtypes.float_data, :flux_mat, [:dofs_per_loop, :?], [:i,IR_operation_node(IRtypes.math_op, [:+,left_index,:j,-1])])
                        #             ])
                        #         ]),
                        #         IR_operation_node(IRtypes.assign_op, [
                        #             IR_data_node(IRtypes.int_data, :mat_ind, [], []),
                        #             IR_operation_node(IRtypes.math_op, [:+,
                        #                 IR_data_node(IRtypes.int_data, :mat_ind, [], []),
                        #                 1
                        #             ])
                        #         ])
                        #     ]))
                        # ]))
                    ])
                )
            ])
            # Add to global matrix
            # for fi = 1:faces_per_element
            push!(result_block.parts, IR_comment_node("Add to global matrix."));
            push!(result_block.parts, IR_loop_node(IRtypes.space_loop, :faces, :fi, 1, :faces_per_element, face_loop_body));
        end
        
    end
    
    return result_block;
end

# Generate the assembly loop structures and insert the content
function generate_assembly_loop_fvm(var, indices, need_parent)
    IRtypes = IR_entry_types();
    
    index_names = [];
    index_updates = [];
    
    # Make names for all of the index variables like INDEX_VAL_something
    elements_included = false;
    ind_shift = 0;
    for i=1:length(indices)
        if indices[i].symbol === :FINCHELEMENTS
            elements_included = true;
            push!(index_names, "elements");
        else
            # INDEX_VAL_something
            push!(index_names, IR_data_node(IRtypes.int_data, Symbol("INDEX_VAL_"*string(indices[i].symbol))));
            # index_values[4] = INDEX_VAL_something
            push!(index_updates, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.int_data, :index_values, [:?], [indices[i].tag]),
                index_names[end]
            ]))
        end
    end
    
    # If elements were not included, make them the outermost loop
    if !elements_included && length(index_names) > 0
        index_names = ["elements"; index_names];
        index_updates = [IR_comment_node(""); index_updates];
        ind_shift = 1;
    end
    
    # Make an offset that is set at the beginning of the loop body
    # We have to assume that unknown variables have the same set of indices.
    # If they don't, we need much more complicated logic.
    index_offset = 0; # offset in variable values for this index
    if length(var[1].indexer) == 0
        # no indexer. 
    else
        # has indices
        index_offset = Symbol("INDEX_VAL_"*string(var[1].indexer[1].symbol));
        prev_ind = length(var[1].indexer[1].range);
        for i=2:length(var[1].indexer)
            index_offset = IR_operation_node(IRtypes.math_op, [:+, index_offset,
                IR_operation_node(IRtypes.math_op, [:*, prev_ind, 
                    IR_operation_node(IRtypes.math_op, [:-, Symbol("INDEX_VAL_"*string(var[1].indexer[i].symbol)), 1])])]);
            prev_ind *= length(var[1].indexer[i].range);
        end
        index_offset = IR_operation_node(IRtypes.math_op, [:-, index_offset, 1]);
    end
    
    el_order = IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.int_data, :elemental_order, [:num_elements], [:ei])]);
    child = IR_operation_node(IRtypes.member_op, [:fv_info, 
                IR_operation_node(IRtypes.member_op, [:parentMaps,
                    IR_data_node(IRtypes.int_data, :parent2child, [:num_children, :num_parents], [:child_id, :parent_id])
                ])
            ]);
    
    # Placeholder computation that will be filled in later
    if need_parent
        placeholder = IR_block_node([IR_operation_node(IRtypes.assign_op, [:eid, child])], "assembly")
    else
        placeholder = IR_block_node([IR_operation_node(IRtypes.assign_op, [:eid, el_order])], "assembly");
    end
    for i=1:length(index_updates)
        push!(placeholder.parts, index_updates[i]);
    end
    push!(placeholder.parts, IR_operation_node(IRtypes.assign_op, [:index_offset, index_offset]));
    push!(placeholder.parts, IR_comment_node("Begin assembly code"));
    
    # generate the loop structures
    if length(index_names) > 0
        # The innermost loop holds placeholder
        if index_names[end] == "elements"
            if need_parent
                assembly_loop = IR_loop_node(IRtypes.space_loop, :parents, :parent_id, 1, :num_parents, 
                                    IR_loop_node(IRtypes.space_loop, :children, :child_id, 1, :num_children, placeholder));
            else
                assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :ei, 1, :num_elements, placeholder);
            end
            
        else
            assembly_loop = IR_loop_node(IRtypes.index_loop, indices[end].symbol, 
                            index_names[end].label, indices[end].range[1], indices[end].range[end], placeholder);
        end
        
        # work outwards nesting assembly_loop
        for i=(length(index_names)-1):-1:1
            if index_names[i] == "elements"
                if need_parent
                    assembly_loop = IR_loop_node(IRtypes.space_loop, :parents, :parent_id, 1, :num_parents, 
                                        IR_loop_node(IRtypes.space_loop, :children, :child_id, 1, :num_children, assembly_loop));
                else
                    assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :ei, 1, :num_elements, assembly_loop);
                end
                
            elseif i > ind_shift
                assembly_loop = IR_loop_node(IRtypes.index_loop, indices[i-ind_shift].symbol, 
                                index_names[i].label, indices[i-ind_shift].range[1], 
                                indices[i-ind_shift].range[end], assembly_loop);
            end
        end
    else # only an element loop
        if need_parent
            assembly_loop = IR_loop_node(IRtypes.space_loop, :parents, :parent_id, 1, :num_parents,
                                IR_loop_node(IRtypes.space_loop, :children, :child_id, 1, :num_children, placeholder));
        else
            assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :ei, 1, :num_elements, placeholder);
        end
    end
    
    return assembly_loop
end

# Generates the time stepping loop using the supplied assembly block and stepper
function generate_time_stepping_loop_fvm(stepper, assembly, prob)
    IRtypes = IR_entry_types();
    tloop_body = IR_block_node([]);
    stepper_nsteps = IR_operation_node(IRtypes.member_op, [:time_stepper, :Nsteps]);
    zero_bdry_done = IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY,
        IR_data_node(IRtypes.int_data, :bdry_done, [:num_faces], []),
        0, :num_faces]);
    zero_flux_done = IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY,
        IR_data_node(IRtypes.int_data, :face_flux_done, [:num_faces], []),
        :false, :num_faces]);
    reset_nonzero_counter = IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.int_data, :next_nonzero_index), 1]);
    
    # Exchange ghosts if needed
    ghost_exchange = IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>), :num_partitions, 1]),
                        IR_block_node([IR_operation_node(IRtypes.named_op, [:GHOST_EXCHANGE_FV])]));
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
    
    # For partitioned meshes
    gsolvec = IR_data_node(IRtypes.float_data, :global_solution, [:fv_dofs_global], []);
    solvec = IR_data_node(IRtypes.float_data, :solution, [:dofs_partition], []);
    gather_system = IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>), :num_partitions, 1]),
        IR_block_node([IR_operation_node(IRtypes.named_op, [:GLOBAL_GATHER_SYSTEM, :FV])]),
        IR_block_node([IR_operation_node(IRtypes.named_op, [:GLOBAL_FORM_MATRIX])])
    );
    distribute_solution = IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>), :num_partitions, 1]),
        IR_block_node([IR_operation_node(IRtypes.named_op, [:GLOBAL_DISTRIBUTE_VECTOR, gsolvec, solvec, :FV])]),
        IR_block_node([IR_operation_node(IRtypes.assign_op, [solvec, gsolvec])])
    );
    gather_solve_distribute = IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>), :num_partitions, 1]),
        IR_block_node([
            IR_operation_node(IRtypes.named_op, [:GLOBAL_GATHER_SYSTEM, :FV]),
            wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :global_solution, :full_global_matrix, :full_global_vector])),
            IR_operation_node(IRtypes.named_op, [:GLOBAL_DISTRIBUTE_VECTOR, gsolvec, solvec, :FV])
        ]),
        IR_block_node([
            IR_operation_node(IRtypes.named_op, [:GLOBAL_FORM_MATRIX]),
            wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector]))
        ])
    );
    
    # A progress meter
    # progress = 100 * ti / time_stepper.Nsteps;
    # last_minor = 0;
    # last_major = 0;
    # if progress >= last_major + 10
    #   last_major += 10
    #   last_minor += 2
    #   print(string(last_major + 10));
    # elseif progress >= last_minor + 2
    #   last_minor += 2
    #   print(".")
    # end
    last_minor = IR_data_node(IRtypes.int_data, :last_minor_progress);
    last_major = IR_data_node(IRtypes.int_data, :last_major_progress);
    root_print = IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :proc_rank, 0]),
                        IR_block_node([]));
    progress_init = IR_block_node([
        IR_operation_node(IRtypes.assign_op, [last_minor, 0]),
        IR_operation_node(IRtypes.assign_op, [last_major, 0]),
        IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :proc_rank, 0]),
            IR_block_node([IR_operation_node(IRtypes.named_op, [:PRINT_STRING, "Time step progress(%) 0"])]))
        
    ])
    progress_update = IR_block_node([
        IR_conditional_node(
            IR_operation_node(IRtypes.math_op, [:(>=), 
                IR_operation_node(IRtypes.math_op, [:*, 100.0, IR_operation_node(IRtypes.math_op, [:/, :ti, stepper_nsteps])]),
                IR_operation_node(IRtypes.math_op, [:+, last_major, 10])]),
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [last_major, IR_operation_node(IRtypes.math_op, [:+, last_major, 10])]),
                IR_operation_node(IRtypes.assign_op, [last_minor, IR_operation_node(IRtypes.math_op, [:+, last_minor, 2])]),
                IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :proc_rank, 0]),
                    IR_block_node([IR_operation_node(IRtypes.named_op, [:PRINT_STRING, last_major])]))
            ]),
            IR_block_node([
                IR_conditional_node(
                    IR_operation_node(IRtypes.math_op, [:(>=), 
                        IR_operation_node(IRtypes.math_op, [:*, 100.0, IR_operation_node(IRtypes.math_op, [:/, :ti, stepper_nsteps])]),
                        IR_operation_node(IRtypes.math_op, [:+, last_minor, 2])]),
                    IR_block_node([
                        IR_operation_node(IRtypes.assign_op, [last_minor, IR_operation_node(IRtypes.math_op, [:+, last_minor, 2])]),
                        IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :proc_rank, 0]),
                            IR_block_node([IR_operation_node(IRtypes.named_op, [:PRINT_STRING, "."])]))
                    ])
                )
            ])
        )
    ])
    
    # Which type of stepper is it?
    if stepper.type == EULER_EXPLICIT
        # This part is used if the solution is updated like sol = sol + dt*(dsol)
        sol_i = IR_data_node(IRtypes.float_data, :solution, [:fv_dofs_partition], [:update_i]);
        dsol_i = IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], [:update_i]);
        update_loop = IR_loop_node(IRtypes.dof_loop, :dof, :update_i, 1, :fv_dofs_partition, IR_block_node([
            IR_operation_node(IRtypes.assign_op, [
                sol_i,
                IR_operation_node(IRtypes.math_op, [:+, sol_i, IR_operation_node(IRtypes.math_op, [:*, :dt, dsol_i])])
            ])
        ]))
        tloop_body.parts = [
            zero_bdry_done,
            zero_flux_done,
            reset_nonzero_counter,
            ghost_exchange,
            pre_step_call,
            wrap_in_timer(:step_assembly, assembly),
            # before updating the bdry vals may need to be put in vars and zeroed in solution
            wrap_in_timer(:update_sol, update_loop),
            IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :solution]),
            IR_operation_node(IRtypes.named_op, [:SCATTER_VARS, :solution]),
            post_step_call,
            IR_operation_node(IRtypes.assign_op, [
                :t,
                IR_operation_node(IRtypes.math_op, [:+, :t, :dt])
            ]),
            progress_update
        ];
        
        # # This part is used if the update is coded into the system sol = A\b (not sol = sol + dt*(A\b))
        # tloop_body.parts = [
        #     zero_bdry_done,
        #     wrap_in_timer(:step_assembly, assembly),
        #     wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector])),
        #     IR_operation_node(IRtypes.named_op, [:SCATTER_VARS, :solution]),
        #     IR_operation_node(IRtypes.assign_op, [
        #         :t,
        #         IR_operation_node(IRtypes.math_op, [:+, :t, stepper.dt])
        #     ])
        # ];
        
        time_loop = IR_block_node([
            progress_init,
            IR_loop_node(IRtypes.time_loop, :time, :ti, 1, stepper_nsteps, tloop_body)
        ]);
       
        # multistage explicit steppers
    elseif stepper.type == LSRK4
        # LSRK4 is a special case, low storage
    
        # Low storage RK4: 
        # p0 = u
        #   ki = ai*k(i-1) + dt*f(p(i-1), t+ci*dt)
        #   pi = p(i-1) + bi*ki
        # u = p5
        
        #=
                tmppi = get_var_vals(var, tmppi);
                tmpki = zeros(size(sol));
                for rki=1:stepper.stages
                    rktime = t + stepper.c[rki]*stepper.dt;
                    # p(i-1) is currently in u
                    
                    if config.num_partitions > 1 exchange_ghosts(var, grid, i); end
                    pre_step_function();
                    
                    sol = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, rktime, stepper.dt, assemble_loops=assemble_func);
                    
                    if rki == 1 # because a1 == 0
                        tmpki = stepper.dt .* sol;
                    else
                        tmpki = stepper.a[rki].*tmpki + stepper.dt.*sol;
                    end
                    tmppi = tmppi + stepper.b[rki].*tmpki
                    
                    FV_copy_bdry_vals_to_vector(var, tmppi, grid, dofs_per_node);
                    place_vector_in_vars(var, tmppi, stepper);
                    
                    post_step_function();
                end
        =#
        tmp_pi = IR_data_node(IRtypes.float_data, :tmppi, [:fv_dofs_partition], []);
        tmp_ki = IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], []);
        piki_loop_one = IR_loop_node(IRtypes.space_loop, :dofs, :piki_i, 1, :fv_dofs_partition, IR_block_node([
            # tmpki .= stepper.dt .* global_vector;
            # tmppi .= tmppi + stepper.b[rki].*tmpki;
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], [:piki_i]),
                IR_operation_node(IRtypes.math_op, [:*, :dt, IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], [:piki_i])])
            ]),
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :tmppi, [:fv_dofs_partition], [:piki_i]),
                IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.float_data, :tmppi, [:fv_dofs_partition], [:piki_i]), 
                    IR_operation_node(IRtypes.math_op, [:*, 
                        IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_data, :b, [:?], [:rki])]), 
                        IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], [:piki_i])])
                ])
            ])
        ]))
        piki_loop_two = IR_loop_node(IRtypes.space_loop, :dofs, :piki_i, 1, :fv_dofs_partition, IR_block_node([
            # tmpki .= stepper.a[rki].*tmpki + stepper.dt .* global_vector;
            # tmppi .= tmppi + stepper.b[rki].*tmpki;
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], [:piki_i]),
                IR_operation_node(IRtypes.math_op, [:+, 
                    IR_operation_node(IRtypes.math_op, [:*, 
                        IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_data, :a, [:?], [:rki])]), 
                        IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], [:piki_i])]),
                    IR_operation_node(IRtypes.math_op, [:*, :dt, IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], [:piki_i])])
                ])
            ]),
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :tmppi, [:fv_dofs_partition], [:piki_i]),
                IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.float_data, :tmppi, [:fv_dofs_partition], [:piki_i]), 
                    IR_operation_node(IRtypes.math_op, [:*, 
                        IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_data, :b, [:?], [:rki])]), 
                        IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], [:piki_i])])
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
                        IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_data, :c, [:?], [:rki])])])])
            ]),
            # assemble
            zero_bdry_done,
            zero_flux_done,
            reset_nonzero_counter,
            ghost_exchange,
            pre_step_call,
            wrap_in_timer(:step_assembly, assembly),
            # update tmppi and tmpki
            piki_condition,
            # copy_bdry_vals_to_vector(var, tmppi, grid_data, dofs_per_node);
            IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :tmppi]),
            IR_operation_node(IRtypes.named_op, [:SCATTER_VARS, :tmppi]),
            post_step_call
        ]);
        stage_loop = IR_loop_node(IRtypes.time_loop, :stages, :rki, 1, stepper.stages, stage_loop_body);
        
        push!(tloop_body.parts, IR_operation_node(IRtypes.named_op, [:GATHER_VARS, :tmppi]));
        push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:last_t, :t]));
        push!(tloop_body.parts, stage_loop);
        push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:t, IR_operation_node(IRtypes.math_op, [:+, :last_t, :dt])]));
        push!(tloop_body.parts, progress_update);
        
        time_loop = IR_block_node([
            IR_operation_node(IRtypes.assign_op, [
                tmp_pi,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :fv_dofs_partition])]),
            IR_operation_node(IRtypes.assign_op, [
                tmp_ki,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :fv_dofs_partition])]),
            progress_init,
            IR_loop_node(IRtypes.time_loop, :time, :ti, 1, stepper_nsteps, tloop_body)
        ]);
    elseif stepper.stages > 1
        # Explicit multi-stage methods: 
        # x = x + dt*sum(bi*ki)
        # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
        
        #=
                last_result = zeros(size(b));
                tmpresult = zeros(size(b));
                tmpki = zeros(length(b), stepper.stages);
                # will hold the final result
                last_result = get_var_vals(var, last_result);
                # will be placed in var.values for each stage
                # Storage for each stage
                # tmpki = zeros(length(b), stepper.stages);
                
                
                for stage=1:stepper.stages
                    stime = t + stepper.c[stage]*stepper.dt;
                    
                    # Update the values in vars to be used in this stage
                    if stage > 1
                        initialized_tmpvals = false;
                        for j=1:stage
                            if stepper.a[stage, j] > 0
                                if !initialized_tmpvals
                                    initialized_tmpvals = true;
                                    for k=1:length(sol)
                                        tmpvals[k] = sol[k] + stepper.dt * stepper.a[stage, j] * tmpki[k,j];
                                    end
                                else
                                    for k=1:length(sol)
                                        tmpvals[k] += stepper.dt * stepper.a[stage, j] * tmpki[k,j];
                                    end
                                end
                            end
                        end
                        
                        if initialized_tmpvals
                            FV_copy_bdry_vals_to_vector(var, tmpvals, grid, dofs_per_node);
                            place_vector_in_vars(var, tmpvals, stepper);
                        end
                        post_step_function(); # seems weird, but imagine this is happening after stage-1
                    end
                    
                    if config.num_partitions > 1 exchange_ghosts(var, grid, i); end
                    pre_step_function();
                    
                    tmpki[:,stage] = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, stime, stepper.dt, assemble_loops=assemble_func);
                    
                end
                # Stages are done. Assemble the final result for this step
                for stage=1:stepper.stages
                    sol += stepper.dt * stepper.b[stage] .* tmpki[:, stage];
                end
                FV_copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node);
                place_vector_in_vars(var, sol, stepper);
                
                post_step_function();
        =#
        
        tmp_last = IR_data_node(IRtypes.float_data, :last_result, [:fv_dofs_partition], []);
        tmp_result = IR_data_node(IRtypes.float_data, :tmpresult, [:fv_dofs_partition], []);
        tmp_ki = IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], []);
        
        update_ki_loop_one = IR_loop_node(IRtypes.space_loop, :dofs, :k, 1, :fv_dofs_partition, IR_block_node([
            # # Update the values in vars to be used in the next stage
            # for k=1:fv_dofs_partition
            #     tmpki[k,stage] = global_vector[k];
            #     tmpresult[k] = last_result[k];
            #     for j=1:stage
            #         tmpresult[k] += stepper.dt * stepper.a[stage+1, j] * tmpki[k,j];
            #     end
            # end
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], [:k, :rki]),
                IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], [:k])
            ]),
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :tmpresult, [:fv_dofs_partition], [:k]),
                IR_data_node(IRtypes.float_data, :last_result, [:fv_dofs_partition], [:k])
            ]),
            IR_loop_node(IRtypes.time_loop, :stage, :j, 1, :rki, IR_block_node([
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_data, :tmpresult, [:fv_dofs_partition], [:k]),
                    IR_operation_node(IRtypes.math_op, [:+,
                            IR_data_node(IRtypes.float_data, :tmpresult, [:fv_dofs_global], [:k]),
                        IR_operation_node(IRtypes.math_op, [:*, :dt, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, 
                                IR_data_node(IRtypes.float_data, :a, [:?], [IR_operation_node(IRtypes.math_op, [:+, :rki, 1]), :j])]),
                            IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], [:k, :j])
                        ])
                    ])
                ])
            ]))
        ]))
        update_ki_loop_two = IR_loop_node(IRtypes.space_loop, :dofs, :k, 1, :fv_dofs_partition, IR_block_node([
            # for k=1:fv_dofs_partition
            #     tmpki[k,stage] = global_vector[k];
            # end
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], [:k, :rki]),
                IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], [:k])
            ])
        ]))
        update_ki_condition = IR_conditional_node(IR_operation_node(IRtypes.math_op, [:<, :rki, stepper.stages]),
            IR_block_node([
                update_ki_loop_one,
                # copy_bdry_vals_to_vector(var, tmpresult, grid_data, dofs_per_node);
                IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :tmpresult]),
                # place_vector_in_vars(var, tmpresult, stepper);
                IR_operation_node(IRtypes.named_op, [:SCATTER_VARS, :tmpresult])
            ]),
            IR_block_node([update_ki_loop_two]));
        
        stage_loop_body = IR_block_node([
            # last_t = t;
            # t = last_t + stepper.c[rki]*stepper.dt;
            IR_operation_node(IRtypes.assign_op, [
                :t,
                IR_operation_node(IRtypes.math_op, [:+, :t, 
                    IR_operation_node(IRtypes.math_op, [:*, :dt, 
                        IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_data, :c, [:?], [:rki])])])])
            ]),
            # assemble
            zero_bdry_done,
            zero_flux_done,
            reset_nonzero_counter,
            ghost_exchange,
            pre_step_call,
            wrap_in_timer(:step_assembly, assembly),
            # update tmpki
            update_ki_condition,
            post_step_call
        ]);
        stage_loop = IR_loop_node(IRtypes.time_loop, :stages, :rki, 1, stepper.stages, stage_loop_body);
        
        # last_result[k] += stepper.dt * stepper.b[rki] * tmpki[k, rki];
        combine_loop_body = IR_block_node([
            IR_loop_node(IRtypes.time_loop, :stages, :rki, 1, stepper.stages, IR_block_node([
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_data, :last_result, [:fv_dofs_partition], [:k]),
                    IR_operation_node(IRtypes.math_op, [:+,
                        IR_data_node(IRtypes.float_data, :last_result, [:fv_dofs_partition], [:k]), 
                        IR_operation_node(IRtypes.math_op, [:*, :dt, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_data, :b, [:?], [:rki])]),
                            IR_data_node(IRtypes.float_data, :tmpki, [:fv_dofs_partition], [:k, :rki])
                        ])
                    ])
                ])
            ]))
        ])
        combine_loop = IR_loop_node(IRtypes.space_loop, :dofs, :k, 1, :fv_dofs_partition, combine_loop_body);
        
        push!(tloop_body.parts, IR_operation_node(IRtypes.named_op, [:GATHER_VARS, :last_result]));
        push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:last_t, :t]));
        push!(tloop_body.parts, stage_loop);
        push!(tloop_body.parts, combine_loop);
        # copy_bdry_vals_to_vector(var, last_result, grid_data, dofs_per_node);
        push!(tloop_body.parts, IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :last_result]))
        # place_vector_in_vars(var, last_result, stepper);
        push!(tloop_body.parts, IR_operation_node(IRtypes.named_op, [:SCATTER_VARS, :last_result]));
        # final post-step fun
        push!(tloop_body.parts, post_step_call);
        # update time
        push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:t, IR_operation_node(IRtypes.math_op, [:+, :last_t, :dt])]));
        push!(tloop_body.parts, progress_update);
        
        time_loop = IR_block_node([
            IR_operation_node(IRtypes.assign_op, [
                tmp_last,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :fv_dofs_partition])]),
            IR_operation_node(IRtypes.assign_op, [
                tmp_result,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :fv_dofs_partition])]),
            IR_operation_node(IRtypes.assign_op, [
                tmp_ki,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_data, :fv_dofs_partition, stepper.stages])]),
            progress_init,
            IR_loop_node(IRtypes.time_loop, :time, :ti, 1, stepper_nsteps, tloop_body)
        ]);
        
    elseif stepper.type == EULER_IMPLICIT || stepper.type == CRANK_NICHOLSON
        # # This part is used if the solution is updated like sol = sol + dt*(dsol)
        # sol_i = IR_data_node(IRtypes.float_data, :solution, [:fv_dofs_partition], [:update_i]);
        # dsol_i = IR_data_node(IRtypes.float_data, :global_vector, [:fv_dofs_partition], [:update_i]);
        # update_loop = IR_loop_node(IRtypes.dof_loop, :dof, :update_i, 1, :fv_dofs_partition, IR_block_node([
        #     IR_operation_node(IRtypes.assign_op, [
        #         sol_i,
        #         IR_operation_node(IRtypes.math_op, [:+, sol_i, IR_operation_node(IRtypes.math_op, [:*, time_stepper.dt, dsol_i])])
        #     ])
        # ]))
        # tloop_body.parts = [
        #     zero_bdry_done,
        #     zero_flux_done,
        #     wrap_in_timer(:step_assembly, assembly),
        #     # before updating the bdry vals may need to be put in vars and zeroed in solution
        #     wrap_in_timer(:update_sol, update_loop),
        #     IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :solution]),
        #     IR_operation_node(IRtypes.named_op, [:SCATTER_VARS, :solution]),
        #     IR_operation_node(IRtypes.assign_op, [
        #         :t,
        #         IR_operation_node(IRtypes.math_op, [:+, :t, stepper.dt])
        #     ])
        # ];
        
        # This part is used if the update is coded into the system sol = A\b (not sol = sol + dt*(A\b))
        tloop_body.parts = [
            zero_bdry_done,
            zero_flux_done,
            reset_nonzero_counter,
            ghost_exchange,
            pre_step_call,
            wrap_in_timer(:step_assembly, assembly),
            # gather_system,
            # # IR_operation_node(IRtypes.named_op, [:GLOBAL_FORM_MATRIX]),
            # wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :global_solution, :global_matrix, :global_vector])),
            # distribute_solution,
            gather_solve_distribute,
            IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :solution]),
            IR_operation_node(IRtypes.named_op, [:SCATTER_VARS, :solution]),
            post_step_call,
            IR_operation_node(IRtypes.assign_op, [
                :t,
                IR_operation_node(IRtypes.math_op, [:+, :t, :dt])
            ]),
            progress_update
        ];
        
        time_loop = IR_block_node([
            progress_init,
            IR_loop_node(IRtypes.time_loop, :time, :ti, 1, stepper_nsteps, tloop_body)
        ]);
    end
    
    return time_loop;
end

# A special function for handling LHS surface parts.
# Swaps variable entities with a number depending on side.
# negate will make side 2 negative.
# Returns a new copy of ex, rather than modifying ex
function replace_lhs_surface_var_entities(ex, var, side, negate=true)
    if typeof(ex) == Expr
        new_ex = copy(ex);
        for i=1:length(new_ex.args)
            new_ex.args[i] = replace_lhs_surface_var_entities(new_ex.args[i], var, side, negate);
        end
        return new_ex;
        
    elseif typeof(ex) == SymEntity
        if is_unknown_var(ex, var)
            which_side = get_face_side_info(ex);
            if side == 0
                ex = 1;
            elseif which_side == 0 || which_side == 3
                ex = 0.5;
            elseif which_side == side
                ex = 1;
            else
                ex = 0;
            end
            if side == 2 && negate
                ex = -ex;
            end
        end
    else
        # number?
    end
    return ex;
end