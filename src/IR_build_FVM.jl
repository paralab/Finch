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
function build_IR_fvm(lhs_vol, lhs_surf, rhs_vol, rhs_surf, var, indices, config, prob, time_stepper, fv_info)
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
    fv_order = fv_info.fluxOrder;
    
    # Options
    average_coefficients = false; # true to integrate and average coefficients over cells and faces.
    
    # A matrix only needs to be created if the stepper is implicit
    need_matrix = time_stepper.implicit;
    
    IRtypes = IR_entry_types();
    
    # These will hold the IR
    allocate_block = IR_block_node([]);
    vol_coef_block = IR_block_node([]);
    surf_coef_block = IR_block_node([]);
    
    source_block = IR_block_node([]);
    flux_block = IR_block_node([]);
    toglobal_block = IR_block_node([]);
    
    # Allocate the global matrix and vector
    if need_matrix
        push!(allocate_block.parts, IR_comment_node("Allocate global matrix(IJV form)"));
        # allocated size of nonzeros is dofs*dofs*(nel + nfaces*4)
        allocatedNZ = IR_data_node(IRtypes.int_64_data, :allocated_nonzeros, [], []);
        push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
            allocatedNZ,
            IR_operation_node(IRtypes.math_op, [:*, :dofs_per_node, :dofs_per_node,
                IR_operation_node(IRtypes.math_op, [:+, :num_elements, IR_operation_node(IRtypes.math_op, [:*, :num_faces,4])])
            ])
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
        # I and J should be initialized to 1 for julia
        # How about c++ ?? do we need a named op for init of IJV?
        push!(allocate_block.parts, IR_comment_node("I and J vectors should be ones"));
        push!(allocate_block.parts, IR_operation_node(IRtypes.named_op, [
            :FILL_ARRAY, globalmat_I, 1, :allocated_nonzeros]));
        push!(allocate_block.parts, IR_operation_node(IRtypes.named_op, [
            :FILL_ARRAY, globalmat_J, 1, allocated_nonzeros]));
        #
        # I and J indices only need to be set once, so init them here
        
    end
    
    # Global vectors
    push!(allocate_block.parts, IR_comment_node("Allocate global vectors."));
    globalvec = IR_data_node(IRtypes.float_64_data, :global_vector, [:fv_dofs_global], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        globalvec,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :fv_dofs_global])
        ]));
    
    solvec = IR_data_node(IRtypes.float_64_data, :solution, [:fv_dofs_global], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        solvec,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :fv_dofs_global])
        ]));
        
    # Allocate the elemental matrix and vector
    if need_matrix
        # push!(allocate_block.parts, IR_comment_node("Allocate elemental matrix."));
        # elementmat = IR_data_node(IRtypes.float_64_data, :element_matrix, [:dofs_per_node, :dofs_per_node], []);
        # push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        #     elementmat,
        #     IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_per_node, :dofs_per_node])
        #     ]));
    end
    
    # elemental flux and source
    push!(allocate_block.parts, IR_comment_node("Allocate elemental source and flux."));
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        IR_data_node(IRtypes.float_64_data, :source, [:dofs_per_loop], []),
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_per_loop])
        ]));
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        IR_data_node(IRtypes.float_64_data, :flux, [:dofs_per_loop], []),
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_per_loop])
        ]));
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        IR_data_node(IRtypes.float_64_data, :flux_tmp, [:dofs_per_loop], []),
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_per_loop])
        ]));
    
    
    # bdry done flag for each face
    push!(allocate_block.parts, IR_comment_node("Boundary done flag for each face."));
    bdry_done = IR_data_node(IRtypes.int_64_data, :bdry_done, [:num_faces], []);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        bdry_done,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.int_64_data, :num_faces])
        ]));
    
    # coefficient prep
    # a list of all entities and rhs only ones
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
    
    # Prepare coefficient and other elemental values
    (vol_coef, surf_coef) = prepare_coefficient_values(all_entities, var, dimension, counts, fv_info);
    # volume coefficients (for source)
    push!(vol_coef_block.parts, IR_comment_node("Evaluate volume coefficients."));
    append!(vol_coef_block.parts, vol_coef.parts);
    push!(vol_coef_block.parts, IR_operation_node(IRtypes.assign_op, [:volume, 
        IR_operation_node(IRtypes.member_op, [:geometric_factors, 
            IR_data_node(IRtypes.float_64_data, :volume, [:num_elements], [:eid])])]))
    
    # surface coefficients (for flux)
    push!(surf_coef_block.parts, IR_comment_node("Evaluate surface coefficients."));
    append!(surf_coef_block.parts, surf_coef.parts);
    push!(surf_coef_block.parts, IR_operation_node(IRtypes.assign_op, [:area, 
        IR_operation_node(IRtypes.member_op, [:geometric_factors, 
            IR_data_node(IRtypes.float_64_data, :area, [:num_faces], [:fid])])]))
    push!(surf_coef_block.parts, IR_operation_node(IRtypes.assign_op, [:area_over_volume, 
        IR_operation_node(IRtypes.math_op, [:(/), :area, :volume])]))
    
    # source block
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
    
    # flux block (includes face loop)
    push!(flux_block.parts, IR_comment_node("Compute flux terms (surface integral) in a face loop"));
    fid = IR_data_node(IRtypes.int_64_data, :fid, [], []); # Face ID
    fbid = IR_data_node(IRtypes.int_64_data, :fbid, [], []); # Face BID
    leftel = IR_data_node(IRtypes.int_64_data, :left_el, [], []); # left element
    rightel = IR_data_node(IRtypes.int_64_data, :right_el, [], []); # right element
    face_loop_body = IR_block_node([
        # fid = fv_grid.element2face[i, eid]; # face ID 
        IR_operation_node(IRtypes.assign_op, [fid,
            IR_operation_node(IRtypes.member_op, [:mesh, 
                IR_data_node(IRtypes.int_64_data, :element2face, [:faces_per_element, :num_elements], [:fi, :eid])])
        ]),
        # fbid = fv_grid.facebid[fid]; # BID of this face
        IR_operation_node(IRtypes.assign_op, [fbid,
            IR_operation_node(IRtypes.member_op, [:mesh, 
                IR_data_node(IRtypes.int_64_data, :facebid, [:num_faces], [:fid])])
        ]),
        # (leftel, rightel) = fv_grid.face2element[:,fid];
        IR_operation_node(IRtypes.assign_op, [leftel,
            IR_operation_node(IRtypes.member_op, [:mesh, 
                IR_data_node(IRtypes.int_64_data, :face2element, [2, :num_faces], [1, :fid])])
        ]),
        IR_operation_node(IRtypes.assign_op, [rightel,
            IR_operation_node(IRtypes.member_op, [:mesh, 
                IR_data_node(IRtypes.int_64_data, :face2element, [2, :num_faces], [2, :fid])])
        ]),
        # if (eid == rightel || rightel == 0) (neighbor = leftel) else (neighbor = rightel)
        IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(||), 
                IR_operation_node(IRtypes.math_op, [:(==), :eid, rightel]),
                IR_operation_node(IRtypes.math_op, [:(==), rightel, 0])]),
            IR_block_node([IR_operation_node(IRtypes.assign_op, [:neighbor, leftel])]),
            IR_block_node([IR_operation_node(IRtypes.assign_op, [:neighbor, rightel])])
        )
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
    
    push!(face_loop_body.parts, IR_comment_node("Apply boundary conditions"));
    if need_matrix && !(lhs_surf === nothing)  && !is_empty_expression(lhs_surf)
        push!(face_loop_body.parts, IR_operation_node(IRtypes.function_op, [
            :apply_boundary_conditions_face, 
            :var, :eid, :fid, :fbid, :mesh, :refel, :geometric_factors, :prob, :t, :flux_matrix, :flux_tmp, :bdry_done, :index_offset]));
            # add flux_tmp to flux after handling BCs
            # TODO
    else
        push!(face_loop_body.parts, IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(>), :fbid, 0]),
            IR_block_node([IR_operation_node(IRtypes.function_op, [
                :apply_boundary_conditions_face_rhs, 
                :var, :eid, :fid, :fbid, :mesh, :refel, :geometric_factors, :prob, :t, :flux_tmp, :bdry_done, :index_offset])
            ])));
        # add flux_tmp to flux after handling BCs
        if dofsper == 1
            push!(face_loop_body.parts, IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_64_data, :flux, [:dofs_per_loop], [1]), 
                IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.float_64_data, :flux, [:dofs_per_loop], [1]), 
                    IR_operation_node(IRtypes.math_op, [:*, IR_data_node(IRtypes.float_64_data, :flux_tmp, [:dofs_per_loop], [1]), :area_over_volume])])]));
        else
            push!(face_loop_body.parts, IR_loop_node(IRtypes.dof_loop, :dofs, :dofi, 1, :dofs_per_loop, IR_block_node([
                IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_64_data, :flux, [:dofs_per_loop], [:dofi]), 
                IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.float_64_data, :flux, [:dofs_per_loop], [:dofi]), 
                    IR_operation_node(IRtypes.math_op, [:*, IR_data_node(IRtypes.float_64_data, :flux_tmp, [:dofs_per_loop], [:dofi]), :area_over_volume])])])
            ])));
        end
    end
    
    push!(flux_block.parts, IR_loop_node(IRtypes.space_loop, :faces, :fi, 1, :faces_per_element, face_loop_body));
    
    # add to global sysem
    push!(toglobal_block.parts, IR_comment_node("Place elemental parts in global system."));
    push!(toglobal_block.parts, generate_local_to_global_fvm(dofsper, offset_ind, vec_only=!need_matrix));
    
    # assembly loop
    assembly_loop = generate_assembly_loop_fvm(var, indices);
    
    # find the innermost assembly loop
    inner_loop = assembly_loop;
    while length(inner_loop.body.parts) > 0 && typeof(inner_loop.body.parts[1]) == IR_loop_node
        inner_loop = inner_loop.body.parts[1];
    end
    
    # zero elemental source and flux
    if dofsper == 1
        zero_source = IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_64_data, :source, [:dofs_per_loop], [1]), 0.0]);
        zero_flux = IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_64_data, :flux, [:dofs_per_loop], [1]), 0.0]);
        zero_flux_t = IR_operation_node(IRtypes.assign_op, [IR_data_node(IRtypes.float_64_data, :flux_tmp, [:dofs_per_loop], [1]), 0.0]);
    else
        zero_source = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY, :source, 0.0, :dofs_per_loop]);
        zero_flux = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY, :flux, 0.0, :dofs_per_loop]);
        zero_flux_t = IR_operation_node(IRtypes.named_op, [:FILL_ARRAY, :flux_tmp, 0.0, :dofs_per_loop]);
    end
    # fill the elemental loop
    append!(inner_loop.body.parts, [
        zero_source,
        zero_flux,
        zero_flux_t,
        source_block,
        flux_block,
        toglobal_block
    ])
    
    # time loop
    need_time_loop = prob.time_dependent; # This should be true for FV
    if need_time_loop
        if need_matrix # implicit steppers
            # TODO
        else # explicit steppers (no matrix)
            step_loop = generate_time_stepping_loop_fvm(time_stepper, assembly_loop);
            time_loop = IR_block_node([
                IR_operation_node(IRtypes.named_op, [:GATHER_SOLUTION, solvec]),
                IR_operation_node(IRtypes.assign_op, [:t, 0]),
                IR_operation_node(IRtypes.assign_op, [:dt, time_stepper.dt]),
                IR_comment_node("###############################################"),
                IR_comment_node("Time stepping loop"),
                wrap_in_timer(:time_steps, step_loop)
            ]);
        end
        
    else # no time stepping
        # Do we need this?
    end
    
    # Put them all together in a master block
    master_block = IR_block_node([
        wrap_in_timer(:allocate, allocate_block),
        time_loop
    ]);
    
    return master_block;
end

# Compute or fetch all needed values
function prepare_coefficient_values(entities, var, dimension, counts, fv_info)
    IRtypes = IR_entry_types();
    row_col_matrix_index = IR_operation_node(IRtypes.named_op, [:ROWCOL_TO_INDEX, :row, :col, :nodes_per_element]);
    col_row_matrix_index = IR_operation_node(IRtypes.named_op, [:ROWCOL_TO_INDEX, :col, :row, :nodes_per_element]);
    
    fv_order = fv_info.fluxOrder;
    
    # These parts will be returned
    vol_block = IR_block_node([]);
    surf_block = IR_block_node([]);
    
    unique_entity_names = []; # avoid duplicate names
    
    # First get x,y,z for cell centers
    cell_coords = IR_block_node([]);
    need_cell_coords = false;
    push!(cell_coords.parts, IR_operation_node(IRtypes.assign_op, [:x, 
        IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :cellCenters, [dimension, :num_elements], [1,:eid])])]));
    if dimension > 1
        push!(cell_coords.parts, IR_operation_node(IRtypes.assign_op, [:y, 
            IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :cellCenters, [dimension, :num_elements], [2,:eid])])]));
    else
        push!(cell_coords.parts, IR_operation_node(IRtypes.assign_op, [:y, 0]));
    end
    if dimension > 2
        push!(cell_coords.parts, IR_operation_node(IRtypes.assign_op, [:z, 
            IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :cellCenters, [dimension, :num_elements], [3,:eid])])]));
    else
        push!(cell_coords.parts, IR_operation_node(IRtypes.assign_op, [:z, 0]));
    end
    
    # And the same for face centers
    face_coords = IR_block_node([]);
    need_face_coords = false;
    push!(face_coords.parts, IR_operation_node(IRtypes.assign_op, [:x, 
        IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :faceCenters, [dimension, :num_faces], [1,:fid])])]));
    if dimension > 1
        push!(face_coords.parts, IR_operation_node(IRtypes.assign_op, [:y, 
            IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :faceCenters, [dimension, :num_faces], [2,:fid])])]));
    else
        push!(face_coords.parts, IR_operation_node(IRtypes.assign_op, [:y, 0]));
    end
    if dimension > 2
        push!(face_coords.parts, IR_operation_node(IRtypes.assign_op, [:z, 
            IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :faceCenters, [dimension, :num_faces], [3,:fid])])]));
    else
        push!(face_coords.parts, IR_operation_node(IRtypes.assign_op, [:z, 0]));
    end
    
    # Are face normals needed?
    need_normals = false;
    normal_part = IR_block_node([]);
    # FACENORMAL1[1] -> value__FACENORMAL1_1 = mesh.facenormals[1,fid];
    push!(normal_part.parts, IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_1, 
        IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.float_64_data, :facenormals, [dimension, :num_faces], [1,:fid])])]));
    if dimension > 1
        push!(normal_part.parts, IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_2, 
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.float_64_data, :facenormals, [dimension, :num_faces], [2,:fid])])]));
    else
        # push!(normal_part.parts, IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_2, 0.0]));
    end
    if dimension > 2
        push!(normal_part.parts, IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_3, 
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.float_64_data, :facenormals, [dimension, :num_faces], [3,:fid])])]));
    else
        # push!(normal_part.parts, IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_3, 0.0]));
    end
    
    if dimension == 1
        push!(normal_part.parts, IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :eid, :left_el]),
            # normal is correct, make FACENORMAL2 = -FACENORMAL1
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_1, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_1])])]),
            # normal is reversed
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_1, :value__FACENORMAL1_1]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_1, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_1])])])
        ))
    elseif dimension == 2
        push!(normal_part.parts, IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :eid, :left_el]),
            # normal is correct, make FACENORMAL2 = -FACENORMAL1
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_1, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_1])]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_2, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_2])])]),
            # normal is reversed
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_1, :value__FACENORMAL1_1]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_2, :value__FACENORMAL1_2]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_1, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_1])]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_2, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_2])])])
        ))
    elseif dimension == 3
        push!(normal_part.parts, IR_conditional_node(IR_operation_node(IRtypes.math_op, [:(==), :eid, :left_el]),
            # normal is correct, make FACENORMAL2 = -FACENORMAL1
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_1, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_1])]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_2, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_2])]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_3, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_3])])]),
            # normal is reversed
            IR_block_node([
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_1, :value__FACENORMAL1_1]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_2, :value__FACENORMAL1_2]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL2_3, :value__FACENORMAL1_3]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_1, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_1])]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_2, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_2])]),
                IR_operation_node(IRtypes.assign_op, [:value__FACENORMAL1_3, IR_operation_node(IRtypes.math_op, [:-, :value__FACENORMAL1_3])])])
        ))
    end
    
    # Is the distance between cell centers needed for derivatives?
    # normal scaled by distance between cell centers"
    # cell_dist = sqrt( (xn-xe)^2 + ()^2 + ()^2 )
    # cell_dx = cell_dist * value__FACENORMAL1_1
    need_deriv_dist = false;
    # get cell_xe, cell_xn, ...
    cell_coords_2 = IR_block_node([
            IR_operation_node(IRtypes.assign_op, [:cell_xe, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :cellCenters, [dimension, :num_elements], [1,:eid])])
            ]),
            IR_operation_node(IRtypes.assign_op, [:cell_xn, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :cellCenters, [dimension, :num_elements], [1,:neighbor])])
            ])
        ]);
    if dimension > 1
        append!(cell_coords_2.parts, [
            IR_operation_node(IRtypes.assign_op, [:cell_ye, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :cellCenters, [dimension, :num_elements], [2,:eid])])
            ]),
            IR_operation_node(IRtypes.assign_op, [:cell_yn, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :cellCenters, [dimension, :num_elements], [2,:neighbor])])
            ])
        ])
    end
    if dimension > 2
        append!(cell_coords_2.parts, [
            IR_operation_node(IRtypes.assign_op, [:cell_ze, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :cellCenters, [dimension, :num_elements], [3,:eid])])
            ]),
            IR_operation_node(IRtypes.assign_op, [:cell_zn, 
                IR_operation_node(IRtypes.member_op, [:fv_info, IR_data_node(IRtypes.float_64_data, :cellCenters, [dimension, :num_elements], [3,:neighbor])])
            ])
        ])
    end
    # find the square distance (xn-xe)^2 + ...
    square_dist = IR_operation_node(IRtypes.math_op, [:*,
        IR_operation_node(IRtypes.math_op, [:-,
            IR_data_node(IRtypes.float_64_data, :cell_xn, [], []),
            IR_data_node(IRtypes.float_64_data, :cell_xe, [], [])
        ]),
        IR_operation_node(IRtypes.math_op, [:-,
            IR_data_node(IRtypes.float_64_data, :cell_xn, [], []),
            IR_data_node(IRtypes.float_64_data, :cell_xe, [], [])
        ])
    ])
    if dimension > 1
        square_dist = IR_operation_node(IRtypes.math_op, [:+, square_dist, 
            IR_operation_node(IRtypes.math_op, [:*,
                IR_operation_node(IRtypes.math_op, [:-,
                    IR_data_node(IRtypes.float_64_data, :cell_yn, [], []),
                    IR_data_node(IRtypes.float_64_data, :cell_ye, [], [])
                ]),
                IR_operation_node(IRtypes.math_op, [:-,
                    IR_data_node(IRtypes.float_64_data, :cell_yn, [], []),
                    IR_data_node(IRtypes.float_64_data, :cell_ye, [], [])
                ])
            ])
        ])
    end
    if dimension > 2
        push!(square_dist.args, 
            IR_operation_node(IRtypes.math_op, [:*,
                IR_operation_node(IRtypes.math_op, [:-,
                    IR_data_node(IRtypes.float_64_data, :cell_zn, [], []),
                    IR_data_node(IRtypes.float_64_data, :cell_ze, [], [])
                ]),
                IR_operation_node(IRtypes.math_op, [:-,
                    IR_data_node(IRtypes.float_64_data, :cell_zn, [], []),
                    IR_data_node(IRtypes.float_64_data, :cell_ze, [], [])
                ])
            ])
        )
    end
    # Put it all together
    deriv_dist = IR_block_node([IR_comment_node("Normal scaled by distance between cell centers"),
        cell_coords_2,
        IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.float_64_data, :cell_dist, [], []),
            IR_operation_node(IRtypes.math_op, [:sqrt, square_dist])
        ]),
        IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.float_64_data, :dxyz_1, [], []),
            IR_operation_node(IRtypes.math_op, [:*, :cell_dist, :value__FACENORMAL1_1])
        ])
    ]);
    if dimension > 1
        push!(deriv_dist.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :dxyz_2, [], []),
                IR_operation_node(IRtypes.math_op, [:*, :cell_dist, :value__FACENORMAL1_2])
            ]))
    end
    if dimension > 2
        push!(deriv_dist.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :dxyz_3, [], []),
                IR_operation_node(IRtypes.math_op, [:*, :cell_dist, :value__FACENORMAL1_3])
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
                        IR_data_node(IRtypes.float_64_data, Symbol(cname), [], []), cval]))
                else
                    push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_64_data, Symbol(cname), [], []), cval]))
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
                    need_cell_coords = true;
                    coef_index = get_coef_index(entities[i]);
                    
                    push!(vol_block.parts, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_64_data, Symbol(cname), [], []),
                        IR_operation_node(IRtypes.named_op, [:COEF_EVAL, coef_index, index_IR, :x, :y, :z, :t, :eid])]));
                    # If derivatives are needed, I'm not sure what to do here
                    
                else # surface
                    need_face_coords = true;
                    coef_index = get_coef_index(entities[i]);
                    
                    push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.float_64_data, Symbol(cname), [], []),
                        IR_operation_node(IRtypes.named_op, [:COEF_EVAL, coef_index, index_IR, :x, :y, :z, :t, :eid])]));
                    # If derivatives are needed, I'm not sure what to do here
                end
                
            elseif ctype == 3 # a known variable value
                cellside = get_face_side_info(entities[i]); # in code_generator_utils.jl
                if vors == "surface"
                    l2gsymbol = :els1
                else
                    l2gsymbol = :el
                end
                if cellside == 1
                    l2gsymbol = :els1
                elseif cellside == 2
                    l2gsymbol = :els2
                end
                
                # Make an index string for indexed variables
                if typeof(entities[i].index) <: Array
                    # It is an indexed variable
                    if length(entities[i].index) == 1
                        indsymbol = symbol("INDEX_VAL_"*entities[i].index[1]);
                    else
                        # There is more than one index. Need to form an expression for it.
                        indsymbol = "(INDEX_VAL_"*entities[i].index[1];
                        indices = variables[cval].indexer;
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
                    if variables[cval].discretization == FV
                        # variables[cval].values[index,eid]
                        push!(vol_block.parts, IR_operation_node(IRtypes.assign_op,[
                            IR_data_node(IRtypes.float_64_data, Symbol(cname), [], []),
                            IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :eid])
                        ]));
                    else
                        # sum(variables[cval].values[index,mesh.loc2glb[:,eid]]) / nodes_per_element
                        # TODO
                    end
                    
                else # surface
                    if variables[cval].discretization == FV
                        # Need to reconstruct it using neighboring cells
                        if fv_order == 1
                            if length(entities[i].derivs) > 0
                                need_normals = true;
                                need_deriv_dist = true;
                                # cname = variables[cval].values[index, els2] - variables[cval].values[index, els1]
                                # if (els1 != els2 && abs(normal[entities[i].derivs[1]]) > 1e-10)
                                #     cname = cname / dxyz[entities[i].derivs[1]]
                                # else
                                #     cname = 0
                                # end
                                if entities[i].derivs[1] == 1
                                    normal_n = :value__FACENORMAL1_1;
                                elseif entities[i].derivs[1] == 2
                                    normal_n = :value__FACENORMAL1_2;
                                else
                                    normal_n = :value__FACENORMAL1_3;
                                end
                                push!(surf_block.parts, IR_conditional_node(
                                IR_operation_node(IRtypes.math_op, [:(&&), 
                                    IR_operation_node(IRtypes.math_op, [:(!=), :eid, :neighbor]),
                                    IR_operation_node(IRtypes.math_op, [:(>), 
                                        IR_operation_node(IRtypes.function_op, [:abs, normal_n]),
                                        1e-10])
                                ]), IR_block_node([
                                
                                IR_operation_node(IRtypes.assign_op,[
                                    IR_data_node(IRtypes.float_64_data, Symbol(cname), [], []),
                                    IR_operation_node(IRtypes.math_op, [:/,
                                        IR_operation_node(IRtypes.math_op, [:-,
                                            IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :neighbor]),
                                            IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :eid])
                                        ]),
                                        IR_data_node(IRtypes.float_64_data, Symbol("dxyz_"*string(entities[i].derivs[1])), [], [])
                                    ])
                                ])
                                ]), IR_block_node([
                                    
                                IR_operation_node(IRtypes.assign_op,[
                                    IR_data_node(IRtypes.float_64_data, Symbol(cname), [], []), 0.0])
                                ])));
                                
                            else
                                if cellside == 0 || cellside == 3
                                    # No side was specified or central approx, so use the average
                                    # cname = 0.5 * (variables[cval].values[index,eid] + variables[cval].values[index,neighbor])
                                    push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                                        IR_data_node(IRtypes.float_64_data, Symbol(cname), [], []),
                                        IR_operation_node(IRtypes.math_op, [:*, 0.5,
                                            IR_operation_node(IRtypes.math_op, [:+,
                                                IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :neighbor]),
                                                IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :eid])
                                            ])
                                        ])
                                    ]))
                                elseif cellside == 1
                                    # cname = variables[cval].values[index, whichone]
                                    push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                                        IR_data_node(IRtypes.float_64_data, Symbol(cname)),
                                        IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :eid])
                                    ]))
                                elseif cellside == 2
                                    # cname = variables[cval].values[index, whichone]
                                    push!(surf_block.parts, IR_operation_node(IRtypes.assign_op,[
                                        IR_data_node(IRtypes.float_64_data, Symbol(cname)),
                                        IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, indsymbol, :neighbor])
                                    ]))
                                # elseif cellside == 3 # central
                                    # same as 0 for this case
                                    # cname = 0.5 * (variables[cval].values[index,eid] + variables[cval].values[index,neighbor])
                                elseif cellside == 4 # neighborhood
                                    # This is a special case only for callback functions.
                                    # Rather than representing a value, it constructs a Neighborhood object to be passed.
                                    # code *= cname * " = Finch.Neighborhood(els, cellx, [Finch.variables["*string(cval)*"].values["*indstr*", els[1]], Finch.variables["*string(cval)*"].values["*indstr*", els[2]]]);\n";
                                    # TODO
                                end
                            end
                            
                        else # higher order
                            
                        end
                        
                    else # defined at nodes
                        
                    end
                end
            end
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
    
    # This will be returned
    compute_block = IR_block_node([]);
    
    # source or flux being computed?
    if vors == "surface"
        if dofsper > 1
            if lorr == LHS
                result_size = [:dofs_per_loop, :?];
            else # RHS
                result_size = [:dofs_per_loop];
            end
            
        else # one dof
            if lorr == LHS
                result_size = [:dofs_per_loop, :faces_per_element];
            else # RHS
                result_size = [:dofs_per_loop];
            end
        end
        result_part = IR_data_node(IRtypes.float_64_data, :flux_tmp, result_size, [])
        
    else # volume
        if dofsper > 1
            if lorr == LHS
                result_size = [:dofs_per_loop, :dofs_per_loop];
            else # RHS
                result_size = [:dofs_per_loop];
            end
            
        else # one dof
            result_size = [:dofs_per_loop];
        end
        result_part = IR_data_node(IRtypes.float_64_data, :source, result_size, [])
    end
    
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
        
        for vi=1:length(var) # variables
            # Process the terms for this variable
            for ci=1:length(terms[vi]) # components
                for i=1:length(terms[vi][ci])
                    (term_IR, what, var_ind) = generate_term_calculation_fvm(terms[vi][ci][i], var, lorr, vors);
                    
                    # Find the appropriate subvector for this term
                    subveci = offset_ind[vi] + ci;
                    subvecj = var_ind;
                    if lorr == LHS
                        subvec_ind = subveci + dofsper * (subvecj-1);
                    else
                        subvec_ind = subveci;
                    end
                    
                    push!(submatrices[subvec_ind], term_IR);
                end
            end
            
        end # vi
        
        # if typeof(var) <: Array
        #     for vi=1:length(var) # variables
        #         # Process the terms for this variable
        #         for ci=1:length(terms[vi]) # components
        #             for i=1:length(terms[vi][ci])
        #                 (term_IR, what, var_ind) = generate_term_calculation_fvm(terms[vi][ci][i], var, lorr, vors);
                        
        #                 # Find the appropriate subvector for this term
        #                 subveci = offset_ind[vi] + ci;
        #                 subvecj = var_ind;
        #                 if lorr == LHS
        #                     subvec_ind = subveci + dofsper * (subvecj-1);
        #                 else
        #                     subvec_ind = subveci;
        #                 end
                        
        #                 push!(submatrices[subvec_ind], term_IR);
        #             end
        #         end
                
        #     end # vi
            
        # else # only one variable
        #     # Process the terms for this variable
        #     for ci=1:length(terms) # components
        #         for i=1:length(terms[ci])
        #             (term_IR, what, var_ind) = generate_term_calculation_fvm(terms[ci][i], var, lorr, vors);
                        
        #             # Find the appropriate subvector for this term
        #             subveci = ci;
        #             subvecj = var_ind;
        #             if lorr == LHS
        #                 subvec_ind = subveci + dofsper * (subvecj-1);
        #             else
        #                 subvec_ind = subveci;
        #             end
                    
        #             push!(submatrices[subvec_ind], term_IR);
        #         end
        #     end
        # end
        
        # Put the submatrices together
        num_nonzero_blocks = 0;
        
        if lorr == LHS
            linalg_matrix_block_args = [];
            push!(linalg_matrix_block_args, :LINALG_MAT_BLOCKS);
            push!(linalg_matrix_block_args, 0);
            push!(linalg_matrix_block_args, :nodes_per_element);
            push!(linalg_matrix_block_args, result_part);
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
        terms = terms[1][1];
        term_vec = Vector{IR_part}(undef,0);
        
        #process each term
        for i=1:length(terms)
            (term_IR, what, var_ind) = generate_term_calculation_fvm(terms[i], var, lorr, vors);
            
            push!(term_vec, term_IR);
        end
        
        if lorr == LHS
            linalg_matrix_block_args = [];
            push!(linalg_matrix_block_args, :LINALG_MAT_BLOCKS);
            push!(linalg_matrix_block_args, 1);
            push!(linalg_matrix_block_args, :nodes_per_element);
            push!(linalg_matrix_block_args, result_part);
            
            # push!(compute_block.parts, IR_operation_node(IRtypes.named_op, linalg_matrix_block_args));
            
        else # RHS
            linalg_vector_block_args = [];
            push!(linalg_vector_block_args, :LINALG_VEC_BLOCKS);
            push!(linalg_vector_block_args, 1);
            push!(linalg_vector_block_args, 1); # number of entries per dof is 1 for FVM
            push!(linalg_vector_block_args, result_part);
            
            # push!(compute_block.parts, IR_operation_node(IRtypes.named_op, linalg_vector_block_args));
        end
        
        if length(term_vec) > 1
            new_term_vec = [];
            push!(new_term_vec, :+);
            append!(new_term_vec, term_vec);
            result_rhs = IR_operation_node(IRtypes.math_op, new_term_vec);
            
        else
            result_rhs = term_vec[1];
        end
        push!(linalg_vector_block_args, 1);
        push!(linalg_vector_block_args, result_rhs);
        
        if lorr == LHS
            push!(compute_block.parts, IR_operation_node(IRtypes.named_op, linalg_matrix_block_args));
            
        else # RHS
            # push!(compute_block.parts, IR_operation_node(IRtypes.named_op, linalg_vector_block_args));
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
        
        if coef_part === nothing # Only an unknown variable
            # handle negative sign
            if typeof(var_part) == Expr && (var_part.args[1] === :.- || var_part.args[1] === :-)
                if vors == "surface"
                    if which_side == 0 || which_side == 3 # average them
                        result1 = -0.5;
                        result2 = -0.5;
                    elseif which_side == 1
                        result1 = -1;
                    elseif which_side == 2
                        result2 = -1;
                    else
                        # neighborhoods are not supported for implicit steppers yet
                    end
                else
                    result1 = -1;
                end
                
            else # no negative sign
                if vors == "surface"
                    if which_side == 0 || which_side == 3 # average them
                        result1 = 0.5;
                        result2 = 0.5;
                    elseif which_side == 1
                        result1 = 1;
                    elseif which_side == 2
                        result2 = 1;
                    else
                        # neighborhoods are not supported for implicit steppers yet
                    end
                else
                    result1 = 1;
                end
            end
            
        else # a coefficient part exists
            if vors == "surface"
                # This is tricky. First make two copies of the coef_part.
                # Since the coef_part could actually contain a variable hidden in an expression like a conditional,
                # Set the correct side to 1 and the other side to 0, or if no side, set both to 0.5
                # Also make side 2 negative... because I need to make upwinding work.
                # This is a bad way to rig this up, but I'm out of ideas here.
                coef_side1 = copy(coef_part);
                coef_side2 = copy(coef_part);
                coef_side1 = replace_lhs_surface_var_entities(coef_side1, var, 1);
                coef_side2 = replace_lhs_surface_var_entities(coef_side2, var, 2);
                if which_side == 0 || which_side == 3 # average them
                    result1 = arithmetic_expr_to_IR(coef_side1);
                    result2 = arithmetic_expr_to_IR(coef_side2);
                elseif which_side == 1
                    result1 = arithmetic_expr_to_IR(coef_side1);
                elseif which_side == 2
                    result2 = arithmetic_expr_to_IR(coef_side2);
                else
                    # neighborhoods are not supported for implicit steppers yet
                end
            else
                result1 = arithmetic_expr_to_IR(coef_part);
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
    if dofs_per_node == 1
        if vec_only
            # global_vector[eid] = source + flux;
            result_block = IR_block_node([
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_64_data, :global_vector, [:fv_dofs_global], [:eid]),
                    IR_operation_node(IRtypes.math_op, [:+,
                        IR_data_node(IRtypes.float_64_data, :source, [:dofs_per_loop], [1]),
                        IR_data_node(IRtypes.float_64_data, :flux, [:dofs_per_loop], [1])
                    ])
                ])
            ])
            
        else # vec and matrix
            # TODO
        end
        
    else # more than one dof
        if vec_only
            # for i=1:dofs_per_node
            #     dofid = i + dofs_per_node * (eid-1)
            #     solution[dofid] = sourcevec[i] + fluxvec[i];
            # end
            result_block = IR_block_node([
                IR_loop_node(IRtypes.space_loop, :dofs, :i, 1, :dofs_per_loop, IR_block_node([
                    IR_operation_node(IRtypes.assign_op, [
                        IR_data_node(IRtypes.int_64_data, :dofid, [], []),
                        IR_operation_node(IRtypes.math_op, [:+,
                            :index_offset,
                            :i,
                            IR_operation_node(IRtypes.math_op, [:*,
                                dofs_per_node,
                                IR_operation_node(IRtypes.math_op, [:-,
                                    IR_data_node(IRtypes.int_64_data, :eid, [], []),
                                    1
                                ])
                            ])
                        ])
                    ]),
                    IR_operation_node(IRtypes.assign_op, [
                        IR_data_node(IRtypes.float_64_data, :global_vector, [:fv_dofs_global], [:dofid]),
                        IR_operation_node(IRtypes.math_op, [:+,
                            IR_data_node(IRtypes.float_64_data, :source, [:dofs_per_element], [:i]),
                            IR_data_node(IRtypes.float_64_data, :flux, [:dpfs_per_element], [:i])
                        ])
                    ])
                ]))
            ])
            
        else # vec and matrix
            
        end
        
    end
    
    return result_block;
end

# Generate the assembly loop structures and insert the content
function generate_assembly_loop_fvm(var, indices)
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
function generate_time_stepping_loop_fvm(stepper, assembly)
    IRtypes = IR_entry_types();
    tloop_body = IR_block_node([]);
    zero_bdry_done = IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY,
        IR_data_node(IRtypes.int_64_data, :bdry_done, [:num_faces], []),
        0, :num_faces]);
    if stepper.type == EULER_EXPLICIT
        # This part is used if the solution is updated like sol = sol + dt*(dsol)
        sol_i = IR_data_node(IRtypes.float_64_data, :solution, [:fv_dofs_global], [:update_i]);
        dsol_i = IR_data_node(IRtypes.float_64_data, :global_vector, [:fv_dofs_global], [:update_i]);
        update_loop = IR_loop_node(IRtypes.dof_loop, :dof, :update_i, 1, :fv_dofs_global, IR_block_node([
            IR_operation_node(IRtypes.assign_op, [
                sol_i,
                IR_operation_node(IRtypes.math_op, [:+, sol_i, IR_operation_node(IRtypes.math_op, [:*, time_stepper.dt, dsol_i])])
            ])
        ]))
        tloop_body.parts = [
            zero_bdry_done,
            wrap_in_timer(:step_assembly, assembly),
            # before updating the bdry vals may need to be put in vars and zeroed in solution
            wrap_in_timer(:update_sol, update_loop),
            IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :solution]),
            wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :solution])),
            IR_operation_node(IRtypes.assign_op, [
                :t,
                IR_operation_node(IRtypes.math_op, [:+, :t, stepper.dt])
            ])
        ];
        
        # # This part is used if the update is coded into the system sol = A\b (not sol = sol + dt*(A\b))
        # tloop_body.parts = [
        #     zero_el_vec,
        #     zero_bdry_done,
        #     wrap_in_timer(:step_assembly, assembly),
        #     wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector])),
        #     wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :solution])),
        #     IR_operation_node(IRtypes.assign_op, [
        #         :t,
        #         IR_operation_node(IRtypes.math_op, [:+, :t, stepper.dt])
        #     ])
        # ];
        
        time_loop = IR_loop_node(IRtypes.time_loop, :time, :ti, 1, stepper.Nsteps, tloop_body);
       
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
                for rki=1:stepper.stages
                    rktime = t + stepper.c[rki]*stepper.dt;
                    # p(i-1) is currently in u
                    
                    b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, rktime, stepper.dt; rhs_only = true, assemble_loops=assemble_func);
                    
                    sol = linear_system_solve(A,b);
                    
                    # At this point sol holds the boundary values
                    # directly write them to the variable values and zero sol.
                    copy_bdry_vals_to_variables(var, sol, grid_data, dofs_per_node, zero_vals=true);
                    
                    if rki == 1 # because a1 == 0
                        tmpki .= stepper.dt .* sol;
                    else
                        tmpki .= stepper.a[rki].*tmpki + stepper.dt.*sol;
                    end
                    tmppi .= tmppi + stepper.b[rki].*tmpki
                    
                    copy_bdry_vals_to_vector(var, tmppi, grid_data, dofs_per_node);
                    place_sol_in_vars(var, tmppi, stepper);
                end
        =#
        tmp_pi = IR_data_node(IRtypes.float_64_data, :tmppi, [:fv_dofs_global], []);
        tmp_ki = IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], []);
        piki_loop_one = IR_loop_node(IRtypes.space_loop, :dofs, :piki_i, 1, :fv_dofs_global, IR_block_node([
            # tmpki .= stepper.dt .* sol;
            # tmppi .= tmppi + stepper.b[rki].*tmpki;
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], [:piki_i]),
                IR_operation_node(IRtypes.math_op, [:*, :dt, IR_data_node(IRtypes.float_64_data, :solution, [:fv_dofs_global], [:piki_i])])
            ]),
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :tmppi, [:fv_dofs_global], [:piki_i]),
                IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.float_64_data, :tmppi, [:fv_dofs_global], [:piki_i]), 
                    IR_operation_node(IRtypes.math_op, [:*, 
                        IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :b, [:?], [:rki])]), 
                        IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], [:piki_i])])
                ])
            ])
        ]))
        piki_loop_two = IR_loop_node(IRtypes.space_loop, :dofs, :piki_i, 1, :fv_dofs_global, IR_block_node([
            # tmpki .= stepper.a[rki].*tmpki + stepper.dt .* sol;
            # tmppi .= tmppi + stepper.b[rki].*tmpki;
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], [:piki_i]),
                IR_operation_node(IRtypes.math_op, [:+, 
                    IR_operation_node(IRtypes.math_op, [:*, 
                        IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :a, [:?], [:rki])]), 
                        IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], [:piki_i])]),
                    IR_operation_node(IRtypes.math_op, [:*, :dt, IR_data_node(IRtypes.float_64_data, :solution, [:fv_dofs_global], [:piki_i])])
                ])
            ]),
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :tmppi, [:fv_dofs_global], [:piki_i]),
                IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.float_64_data, :tmppi, [:fv_dofs_global], [:piki_i]), 
                    IR_operation_node(IRtypes.math_op, [:*, 
                        IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :b, [:?], [:rki])]), 
                        IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], [:piki_i])])
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
                        IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :c, [:?], [:rki])])])])
            ]),
            # assemble
            zero_el_vec,
            zero_bdry_done,
            wrap_in_timer(:step_assembly, assembly),
            # solve
            wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector])),
            # copy_bdry_vals_to_variables(var, sol, grid_data, dofs_per_node, true);
            IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :solution, :var, :true]),
            # update tmppi and tmpki
            piki_condition,
            # copy_bdry_vals_to_vector(var, tmppi, grid_data, dofs_per_node);
            IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :tmppi]),
            wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :tmppi]))
        ]);
        stage_loop = IR_loop_node(IRtypes.time_loop, :stages, :rki, 1, stepper.stages, stage_loop_body);
        
        push!(tloop_body.parts, IR_operation_node(IRtypes.named_op, [:GATHER_SOLUTION, :tmppi]));
        push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:last_t, :t]));
        push!(tloop_body.parts, stage_loop);
        push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:t, IR_operation_node(IRtypes.math_op, [:+, :last_t, :dt])]));
        
        time_loop = IR_block_node([
            IR_operation_node(IRtypes.assign_op, [
                tmp_pi,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :fv_dofs_global])]),
            IR_operation_node(IRtypes.assign_op, [
                tmp_ki,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :fv_dofs_global])]),
                
            IR_loop_node(IRtypes.time_loop, :time, :ti, 1, stepper.Nsteps, tloop_body)
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
                    
                    b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, stime, stepper.dt; rhs_only = true, assemble_loops=assemble_func)
                    
                    solution .= linear_system_solve(A,b);
                    
                    # At this point solution holds the boundary values
                    # directly write them to the variable values and zero sol.
                    copy_bdry_vals_to_variables(var, solution, grid_data, dofs_per_node, zero_vals=true);
                    
                    if stage < time_stepper.stages
                        # Update the values in vars to be used in the next stage
                        for k=1:fv_dofs_global
                            tmpki[k,stage] = solution[k];
                            tmpresult[k] = last_result[k];
                            for j=1:stage
                                tmpresult[k] += stepper.dt * stepper.a[stage+1, j] * tmpki[k,j];
                            end
                        end
                        
                        copy_bdry_vals_to_vector(var, tmpresult, grid_data, dofs_per_node);
                        place_sol_in_vars(var, tmpresult, stepper);
                        
                    else
                        for k=1:fv_dofs_global
                            tmpki[k,stage] = solution[k];
                        end
                    end
                    
                end
                for k=1:fv_dofs_global
                    for stage=1:stepper.stages
                        last_result[k] += stepper.dt * stepper.b[stage] .* tmpki[k, stage];
                    end
                end
                copy_bdry_vals_to_vector(var, last_result, grid_data, dofs_per_node);
                place_sol_in_vars(var, last_result, stepper);
        =#
        
        tmp_last = IR_data_node(IRtypes.float_64_data, :last_result, [:fv_dofs_global], []);
        tmp_result = IR_data_node(IRtypes.float_64_data, :tmpresult, [:fv_dofs_global], []);
        tmp_ki = IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], []);
        
        update_ki_loop_one = IR_loop_node(IRtypes.space_loop, :dofs, :k, 1, :fv_dofs_global, IR_block_node([
            # # Update the values in vars to be used in the next stage
            # for k=1:fv_dofs_global
            #     tmpki[k,stage] = solution[k];
            #     tmpresult[k] = last_result[k];
            #     for j=1:stage
            #         tmpresult[k] += stepper.dt * stepper.a[stage+1, j] * tmpki[k,j];
            #     end
            # end
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], [:k, :rki]),
                IR_data_node(IRtypes.float_64_data, :solution, [:fv_dofs_global], [:k])
            ]),
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :tmpresult, [:fv_dofs_global], [:k]),
                IR_data_node(IRtypes.float_64_data, :last_result, [:fv_dofs_global], [:k])
            ]),
            IR_loop_node(IRtypes.time_loop, :stage, :j, 1, :rki, IR_block_node([
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_64_data, :tmpresult, [:fv_dofs_global], [:k]),
                    IR_operation_node(IRtypes.math_op, [:+,
                            IR_data_node(IRtypes.float_64_data, :tmpresult, [:dofs_global], [:k]),
                        IR_operation_node(IRtypes.math_op, [:*, :dt, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, 
                                IR_data_node(IRtypes.float_64_data, :a, [:?], [IR_operation_node(IRtypes.math_op, [:+, :rki, 1]), :j])]),
                            IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], [:k, :j])
                        ])
                    ])
                ])
            ]))
        ]))
        update_ki_loop_two = IR_loop_node(IRtypes.space_loop, :dofs, :k, 1, :fv_dofs_global, IR_block_node([
            # for k=1:fv_dofs_global
            #     tmpki[k,stage] = solution[k];
            # end
            IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], [:k, :rki]),
                IR_data_node(IRtypes.float_64_data, :solution, [:fv_dofs_global], [:k])
            ])
        ]))
        update_ki_condition = IR_conditional_node(IR_operation_node(IRtypes.math_op, [:<, :rki, stepper.stages]),
            IR_block_node([
                update_ki_loop_one,
                # copy_bdry_vals_to_vector(var, tmpresult, grid_data, dofs_per_node);
                IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :tmpresult]),
                # place_sol_in_vars(var, tmpresult, stepper);
                wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :tmpresult]))
            ]),
            IR_block_node([update_ki_loop_two]));
        
        stage_loop_body = IR_block_node([
            # last_t = t;
            # t = last_t + stepper.c[rki]*stepper.dt;
            IR_operation_node(IRtypes.assign_op, [
                :t,
                IR_operation_node(IRtypes.math_op, [:+, :t, 
                    IR_operation_node(IRtypes.math_op, [:*, :dt, 
                        IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :c, [:?], [:rki])])])])
            ]),
            # assemble
            zero_el_vec,
            zero_bdry_done,
            wrap_in_timer(:step_assembly, assembly),
            # solve
            wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector])),
            # copy_bdry_vals_to_variables(var, sol, grid_data, dofs_per_node, true);
            IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :solution, :var, :true]),
            # update tmpki
            update_ki_condition
        ]);
        stage_loop = IR_loop_node(IRtypes.time_loop, :stages, :rki, 1, stepper.stages, stage_loop_body);
        
        # last_result[k] += stepper.dt * stepper.b[rki] * tmpki[k, rki];
        combine_loop_body = IR_block_node([
            IR_loop_node(IRtypes.time_loop, :stages, :rki, 1, stepper.stages, IR_block_node([
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.float_64_data, :last_result, [:fv_dofs_global], [:k]),
                    IR_operation_node(IRtypes.math_op, [:+,
                        IR_data_node(IRtypes.float_64_data, :last_result, [:fv_dofs_global], [:k]), 
                        IR_operation_node(IRtypes.math_op, [:*, :dt, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.float_64_data, :b, [:?], [:rki])]),
                            IR_data_node(IRtypes.float_64_data, :tmpki, [:fv_dofs_global], [:k, :rki])
                        ])
                    ])
                ])
            ]))
        ])
        combine_loop = IR_loop_node(IRtypes.space_loop, :dofs, :k, 1, :fv_dofs_global, combine_loop_body);
        
        push!(tloop_body.parts, IR_operation_node(IRtypes.named_op, [:GATHER_SOLUTION, :last_result]));
        push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:last_t, :t]));
        push!(tloop_body.parts, stage_loop);
        push!(tloop_body.parts, combine_loop);
        # copy_bdry_vals_to_vector(var, last_result, grid_data, dofs_per_node);
        push!(tloop_body.parts, IR_operation_node(IRtypes.named_op, [:BDRY_TO_VECTOR, :last_result]))
        # place_sol_in_vars(var, last_result, stepper);
        push!(tloop_body.parts, wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :last_result])));
        # update time
        push!(tloop_body.parts, IR_operation_node(IRtypes.assign_op, [:t, IR_operation_node(IRtypes.math_op, [:+, :last_t, :dt])]));
        
        time_loop = IR_block_node([
            IR_operation_node(IRtypes.assign_op, [
                tmp_last,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :fv_dofs_global])]),
            IR_operation_node(IRtypes.assign_op, [
                tmp_result,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :fv_dofs_global])]),
            IR_operation_node(IRtypes.assign_op, [
                tmp_ki,
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :fv_dofs_global, stepper.stages])]),
                
            IR_loop_node(IRtypes.time_loop, :time, :ti, 1, stepper.Nsteps, tloop_body)
        ]);
        
    else
        # implicit steppers TODO
    end
    
    return time_loop;
end

# TDM(A,b,C) = transpose(A) * diagm(b) * C
# P_ij =  A'_ik * b_k * C_kj = A_ki * b_k * C_kj
# return the IR for A[k,i] * b[k] * C[k,j]
# A,b,C are symbols or IR_data_node
# i,j,k are symbols, numbers, IR_part
function generate_linalg_TDM_product(A, b, C, i, j, k)
    IRtypes = IR_entry_types();
    if typeof(A) == IR_data_node
        A_part = IR_data_node(IRtypes.float_64_data, A.label, A.size, [k,i]);
    elseif typeof(A) <: IR_part
        A_part = apply_indexed_access(A, [k,i], IRtypes);
    else
        A_part = IR_data_node(IRtypes.float_64_data, A, [:?,:?], [k,i]);
    end
    if typeof(b) == IR_data_node
        b_part = IR_data_node(IRtypes.float_64_data, b.label, b.size, [k]);
    elseif typeof(b) <: IR_part
        b_part = apply_indexed_access(b, [k], IRtypes);
    else
        b_part = IR_data_node(IRtypes.float_64_data, b, [:?], [k]);
    end
    if typeof(C) == IR_data_node
        C_part = IR_data_node(IRtypes.float_64_data, C.label, C.size, [k,j]);
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
        A_part = IR_data_node(IRtypes.float_64_data, A.label, A.size, [j,i]);
    elseif typeof(A) <: IR_part
        A_part = apply_indexed_access(A, [j,i], IRtypes);
    else
        A_part = IR_data_node(IRtypes.float_64_data, A, [:?,:?], [j,i]);
    end
    if typeof(b) == IR_data_node
        b_part = IR_data_node(IRtypes.float_64_data, b.label, [j]);
    elseif typeof(b) <: IR_part
        b_part = apply_indexed_access(b, [j], IRtypes);
    else
        b_part = IR_data_node(IRtypes.float_64_data, b, [:?], [j]);
    end
    
    return IR_operation_node(IRtypes.math_op, [:*, A_part, b_part]);
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