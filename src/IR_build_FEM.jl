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
function build_IR_fem(lhs_vol, lhs_surf, rhs_vol, rhs_surf, var, config, prob, time_stepper)
    dimension = config.dimension;
    # Count variables, dofs, and store offsets
    multivar = typeof(var) <:Array;
    varcount = 1;
    dofsper = 0;
    offset_ind = [0];
    if multivar
        varcount = length(var);
        offset_ind = zeros(Int, varcount);
        dofsper = length(var[1].symvar);
        for i=2:length(var)
            offset_ind[i] = dofsper;
            dofsper = dofsper + length(var[i].symvar);
        end
    else
        var = [var]; # put it in an array for consistency
    end
    
    IRtypes = IR_entry_types();
    
    # These will hold the IR
    allocate_block = IR_block_node([]);
    prepare_block = IR_block_node([]);
    derivmat_block = IR_block_node([]);
    time_prepare_block = IR_block_node([]);
    time_derivmat_block = IR_block_node([]);
    matrix_block = IR_block_node([]);
    vector_block = IR_block_node([]);
    bdry_block = IR_block_node([]);
    toglobal_block = IR_block_node([]);
    time_bdry_block = IR_block_node([]);
    time_toglobal_block = IR_block_node([]);
    
    # Allocate the global matrix and vector
    push!(allocate_block.parts, IR_comment_node("Allocate global matrix(IJV form) and vector."));
    allocatedNZ = IR_data_node(IRtypes.scalar_data, :allocated_nonzeros);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        allocatedNZ,
        IR_operation_node(IRtypes.math_op, [:*, :num_elements, :dofs_per_element, :dofs_per_element])
        ]));
    globalmat_I = IR_data_node(IRtypes.array_data, :global_matrix_I);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        globalmat_I,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.int_64_data, :allocated_nonzeros])
        ]));
    globalmat_J = IR_data_node(IRtypes.array_data, :global_matrix_J);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        globalmat_J,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.int_64_data, :allocated_nonzeros])
        ]));
    globalmat_V = IR_data_node(IRtypes.array_data, :global_matrix_V);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        globalmat_V,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :allocated_nonzeros])
        ]));
        
    globalvec = IR_data_node(IRtypes.array_data, :global_vector);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        globalvec,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])
        ]));
    # I and J should be initialized to 1
    push!(allocate_block.parts, IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY,
        globalmat_I,
        1
        ]));
    push!(allocate_block.parts, IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY,
        globalmat_J,
        1
        ]));
    
    # Allocate the elemental matrix and vector
    push!(allocate_block.parts, IR_comment_node("Allocate elemental matrix and vector."));
    elementmat = IR_data_node(IRtypes.array_data, :element_matrix);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        elementmat,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_per_element, :dofs_per_element])
        ]));
    
    elementvec = IR_data_node(IRtypes.array_data, :element_vector);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        elementvec,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_per_element])
        ]));
    
    # bdry done flag for each node
    push!(allocate_block.parts, IR_comment_node("Boundary done flag for each node."));
    bdry_done = IR_data_node(IRtypes.array_data, :bdry_done);
    push!(allocate_block.parts, IR_operation_node(IRtypes.assign_op, [
        bdry_done,
        IR_operation_node(IRtypes.allocate_op, [IRtypes.boolean_data, :nnodes_global])
        ]));
    
    # coefficient prep
    # a list of all entities and rhs only ones
    all_entities = [];
    rhs_entities = [];
    counts = zeros(Int,4); # how many entities for each piece 
    
    # LHS volume
    if !(lhs_vol === nothing)
        entities = extract_entities(lhs_vol, multivar);
        append!(all_entities, entities);
        counts[1] = length(entities);
    end
    
    # RHS volume
    if !(rhs_vol === nothing)
        entities = extract_entities(rhs_vol, multivar);
        append!(all_entities, entities);
        append!(rhs_entities, entities);
        counts[2] = length(entities);
    end
    
    # LHS surface
    if !(lhs_surf === nothing)
        entities = extract_entities(lhs_surf, multivar);
        append!(all_entities, entities);
        counts[3] = length(entities);
    end
    
    # RHS surface
    if !(rhs_surf === nothing)
        entities = extract_entities(rhs_surf, multivar);
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
            IR_data_node(IRtypes.array_data, :detJ, [1, :eid])])]))
    
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
            IR_data_node(IRtypes.array_data, :detJ, [1, :eid])])]))
    
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
    # apply_boundary_conditions_elemental(var, eid, grid, refel, geo_facs, prob, t, elmat, elvec)
    push!(bdry_block.parts, IR_comment_node("Apply boundary conditions."));
    push!(bdry_block.parts, IR_operation_node(IRtypes.function_op, [
        IR_operation_node(IRtypes.member_op, [:CGSolver, :apply_boundary_conditions_elemental]), 
        :var, :eid, :mesh, :refel, :geometric_factors, :prob, :t, :element_matrix, :element_vector, :bdry_done]));
    push!(time_bdry_block.parts, IR_comment_node("Apply boundary conditions."));
    push!(time_bdry_block.parts, IR_operation_node(IRtypes.function_op, [
        IR_operation_node(IRtypes.member_op, [:CGSolver, :apply_boundary_conditions_elemental_rhs]), 
        :var, :eid, :mesh, :refel, :geometric_factors, :prob, :t, :element_vector, :bdry_done]));
    
    # add to global sysem
    push!(toglobal_block.parts, IR_comment_node("Place elemental parts in global system."));
    push!(toglobal_block.parts, generate_local_to_global_fem(dofsper, offset_ind));
    push!(time_toglobal_block.parts, IR_comment_node("Place elemental parts in global system."));
    push!(time_toglobal_block.parts, generate_local_to_global_fem(dofsper, offset_ind, vec_only=true));
    
    # assembly loop
    indices = [];
    # check for indexed variables TODO
    assembly_loop = generate_assembly_loop_fem(var, indices);
    rhs_assembly_loop = generate_assembly_loop_fem(var, indices);
    
    # fill the loop
    assembly_loop.body = IR_block_node([
        prepare_block,
        derivmat_block,
        wrap_in_timer(:el_matrix, matrix_block),
        wrap_in_timer(:el_vector, vector_block),
        wrap_in_timer(:bdry_cond, bdry_block),
        toglobal_block
    ])
    rhs_assembly_loop.body = IR_block_node([
        time_prepare_block,
        time_derivmat_block,
        wrap_in_timer(:el_vector, vector_block),
        wrap_in_timer(:bdry_cond, time_bdry_block),
        time_toglobal_block
    ])
    
    # time loop
    need_time_loop = prob.time_dependent;
    if need_time_loop
        step_loop = generate_time_stepping_loop_fem(time_stepper, rhs_assembly_loop);
        # step_loop.body = IR_block_node([rhs_assembly_loop]); # TODO add solve and scatter
        time_loop = IR_block_node([
            IR_operation_node(IRtypes.assign_op, [:t, 0]),
            IR_operation_node(IRtypes.assign_op, [:dt, time_stepper.dt]),
            IR_comment_node("Initial loop to build matrix"),
            wrap_in_timer(:first_assembly, assembly_loop),
            IR_operation_node(IRtypes.named_op, [:GLOBAL_FINALIZE]),
            IR_operation_node(IRtypes.named_op, [:GATHER_SOLUTION, :solution]),
            
            IR_comment_node("###############################################"),
            IR_comment_node("Time stepping loop"),
            step_loop
        ]);
    else
        time_loop = IR_block_node([
            IR_operation_node(IRtypes.assign_op, [:t, 0]),
            wrap_in_timer(:first_assembly, assembly_loop),
            IR_operation_node(IRtypes.named_op, [:GLOBAL_FINALIZE]),
            wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector])),
            wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :solution])),
        ])
    end
    
    # Put them all together in a master block
    master_block = IR_block_node([
        allocate_block,
        time_loop
    ]);
    
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
                IR_data_node(IRtypes.array_data, Symbol("RQ"*string(i))),
                IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :qnodes_per_element, :nodes_per_element])
                ]));
        end
        if needed_derivative_matrices[i+4]
            push!(allocate_part, IR_operation_node(
                IRtypes.assign_op,[
                IR_data_node(IRtypes.array_data, Symbol("RD"*string(i))),
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
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.array_data, :loc2glb, [:ni, :eid])])]),
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
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.array_data, :allnodes, [1,:nodeID])])])
    ];
    if dimension > 1
        # y = mesh.allnodes[nodeID*dimension+1]
        push!(coef_loop_body, IR_operation_node(IRtypes.assign_op, [
            :y, 
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.array_data, :allnodes, [2,:nodeID])])]))
    else
        # y = 0
        push!(coef_loop_body, IR_operation_node(IRtypes.assign_op, [:y, 0]));
    end
    if dimension > 2
        # z = mesh.allnodes[nodeID*dimension+2]
        push!(coef_loop_body, IR_operation_node(IRtypes.assign_op, [
            :z, 
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.array_data, :allnodes, [3,:nodeID])])]))
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
               # TODO
                
            elseif ctype == 2 # a coefficient function
                if vors == "volume"
                    need_vol_interp_loop = true;
                    need_vol_coef_loop = true;
                    # Need to allocate NP 
                    np_allocate = IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :nodes_per_element]);
                    nqp_allocate = IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :qnodes_per_element]);
                    nodal_coef_name = "NODAL"*cname;
                    # For the nodal values
                    push!(allocate_part, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name)),
                        np_allocate]));
                    # For the interpolated/differentiated quadrature values
                    push!(allocate_part, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.array_data, Symbol(cname)),
                        nqp_allocate]));
                    
                    coef_index = get_coef_index(entities[i]);
                    
                    push!(coef_loop_body, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name), [:ni]),
                        IR_operation_node(IRtypes.named_op, [:COEF_EVAL, coef_index, entities[i].index, :x, :y, :z, :t, :nodeID])]));
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        quad_coef_node = IR_data_node(IRtypes.array_data, Symbol(cname), [:col]);
                        nodal_coef_node = IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name), [:row]);
                        push!(coef_init_loop_body, IR_operation_node(IRtypes.assign_op, [quad_coef_node, 0.0]));
                        for di=1:length(entities[i].derivs)
                            deriv_quad_mat = IR_data_node(IRtypes.array_data, Symbol("RQ"*string(entities[i].derivs[di])), [col_row_matrix_index]);
                            push!(coef_interp_loop_body, IR_operation_node(IRtypes.assign_op,[
                                quad_coef_node,
                                IR_operation_node(IRtypes.math_op, [:(+), quad_coef_node, IR_operation_node(IRtypes.math_op, [:(*), deriv_quad_mat, nodal_coef_node])])
                                ]));
                        end
                    else
                        quad_coef_node = IR_data_node(IRtypes.array_data, Symbol(cname), [:col]);
                        nodal_coef_node = IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name), [:row]);
                        # refelQ = IR_operation_node(IRtypes.member_op, [:refel, IR_data_node(IRtypes.array_data, :Q, [row_col_matrix_index])]);
                        refelQ = IR_data_node(IRtypes.array_data, :Q, [col_row_matrix_index]);
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
                # TODO
                # 
                if vors == "volume"
                    need_vol_interp_loop = true;
                    need_vol_coef_loop = true;
                    # Need to allocate NP 
                    np_allocate = IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :nodes_per_element]);
                    nqp_allocate = IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :qnodes_per_element]);
                    nodal_coef_name = "NODAL"*cname;
                    # For the nodal values
                    push!(allocate_part, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name)),
                        np_allocate]));
                    # For the interpolated/differentiated quadrature values
                    push!(allocate_part, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.array_data, Symbol(cname)),
                        nqp_allocate]));
                    
                    push!(coef_loop_body, IR_operation_node(IRtypes.assign_op,[
                        IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name), [:ni]),
                        IR_operation_node(IRtypes.named_op, [:KNOWN_VAR, cval, entities[i].index, :nodeID])]));
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        quad_coef_node = IR_data_node(IRtypes.array_data, Symbol(cname), [:col]);
                        nodal_coef_node = IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name), [:row]);
                        push!(coef_init_loop_body, IR_operation_node(IRtypes.assign_op, [quad_coef_node, 0.0]));
                        for di=1:length(entities[i].derivs)
                            deriv_quad_mat = IR_data_node(IRtypes.array_data, Symbol("RQ"*string(entities[i].derivs[di])), [col_row_matrix_index]);
                            push!(coef_interp_loop_body, IR_operation_node(IRtypes.assign_op,[
                                quad_coef_node,
                                IR_operation_node(IRtypes.math_op, [:(+), quad_coef_node, IR_operation_node(IRtypes.math_op, [:(*), deriv_quad_mat, nodal_coef_node])])
                                ]));
                        end
                    else
                        quad_coef_node = IR_data_node(IRtypes.array_data, Symbol(cname), [:col]);
                        nodal_coef_node = IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name), [:row]);
                        # refelQ = IR_operation_node(IRtypes.member_op, [:refel, IR_data_node(IRtypes.array_data, :Q, [row_col_matrix_index])]);
                        refelQ = IR_data_node(IRtypes.array_data, :Q, [col_row_matrix_index]);
                        push!(coef_init_loop_body, IR_operation_node(IRtypes.assign_op, [quad_coef_node, 0.0]));
                        push!(coef_interp_loop_body, IR_operation_node(IRtypes.assign_op,[
                            quad_coef_node,
                            IR_operation_node(IRtypes.math_op, [:(+), quad_coef_node, IR_operation_node(IRtypes.math_op, [:(*), refelQ, nodal_coef_node])])
                            ]));
                    end
                else
                    # TODO
                end
            elseif ctype == 4 # an indexed coefficient
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
    compute_block = IR_block_node([]);
    
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
            push!(linalg_matrix_block_args, :LINALG_MATRIX_BLOCK);
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
                        #                             :LINALG_MATRIX_BLOCK, 
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
        else
            linalg_vector_block_args = [];
            push!(linalg_vector_block_args, :LINALG_VECTOR_BLOCK);
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
                    #                             :LINALG_VECTOR_BLOCK, 
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
                    :LINALG_MATRIX_BLOCK, 1, :nodes_per_element, :element_matrix, 1, 1, combined_rhs]));
            else
                push!(compute_block.parts, IR_operation_node(IRtypes.named_op, [ 
                    :LINALG_VECTOR_BLOCK, 1, :nodes_per_element, :element_vector, 1, combined_rhs]));
            end
        else
            if lorr == LHS
                push!(compute_block.parts, IR_operation_node(IRtypes.named_op, [
                    :LINALG_MATRIX_BLOCK, 1, :nodes_per_element, :element_matrix, 1, 1, term_vec[1]]));
            else
                push!(compute_block.parts, IR_operation_node(IRtypes.named_op, [ 
                    :LINALG_VECTOR_BLOCK, 1, :nodes_per_element, :element_vector, 1, term_vec[1]]));
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
            test_part = IR_data_node(IRtypes.matrix_data, Symbol("RQ"*string(deriv_index)), [:qnodes_per_element, :nodes_per_element]);
        else
            # test_part = IR_operation_node(IRtypes.member_op, [:refel, IR_data_node(IRtypes.matrix_data, :Q, [:qnodes_per_element, :nodes_per_element])]);
            test_part = IR_data_node(IRtypes.matrix_data, :Q, [:qnodes_per_element, :nodes_per_element]);
        end
    else
        test_part = nothing;
    end
    if typeof(trial_ex) == SymEntity
        if length(trial_ex.derivs) > 0
            deriv_index = trial_ex.derivs[1];
            trial_part = IR_data_node(IRtypes.matrix_data, Symbol("RQ"*string(deriv_index)), [:qnodes_per_element, :nodes_per_element]);
        else
            # trial_part = IR_operation_node(IRtypes.member_op, [:refel, IR_data_node(IRtypes.matrix_data, :Q, [:qnodes_per_element, :nodes_per_element])]);
            trial_part = IR_data_node(IRtypes.matrix_data, :Q, [:qnodes_per_element, :nodes_per_element]);
        end
    else
        trial_part = nothing;
    end
    
    # Turn the coefficient part into IR
    wg_part = IR_data_node(IRtypes.array_data, :wg);
    if false # non-constant J
        detj_part = IR_data_node(IRtypes.array_data, :detj);
    else
        detj_part = :detj;
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

# This part inserts the elemental matrix and vector into the global one
function generate_local_to_global_fem(dofs_per_node, offset_ind; vec_only=false)
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
                        IR_operation_node(IRtypes.math_op, [:*, 
                            IR_operation_node(IRtypes.math_op, [:-, :eid, 1]), :nodes_per_element]),
                            :nodes_per_element])])
            ]));
        end
        
        push!(vector_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
            :glb_i,
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.array_data, :loc2glb, [:ni, :eid])])
        ]));
        push!(vector_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.array_data, :global_vector, [:glb_i]),
            IR_operation_node(IRtypes.math_op, [:+, 
                IR_data_node(IRtypes.array_data, :global_vector, [:glb_i]),
                IR_data_node(IRtypes.array_data, :element_vector, [:ni])])
        ]));
        
        # only to the matrix part if not vec_only
        if !vec_only
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                :glb_j,
                IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.array_data, :loc2glb, [:nj, :eid])])
            ]));
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.array_data, :global_matrix_I, [:next_ind]),
                :glb_i
            ]));
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.array_data, :global_matrix_J, [:next_ind]),
                :glb_j
            ]));
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.array_data, :global_matrix_V, [:next_ind]),
                IR_data_node(IRtypes.array_data, :element_matrix, [:ni, :nj])
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
        # next_ind = +(1, *(*(-(eid, 1), dofs_per_element), dofs_per_element));
        # for ni = 1:nodes_per_element
        #     glb_i = mesh.loc2glb[ni, eid];
        #     for di = 1:dofs_per_node
        #         glb_dofi = (glb_i - 1)*dofs_per_node + di;
        #         global_vector[glb_dofi] = +(global_vector[glb_dofi], element_vector[(di-1)*nodes_per_element + ni]);
        #         for nj = 1:nodes_per_element
        #             glb_j = mesh.loc2glb[nj, eid];
        #             for dj = 1:dofs_per_node
        #                 glb_dofj = (glb_j - 1)*dofs_per_node + dj;
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
                        IR_operation_node(IRtypes.math_op, [:*, 
                            IR_operation_node(IRtypes.math_op, [:-, :eid, 1]), :dofs_per_element]),
                            :dofs_per_element])])
            ]));
        end
        
        push!(vector_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
            :glb_i,
            IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.array_data, :loc2glb, [:ni, :eid])])
        ]));
        push!(vector_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
            :glb_dofi,
            generate_the_one_pattern(:glb_i, 1, :dofs_per_node, :di)
        ]));
        push!(vector_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
            IR_data_node(IRtypes.array_data, :global_vector, [:glb_dofi]),
            IR_operation_node(IRtypes.math_op, [:+, 
                IR_data_node(IRtypes.array_data, :global_vector, [:glb_dofi]),
                IR_data_node(IRtypes.array_data, :element_vector, [generate_the_one_pattern(:di, 1, :nodes_per_element, :ni)])])
        ]));
        
        # only needed if not vec_only
        if !vec_only
            push!(matrix_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                :glb_j,
                IR_operation_node(IRtypes.member_op, [:mesh, IR_data_node(IRtypes.array_data, :loc2glb, [:nj, :eid])])
            ]));
            push!(matrix_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                :glb_dofj,
                generate_the_one_pattern(:glb_j, 1, :dofs_per_node, :dj)
            ]));
            push!(matrix_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.array_data, :global_matrix_I, [:next_ind]),
                :glb_dofi
            ]));
            push!(matrix_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.array_data, :global_matrix_J, [:next_ind]),
                :glb_dofj
            ]));
            push!(matrix_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                IR_data_node(IRtypes.array_data, :global_matrix_V, [:next_ind]),
                IR_data_node(IRtypes.array_data, :element_matrix, [
                    generate_the_one_pattern(:di, 1, :nodes_per_element, :ni),
                    generate_the_one_pattern(:dj, 1, :nodes_per_element, :nj)])
            ]));
            push!(matrix_dof_loop_body.parts, IR_operation_node(IRtypes.assign_op, [
                :next_ind,
                IR_operation_node(IRtypes.math_op, [:+, :next_ind, 1])
            ]));
            
            push!(matrix_loop_body.parts, IR_loop_node(IRtypes.space_loop, :nodes, :dj, 1, :dofs_per_node, matrix_dof_loop_body));
            push!(vector_dof_loop_body.parts, IR_loop_node(IRtypes.space_loop, :nodes, :nj, 1, :nodes_per_element, matrix_loop_body));
        end
        
        push!(vector_loop_body.parts, IR_loop_node(IRtypes.space_loop, :nodes, :di, 1, :dofs_per_node, vector_dof_loop_body));
        toglobal_loop = IR_loop_node(IRtypes.space_loop, :nodes, :ni, 1, :nodes_per_element, vector_loop_body);
        push!(result_block.parts, toglobal_loop)
        
        toglobal_loop = IR_loop_node(IRtypes.space_loop, :nodes, :ni, 1, :nodes_per_element, vector_loop_body);
    end
    
    return result_block;
end

# Generate the assembly loop structures and insert the content
function generate_assembly_loop_fem(var, indices=[])
    IRtypes = IR_entry_types();
    # Each of the indices must be passed to the functions in a named tuple.
    # Pass all defined indexers.
    index_names = [];
    elements_included = false;
    ind_shift = 0;
    for i=1:length(indices)
        if indices[i] == "elements" || indices[i] == "cells"
            elements_included = true;
            push!(index_names, "elements");
        else
            push!(index_names, IR_data_node(IRtypes.scalar_data, Symbol("index_val_"*string(indices[i].symbol))));
        end
    end
    
    # If elements were not included, make them the outermost loop
    if !elements_included && length(index_names) > 0
        index_names = ["elements"; index_names];
        ind_shift = 1;
    end
    
    # Placeholder computation that will be replaced by the actual one
    placeholder = IR_block_node([]);
    
    # generate the loop structures
    if length(index_names) > 0
        # The innermost loop holds placeholder
        if index_names[end] == "elements"
            assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :eid, 1, :num_elements, placeholder);
        else
            assembly_loop = IR_loop_node(IRtypes.index_loop, indices[end].symbol, 
                            index_names[end].var, indices[end].range[1], indices[end].range[end], placeholder);
        end
        
        # work outwards nesting assembly_loop
        for i=(length(index_names)-1):-1:1
            if index_names[end] == "elements"
                assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :eid, 1, :num_elements, assembly_loop);
            else
                assembly_loop = IR_loop_node(IRtypes.index_loop, indices[i-ind_shift].symbol, 
                                index_names[i].var, indices[i-ind_shift].range[1], 
                                indices[i-ind_shift].range[i-ind_shift], assembly_loop);
            end
        end
    else # only an element loop
        assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :eid, 1, :num_elements, placeholder);
    end
    
    return assembly_loop
end

# Generates the time stepping loop using the supplied assembly block and stepper
function generate_time_stepping_loop_fem(stepper, assembly)
    IRtypes = IR_entry_types();
    tloop_body = IR_block_node([]);
    zero_el_vec = IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY,
        IR_data_node(IRtypes.array_data, :global_vector),
        0
        ]);
    zero_bdry_done = IR_operation_node(IRtypes.named_op, [
        :FILL_ARRAY,
        IR_data_node(IRtypes.array_data, :bdry_done),
        :false
        ]);
    if stepper.stages < 2
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
            wrap_in_timer(:step_assembly, assembly),
            wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector])),
            wrap_in_timer(:scatter, IR_operation_node(IRtypes.named_op, [:SCATTER_SOLUTION, :solution])),
            IR_operation_node(IRtypes.assign_op, [
                :t,
                IR_operation_node(IRtypes.math_op, [:+, :t, stepper.dt])
            ])
        ];
        
        time_loop = IR_loop_node(IRtypes.time_loop, :time, :t, 1, stepper.Nsteps, tloop_body);
        
    else # multistage explicit steppers
        # LSRK4 is a special case, low storage
        if stepper.type == LSRK4
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
            tmp_pi = IR_data_node(IRtypes.array_data, :tmppi);
            tmp_ki = IR_data_node(IRtypes.array_data, :tmpki);
            piki_loop_one = IR_loop_node(IRtypes.space_loop, :dofs, :piki_i, 1, :dofs_global, IR_block_node([
                # tmpki .= stepper.dt .* sol;
                # tmppi .= tmppi + stepper.b[rki].*tmpki;
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.array_data, :tmpki, [:piki_i]),
                    IR_operation_node(IRtypes.math_op, [:*, :dt, IR_data_node(IRtypes.array_data, :solution, [:piki_i])])
                ]),
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.array_data, :tmppi, [:piki_i]),
                    IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.array_data, :tmppi, [:piki_i]), 
                        IR_operation_node(IRtypes.math_op, [:*, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.array_data, :b, [:rki])]), 
                            IR_data_node(IRtypes.array_data, :tmpki, [:piki_i])])
                    ])
                ])
            ]))
            piki_loop_two = IR_loop_node(IRtypes.space_loop, :dofs, :piki_i, 1, :dofs_global, IR_block_node([
                # tmpki .= stepper.a[rki].*tmpki + stepper.dt .* sol;
                # tmppi .= tmppi + stepper.b[rki].*tmpki;
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.array_data, :tmpki, [:piki_i]),
                    IR_operation_node(IRtypes.math_op, [:+, 
                        IR_operation_node(IRtypes.math_op, [:*, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.array_data, :a, [:rki])]), 
                            IR_data_node(IRtypes.array_data, :tmpki, [:piki_i])]),
                        IR_operation_node(IRtypes.math_op, [:*, :dt, IR_data_node(IRtypes.array_data, :solution, [:piki_i])])
                    ])
                ]),
                IR_operation_node(IRtypes.assign_op, [
                    IR_data_node(IRtypes.array_data, :tmppi, [:piki_i]),
                    IR_operation_node(IRtypes.math_op, [:+, IR_data_node(IRtypes.array_data, :tmppi, [:piki_i]), 
                        IR_operation_node(IRtypes.math_op, [:*, 
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.array_data, :b, [:rki])]), 
                            IR_data_node(IRtypes.array_data, :tmpki, [:piki_i])])
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
                            IR_operation_node(IRtypes.member_op, [:time_stepper, IR_data_node(IRtypes.array_data, :c, [:rki])])])])
                ]),
                # assemble
                zero_el_vec,
                zero_bdry_done,
                wrap_in_timer(:step_assembly, assembly),
                # solve
                wrap_in_timer(:lin_solve, IR_operation_node(IRtypes.named_op, [:GLOBAL_SOLVE, :solution, :global_matrix, :global_vector])),
                # copy_bdry_vals_to_variables(var, sol, grid_data, dofs_per_node, zero_vals=true);
                IR_operation_node(IRtypes.function_op, [
                    IR_operation_node(IRtypes.member_op, [:CGSolver, :copy_bdry_vals_to_variables]), 
                    :var, :solution, :mesh, :dofs_per_node, IR_operation_node(IRtypes.assign_op, [:zero_vals, :true])]),
                # update tmppi and tmpki
                piki_condition,
                # copy_bdry_vals_to_vector(var, tmppi, grid_data, dofs_per_node);
                IR_operation_node(IRtypes.function_op, [
                    IR_operation_node(IRtypes.member_op, [:CGSolver, :copy_bdry_vals_to_vector]), 
                    :var, :tmppi, :mesh, :dofs_per_node]),
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
                    IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])]),
                IR_operation_node(IRtypes.assign_op, [
                    tmp_ki,
                    IR_operation_node(IRtypes.allocate_op, [IRtypes.float_64_data, :dofs_global])]),
                    
                IR_loop_node(IRtypes.time_loop, :time, :t, 1, stepper.Nsteps, tloop_body)
            ]);
        else
            # Explicit multi-stage methods: 
            # x = x + dt*sum(bi*ki)
            # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
            
            #=
                    # will hold the final result
                    last_result = get_var_vals(var, last_result);
                    # will be placed in var.values for each stage
                    # Storage for each stage
                    # tmpki = zeros(length(b), stepper.stages);
                    for stage=1:stepper.stages
                        stime = t + stepper.c[stage]*stepper.dt;
                        
                        # Update the values in vars to be used in this stage
                        if stage > 1
                            initialized_tmpresult = false;
                            for j=1:stage
                                if stepper.a[stage, j] > 0
                                    if !initialized_tmpresult
                                        initialized_tmpresult = true;
                                        for k=1:length(last_result)
                                            tmpresult[k] = last_result[k] + stepper.dt * stepper.a[stage, j] * tmpki[k,j];
                                        end
                                    else
                                        for k=1:length(last_result)
                                            tmpresult[k] += stepper.dt * stepper.a[stage, j] * tmpki[k,j];
                                        end
                                    end
                                end
                            end
                            
                            if initialized_tmpresult
                                copy_bdry_vals_to_vector(var, tmpresult, grid_data, dofs_per_node);
                                place_sol_in_vars(var, tmpresult, stepper);
                            end
                        end
                        
                        b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, stime, stepper.dt; rhs_only = true, assemble_loops=assemble_func)
                        
                        tmpki[:,stage] = linear_system_solve(A,b);
                        
                        # At this point tmpki[:,stage] holds the boundary values
                        # directly write them to the variable values and zero sol.
                        copy_bdry_vals_to_variables(var, tmpki[:,stage], grid_data, dofs_per_node, zero_vals=true);
                        
                    end
                    for stage=1:stepper.stages
                        last_result += stepper.dt * stepper.b[stage] .* tmpki[:, stage];
                    end
                    copy_bdry_vals_to_vector(var, last_result, grid_data, dofs_per_node);
                    place_sol_in_vars(var, last_result, stepper);
            =#
            
            time_loop = IR_block_node([
                IR_loop_node(IRtypes.time_loop, :time, :t, 1, stepper.Nsteps, tloop_body);
            ]);
        end
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
        A_part = IR_data_node(IRtypes.array_data, A.var, [k,i]);
    elseif typeof(A) <: IR_part
        A_part = apply_indexed_access(A, [k,i], IRtypes);
    else
        A_part = IR_data_node(IRtypes.array_data, A, [k,i]);
    end
    if typeof(b) == IR_data_node
        b_part = IR_data_node(IRtypes.array_data, b.var, [k]);
    elseif typeof(b) <: IR_part
        b_part = apply_indexed_access(b, [k], IRtypes);
    else
        b_part = IR_data_node(IRtypes.array_data, b, [k]);
    end
    if typeof(C) == IR_data_node
        C_part = IR_data_node(IRtypes.array_data, C.var, [k,j]);
    elseif typeof(C) <: IR_part
        C_part = apply_indexed_access(C, [k,j], IRtypes);
    else
        C_part = IR_data_node(IRtypes.array_data, C, [k,j]);
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
        A_part = IR_data_node(IRtypes.array_data, A.var, [j,i]);
    elseif typeof(A) <: IR_part
        A_part = apply_indexed_access(A, [j,i], IRtypes);
    else
        A_part = IR_data_node(IRtypes.array_data, A, [j,i]);
    end
    if typeof(b) == IR_data_node
        b_part = IR_data_node(IRtypes.array_data, b.var, [j]);
    elseif typeof(b) <: IR_part
        b_part = apply_indexed_access(b, [j], IRtypes);
    else
        b_part = IR_data_node(IRtypes.array_data, b, [j]);
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