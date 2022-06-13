#=
Functions for building an IR from the DSL for FEM type problems.
These consist of assembling a global matrix/vector from local ones.
Quadrature is done with matrices from the refel.
=#

# builds the whole IR
#  allocate - prepare_coefficient_values
#  time loop - TODO
#    pre step - TODO
#    element/index loops - generate_assembly_loop
#      prepare values - prepare_coefficient_values
#      elemental mat - make_elemental_computation
#      elemental vec - make_elemental_computation
#      bdry - TODO (make_elemental_computation)
#    solve - handled by codegen
#    post step - TODO
#  finalize - handled by codegen
function build_IR_fem(lhs_vol, lhs_surf, rhs_vol, rhs_surf, var, dimension, solver)
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
    
    # Gather some needed numbers
    # dimension = size(grid_data.allnodes,2);
    # dofs_per_node = 0;
    # dofs_per_loop = 0;
    # for vi=1:length(var)
    #     dofs_per_loop += length(var[vi].symvar);
    #     dofs_per_node += var[vi].total_components;
    # end
    # nnodes_partition = size(grid_data.allnodes,2);
    # nnodes_global = grid_data.nnodes_global;
    # dofs_global = dofs_per_node * nnodes_global;
    # dofs_partition = dofs_per_node * nnodes_partition;
    # nodes_per_element = refel.Np;
    # dofs_per_element = dofs_per_node * nodes_per_element;
    # num_elements = size(grid_data.loc2glb,2);
    
    IRtypes = IR_entry_types();
    
    # These will hold the IR
    allocate_block = IR_block_node([]);
    prepare_block = IR_block_node([]);
    derivmat_block = IR_block_node([]);
    matrix_block = IR_block_node([]);
    vector_block = IR_block_node([]);
    bdry_block = IR_block_node([]);
    
    # separate vol and surf entities
    vol_entities = [];
    surf_entities = [];
    
    # Allocate the elemental matrix and vector
    elementmat = IR_data_node(IRtypes.array_data, :element_matrix);
    push!(allocate_block.parts, IR_statement_node(
        IRtypes.allocate_statement,
        elementmat,
        IR_allocate_node(IRtypes.float_64_data, [:dofs_per_element, :dofs_per_element])));
    
    elementvec = IR_data_node(IRtypes.array_data, :element_vector);
    push!(allocate_block.parts, IR_statement_node(
        IRtypes.allocate_statement,
        elementvec,
        IR_allocate_node(IRtypes.float_64_data, [:dofs_per_element])));
    
    # coefficient prep
    # LHS volume
    if !(lhs_vol === nothing)
        entities = extract_entities(lhs_vol, multivar);
        append!(vol_entities, entities);
        (allocate, coef) = prepare_coefficient_values(entities, var, dimension, LHS, "volume");
        append!(allocate_block.parts, allocate);
        append!(prepare_block.parts, coef);
    end
    
    # LHS surface
    if !(lhs_surf === nothing)
        entities = extract_entities(lhs_surf, multivar);
        append!(surf_entities, entities);
        (allocate, coef) = prepare_coefficient_values(entities, var, dimension, LHS, "surface");
        append!(allocate_block.parts, allocate);
        append!(prepare_block.parts, coef);
    end
    
    # RHS volume
    if !(rhs_vol === nothing)
        entities = extract_entities(rhs_vol, multivar);
        append!(vol_entities, entities);
        (allocate, coef) = prepare_coefficient_values(entities, var, dimension, RHS, "volume");
        append!(allocate_block.parts, allocate);
        append!(prepare_block.parts, coef);
    end
    
    # RHS surface
    if !(rhs_surf === nothing)
        entities = extract_entities(rhs_surf, multivar);
        append!(surf_entities, entities);
        (allocate, coef) = prepare_coefficient_values(entities, var, dimension, RHS, "surface");
        append!(allocate_block.parts, allocate);
        append!(prepare_block.parts, coef);
    end
    
    # derivative matrix prep
    (allocate, deriv) = prepare_derivative_matrices(vol_entities, surf_entities, var);
    append!(allocate_block.parts, allocate);
    derivmat_block.parts = deriv;
    
    # computation
    if !(lhs_vol === nothing)
        lhsvol_terms = process_terms(lhs_vol);
        push!(matrix_block.parts, make_elemental_computation_fem(lhsvol_terms, var, dofsper, offset_ind, LHS, "volume"));
    end
    if !(lhs_surf === nothing) 
        lhssurf_terms = process_terms(lhs_surf);
        push!(matrix_block.parts, make_elemental_computation_fem(lhssurf_terms, var, dofsper, offset_ind, LHS, "surface"));
    end
    if !(rhs_vol === nothing)
        rhsvol_terms = process_terms(rhs_vol);
        push!(vector_block.parts, make_elemental_computation_fem(rhsvol_terms, var, dofsper, offset_ind, RHS, "volume"));
    end
    if !(rhs_surf === nothing)
        rhssurf_terms = process_terms(rhs_surf);
        push!(vector_block.parts, make_elemental_computation_fem(rhssurf_terms, var, dofsper, offset_ind, RHS, "surface"));
    end
    
    # assembly loop
    indices = [];
    # check for indexed variables TODO
    assembly_loop = generate_assembly_loop_fem(var, indices);
    
    # fill the loop
    assembly_loop.body = IR_block_node([
        prepare_block,
        derivmat_block,
        matrix_block,
        vector_block,
        bdry_block
    ])
    
    # time loop
    need_time_loop = false; # TODO
    if need_time_loop
        #TODO
    else
        time_loop = IR_block_node([
            assembly_loop, # assembly loop
            IR_eval_node(IRtypes.special_eval, :SOLVE, []) # solve
        ])
    end
    
    # Is this needed?
    finalize_block = IR_block_node([]);
    
    # Put them all together in a master block
    master_block = IR_block_node([
        allocate_block,
        time_loop,
        finalize_block
    ]);
    
    return master_block;
end

function prepare_derivative_matrices(entities_vol, entities_surf, var)
    IRtypes = IR_entry_types();
    deriv_part = Vector{IR_part}(undef,0); # Derivative matrix building inside elemental loop
    allocate_part = Vector{IR_part}(undef,0); # Allocations that are done outside of the loops
    
    needed_derivative_matrices = fill(false, 8); # 1,2,3 = x,y,z quadrature points, 5,6,7 = nodes
    
    # Loop over entities to check for derivatives
    for i=1:length(entities_vol)
        # Only check for derivatives
        if length(entities_vol[i].derivs) > 0
            for di=1:length(entities_vol[i].derivs)
                needed_derivative_matrices[entities_vol[i].derivs[di]] = true;
            end
        end
    end
    
    for i=1:length(entities_surf)
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
            push!(allocate_part, IR_statement_node(
                IRtypes.allocate_statement,
                IR_data_node(IRtypes.array_data, Symbol("RQ"*string(i))),
                IR_allocate_node(IRtypes.float_64_data, [:qnodes_per_element, :nodes_per_element])));
        end
        if needed_derivative_matrices[i+4]
            push!(allocate_part, IR_statement_node(
                IRtypes.allocate_statement,
                IR_data_node(IRtypes.array_data, Symbol("RD"*string(i))),
                IR_allocate_node(IRtypes.float_64_data, [:nodes_per_element, :nodes_per_element])));
        end
        
        # Build derivative matrices
        if needed_derivative_matrices[i]
            push!(deriv_part, 
                IR_statement_node(IRtypes.call_statement, nothing, IR_eval_node(IRtypes.func_eval, :build_derivative_matrix, [i, :eid, 0, Symbol("RQ"*string(i))])));
        end
        if needed_derivative_matrices[i+4]
            push!(deriv_part, 
                IR_statement_node(IRtypes.call_statement, nothing, IR_eval_node(IRtypes.func_eval, :build_derivative_matrix, [i, :eid, 1, Symbol("RD"*string(i))])));
        end
    end
    
    return (allocate_part, deriv_part);
end

# Allocate, compute, or fetch all needed values
function prepare_coefficient_values(entities, var, dimension, lorr, vors)
    # # First, gather some needed numbers
    # dimension = size(grid_data.allnodes,2);
    # dofs_per_node = 0;
    # dofs_per_loop = 0;
    # for vi=1:length(var)
    #     dofs_per_loop += length(var[vi].symvar);
    #     dofs_per_node += var[vi].total_components;
    # end
    # nnodes_partition = size(grid_data.allnodes,2);
    # nnodes_global = grid_data.nnodes_global;
    # dofs_global = dofs_per_node * nnodes_global;
    # dofs_partition = dofs_per_node * nnodes_partition;
    # nodes_per_element = refel.Np;
    # dofs_per_element = dofs_per_node * nodes_per_element;
    # num_elements = size(grid_data.loc2glb,2);
    
    IRtypes = IR_entry_types();
    row_col_matrix_index = IR_eval_node(IRtypes.special_eval, :ROWCOL_TO_INDEX, [:row, :col, :nnodes]);
    
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
    
    unique_entity_names = []; # avoid duplicate names
    
    # This loop evaluates coefficient functions INSIDE THE ELEMENTAL LOOP
    # If these are to be precomputed, this must change.
    # First get x,y,z,t,nodeID
    coef_loop_body = [
        # nodeID = mesh_loc2glb[eid, ni]
        IR_statement_node(IRtypes.assign_statement, :nodeID, 
            IR_data_node(IRtypes.array_data, :mesh_loc2glb, [:eid, :ni])),
        # t = the current step time
        # IR_statement_node(IRtypes.assign_statement, :t, 0.0),
        # x = mesh_allnodes[nodeID*dimension]
        IR_statement_node(IRtypes.assign_statement, :x, 
            IR_data_node(IRtypes.array_data, :mesh_allnodes, [
                IR_eval_node(IRtypes.math_eval, :(*), [:nodeID, :dimension])
            ]))
    ];
    if dimension > 1
        # y = mesh_allnodes[nodeID*dimension+1]
        push!(coef_loop_body, IR_statement_node(IRtypes.assign_statement, :y, 
            IR_data_node(IRtypes.array_data, :mesh_allnodes, [
                IR_eval_node(IRtypes.math_eval, :(+), [IR_eval_node(IRtypes.math_eval, :(*), [:nodeID, :dimension]), 1])
            ])))
    else
        # y = 0
        push!(coef_loop_body, IR_statement_node(IRtypes.assign_statement, :y, 0));
    end
    if dimension > 2
        # z = mesh_allnodes[nodeID*dimension+2]
        push!(coef_loop_body, IR_statement_node(IRtypes.assign_statement, :z, 
            IR_data_node(IRtypes.array_data, :mesh_allnodes, [
                IR_eval_node(IRtypes.math_eval, :(+), [IR_eval_node(IRtypes.math_eval, :(*), [:nodeID, :dimension]), 2])
            ])))
    else
        # z = 0
        push!(coef_loop_body, IR_statement_node(IRtypes.assign_statement, :z, 0));
    end
    
    # Loop over entities to perpare for each one
    for i=1:length(entities)
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
                
        if is_test_function(entities[i])
            # ?
        elseif is_unknown_var(entities[i], var) && lorr == LHS
            # ?
        else  # Is a coefficient(number or function) or variable(array)?
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
                # Need to allocate NP 
                np_allocate = IR_allocate_node(IRtypes.float_64_data, [:nodes_per_element]);
                nqp_allocate = IR_allocate_node(IRtypes.float_64_data, [:qnodes_per_element]);
                nodal_coef_name = "NODAL"*cname;
                # For the nodal values
                push!(allocate_part, IR_statement_node(IRtypes.allocate_statement,
                    IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name)),
                    np_allocate));
                # For the interpolated/differentiated quadrature values
                push!(allocate_part, IR_statement_node(IRtypes.allocate_statement,
                    IR_data_node(IRtypes.array_data, Symbol(cname)),
                    nqp_allocate));
                
                coef_index = get_coef_index(entities[i]);
                
                push!(coef_loop_body, IR_statement_node(IRtypes.assign_statement,
                    IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name), [:ni]),
                    IR_eval_node(IRtypes.special_eval, :COEF_EVAL, [coef_index, entities[i].index, :x, :y, :z, :t, :nodeID])));
                
                if vors == "volume"
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        for di=1:length(entities[i].derivs)
                            # TODO
                        end
                    else
                        quad_coef_node = IR_data_node(IRtypes.array_data, Symbol(cname), [:row]);
                        nodal_coef_node = IR_data_node(IRtypes.array_data, Symbol(nodal_coef_name), [:col]);
                        refelQ = IR_data_node(IRtypes.array_data, :refel_Q, [row_col_matrix_index]);
                        push!(coef_init_loop_body, IR_statement_node(IRtypes.assign_statement, quad_coef_node, 0.0));
                        push!(coef_interp_loop_body, IR_statement_node(IRtypes.assign_statement,
                            quad_coef_node,
                            IR_eval_node(IRtypes.math_eval, :(+), [quad_coef_node, IR_eval_node(IRtypes.math_eval, :(*), [refelQ, nodal_coef_node])])));
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
            elseif ctype == 4 # an indexed coefficient
                # TODO
            end
        end # if coefficient
    end # entity loop
    
    # Build the coef eval loop
    coef_eval_loop = IR_loop_node(IRtypes.space_loop, :nodes, :ni, 0, 0, coef_loop_body);
    
    push!(coef_init_loop_body, IR_loop_node(IRtypes.space_loop, :nodes, :row, 0, 0, coef_interp_loop_body))
    coef_interp_loop = IR_loop_node(IRtypes.space_loop, :qnodes, :col, 0, 0, coef_init_loop_body);
    
    # If there are no needed coefficients, return empty lists
    if length(coef_loop_body) > 0
        push!(coef_part, coef_eval_loop);
    end
    if length(coef_init_loop_body) > 0
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
                            term_IR = IR_eval_node(IRtypes.special_eval, :LINALG_MDM, [test_part, coef_part, trial_part]);
                        else
                            term_IR = IR_eval_node(IRtypes.special_eval, :LINALG_Mv, [test_part, coef_part]);
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
                        term_IR = IR_eval_node(IRtypes.special_eval, :LINALG_MDM, [test_part, coef_part, trial_part]);
                    else
                        term_IR = IR_eval_node(IRtypes.special_eval, :LINALG_Mv, [test_part, coef_part]);
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
        # Let's let the code generator decide how to do this
        if lorr == LHS
            for smi=1:dofsper
                for smj=1:dofsper
                    submat_ind = smj + (smi-1)*dofsper;
                    if length(submatrices[i,j]) > 0
                        submat_rhs = IR_eval_node(IRtypes.math_eval, :+, submatrices[i,j]);
                        push!(compute_block.parts, IR_statement_node(IRtypes.special_statement, :LINALG_MATRIX_BLOCK, 
                                                    [smi, smj, :nodes_per_element, :element_matrix, submat_rhs]));
                    end
                end
            end
        else
            for smj=1:dofsper
                if length(submatrices[smj]) > 0
                    submat_rhs = IR_eval_node(IRtypes.math_eval, :+, submatrices[smj]);
                    push!(compute_block.parts, IR_statement_node(IRtypes.special_statement, :LINALG_VECTOR_BLOCK, 
                                                [smj, :nodes_per_element, :element_vector, submat_rhs]));
                end
            end
        end
        
        
    else # one dof
        terms = terms[1];
        term_vec = Vector{IR_part}(undef,0);
        
        #process each term
        for i=1:length(terms)
            (test_part, trial_part, coef_part, test_ind, trial_ind) = generate_term_calculation_fem(terms[i], var, lorr, vors);
            
            # Turn these three parts into an expression like A'DB or A'Dv = A'd
            # Where D is a diagonal matrix specified by a vector in the IR
            if lorr == LHS
                term_IR = IR_eval_node(IRtypes.special_eval, :LINALG_MDM, [test_part, coef_part, trial_part]);
            else
                term_IR = IR_eval_node(IRtypes.special_eval, :LINALG_Mv, [test_part, coef_part]);
            end
            
            push!(term_vec, term_IR);
            
        end
        if length(term_vec) > 0
            combined_rhs = IR_eval_node(IRtypes.math_eval, :+, term_vec);
            if lorr == LHS
                push!(compute_block.parts, IR_statement_node(IRtypes.call_statement, nothing, 
                        IR_eval_node(IRtypes.special_eval, :LINALG_MATRIX_BLOCK, [1, 1, :nodes_per_element, :element_matrix, combined_rhs])));
            else
                push!(compute_block.parts, IR_statement_node(IRtypes.call_statement, nothing, 
                        IR_eval_node(IRtypes.special_eval, :LINALG_VECTOR_BLOCK, [1, :nodes_per_element, :element_vector, combined_rhs])));
            end
        else
            # This shouldn't be possible?
        end
        
    end
    
    return compute_block;
end

# This takes a term expression that should have a form like test_part * (coef_parts) * trial_part
# The parts are separated, test and trial parts are translated into quadrature matrices (refel_Q or RQn)
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
            test_part = IR_data_node(IRtypes.matrix_data, Symbol("RQ"*string(deriv_index)), [:nqnodes, :nnodes]);
        else
            test_part = IR_data_node(IRtypes.matrix_data, Symbol("refel_Q"), [:nqnodes, :nnodes]);
        end
    else
        test_part = nothing;
    end
    if typeof(trial_ex) == SymEntity
        if length(trial_ex.derivs) > 0
            deriv_index = trial_ex.derivs[1];
            trial_part = IR_data_node(IRtypes.matrix_data, Symbol("RQ"*string(deriv_index)), [:nqnodes, :nnodes]);
        else
            trial_part = IR_data_node(IRtypes.matrix_data, Symbol("refel_Q"), [:nqnodes, :nnodes]);
        end
    else
        trial_part = nothing;
    end
    
    # Turn the coefficient part into IR
    if !(coef_ex === nothing)
        coef_part = arithmetic_expr_to_IR(coef_ex);
        if trial_negative || test_negative
            coef_part = IR_eval_node(IRtypes.math_eval, :*, [:refel_wg, :detj, coef_part, -1]);
        else
            coef_part = IR_eval_node(IRtypes.math_eval, :*, [:refel_wg, :detj, coef_part]);
        end
        
    else
        if trial_negative || test_negative
            coef_part = IR_eval_node(IRtypes.math_eval, :*, [:refel_wg, :detj, -1]);
        else
            coef_part = IR_eval_node(IRtypes.math_eval, :*, [:refel_wg, :detj]);
        end
        
    end
    
    return (test_part, trial_part, coef_part, test_ind, trial_ind);
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
            assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :eid, 0, 0, placeholder);
        else
            assembly_loop = IR_loop_node(IRtypes.index_loop, indices[end].symbol, 
                            index_names[end].var, indices[end].range[1], indices[end].range[end], placeholder);
        end
        
        # work outwards nesting assembly_loop
        for i=(length(index_names)-1):-1:1
            if index_names[end] == "elements"
                assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :eid, 0, 0, assembly_loop);
            else
                assembly_loop = IR_loop_node(IRtypes.index_loop, indices[i-ind_shift].symbol, 
                                index_names[i].var, indices[i-ind_shift].range[1], 
                                indices[i-ind_shift].range[i-ind_shift], assembly_loop);
            end
        end
    else # only an element loop
        assembly_loop = IR_loop_node(IRtypes.space_loop, :elements, :eid, 0, 0, placeholder);
    end
    
    return assembly_loop
end
