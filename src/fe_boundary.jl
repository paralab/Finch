
# Apply boundary conditions to one element
function apply_boundary_conditions_elemental(var::Union{Variable, Vector{Variable}}, eid::Int, grid::Grid, refel::Refel,
                                            geo_facs::GeometricFactors, prob::Finch_prob, t::Union{Int,Float64},
                                            elmat::Matrix{Float64}, elvec::Vector{Float64}, bdry_done::Vector{Bool})
    # Check each node to see if the bid is > 0 (on boundary)
    nnodes = refel.Np;
    norm_dot_grad = nothing;
    for ni=1:nnodes
        node_id = grid.loc2glb[ni,eid];
        node_bid = grid.nodebid[node_id];
        face_id = -1; # This may need to be figured out, but it is not clear which face for vertices
        if node_bid > 0
            # This is a boundary node in node_bid
            # Handle the BC for each variable
            if typeof(var) <: Array;
                row_index = ni;
                col_offset = 0;
                for vi=1:length(var)
                    for compo=1:var[vi].total_components
                        bc_type = prob.bc_type[var[vi].index, node_bid];
                        if bc_type == NO_BC
                            # do nothing
                        elseif bc_type == DIRICHLET
                            # zero the row
                            for nj=1:size(elmat,2)
                                elmat[row_index, nj] = 0;
                            end
                            elvec[row_index] = 0;
                            if !bdry_done[node_id]
                                # elmat row is identity
                                elmat[row_index, row_index] = 1;
                                # elvec row is value
                                elvec[row_index] = evaluate_at_node(prob.bc_func[var[vi].index, node_bid][compo], node_id, face_id, t, grid);
                            end
                            
                        elseif bc_type == NEUMANN
                            # zero the row
                            for nj=1:size(elmat,2)
                                elmat[row_index, nj] = 0;
                            end
                            elvec[row_index] = 0;
                            if !bdry_done[node_id]
                                # elmat row is grad dot norm
                                if norm_dot_grad === nothing
                                    # This only needs to be made once per element
                                    norm_dot_grad = get_norm_dot_grad(eid, grid, refel, geo_facs);
                                end
                                for nj = 1:nnodes
                                    elmat[row_index, col_offset+nj] = norm_dot_grad[ni, nj];
                                end
                                # elvec row is value
                                elvec[row_index] = evaluate_at_node(prob.bc_func[var[vi].index, node_bid][compo], node_id, face_id, t, grid);
                            end
                        elseif bc_type == ROBIN
                            printerr("Robin BCs not ready.");
                        else
                            printerr("Unsupported boundary condition type: "*bc_type);
                        end
                        
                        row_index += nnodes;
                        col_offset += nnodes;
                    end
                end
            else
                row_index = ni;
                col_offset = 0;
                for compo=1:var.total_components
                    bc_type = prob.bc_type[var.index, node_bid];
                    if bc_type == NO_BC
                        # do nothing
                    elseif bc_type == DIRICHLET
                        # zero the row
                        for nj=1:size(elmat,2)
                            elmat[row_index, nj] = 0;
                        end
                        elvec[row_index] = 0;
                        if !bdry_done[node_id]
                            # elmat row is identity
                            elmat[row_index, row_index] = 1;
                            # elvec row is value
                            elvec[row_index] = evaluate_at_node(prob.bc_func[var.index, node_bid][compo], node_id, face_id, t, grid);
                        end
                    elseif bc_type == NEUMANN
                        # zero the row
                        for nj=1:size(elmat,2)
                            elmat[row_index, nj] = 0;
                        end
                        elvec[row_index] = 0;
                        if !bdry_done[node_id]
                            # elmat row is grad dot norm
                            if norm_dot_grad === nothing
                                # This only needs to be made once per element
                                norm_dot_grad = get_norm_dot_grad(eid, grid, refel, geo_facs);
                            end
                            for nj = 1:nnodes
                                elmat[row_index, col_offset+nj] = norm_dot_grad[ni, nj];
                            end
                            # elvec row is value
                            elvec[row_index] = evaluate_at_node(prob.bc_func[var.index, node_bid][compo], node_id, face_id, t, grid);
                        end
                    elseif bc_type == ROBIN
                        printerr("Robin BCs not ready.");
                    else
                        printerr("Unsupported boundary condition type: "*bc_type);
                    end
                    
                    row_index += nnodes;
                    col_offset += nnodes;
                end
            end
            
            # Set the flag to done
            bdry_done[node_id] = true;
        end
    end
end

# Apply boundary conditions to one element
function apply_boundary_conditions_elemental_rhs(var::Union{Variable, Vector{Variable}}, eid::Int, grid::Grid, refel::Refel,
                                            geo_facs::GeometricFactors, prob::Finch_prob, t::Union{Int,Float64},
                                            elvec::Vector{Float64}, bdry_done::Vector{Bool})
    # Check each node to see if the bid is > 0 (on boundary)
    nnodes = refel.Np;
    for ni=1:nnodes
        node_id = grid.loc2glb[ni,eid];
        node_bid = grid.nodebid[node_id];
        face_id = -1; # This may need to be figured out, but it is not clear which face for vertices
        if node_bid > 0
            # This is a boundary node in node_bid
            # Handle the BC for each variable
            if typeof(var) <: Array;
                row_index = ni;
                for vi=1:length(var)
                    for compo=1:var[vi].total_components
                        bc_type = prob.bc_type[var[vi].index, node_bid];
                        if bc_type == NO_BC
                            # do nothing
                        elseif bc_type == DIRICHLET
                            elvec[row_index] = 0;
                            if !bdry_done[node_id]
                                # elvec row is value
                                elvec[row_index] = evaluate_at_node(prob.bc_func[var[vi].index, node_bid][compo], node_id, face_id, t, grid);
                            end
                            
                        elseif bc_type == NEUMANN
                            elvec[row_index] = 0;
                            if !bdry_done[node_id]
                                # elvec row is value
                                elvec[row_index] = evaluate_at_node(prob.bc_func[var[vi].index, node_bid][compo], node_id, face_id, t, grid);
                            end
                        elseif bc_type == ROBIN
                            printerr("Robin BCs not ready.");
                        else
                            printerr("Unsupported boundary condition type: "*bc_type);
                        end
                        
                        row_index += nnodes;
                    end
                end
            else
                row_index = ni;
                for compo=1:var.total_components
                    bc_type = prob.bc_type[var.index, node_bid];
                    if bc_type == NO_BC
                        # do nothing
                    elseif bc_type == DIRICHLET
                        elvec[row_index] = 0;
                        if !bdry_done[node_id]
                            # elvec row is value
                            elvec[row_index] = evaluate_at_node(prob.bc_func[var.index, node_bid][compo], node_id, face_id, t, grid);
                        end
                    elseif bc_type == NEUMANN
                        elvec[row_index] = 0;
                        if !bdry_done[node_id]
                            # elvec row is value
                            elvec[row_index] = evaluate_at_node(prob.bc_func[var.index, node_bid][compo], node_id, face_id, t, grid);
                        end
                    elseif bc_type == ROBIN
                        printerr("Robin BCs not ready.");
                    else
                        printerr("Unsupported boundary condition type: "*bc_type);
                    end
                    
                    row_index += nnodes;
                end
            end
            
            # Set the flag to done
            bdry_done[node_id] = true;
        end
    end
end

#############################################################################################################
# Utilities

# Returns the matrix operator for norm dot grad for each face node (other rows are zero)
# For edge/vertex nodes, the bid assigned to the node determines the face, but if
# the bid matches multiple faces, the choice is not guaranteed.
function get_norm_dot_grad(eid::Int, grid::Grid, refel::Refel, geometric_factors::GeometricFactors)
    nnodes = refel.Np;
    dim = refel.dim;
    # Will need all the the Ddr derivative matrices for this element
    RD1 = zeros(Float64, nnodes, nnodes);
    build_derivative_matrix(refel, geometric_factors, 1, eid, 1, RD1);
    if dim > 1
        RD2 = zeros(Float64, nnodes, nnodes);
        build_derivative_matrix(refel, geometric_factors, 2, eid, 1, RD2);
    end
    if dim > 2
        RD3 = zeros(Float64, nnodes, nnodes);
        build_derivative_matrix(refel, geometric_factors, 3, eid, 1, RD3);
    end
    
    # The result will be this
    ndotgrad = zeros(nnodes,nnodes);
    
    # Match each node to a face
    nid = zeros(Int, nnodes);
    nbid = zeros(Int, nnodes);
    # First fill this with the bid of each node
    for ni=1:nnodes
        nid[ni] = grid.loc2glb[ni, eid];
        nbid[ni] = grid.nodebid[nid[ni]];
    end
    # Loop over faces of this element
    for fi=1:refel.Nfaces
        fid = grid.element2face[fi, eid];
        fbid = grid.facebid[fid];
        if fbid > 0
            fnorm = grid.facenormals[:,fid];
            for fni = 1:refel.Nfp[fi] # for nodes on this face
                for ni=1:nnodes # for nodes in this element
                    if nid[ni] == grid.face2glb[fni, 1, fid]
                        # This ni corresponds to this fni
                        # make sure bid is the same
                        if fbid == nbid[ni]
                            # It's a match
                            # row ni in the matrix will be set
                            for col = 1:nnodes
                                if dim == 1
                                    ndotgrad[ni, col] = RD1[ni, col] * fnorm[1];
                                elseif dim == 2
                                    ndotgrad[ni, col] = RD1[ni, col] * fnorm[1] + RD2[ni, col] * fnorm[2];
                                elseif dim == 3
                                    ndotgrad[ni, col] = RD1[ni, col] * fnorm[1] + RD2[ni, col] * fnorm[2] + RD3[ni, col] * fnorm[3];
                                end
                            end
                        end
                        break;
                    end
                end
            end
        end
    end
    
    return ndotgrad;
end

function copy_bdry_vals_to_variables(var, vec, grid, dofs_per_node; zero_vals=true)
    if typeof(var) <: Array
        dofind = 0;
        for vi=1:length(var)
            for compo=1:length(var[vi].symvar)
                dofind = dofind + 1;
                for bid=1:size(prob.bc_type,2)
                    if prob.bc_type[var[vi].index, bid] == DIRICHLET
                        for i = 1:length(grid.bdryface[bid]) # loop over faces with this BID
                            fid = grid.bdryface[bid][i];
                            face_nodes = grid.face2glb[:,1,fid];
                            for ni=1:length(face_nodes)
                                node = face_nodes[ni];
                                var[vi].values[compo, node] = vec[(node-1)*dofs_per_node + dofind];
                                if zero_vals
                                    vec[(node-1)*dofs_per_node + dofind] = 0;
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        for d=1:dofs_per_node
            dofind = d;
            for bid=1:size(prob.bc_type,2)
                if prob.bc_type[var.index, bid] == DIRICHLET
                    for i = 1:length(grid.bdryface[bid]) # loop over faces with this BID
                        fid = grid.bdryface[bid][i];
                        face_nodes = grid.face2glb[:,1,fid];
                        for ni=1:length(face_nodes)
                            node = face_nodes[ni];
                            var.values[d, node] = vec[(node-1)*dofs_per_node + dofind];
                            if zero_vals
                                vec[(node-1)*dofs_per_node + dofind] = 0;
                            end
                        end
                    end
                end
            end
        end
    end
end

function copy_bdry_vals_to_vector(var, vec, grid, dofs_per_node)
    if typeof(var) <: Array
        dofind = 0;
        for vi=1:length(var)
            for compo=1:length(var[vi].symvar)
                dofind = dofind + 1;
                for bid=1:size(prob.bc_type,2)
                    if prob.bc_type[var[vi].index, bid] == DIRICHLET
                        for i = 1:length(grid.bdryface[bid]) # loop over faces with this BID
                            fid = grid.bdryface[bid][i];
                            face_nodes = grid.face2glb[:,1,fid];
                            for ni=1:length(face_nodes)
                                node = face_nodes[ni];
                                vec[(node-1)*dofs_per_node + dofind] = var[vi].values[compo, node];
                            end
                        end
                    end
                end
            end
        end
    else
        for d=1:dofs_per_node
            dofind = d;
            for bid=1:size(prob.bc_type,2)
                if prob.bc_type[var.index, bid] == DIRICHLET
                    for i = 1:length(grid.bdryface[bid]) # loop over faces with this BID
                        fid = grid.bdryface[bid][i];
                        face_nodes = grid.face2glb[:,1,fid];
                        for ni=1:length(face_nodes)
                            node = face_nodes[ni];
                            vec[(node-1)*dofs_per_node + dofind] = var.values[d, node];
                        end
                    end
                end
            end
        end
    end
end

# This evaluates the BC at a specific node.
# That could mean:
# - the value of constant BCs
# - evaluate a genfunction.func for BCs defined by strings->genfunctions
# - evaluate a function for BCs defined by callback functions
function evaluate_at_nodes(val, nodes, face, t)
    N = length(nodes);
    dim = config.dimension;
    if typeof(val) <: Number
        result = fill(val, N);
        
    elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
        result = zeros(N);
        if dim == 1
            for i=1:N
                result[i]=val.value[1].func(grid_data.allnodes[1,nodes[i]],0,0,t,nodes[i],face);
            end
        elseif dim == 2
            for i=1:N
                result[i]=val.value[1].func(grid_data.allnodes[1,nodes[i]],grid_data.allnodes[2,nodes[i]],0,t,nodes[i],face);
            end
        else
            for i=1:N
                result[i]=val.value[1].func(grid_data.allnodes[1,nodes[i]],grid_data.allnodes[2,nodes[i]],grid_data.allnodes[3,nodes[i]],t,nodes[i],face);
            end
        end
        
    elseif typeof(val) == GenFunction
        result = zeros(N);
        if dim == 1
            for i=1:N
                result[i]=val.func(grid_data.allnodes[1,nodes[i]],0,0,t,nodes[i],face);
            end
        elseif dim == 2
            for i=1:N
                result[i]=val.func(grid_data.allnodes[1,nodes[i]],grid_data.allnodes[2,nodes[i]],0,t,nodes[i],face);
            end
        else
            for i=1:N
                result[i]=val.func(grid_data.allnodes[1,nodes[i]],grid_data.allnodes[2,nodes[i]],grid_data.allnodes[3,nodes[i]],t,nodes[i],face);
            end
        end
        
    elseif typeof(val) == CallbackFunction
        result = zeros(N);
        for i=1:N
            #form a dict for the arguments. x,y,z,t are always included
            arg_list = [];
            append!(arg_list, [("x", grid_data.allnodes[1,nodes[i]]),
                        ("y", dim>1 ? grid_data.allnodes[2,nodes[i]] : 0),
                        ("z", dim>2 ? grid_data.allnodes[3,nodes[i]] : 0),
                        ("t", t)]);
            # Add in the requires list
            for r in val.args
                foundit = false;
                # is it an entity?
                if typeof(r) == Variable
                    foundit = true;
                    if size(r.values,1) == 1
                        push!(arg_list, (string(r.symbol), r.values[1,nodes[i]]));
                    else
                        push!(arg_list, (string(r.symbol), r.values[:,nodes[i]]));
                    end
                elseif typeof(r) == Coefficient
                    # TODO: evaluate coefficient at node
                    foundit = true;
                    push!(arg_list, (string(r.symbol), evaluate_at_nodes(r, nodes[i], face, t)));
                elseif typeof(r) == String
                    # This could also be a variable, coefficient, or special thing like normal, 
                    if r in ["x","y","z","t"]
                        foundit = true;
                        # These are already included
                    elseif r == "normal"
                        foundit = true;
                        push!(arg_list, ("normal", grid_data.facenormals[:,face]));
                    else
                        for v in variables
                            if string(v.symbol) == r
                                foundit = true;
                                if size(v.values,1) == 1
                                    push!(arg_list, (r, v.values[1,nodes[i]]));
                                else
                                    push!(arg_list, (r, v.values[:,nodes[i]]));
                                end
                                break;
                            end
                        end
                        for c in coefficients
                            if string(c.symbol) == r
                                foundit = true;
                                push!(arg_list, (r, evaluate_at_nodes(c, nodes[i], face, t)));
                                break;
                            end
                        end
                    end
                else
                    # What else could it be?
                end

                if !foundit
                    # Didn't figure this thing out.
                    printerr("Unknown requirement for callback function "*string(fun)*" : "*string(requires[i]));
                end
            end
            
            # Build the dict
            args = Dict(arg_list);
            
            # call the function
            result[i] = val.func(args);
        end
        
    end
    return result;
end

# This evaluates the BC at one specific node.
# That could mean:
# - the value of constant BCs
# - evaluate a genfunction.func for BCs defined by strings->genfunctions
# - evaluate a function for BCs defined by callback functions
function evaluate_at_node(val, node::Int, face::Int, t::Union{Int, Float64}, grid::Grid)
    if typeof(val) <: Number
        return val;
    end
    
    dim = size(grid.allnodes,1);
    result = 0;
    x = grid.allnodes[1,node];
    y = dim > 1 ? grid.allnodes[2,node] : 0;
    z = dim > 2 ? grid.allnodes[3,node] : 0;
        
    if typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
        result = val.value[1].func(x,y,z,t,node,face);
        
    elseif typeof(val) == GenFunction
        result = val.func(x,y,z,t,node,face);
        
    elseif typeof(val) == CallbackFunction
        #form a dict for the arguments. x,y,z,t are always included
        arg_list = [];
        append!(arg_list, [("x", x), ("y", y), ("z", z), ("t", t)]);
        # Add in the requires list
        for r in val.args
            foundit = false;
            # is it an entity?
            if typeof(r) == Variable
                foundit = true;
                if size(r.values,1) == 1
                    push!(arg_list, (string(r.symbol), r.values[1,node]));
                else
                    push!(arg_list, (string(r.symbol), r.values[:,node]));
                end
            elseif typeof(r) == Coefficient
                # TODO: evaluate coefficient at node
                foundit = true;
                push!(arg_list, (string(r.symbol), evaluate_at_node(r, node, face, t, grid)));
            elseif typeof(r) == String
                # This could also be a variable, coefficient, or special thing like normal, 
                if r in ["x","y","z","t"]
                    foundit = true;
                    # These are already included
                elseif r == "normal"
                    foundit = true;
                    push!(arg_list, ("normal", grid.facenormals[:,face]));
                else
                    for v in variables
                        if string(v.symbol) == r
                            foundit = true;
                            if size(v.values,1) == 1
                                push!(arg_list, (r, v.values[1,node]));
                            else
                                push!(arg_list, (r, v.values[:,node]));
                            end
                            break;
                        end
                    end
                    for c in coefficients
                        if string(c.symbol) == r
                            foundit = true;
                            push!(arg_list, (r, evaluate_at_nodes(c, node, face, t)));
                            break;
                        end
                    end
                end
            else
                # What else could it be?
            end

            if !foundit
                # Didn't figure this thing out.
                printerr("Unknown requirement for callback function "*string(fun)*" : "*string(r));
            end
        end
        
        # Build the dict
        args = Dict(arg_list);
        
        # call the function
        result = val.func(args);
        
    end
    return result;
end


#############################################################################################################
# These are for applying to the global system
# May be removed eventually

# Apply boundary conditions to the system
function apply_boundary_conditions_rhs_only(var, b, t)
    return apply_boundary_conditions_lhs_rhs(var, nothing, b, t)
end

function apply_boundary_conditions_lhs_rhs(var, A, b, t)
    if A === nothing
        rhs_only = true
    else
        rhs_only = false
    end
    multivar = typeof(var) <: Array;
    maxvarindex = 0;
    if multivar
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        var_to_dofs = [];
        for vi=1:length(var)
            tmp = dofs_per_node;
            dofs_per_node += var[vi].total_components;
            push!(var_to_dofs, (tmp+1):dofs_per_node);
            maxvarindex = max(maxvarindex,var[vi].index);
        end
    else
        # one variable
        dofs_per_node = var.total_components;
        maxvarindex = var.index;
    end
    
    # If there is a reason to zero the rhs bdry entries, uncomment this.
    # But this destroys the case NO_BC. 
    #b = zero_rhs_bdry_vals(b, dofs_per_node);
    
    bidcount = length(grid_data.bids); # the number of BIDs
    dirichlet_nodes = zeros(Int, 0); # Nodes corresponding to rows
    dirichlet_rows = zeros(Int, 0); # Dirichlet rows will be identity
    dirichlet_vals = zeros(0); # The values that go in RHS (b)
    neumann_nodes = zeros(Int, 0); # Nodes corresponding to rows
    neumann_rows = zeros(Int, 0); # Neumann rows will be set accordingly
    neumann_vals = zeros(0); # The values that go in RHS (b)
    if !rhs_only
        neumann_Is = zeros(Int,0);
        neumann_Js = zeros(Int,0);
        neumann_Vs = zeros(0);
    end
    if multivar
        dofind = 0;
        for vi=1:length(var)
            for compo=1:var[vi].total_components
                dofind = dofind + 1;
                for bid=1:bidcount
                    if prob.bc_type[var[vi].index, bid] == NO_BC
                        # do nothing
                    elseif prob.bc_type[var[vi].index, bid] == DIRICHLET
                        (tmprows, tmpnodes, tmpvals) = dirichlet_bc(prob.bc_func[var[vi].index, bid][compo], grid_data.bdryface[bid], t, dofind, dofs_per_node);
                        append!(dirichlet_rows, tmprows);
                        append!(dirichlet_vals, tmpvals);
                        append!(dirichlet_nodes, tmpnodes);
                    elseif prob.bc_type[var[vi].index, bid] == NEUMANN
                        if rhs_only
                            (tmprows, tmpnodes, tmpvals) = neumann_bc(prob.bc_func[var[vi].index, bid][compo], grid_data.bdryface[bid], bid, t, dofind, dofs_per_node, rhs_only=true);
                            append!(neumann_rows, tmprows);
                            append!(neumann_vals, tmpvals);
                            append!(neumann_nodes, tmpnodes);
                        else
                            (tmprows, tmpnodes, tmpvals, tmpIs, tmpJs, tmpVs) = neumann_bc(prob.bc_func[var[vi].index, bid][compo], grid_data.bdryface[bid], bid, t, dofind, dofs_per_node);
                            append!(neumann_rows, tmprows);
                            append!(neumann_vals, tmpvals);
                            append!(neumann_nodes, tmpnodes);
                            append!(neumann_Is, tmpIs);
                            append!(neumann_Js, tmpJs);
                            append!(neumann_Vs, tmpVs);
                        end
                        
                    elseif prob.bc_type[var[vi].index, bid] == ROBIN
                        printerr("Robin BCs not ready.");
                    else
                        printerr("Unsupported boundary condition type: "*prob.bc_type[var[vi].index, bid]);
                    end
                end
            end
        end
    else
        for d=1:dofs_per_node
            dofind = d;
            for bid=1:bidcount
                if prob.bc_type[var.index, bid] == NO_BC
                    # do nothing
                elseif prob.bc_type[var.index, bid] == DIRICHLET
                    (tmprows, tmpnodes, tmpvals) = dirichlet_bc(prob.bc_func[var.index, bid][d], grid_data.bdryface[bid], t, dofind, dofs_per_node);
                    append!(dirichlet_rows, tmprows);
                    append!(dirichlet_vals, tmpvals);
                    append!(dirichlet_nodes, tmpnodes);
                elseif prob.bc_type[var.index, bid] == NEUMANN
                    if rhs_only
                        (tmprows, tmpnodes, tmpvals) = neumann_bc(prob.bc_func[var.index, bid][d], grid_data.bdryface[bid], bid, t, dofind, dofs_per_node, rhs_only=true);
                        append!(neumann_rows, tmprows);
                        append!(neumann_vals, tmpvals);
                        append!(neumann_nodes, tmpnodes);
                    else
                        (tmprows, tmpnodes, tmpvals, tmpIs, tmpJs, tmpVs) = neumann_bc(prob.bc_func[var.index, bid][d], grid_data.bdryface[bid], bid, t, dofind, dofs_per_node);
                        append!(neumann_rows, tmprows);
                        append!(neumann_vals, tmpvals);
                        append!(neumann_nodes, tmpnodes);
                        append!(neumann_Is, tmpIs);
                        append!(neumann_Js, tmpJs);
                        append!(neumann_Vs, tmpVs);
                    end
                elseif prob.bc_type[var.index, bid] == ROBIN
                    printerr("Robin BCs not ready.");
                else
                    printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, bid]);
                end
            end
        end
    end
    
    ## Now I have all of the rows, nodes and vals
    # First apply Neumann, then Dirichlet(can overwrite Neumann), then zero parts that are for non-owned nodes if needed.
    if !rhs_only
        if length(neumann_Is)>0
            A = insert_sparse_rows(A, neumann_Is, neumann_Js, neumann_Vs);
        end
        if length(dirichlet_rows)>0
            A = identity_rows(A, dirichlet_rows, length(b));
        end
    end
    for i=1:length(neumann_rows)
        b[neumann_rows[i]] = neumann_vals[i];
    end
    for i=1:length(dirichlet_rows)
        b[dirichlet_rows[i]] = dirichlet_vals[i];
    end
    # For partitioned meshes only
    if config.num_partitions > 1
        rows_to_zero = zeros(Int,0);
        # Zero all entries for boundary nodes that are not owned
        # Base this on the global boundary flags, not local because of annoying cases
        for ni=1:length(grid_data.partition2global)
            if grid_data.global_bdry_index[grid_data.partition2global[ni]] > 0 && grid_data.node_owner[ni] != config.partition_index
                if dofs_per_node == 1
                    b[ni] = 0;
                    push!(rows_to_zero, ni);
                else
                    for di=1:dofs_per_node
                        rowid = (ni-1)*dofs_per_node+di;
                        b[rowid] = 0;
                        push!(rows_to_zero, rowid);
                    end
                end
            end
        end
        
        if !rhs_only && length(rows_to_zero) > 0
            A = identity_rows(A, rows_to_zero, length(b), diagonal_val=0.0);
        end
    end
    
    # Reference points
    if size(prob.ref_point,1) >= maxvarindex
        if multivar
            posind = zeros(Int,0);
            vals = zeros(0);
            for vi=1:length(var)
                if prob.ref_point[var[vi].index,1]
                    eii = prob.ref_point[var[vi].index, 2];
                    tmp = (grid_data.glbvertex[eii[1], eii[2]] - 1)*dofs_per_node + var_to_dofs[vi][1];
                    if length(prob.ref_point[var[vi].index, 3]) > 1
                        tmp = tmp:(tmp+length(prob.ref_point[var[vi].index, 3])-1);
                    end
                    posind = [posind; tmp];
                    vals = [vals; prob.ref_point[var[vi].index, 3]];
                end
            end
            if length(vals) > 0
                if !rhs_only A = identity_rows(A, posind, length(b)); end
                b[posind] = vals;
            end
            
        else
            if prob.ref_point[var.index,1]
                eii = prob.ref_point[var.index, 2];
                posind = (grid_data.glbvertex[eii[1], eii[2]] - 1)*dofs_per_node + 1;
                if length(prob.ref_point[var.index, 3]) > 1
                    posind = posind:(posind+length(prob.ref_point[var[vi].index, 3])-1);
                else
                    posind = [posind];
                end
                if !rhs_only A = identity_rows(A, posind, length(b)); end
                b[posind] = prob.ref_point[var.index, 3];
            end
        end
    end
    
    if rhs_only
        return b;
    else
        return (A, b);
    end
end

# This is just to speed things up (huge improvement)
function is_in_rows(k, rows)
    # rows are sorted
    b = 1;
    e = length(rows);
    c = Int(round((e+1)/2));
    while e-b > 1
        if k == rows[c]
            return true;
        elseif k > rows[c]
            b=c;
        else
            e=c;
        end
        c = Int(round((b+e)/2));
    end
    return k==rows[b] || k==rows[e];
end

# zeros the rows and puts a 1 on the diagonal
function identity_rows(A, urows, N; diagonal_val=1)
    if length(urows) < 1
        return A;
    end
    if issparse(A)
        (I, J, V) = findnz(A);
        eN = length(I);
        newI = zeros(Int,eN);
        newJ = zeros(Int,eN);
        newV = zeros(eN);
        newind = 0;
        (M,N) = size(A);
        
        rows = sort(urows);
        
        # remove duplicates
        newrows = zeros(Int, length(rows));
        newcount = 1;
        newrows[1] = rows[1];
        for i=2:length(rows)
            if rows[i] != rows[i-1] # because sorted
                newcount += 1;
                newrows[newcount] = rows[i];
            end
        end
        rows = newrows[1:newcount];
        
        # Remove bdry row elements
        for k in 1:length(I)
            isbdry = is_in_rows(I[k], rows);
            if !isbdry
                newind = newind+1;
                newI[newind] = I[k];
                newJ[newind] = J[k];
                newV[newind] = V[k];
            end
        end
        
        # Add diogonals
        rn = length(rows);
        if newind + rn <= eN
            for i=1:rn
                newind = newind+1;
                newI[newind] = rows[i];
                newJ[newind] = rows[i];
                newV[newind] = diagonal_val;
            end
            newI = newI[1:newind];
            newJ = newJ[1:newind];
            newV = newV[1:newind];
        else
            toadd = rn - (eN-newind);
            for i=1:(eN-newind)
                newind = newind+1;
                newI[newind] = rows[i];
                newJ[newind] = rows[i];
                newV[newind] = diagonal_val;
            end
            append!(newI, rows[(rn-toadd+1):rn]);
            append!(newJ, rows[(rn-toadd+1):rn]);
            append!(newV, diagonal_val * ones(toadd));
        end
        
        return sparse(newI,newJ,newV, M,N);
        
    else
        for i=1:length(rows)
            A[rows[i],:] .= 0;
            A[rows[i],rows[i]] = diagonal_val;
        end
        return A;
    end
end

# Inserts the rows of S in A
function insert_rows(A, S, rows, N)
    if issparse(A)
        (I, J, V) = findnz(A);
        (M,N) = size(A);
        # determine which elements to remove
        toskip = zeros(Bool, length(I));
        includen = 0;
        for k = 1:length(I)
            for i=1:length(rows)
                if I[k] == rows[i]
                    toskip[k] = 1;
                end
            end
            if toskip[k] == 0
                includen = includen + 1;
            end
        end
        
        # put the remaining elements in a new matrix
        newI = zeros(Int, includen);
        newJ = zeros(Int, includen);
        newV = zeros(includen);
        ind = 1;
        for k = 1:length(I)
            if toskip[k] == 0
                newI[ind] = I[k];
                newJ[ind] = J[k];
                newV[ind] = V[k];
                ind = ind + 1;
            end
        end
        
        # Do something similar to extract elements of S
        (I, J, V) = findnz(S);
        
        toskip = ones(Int, length(I));
        includen = 0;
        for k = 1:length(I)
            for i=1:length(rows)
                if I[k] == rows[i]
                    toskip[k] = 0;
                end
            end
            if toskip[k] == 0
                includen = includen + 1;
            end
        end
        newI2 = zeros(Int, includen);
        newJ2 = zeros(Int, includen);
        newV2 = zeros(includen);
        ind = 1;
        for k = 1:length(I)
            if toskip[k] == 0
                newI2[ind] = I[k];
                newJ2[ind] = J[k];
                newV2[ind] = V[k];
                ind = ind + 1;
            end
        end
        
        append!(newI, newI2);
        append!(newJ, newJ2);
        append!(newV, newV2);
        return sparse(newI,newJ,newV, M, N);
    else
        for i=1:length(rows)
            A[rows[i],:] = S[rows[i],:];
        end
        return A;
    end
end

function insert_sparse_rows(A, SI, SJ, SV)
    (I, J, V) = findnz(A);
    (M,N) = size(A);
    
    # since we need to sort the three arrays in sync, just get the sorted permutation
    sorted_order = sortperm(SI);
    
    SI = SI[sorted_order];
    SJ = SJ[sorted_order];
    SV = SV[sorted_order];
    
    # Zero existing elements in SI rows
    for k in 1:length(I)
        if is_in_rows(I[k], SI)
            V[k] = 0;
        end
    end
    
    # append S values
    append!(I, SI);
    append!(J, SJ);
    append!(V, SV);
    
    # Form sparse matrix
    return sparse(I, J, V, M, N);
end

function zero_rhs_bdry_vals(rhs, dofs)
    for bid=1:length(grid_data.bids)
        for fi=1:length(grid_data.bdryface[bid])
            face_nodes = grid_data.face2glb[:, 1, grid_data.bdryface[bid][fi]];
            for ni=1:length(face_nodes)
                rhs_ind = ((face_nodes[ni]-1)*dofs + 1):(face_nodes[ni]*dofs);
                rhs[rhs_ind] .= 0;
            end
        end
    end
    return rhs;
end

function dirichlet_bc(val, bdryface, t=0, dofind = 1, totaldofs = 1)
    # Assuming every face has the same number of nodes
    nodes_per_face = size(grid_data.face2glb,1);
    nnodes = length(bdryface)*nodes_per_face;
    bdry_rows = zeros(Int, nnodes);
    bdry_nodes = zeros(Int, nnodes);
    bdry_vals = zeros(nnodes);
    
    for fi = 1:length(bdryface)
        if bdryface[fi] < 1
            println("found bad bdryface: "*string(bdryface));
        end
        # global indices for the nodes on the face
        face_nodes = grid_data.face2glb[:, 1, bdryface[fi]];
        
        # offset for multi dof
        if totaldofs > 1
            offsetglb = (face_nodes.-1) .* totaldofs .+ dofind;
        else
            offsetglb = face_nodes;
        end
        
        bdry_rows[((fi-1)*nodes_per_face+1):(fi*nodes_per_face)] = offsetglb;
        bdry_nodes[((fi-1)*nodes_per_face+1):(fi*nodes_per_face)] = face_nodes;
        bdry_vals[((fi-1)*nodes_per_face+1):(fi*nodes_per_face)] = evaluate_at_nodes(val, face_nodes, bdryface[fi], t);
    end
    
    return (bdry_rows, bdry_nodes, bdry_vals);
end

function neumann_bc(val, bdryface, bid, t=0, dofind = 1, totaldofs = 1; rhs_only=false)
    # Assuming every element has the same number of nodes
    nodes_per_element = size(grid_data.loc2glb,1);
    nodes_per_face = size(grid_data.face2glb,1);
    nnodes = length(bdryface)*nodes_per_face;
    bdry_rows = zeros(Int, nnodes);
    bdry_nodes = zeros(Int, nnodes);
    bdry_vals = zeros(nnodes);
    
    # To deal with the possibility of overlapping nodes between elements or other bids
    # limit it to nodes assigned to this BID and assign each node to one of the bdryfaces.
    # Need to form a list of valid nodes for each face in bdryface.
    bidnodes = grid_data.bdry[bid];
    node_available = fill(true, length(bidnodes));
    valid_face_mask = zeros(Bool, nodes_per_face, length(bdryface));
    for fi=1:length(bdryface)
        face_nodes = grid_data.face2glb[:, 1, bdryface[fi]];
        for ni=1:length(face_nodes)
            bnode_ind = indexin([face_nodes[ni]], bidnodes);
            if !(bnode_ind[1] === nothing)
                if node_available[bnode_ind[1]]
                    # This is a valid node for this face. Reserve it
                    node_available[bnode_ind[1]] = false;
                    valid_face_mask[ni,fi] = true;
                end
            end
        end
    end
    total_entries = length(bidnodes) * nodes_per_element;
    
    if !rhs_only
        # These will be sparse rows to insert into A
        S_I = zeros(Int, total_entries);
        S_J = zeros(Int, total_entries);
        S_V = zeros(total_entries);
        next_S_ind = 1;
        
        # This can be precomputed for uniform grid meshes
        if config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
            precomputed = true;
            glb = grid_data.loc2glb[:,1];
            detJ = geo_factors.detJ[1];
            J = geo_factors.J[1]
            if config.dimension == 1
                RD1 = zeros(size(refel.Ddr));
                # Multiply rows of Qr and Ddr by J.rx
                for i=1:size(RD1,2) # loop over columns
                    RD1[:,i] = J.rx .* refel.Ddr[:,i];
                end
                
            elseif config.dimension == 2
                RD1 = zeros(size(refel.Ddr));
                RD2 = zeros(size(refel.Ddr));
                for i=1:size(RD1,2)
                    RD1[:,i] = J.rx .* refel.Ddr[:,i] + J.sx .* refel.Dds[:,i];
                    RD2[:,i] = J.ry .* refel.Ddr[:,i] + J.sy .* refel.Dds[:,i];
                end
                
            elseif config.dimension == 3
                RD1 = zeros(size(refel.Ddr));
                RD2 = zeros(size(refel.Ddr));
                RD3 = zeros(size(refel.Ddr));
                for i=1:size(RD1,2)
                    RD1[:,i] = J.rx .* refel.Ddr[:,i] + J.sx .* refel.Dds[:,i] + J.tx .* refel.Ddt[:,i];
                    RD2[:,i] = J.ry .* refel.Ddr[:,i] + J.sy .* refel.Dds[:,i] + J.ty .* refel.Ddt[:,i];
                    RD3[:,i] = J.rz .* refel.Ddr[:,i] + J.sz .* refel.Dds[:,i] + J.tz .* refel.Ddt[:,i];
                end
            end
        else
            precomputed = false;
        end
    end
    
    next_row_index = 1;
    for fi = 1:length(bdryface)
        # Local indices for the nodes on the face
        flocal = refel.face2local[grid_data.faceRefelInd[1,bdryface[fi]]];
        
        # global indices for the nodes on the face
        face_nodes = grid_data.face2glb[:, 1, bdryface[fi]];
        
        # remove nodes that do not belong to this face and bid
        newflocal = zeros(Int, length(flocal));
        newface_nodes = zeros(Int, length(face_nodes));
        node_count = 0;
        for ni=1:length(face_nodes)
            if valid_face_mask[ni, fi]
                node_count += 1;
                newflocal[node_count] = flocal[ni];
                newface_nodes[node_count] = face_nodes[ni];
            end
        end
        face_nodes = newface_nodes[1:node_count];
        flocal = newflocal[1:node_count];
        
        # offset for multi dof
        if totaldofs > 1
            offsetface = (face_nodes.-1) .* totaldofs .+ dofind;
        else
            offsetface = face_nodes;
        end
        
        if !rhs_only
            # find the relevant element
            e = grid_data.face2element[1,bdryface[fi]] # >0 ? grid_data.face2element[1,bdryface[fi]] : grid_data.face2element[2,bdryface[fi]];
            glb = grid_data.loc2glb[:,e];       # global indices of this element's nodes
            # offset for multi dof
            if totaldofs > 1
                elglb = (glb.-1) .* totaldofs .+ dofind;
            else
                elglb = glb;
            end
            
            # Build diff matrices if needed
            if !precomputed
                detJ = geo_factors.detJ[e];
                J = geo_factors.J[e]
                if config.dimension == 1
                    RD1 = zeros(size(refel.Ddr));
                    # Multiply rows of Qr and Ddr by J.rx
                    for i=1:size(RD1,2) # loop over columns
                        RD1[:,i] = J.rx .* refel.Ddr[:,i];
                    end
                    
                elseif config.dimension == 2
                    RD1 = zeros(size(refel.Ddr));
                    RD2 = zeros(size(refel.Ddr));
                    for i=1:size(RD1,2)
                        RD1[:,i] = J.rx .* refel.Ddr[:,i] + J.sx .* refel.Dds[:,i];
                        RD2[:,i] = J.ry .* refel.Ddr[:,i] + J.sy .* refel.Dds[:,i];
                    end
                    
                elseif config.dimension == 3
                    RD1 = zeros(size(refel.Ddr));
                    RD2 = zeros(size(refel.Ddr));
                    RD3 = zeros(size(refel.Ddr));
                    for i=1:size(RD1,2)
                        RD1[:,i] = J.rx .* refel.Ddr[:,i] + J.sx .* refel.Dds[:,i] + J.tx .* refel.Ddt[:,i];
                        RD2[:,i] = J.ry .* refel.Ddr[:,i] + J.sy .* refel.Dds[:,i] + J.ty .* refel.Ddt[:,i];
                        RD3[:,i] = J.rz .* refel.Ddr[:,i] + J.sz .* refel.Dds[:,i] + J.tz .* refel.Ddt[:,i];
                    end
                end
            end
            
            # Add to S matrix
            for i=1:length(elglb)
                irange = (next_S_ind + (i-1)*node_count):(next_S_ind + i*node_count - 1);
                S_I[irange] = elglb[flocal];
                S_J[irange] .= elglb[i];
            end
            if config.dimension == 1
                S1_V= RD1[flocal,:][:];
                
            elseif config.dimension == 2
                S1_V= RD1[flocal,:][:];
                S2_V= RD2[flocal,:][:];
                
            elseif config.dimension == 3
                S1_V= RD1[flocal,:][:];
                S2_V= RD2[flocal,:][:];
                S3_V= RD3[flocal,:][:];
            end
            
            # Add the right components of S1,S2,S3 according to normal vector
            ind_range = next_S_ind:(next_S_ind + node_count * nodes_per_element - 1);
            norm = grid_data.facenormals[:, bdryface[fi]];
            # println("bid "*string(bid)*" normal "*string(norm));
            if config.dimension == 1
                S_V[ind_range] = norm[1] .* S1_V;
            elseif config.dimension == 2
                S_V[ind_range] = norm[1] .* S1_V + norm[2] .* S2_V;
            elseif config.dimension == 3
                S_V[ind_range] = norm[1] .* S1_V + norm[2] .* S2_V + norm[3] .* S3_V;
            end
            
            # Increase S index
            next_S_ind += node_count * nodes_per_element;
        end
        
        bdry_nodes[next_row_index:(next_row_index+length(face_nodes)-1)] = face_nodes;
        bdry_rows[next_row_index:(next_row_index+length(face_nodes)-1)] = offsetface;
        bdry_vals[next_row_index:(next_row_index+length(face_nodes)-1)] = evaluate_at_nodes(val, face_nodes, bdryface[fi], t);
        next_row_index += length(face_nodes);
    end
    
    bdry_nodes = bdry_nodes[1:(next_row_index-1)];
    bdry_rows = bdry_nodes[1:(next_row_index-1)];
    bdry_vals = bdry_nodes[1:(next_row_index-1)];
    
    if rhs_only
        return (bdry_rows, bdry_nodes, bdry_vals);
    else
        if next_S_ind <= length(S_I)
            # This should not happen
            log_entry("Unexpected entry count when doing Neumann BC. expected "*string(length(S_I))*" found "*string(next_S_ind-1));
            S_I = S_I[1:next_S_ind-1];
            S_J = S_J[1:next_S_ind-1];
            S_V = S_V[1:next_S_ind-1];
        end
        return (bdry_rows, bdry_nodes, bdry_vals, S_I, S_J, S_V);
    end
    
end
