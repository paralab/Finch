#=
Utilities used by the generated solve() functions.
These will be called by the generated solve code
so they should be made as efficient as possible.
=#

# place the values from sol into the variable value arrays
function place_sol_in_vars(var::Vector{Variable}, sol::Vector{Float64})
    tmp = 0;
    totalcomponents = 0;
    for vi=1:length(var)
        totalcomponents = totalcomponents + var[vi].total_components;
    end
    for vi=1:length(var)
        components = var[vi].total_components;
        for compi=1:components
            var[vi].values[compi,:] .= sol[(compi+tmp):totalcomponents:end];
        end
        tmp = tmp + components;
    end
    
    return nothing;
end

# Set the values from variable arrays in a global vector
function get_var_vals(var::Vector{Variable}, vect::Union{Nothing, Vector{Float64}}=nothing)
    tmp = 0;
    totalcomponents = 0;
    for vi=1:length(var)
        totalcomponents = totalcomponents + var[vi].total_components;
    end
    if vect === nothing
        vect = zeros(totalcomponents * size(var[1].values, 2));
    end
    for vi=1:length(var)
        components = var[vi].total_components;
        for compi=1:components
            vect[(compi+tmp):totalcomponents:end] = var[vi].values[compi,:];
        end
        tmp = tmp + components;
    end
    
    return vect;
end

# Only for Dirichlet Boundaries!
# Copy the dirichlet values from vec into var.values
# If zero_vals, the values in vec will be zero after copying
function copy_bdry_vals_to_variables(var::Vector{Variable}, vec::Vector{Float64}, grid::Grid, 
                                    dofs_per_node::Int, zero_vals::Bool=true)
    dofind = 0;
    for vi=1:length(var)
        for compo=1:length(var[vi].symvar)
            dofind = dofind + 1;
            for bid=1:size(prob.bc_type,2)
                if prob.bc_type[var[vi].index, bid] == DIRICHLET
                    for i = 1:length(grid.bdryface[bid]) # loop over faces with this BID
                        fid = grid.bdryface[bid][i];
                        # Handle nodal and cell variables separately
                        if var[vi].location == CELL
                            eid = grid.face2element[1,fid];
                            var[vi].values[compo, eid] = vec[(eid-1)*dofs_per_node + dofind];
                            if zero_vals
                                vec[(eid-1)*dofs_per_node + dofind] = 0;
                            end
                            
                        else # NODAL
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
    end
end

# Only for Dirichlet Boundaries!
# Copy the dirichlet values from var.values into vec_b
function copy_bdry_vals_to_vector(var::Vector{Variable}, vec::Vector{Float64}, grid::Grid, dofs_per_node::Int)
    dofind = 0;
    for vi=1:length(var)
        for compo=1:length(var[vi].symvar)
            dofind = dofind + 1;
            for bid=1:size(prob.bc_type,2)
                if prob.bc_type[var[vi].index, bid] == DIRICHLET
                    for i = 1:length(grid.bdryface[bid]) # loop over faces with this BID
                        fid = grid.bdryface[bid][i];
                        # Handle nodal and cell variables separately
                        if var[vi].location == CELL
                            eid = grid.face2element[1,fid];
                            vec[(eid-1)*dofs_per_node + dofind] = var[vi].values[compo, eid];
                            
                        else # NODAL
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
    end
end
