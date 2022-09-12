#=
Utilities used by the generated solve() functions
These used to be in the CGSolve, FVSolve, etc. modules
=#

# place the values in the variable value arrays
function place_sol_in_vars(var, sol)
    if typeof(var) <: Array
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
    else
        components = var.total_components;
        for compi=1:components
            var.values[compi,:] .= sol[compi:components:end];
        end
    end
end

# Set the values from variable arrays in a global vector
function get_var_vals(var, vect=nothing)
    # place the variable values in a vector
    if typeof(var) <: Array
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
                tmp = tmp + 1;
            end
        end
    else
        components = var.total_components;
        if vect === nothing
            vect = zeros(components * size(var.values, 2));
        end
        for compi=1:components
            vect[compi:components:end] = var.values[compi,:];
        end
    end
    
    return vect;
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
    else
        printerr("Var wasn't in an array? (this is in solver_utils)")
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
    else
        printerr("Var wasn't in an array? (this is in solver_utils)")
    end
end
