#=
Utilities used by the generated solve() functions.
These will be called by the generated solve code
so they should be made as efficient as possible.
=#

# Solves a sparse system IN PLACE
function linear_system_solve!(A::Union{SparseMatrixCSC, LinearMap}, b::Vector, x::Vector, config::Finch_config)
    if config.num_procs > 1 && config.proc_rank > 0
        # Other procs don't have the full A, just return b
        return b;
    end
    
    if config.linalg_usePetsc == false
        if config.linalg_matrixfree
            # How should we handle preconditioners for matrix-free?
            preconditioner = IterativeSolvers.Identity();
            
            if config.linalg_iterative_method == "GMRES"
                if config.linalg_iterative_maxiter == 0
                    config.linalg_iterative_maxiter = 500;
                end
                if config.linalg_iterative_gmresRestart == 0
                    config.linalg_iterative_gmresRestart = 20;
                end
                IterativeSolvers.gmres!(x, A, b, abstol=config.linalg_iterative_abstol, reltol=config.linalg_iterative_reltol,
                            maxiter=config.linalg_iterative_maxiter, restart=config.linalg_iterative_gmresRestart, 
                            Pl=preconditioner, verbose=config.linalg_iterative_verbose);
                            
            else # "CG"
                if config.linalg_iterative_maxiter == 0
                    config.linalg_iterative_maxiter = 500;
                end
                IterativeSolvers.cg!(x, A, b, abstol=config.linalg_iterative_abstol, reltol=config.linalg_iterative_reltol,
                            maxiter=config.linalg_iterative_maxiter, 
                            Pl=preconditioner, verbose=config.linalg_iterative_verbose);
            end
            
        elseif config.linalg_iterative
            if config.linalg_iterative_pc == "ILU"
                preconditioner = IncompleteLU.ilu(A);
            elseif config.linalg_iterative_pc == "AMG"
                preconditioner = AlgebraicMultigrid.aspreconditioner(ruge_stuben(A));
            else
                preconditioner = IterativeSolvers.Identity();
            end
            
            if config.linalg_iterative_method == "GMRES"
                if config.linalg_iterative_maxiter == 0
                    config.linalg_iterative_maxiter = size(A,2);
                end
                if config.linalg_iterative_gmresRestart == 0
                    config.linalg_iterative_gmresRestart = min(20,size(A,2));
                end
                IterativeSolvers.gmres!(x, A, b, abstol=config.linalg_iterative_abstol, reltol=config.linalg_iterative_reltol,
                            maxiter=config.linalg_iterative_maxiter, restart=config.linalg_iterative_gmresRestart, 
                            Pl=preconditioner, verbose=config.linalg_iterative_verbose);
                            
            else # "CG"
                if config.linalg_iterative_maxiter == 0
                    config.linalg_iterative_maxiter = size(A,2);
                end
                IterativeSolvers.cg!(x, A, b, abstol=config.linalg_iterative_abstol, reltol=config.linalg_iterative_reltol,
                            maxiter=config.linalg_iterative_maxiter, 
                            Pl=preconditioner, verbose=config.linalg_iterative_verbose);
            end
            
        else # default
            # Yes, this is not truly in-place. Sparse factorizations are not available in-place.
            x .= A\b;
        end
        
    else
        # For now try this simple setup.
        # There will be many options that need to be available to the user. TODO
        # The matrix should really be constructed from the begninning as a PETSc one. TODO
        
        petsclib = PETSc.petsclibs[1]           #
        inttype = PETSc.inttype(petsclib);      # These should maybe be kept in config
        scalartype = PETSc.scalartype(petsclib);#
        
        (I, J, V) = findnz(A);
        (n,n) = size(A);
        nnz = zeros(inttype, n); # number of non-zeros per row
        for i=1:length(I)
            nnz[I[i]] += 1;
        end
        
        petscA = PETSc.MatSeqAIJ{scalartype}(n,n,nnz);
        for i=1:length(I)
            petscA[I[i],J[i]] = V[i];
        end
        PETSc.assemble(petscA);
        
        if config.linalg_iterative_pc == "AMG"
            pcType = "mg";
        else
            pcType = "ilu";
        end
        
        ksp = PETSc.KSP(petscA; ksp_rtol=config.linalg_iterative_reltol, ksp_atol=config.linalg_iterative_abstol, 
                        ksp_max_it=config.linalg_iterative_maxiter, pc_type=pcType, 
                        ksp_monitor=config.linalg_iterative_verbose); # Options should be available to user
        
        x .= ksp\b;
    end
    
    return x;
end

# Solves a sparse system
function linear_system_solve(A::Union{SparseMatrixCSC, LinearMap}, b::Vector, config::Finch_config)
    if config.num_procs > 1 && config.proc_rank > 0
        # Other procs don't have the full A, just return b
        return b;
    end
    
    if config.linalg_usePetsc == false
        if config.linalg_matrixfree
            # How should we handle preconditioners for matrix-free?
            preconditioner = IterativeSolvers.Identity();
            
            if config.linalg_iterative_method == "GMRES"
                if config.linalg_iterative_maxiter == 0
                    config.linalg_iterative_maxiter = 500;
                end
                if config.linalg_iterative_gmresRestart == 0
                    config.linalg_iterative_gmresRestart = 20;
                end
                return IterativeSolvers.gmres(A, b, abstol=config.linalg_iterative_abstol, reltol=config.linalg_iterative_reltol,
                            maxiter=config.linalg_iterative_maxiter, restart=config.linalg_iterative_gmresRestart, 
                            Pl=preconditioner, verbose=config.linalg_iterative_verbose);
                            
            else # "CG"
                if config.linalg_iterative_maxiter == 0
                    config.linalg_iterative_maxiter = 500;
                end
                return IterativeSolvers.cg(A, b, abstol=config.linalg_iterative_abstol, reltol=config.linalg_iterative_reltol,
                            maxiter=config.linalg_iterative_maxiter, 
                            Pl=preconditioner, verbose=config.linalg_iterative_verbose);
            end
            
        elseif config.linalg_iterative
            if config.linalg_iterative_pc == "ILU"
                preconditioner = IncompleteLU.ilu(A);
            elseif config.linalg_iterative_pc == "AMG"
                preconditioner = AlgebraicMultigrid.aspreconditioner(ruge_stuben(A));
            else
                preconditioner = IterativeSolvers.Identity();
            end
            
            if config.linalg_iterative_method == "GMRES"
                if config.linalg_iterative_maxiter == 0
                    config.linalg_iterative_maxiter = size(A,2);
                end
                if config.linalg_iterative_gmresRestart == 0
                    config.linalg_iterative_gmresRestart = min(20,size(A,2));
                end
                return IterativeSolvers.gmres(A, b, abstol=config.linalg_iterative_abstol, reltol=config.linalg_iterative_reltol,
                            maxiter=config.linalg_iterative_maxiter, restart=config.linalg_iterative_gmresRestart, 
                            Pl=preconditioner, verbose=config.linalg_iterative_verbose);
                            
            else # "CG"
                if config.linalg_iterative_maxiter == 0
                    config.linalg_iterative_maxiter = size(A,2);
                end
                return IterativeSolvers.cg(A, b, abstol=config.linalg_iterative_abstol, reltol=config.linalg_iterative_reltol,
                            maxiter=config.linalg_iterative_maxiter, 
                            Pl=preconditioner, verbose=config.linalg_iterative_verbose);
            end
            
        else # default
            return A\b;
        end
        
    else
        # For now try this simple setup.
        # There will be many options that need to be available to the user. TODO
        # The matrix should really be constructed from the begninning as a PETSc one. TODO
        
        petsclib = PETSc.petsclibs[1]           #
        inttype = PETSc.inttype(petsclib);      # These should maybe be kept in config
        scalartype = PETSc.scalartype(petsclib);#
        
        (I, J, V) = findnz(A);
        (n,n) = size(A);
        nnz = zeros(inttype, n); # number of non-zeros per row
        for i=1:length(I)
            nnz[I[i]] += 1;
        end
        
        petscA = PETSc.MatSeqAIJ{scalartype}(n,n,nnz);
        for i=1:length(I)
            petscA[I[i],J[i]] = V[i];
        end
        PETSc.assemble(petscA);
        
        if config.linalg_iterative_pc == "AMG"
            pcType = "mg";
        else
            pcType = "ilu";
        end
        
        ksp = PETSc.KSP(petscA; ksp_rtol=config.linalg_iterative_reltol, ksp_atol=config.linalg_iterative_abstol, 
                        ksp_max_it=config.linalg_iterative_maxiter, pc_type=pcType, 
                        ksp_monitor=config.linalg_iterative_verbose); # Options should be available to user
        
        return ksp\b;
    end
end

# place the values from sol into the variable value arrays
function place_vector_in_vars(var::Vector{Variable}, vect::Vector)
    tmp = 0;
    totalcomponents = 0;
    for vi=1:length(var)
        totalcomponents = totalcomponents + var[vi].total_components;
    end
    
    for vi=1:length(var)
        components = var[vi].total_components;
        for compi=1:components
            var[vi].values[compi,:] .= @view(vect[(compi+tmp):totalcomponents:end]);
        end
        tmp = tmp + components;
    end
    
    return nothing;
end

# Set the values from variable arrays in a global vector
function get_var_vals(var::Vector{Variable}, vect::Union{Nothing, Vector}=nothing)
    tmp = 0;
    totalcomponents = 0;
    for vi=1:length(var)
        totalcomponents = totalcomponents + var[vi].total_components;
    end
    if vect === nothing
        vect = zeros(config.float_type, totalcomponents * size(var[1].values, 2));
    end
    
    for vi=1:length(var)
        components = var[vi].total_components;
        for compi=1:components
            vect[(compi+tmp):totalcomponents:end] .= @view(var[vi].values[compi,:]);
        end
        tmp = tmp + components;
    end
    
    return vect;
end

# Only for Dirichlet Boundaries!
# Copy the dirichlet values from vec into var.values
# If zero_vals, the values in vec will be zero after copying
function copy_bdry_vals_to_variables(var::Vector{Variable}, vec::Vector, grid::Grid, 
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
function copy_bdry_vals_to_vector(var::Vector{Variable}, vec::Vector, grid::Grid, dofs_per_node::Int)
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

######################################################################################################
## Specific to FVM
######################################################################################################

# Reconstructs var_component at face center x based on a neighborhood of cells around it.
# childid and faceid are relative to the parent
# side could be 0,3=centered, 1=left, 2=right
function FV_reconstruct_face_value(var::Array, component::Int, faceid::Int, side::Int, mesh::Grid, fv_info::FVInfo, derivs...)
    dim = size(mesh.allnodes,1);
    # Is it a centered approximation, or one side?
    if (side == 1 || side == 2) && length(derivs)==0
        left_cells = fv_info.parentMaps.face_neighborhoods[1,faceid];
        right_cells = fv_info.parentMaps.face_neighborhoods[2,faceid];
        # get coordinates and values of cells
        left_cellx = fv_info.cellCenters[:,left_cells];
        left_cellu = var[component, left_cells];
        right_cellx = fv_info.cellCenters[:,right_cells];
        right_cellu = var[component, right_cells];
        # coordinates ofthe face
        x = fv_info.faceCenters[:, faceid];
        
        (left_val, right_val) = FV_reconstruct_value_left_right(left_cellx, right_cellx, left_cellu, right_cellu, x, limiter="vanleer");
        if side == 1
            return left_val;
        else
            return right_val;
        end
        
    else # centered or derivatives
        left_cells = fv_info.parentMaps.face_neighborhoods[1,faceid];
        right_cells = fv_info.parentMaps.face_neighborhoods[2,faceid];
        count = Int(floor((order+1)/2));
        cells = zeros(Int, count*2);
        nnz = 0;
        ind = 0;
        while ind < count
            if length(left_cells) >= ind
                nnz += 1;
                cells[nnz] = left_cells[ind];
            end
            if length(right_cells) >= ind
                nnz += 1;
                cells[nnz] = right_cells[ind];
            end
            ind += 1;
        end
        if nnz < count*2
            cells = cells[1:nnz];
        end
        
        # extract the cell coordinates and values and face center coordinates
        cellx = fv_info.cellCenters[:,cells];
        cellu = var[component, cells];
        x = fv_info.faceCenters[:, faceid];
        
        return polyharmonic_interp(x, cellx, cellu)[1];
    end
end

# # Reconstructs var_component at face center x based on a neighborhood of cells around it.
# # childid and faceid are relative to the parent
# # side could be 0,3=centered, 1=left, 2=right
# function FV_reconstruct_value(var::Array, component::Int, childid::Int, faceid::Int, parentid::Int, side::Int, 
#                                 patch_cells::Vector{Int}, mesh::Grid, fv_info::FVInfo, derivs...)
#     dim = size(mesh.allnodes,1);
#     order = fv_info.fluxOrder;
#     min_needed = max(order, length(derivs)+1);
#     # Reconstructing a cell value or a face value?
#     # only one of these can be positive
#     if faceid > 0 # face value
#         # Is it a centered approximation, or one side?
#         if (side == 1 || side == 2) && length(derivs)==0
#             # collect cell indices
#             left_cells = fv_info.parentMaps.leftCells[faceid];
#             if length(left_cells) > min_needed
#                 left_cells = patch_cells[left_cells[1:min_needed]];
#             else
#                 left_cells = patch_cells[left_cells];
#             end
            
#             right_cells = fv_info.parentMaps.rightCells[faceid];
#             if length(right_cells) > min_needed
#                 right_cells = patch_cells[right_cells[1:min_needed]];
#             else
#                 right_cells = patch_cells[right_cells];
#             end
            
#             # remove zeros
#             first_zero = 1;
#             for i=1:length(left_cells)
#                 if left_cells[i] == 0
#                     break;
#                 end
#                 first_zero += 1;
#             end
#             left_cells = left_cells[1:(first_zero-1)];
            
#             first_zero = 1;
#             for i=1:length(right_cells)
#                 if right_cells[i] == 0
#                     break;
#                 end
#                 first_zero += 1;
#             end
#             right_cells = right_cells[1:(first_zero-1)];
            
#             left_cellx = fv_info.cellCenters[:,left_cells];
#             left_cellu = var[component, left_cells];
#             right_cellx = fv_info.cellCenters[:,right_cells];
#             right_cellu = var[component, right_cells];
#             x = fv_info.faceCenters[:, fv_info.parentMaps.parent2face[faceid, parentid]];
            
#             (left_val, right_val) = FV_reconstruct_value_left_right(left_cellx, right_cellx, left_cellu, right_cellu, x, limiter="vanleer");
#             if side == 1
#                 return left_val;
#             else
#                 return right_val;
#             end
            
#         else # centered or derivatives
#             cellsL = fv_info.parentMaps.leftCells[faceid];
#             cellsR = fv_info.parentMaps.rightCells[faceid];
#             cellsL = patch_cells[cellsL];
#             cellsR = patch_cells[cellsR];
#             cells = zeros(Int, min_needed);
#             nnz = 0;
#             ind = 1;
#             while nnz < min_needed
#                 if cellsL[ind] > 0
#                     nnz += 1;
#                     cells[nnz] = cellsL[ind];
#                 end
#                 if cellsR[ind] > 0 && nnz <= min_needed
#                     nnz += 1;
#                     cells[nnz] = cellsR[ind];
#                 end
#                 ind += 1;
#             end
            
#             # extract the cell coordinates and values and face center coordinates
#             cellx = fv_info.cellCenters[:,cells];
#             cellu = var[component, cells];
#             x = fv_info.faceCenters[:, fv_info.parentMaps.parent2face[faceid, parentid]];
            
#             return polyharmonic_interp(x, cellx, cellu)[1];
#         end
        
#     else # cell value
#         #TODO
#     end
    
# end

# reconstructs u at x based on a cell set
function FV_reconstruct_value(cellx::Matrix, cellu::Vector, x::Vector)
    # If only one cell given, return that value
    if length(cellu) == 1
        return cellu[1];
    end
    
    # Interpolate value at x
    # Note: if extrapolating from only one side, 
    # consider using FV_reconstruct_value_left_right to include slope limiting
    return polyharmonic_interp(x, cellx, cellu)[1];
end

# reconstructs u at x based on left and right cell sets
function FV_reconstruct_value_left_right(leftx::Matrix, rightx::Matrix, 
                                        leftu::Vector, rightu::Vector, 
                                        x::Vector; limiter::String="none")
    left = 0.0;
    right = 0.0;
    leftslope = 0.0;
    rightslope = 0.0;
    centerslope = 0.0;
    
    # if only one cell on either side, just return those values
    # If more, extrapolate
    if length(leftu) == 1
        left = leftu[1];
        
    elseif length(leftu) > 1
        leftslope = (leftu[1]-leftu[2]) / norm(leftx[:,1] - leftx[:,2]);
        left = polyharmonic_interp(x, leftx, leftu)[1];
    end
    
    if length(rightu) == 1
        right = rightu[1];
        
    elseif length(rightu) > 1
        rightslope = (rightu[2]-rightu[1]) / norm(rightx[:,1] - rightx[:,2]);
        right = polyharmonic_interp(x, rightx, rightu)[1];
    end
    
    # if one side is empty(boundary?), just make them equal
    if length(leftu) == 0 && length(rightu) > 0
        left = right;
    elseif length(rightu) == 0 && length(leftu) > 0
        right = left;
        
    else
        center_dist = norm(rightx[1] - leftx[1]);
        if center_dist < 1e-16
            # This could also be a boundary.
            centerslope = 0.0;
        else
            centerslope = (rightu[1] - leftu[1]) / center_dist;
        end
        
    end
    
    # If a limiter is specified, limit the slope
    if !(limiter == "none")
        left_r = abs(leftslope) > 1e-10 ? centerslope / leftslope : 1;
        right_r = abs(rightslope) > 1e-10 ? centerslope / rightslope : 1;
        
        left_phi = 2;
        right_phi = 2;
        if limiter == "vanleer"
            left_phi = (left_r + abs(left_r)) / (1+left_r);
            right_phi = (right_r + abs(right_r)) / (1+right_r);
        end
        
        if length(leftu) > 1
            left = (1-left_phi*0.5)*leftu[1] + (left_phi*0.5)*left;
        end
        if length(rightu) > 1
            right = (1-right_phi*0.5)*rightu[1] + (right_phi*0.5)*right;
        end
    end
    
    return (left, right);
end

# Returns arrays of left and right cells
function FV_get_left_right_cells_old(patch::Vector{Int}, face::Int, maps::ParentMaps, dim::Int, order::Int)
    # Provide all of the cells that can be used for this face in a particular order.
    # Make two arrays, one for left one for right, starting from the nearest neighbor.
    left_cells = [];
    right_cells = [];
    if dim == 1
        left_cell_table_1d = [
            [ # children=1
                [2],
                [1, 2]
            ],
            [ # children=2
                [4, 3],
                [1, 4, 3],
                [2, 1, 4, 3]
            ],
            [ # children=3
                [6, 5, 4],
                [1, 6, 5, 4],
                [2, 1, 6, 5, 4],
                [3, 2, 1, 6, 5, 4]
            ],
            [ # children=4
                [8, 7, 6, 5],
                [1, 8, 7, 6, 5],
                [2, 1, 8, 7, 6, 5],
                [3, 2, 1, 8, 7, 6, 5],
                [4, 3, 2, 1, 8, 7, 6, 5],
            ]
        ]
        right_cell_table_1d = [
            [ # children=1
                [1, 3],
                [3]
            ],
            [ # children=2
                [1, 2, 5, 6],
                [2, 5, 6],
                [5, 6]
            ],
            [ # children=3
                [1, 2, 3, 7, 8, 9],
                [2, 3, 7, 8, 9],
                [3, 7, 8, 9],
                [7, 8, 9]
            ],
            [ # children=4
                [1, 2, 3, 4, 9, 10, 11, 12],
                [2, 3, 4, 9, 10, 11, 12],
                [3, 4, 9, 10, 11, 12],
                [4, 9, 10, 11, 12],
                [9, 10, 11, 12]
            ]
        ]
        nchildren = size(maps.parent2child,1);
        left_cells = left_cell_table_1d[nchildren][face];
        right_cells = right_cell_table_1d[nchildren][face];
        
    elseif dim == 2
        if size(maps.parent2neighbor,1) == 3 # triangles
            # Triangle parents have 9 faces, patches have 16 cells
            # Here Left means toward the center of the central parent
            left_cell_table_triangle = [
                [1, 4, 3, 15, 16, 13, 2],
                [2, 4, 3, 9, 12, 11, 1],
                [2, 4, 1, 7, 8, 5, 3],
                [3, 4, 1, 13, 16, 15, 2],
                [3, 4, 2, 11, 12, 9, 1],
                [1, 4, 2, 5, 8, 7, 3],
                [4, 2, 3, 9, 11, 12],
                [4, 1, 3, 15, 13, 16],
                [4, 1, 2, 5, 7, 8]
            ]
            right_cell_table_triangle = [
                [5, 8, 6, 7],
                [7, 8, 6, 5],
                [9, 12, 10, 11],
                [11, 12, 10, 9],
                [13, 16, 14, 15],
                [15, 16, 14, 13],
                [1, 5, 15, 8, 16],
                [2, 7, 9, 8, 12],
                [3, 11, 13, 12, 16]
            ]
            left_cells = left_cell_table_triangle[face];
            right_cells = right_cell_table_triangle[face];
        else # quads
            # Quad parents have 12 faces, patches have 20 cells
            # Here Left means toward the center of the central parent
            left_cell_table_quad = [
                [1, 4, 2, 16, 15, 3, 13, 14],
                [2, 3, 1, 13, 14, 4, 16, 15],
                [2, 1, 3, 20, 19, 4, 17, 18],
                [3, 4, 2, 17, 18, 1, 20, 19],
                [3, 2, 4, 8, 7, 1, 5, 6],
                [4, 1, 3, 5, 6, 2, 8, 7],
                [4, 3, 1, 12, 11, 2, 9, 10],
                [1, 2, 4, 9, 10, 3, 12, 11],
                [1, 20, 4, 19, 17, 18, 5],
                [2, 8, 1, 7, 5, 6, 9],
                [3, 12, 2, 11, 9, 10, 13],
                [4, 16, 3, 15, 13, 14, 17]
            ]
            right_cell_table_quad = [
                [5, 6, 8, 7],
                [8, 7, 5, 6],
                [9, 10, 12, 11],
                [12, 11, 9, 10],
                [13, 14, 16, 15],
                [16, 15, 13, 14],
                [17, 18, 20, 19],
                [20, 19, 17, 18],
                [2, 9, 3, 10, 12, 11],
                [3, 13, 4, 14, 16, 15],
                [4, 17, 1, 18, 20, 19],
                [1, 5, 2, 6, 8, 7]
            ]
            left_cells = left_cell_table_quad[face];
            right_cells = right_cell_table_quad[face];
        end
        
    elseif dim == 3
        #TODO
    end
    
    if length(left_cells) > order
        left_cells = left_cells[1:order];
    end
    if length(right_cells) > order
        right_cells = right_cells[1:order];
    end
    
    return (patch[left_cells], patch[right_cells]);
end

# Set the i,j indices in the ai and aj vectors for building the sparse matrix
function set_matrix_indices!(ai::Vector{Int}, aj::Vector{Int}, dofs_per_node::Int, grid::Grid)
    dofs_squared = dofs_per_node*dofs_per_node;
    nel = grid.nel_owned;
    nfaces = size(grid.face2element, 2);
    faces_per_element = size(grid.element2face, 1);
    face_done = fill(false, nfaces);
    
    # Elemental loop
    for ei=1:nel
        eid = grid.elemental_order[ei]; # The index of this element
        first_ind = (eid-1)*dofs_per_node + 1; # First global index for this element
        last_ind = eid*dofs_per_node; # last global index
        
        # Diagonal blocks for each cell
        # 1 2
        # 3 4
        for di = 1:dofs_per_node
            first = (eid-1)*dofs_squared + (di-1)*dofs_per_node + 1;
            last = first + dofs_per_node - 1;
            ai[first:last] .= first_ind + di - 1;
            aj[first:last] .= first_ind:last_ind;
        end
    end# element loop
    
    # Diagonal and off diagonal blocks for each face
    # Loop over faces
    for fid = 1:nfaces
        (leftel, rightel) = grid.face2element[:,fid];
        left_first = (leftel-1)*dofs_per_node + 1; # First global index for left element
        left_last = left_first + dofs_per_node - 1; # last global index
        right_first = (rightel-1)*dofs_per_node + 1; # First global index for right element
        right_last = right_first + dofs_per_node - 1; # last global index
        
        # left diagonal dofs
        for di = 1:dofs_per_node
            first = nel*dofs_squared + (fid-1)*dofs_squared*4 + (di-1)*dofs_per_node + 1;
            last = first + dofs_per_node - 1;
            ai[first:last] .= left_first + di - 1;
            aj[first:last] .= left_first:left_last;
        end
        
        # Everything else is only set if rightel is a valid element (not 0)
        if rightel > 0
            # left off-diagonal dofs
            for di = 1:dofs_per_node
                first = nel*dofs_squared + (fid-1)*dofs_squared*4 + dofs_squared + (di-1)*dofs_per_node + 1;
                last = first + dofs_per_node - 1;
                ai[first:last] .= left_first + di - 1;
                aj[first:last] .= right_first:right_last;
            end
            
            # right diagonal dofs
            for di = 1:dofs_per_node
                first = nel*dofs_squared + (fid-1)*dofs_squared*4 + dofs_squared*2 + (di-1)*dofs_per_node + 1;
                last = first + dofs_per_node - 1;
                ai[first:last] .= right_first + di - 1;
                aj[first:last] .= right_first:right_last;
            end
            
            # right off-diagonal dofs
            for di = 1:dofs_per_node
                first = nel*dofs_squared + (fid-1)*dofs_squared*4 + dofs_squared*3 + (di-1)*dofs_per_node + 1;
                last = first + dofs_per_node - 1;
                ai[first:last] .= right_first + di - 1;
                aj[first:last] .= left_first:left_last;
            end
        end
        
    end# face loop
    
    # Some indices may be zero due to boundaries. set those all to [1,1] (not the most efficient, but it's consistent)
    # Since the av part will be zero, this won't change the matrix.
    for i=1:length(ai)
        if ai[i]<1 || aj[i]<1
            ai[i] = 1;
            aj[i] = 1;
        end
    end
end

# Returns a copy of a with zeros removed.
function remove_zero_cells(a::Vector{Int})
    b = similar(a);
    nnz = 0;
    for ai in a
        if !(ai == 0)
            nnz += 1;
            b[nnz] = ai;
        end
    end
    b = b[1:nnz];
    
    return b;
end