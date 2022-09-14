#=
Contains info specific to finite volume such as
- cell enters
- face centers
- interpolation weights
- 

Also has functions useful for FV such as
- cell to node / node to cell interpolation
- 
=#

struct FVInfo
    fluxOrder::Int                  # Order for flux reconstruction
    cellCenters::Array{Float64}     # Coordinates of cell centers
    faceCenters::Array{Float64}     # Coordinates of face centers
    
    cell2node::Array{Array{Int64,1}}            # Neighboring cell indices for each node
    cell2nodeWeight::Array{Array{Float64,1}}    # Cell to node interpolation weights
end

etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5]; # number of vertices for each type
etypetodim= [1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3]; # dimension of each type
etypetonf = [2, 3, 4, 4, 6, 5, 5, 2, 3, 4, 4, 6, 5, 5, 1, 4, 6, 5, 5]; # number of faces for element types
etypetoftype=[0,1, 1, 2, 3, 3, 3, 0, 1, 1, 2, 3, 3, 3, 0, 0, 0, 0, 0]; # type of faces for this element type

# Build the FVInfo from a given grid
function build_FV_info(grid, order=1)
    nel = size(grid.loc2glb, 2);
    nface = size(grid.face2glb, 3);
    nnode = size(grid.allnodes, 2);
    dim = size(grid.allnodes, 1);
    
    cell_centers = zeros(dim, nel);
    face_centers = zeros(dim, nface);
    cell2node = Array{Array{Int64,1},1}(undef,nnode);
    cell2nodeWeight = Array{Array{Float64,1},1}(undef,nnode);
    
    # Find cell centers
    for ei=1:nel
        center = zeros(dim);
        nvtx = size(grid.loc2glb, 1);
        for ni=1:nvtx
            center .+= grid.allnodes[:, grid.loc2glb[ni, ei]];
        end
        cell_centers[:, ei] = center ./ nvtx;
    end
    
    # Find face centers
    for fi=1:nface
        center = zeros(dim);
        nvtx = size(grid.face2glb, 1);
        for ni=1:nvtx
            center .+= grid.allnodes[:, grid.face2glb[ni, 1, fi]];
        end
        face_centers[:, fi] = center ./ nvtx;
    end
    
    # Build cell2node map
    # Check all elements to find ones touching this node.
    nvtx = size(grid.loc2glb, 1);
    for ni=1:nnode
        cell2node[ni] = zeros(Int,0);
        for ei=1:nel
            for nj=1:nvtx
                if grid.loc2glb[nj, ei] == ni
                    # This cell touches this node
                    push!(cell2node[ni], ei);
                    break;
                end
            end
        end
    end
    
    # Compute cell2node interpolation weights
    for ni=1:nnode
        cell2nodeWeight[ni] = zeros(length(cell2node[ni]));
        sumweight = 0;
        for ei=1:length(cell2node[ni])
            dist = norm(grid.allnodes[:,ni] - cell_centers[:, cell2node[ni][ei]]);
            sumweight += 1 / dist;
            cell2nodeWeight[ni][ei] = 1 / dist;
        end
        for ei=1:length(cell2node[ni])
            cell2nodeWeight[ni][ei] /= sumweight;
        end
    end
    
    return FVInfo(order, cell_centers, face_centers, cell2node, cell2nodeWeight);
end

# Interpolate nodal values for one cell.
# For now this just assigns all nodes the cell average.
# This shoould be improved eventually, but needs care to be efficient.
function FV_cell_to_node(cell_value, grid)
    return fill(cell_value, size(grid.loc2glb,1));
end

# Approximates a cell average by simply averaging nodal values.
function FV_node_to_cell(node_values)
    return sum(node_values)/length(node_values);
end

# Interpolate nodal values from neighboring cells.
# Returns an array of nodal values.
function FV_cell_to_node_all(cell_values, node_values; dofs = 1)
    Nnodes = size(fv_grid.allnodes,2);
    if node_values === nothing
        node_values = zeros(Nnodes * dofs);
    end
    
    for ni=1:Nnodes
        if dofs > 1
            st = (ni-1)*dofs + 1;
            en = ni*dofs;
            node_values[st:en] .= 0;
            ncell = length(fv_info.cell2nodeWeight[ni]);
            for ci=1:ncell
                cind = fv_info.cell2node[ni][ci];
                cst = (cind-1)*dofs + 1;
                cen = cind*dofs;
                node_values[st:en] .+= fv_info.cell2nodeWeight[ni][ci] .* cell_values[cst:cen];
            end
            
        else
            node_values[ni] = 0;
            ncell = length(fv_info.cell2nodeWeight[ni]);
            for ci=1:ncell
                node_values[ni] += fv_info.cell2nodeWeight[ni][ci] * cell_values[fv_info.cell2node[ni][ci]];
            end
        end
    end
    
    return node_values;
end

# Find cell averages from nodal values
# Returns an array of cell values.
function FV_node_to_cell_all(node_values, cell_values = nothing)
    Ncells = size(fv_grid.loc2glb,2);
    if cell_values === nothing
        cell_values = zeros(Ncells);
    end
    
    for ci=1:Ncells
        cell_values[ci] = 0;
        nnode = length(fv_grid.loc2glb[:,ci]);
        for ni=1:nnode
            cell_values[ci] += node_values[fv_grid.loc2glb[ni, ci]];
        end
        cell_values[ci] /= nnode;
    end
    
    return cell_values;
end

# reconstructs u at x based on a cell set
function FV_reconstruct_value(cellx, cellu, x)
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
function FV_reconstruct_value_left_right(leftx, rightx, leftu, rightu, x; limiter="none")
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
function get_left_right_cells(patch, face, maps, dim, order)
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
function set_matrix_indices!(ai, aj, dofs_per_node)
    dofs_squared = dofs_per_node*dofs_per_node;
    nel = fv_grid.nel_owned;
    nfaces = size(fv_grid.face2element, 2);
    face_done = fill(false, nfaces);
    
    # Elemental loop
    for ei=1:nel
        eid = elemental_order[ei]; # The index of this element
        first_ind = (eid-1)*dofs_per_node + 1; # First global index for this element
        last_ind = eid*dofs_per_node; # last global index
        
        # Diagonal blocks for each cell
        for di = 1:dofs_per_node
            first = (eid-1)*dofs_squared + (di-1)*dofs_per_node + 1;
            last = first + dofs_per_node - 1;
            ai[first:last] .= first_ind + di - 1;
            aj[first:last] = first_ind:last_ind;
        end
        
        # Diagonal and off diagonal blocks for each face
        # Loop over this element's faces.
        for i=1:refel.Nfaces
            fid = fv_grid.element2face[i, eid];
            
            if !face_done[fid]
                face_done[fid] = true;
                (leftel, rightel) = fv_grid.face2element[:,fid];
                neighborID = (rightel==eid ? leftel : rightel);
                nfirst_ind = (neighborID-1)*dofs_per_node + 1;
                nlast_ind = neighborID*dofs_per_node;
                
                for di = 1:dofs_per_node
                    # The eid components
                    first = nel * dofs_squared + (fid-1)*dofs_squared*4 + (di-1)*dofs_per_node + 1;
                    last = first + dofs_per_node - 1;
                    ai[first:last] .= first_ind + di - 1;
                    aj[first:last] = first_ind:last_ind;
                    
                    # The neighborID components
                    if neighborID > 0
                        first += dofs_squared;
                        last += dofs_squared;
                        ai[first:last] .= first_ind + di - 1;
                        aj[first:last] = nfirst_ind:nlast_ind;
                    end
                end
                
                # While we're here, set the components for the neighborID rows as well.
                if neighborID > 0
                    for di = 1:dofs_per_node
                        # The eid components
                        first = nel * dofs_squared + (fid-1)*dofs_squared*4 + dofs_squared*2 + (di-1)*dofs_per_node + 1;
                        last = first + dofs_per_node - 1;
                        ai[first:last] .= nfirst_ind + di - 1;
                        aj[first:last] = first_ind:last_ind;
                        
                        # The neighborID components
                        first += dofs_squared;
                        last += dofs_squared;
                        ai[first:last] .= nfirst_ind + di - 1;
                        aj[first:last] = nfirst_ind:nlast_ind;
                    end
                end
                
            else
                # already done
            end
        end# face loop
    end# element loop
    
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
function remove_zero_cells(a)
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