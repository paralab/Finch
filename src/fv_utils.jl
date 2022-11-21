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

# struct FVInfo
#     fluxOrder::Int                  # Order for flux reconstruction
#     cellCenters::Matrix             # Coordinates of cell centers
#     faceCenters::Matrix             # Coordinates of face centers
    
#     cell2node::Vector{Vector}       # Neighboring cell indices for each node
#     cell2nodeWeight::Vector{Vector} # Cell to node interpolation weights
    
#     parentMaps::Union{Nothing, ParentMaps} # For when there is a parent-child mesh
# end

# etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5]; # number of vertices for each type
# etypetodim= [1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3]; # dimension of each type
# etypetonf = [2, 3, 4, 4, 6, 5, 5, 2, 3, 4, 4, 6, 5, 5, 1, 4, 6, 5, 5]; # number of faces for element types
# etypetoftype=[0,1, 1, 2, 3, 3, 3, 0, 1, 1, 2, 3, 3, 3, 0, 0, 0, 0, 0]; # type of faces for this element type

# Build the FVInfo from a given grid
function build_FV_info(grid, order=1, p_maps=nothing)
    int_type = finch_state.config.index_type;
    float_type = finch_state.config.float_type;
    
    nel = size(grid.loc2glb, 2);
    nface = size(grid.face2glb, 3);
    nnode = size(grid.allnodes, 2);
    dim = size(grid.allnodes, 1);
    
    cell_centers = zeros(float_type, dim, nel);
    face_centers = zeros(float_type, dim, nface);
    cell2node = fill(zeros(int_type,0),nnode);
    cell2nodeWeight = fill(zeros(float_type,0),nnode);
    
    # Find cell centers
    center = zeros(float_type, dim);
    nvtx = size(grid.loc2glb, 1);
    for ei=1:nel
        center .= 0;
        for ni=1:nvtx
            center .+= grid.allnodes[:, grid.loc2glb[ni, ei]];
        end
        cell_centers[:, ei] .= center ./ nvtx;
    end
    
    # Find face centers
    nvtx = size(grid.face2glb, 1);
    for fi=1:nface
        center .= 0;
        for ni=1:nvtx
            center .+= grid.allnodes[:, grid.face2glb[ni, 1, fi]];
        end
        face_centers[:, fi] .= center ./ nvtx;
    end
    
    # Build cell2node map
    # Check all elements to find ones touching this node.
    nvtx = size(grid.loc2glb, 1);
    for ei=1:nel
        for ni=1:nvtx
            push!(cell2node[grid.loc2glb[ni, ei]], ei);
        end
    end
    
    # Compute cell2node interpolation weights
    for ni=1:nnode
        cell2nodeWeight[ni] = zeros(float_type, length(cell2node[ni]));
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
    
    if p_maps === nothing
        p_maps = ParentMaps();
    end
    
    return FVInfo(order, cell_centers, face_centers, cell2node, cell2nodeWeight, p_maps);
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
    Nnodes = size(finch_state.fv_grid.allnodes,2);
    if node_values === nothing
        node_values = zeros(Nnodes * dofs);
    end
    
    for ni=1:Nnodes
        if dofs > 1
            st = (ni-1)*dofs + 1;
            en = ni*dofs;
            node_values[st:en] .= 0;
            ncell = length(finch_state.fv_info.cell2nodeWeight[ni]);
            for ci=1:ncell
                cind = finch_state.fv_info.cell2node[ni][ci];
                cst = (cind-1)*dofs + 1;
                cen = cind*dofs;
                node_values[st:en] .+= finch_state.fv_info.cell2nodeWeight[ni][ci] .* cell_values[cst:cen];
            end
            
        else
            node_values[ni] = 0;
            ncell = length(finch_state.fv_info.cell2nodeWeight[ni]);
            for ci=1:ncell
                node_values[ni] += finch_state.fv_info.cell2nodeWeight[ni][ci] * cell_values[finch_state.fv_info.cell2node[ni][ci]];
            end
        end
    end
    
    return node_values;
end

# Find cell averages from nodal values
# Returns an array of cell values.
function FV_node_to_cell_all(node_values, cell_values = nothing)
    Ncells = size(finch_state.fv_grid.loc2glb,2);
    if cell_values === nothing
        cell_values = zeros(Ncells);
    end
    
    for ci=1:Ncells
        cell_values[ci] = 0;
        nnode = length(finch_state.fv_grid.loc2glb[:,ci]);
        for ni=1:nnode
            cell_values[ci] += node_values[finch_state.fv_grid.loc2glb[ni, ci]];
        end
        cell_values[ci] /= nnode;
    end
    
    return cell_values;
end
