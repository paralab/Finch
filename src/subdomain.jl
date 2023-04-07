#=
Utilities related to carving out and handling sub-domains.
In particular, this is for shifted boundary method use.
=#

# Takes a Grid and returns a new Grid with carved out subdomain
function carve_grid(grid::Grid, subdomain::Subdomain, keep_cuts::Bool)
    # useful numbers
    dim = size(grid.allnodes, 1);
    nnodes = size(grid.allnodes, 2);
    nel = size(grid.loc2glb, 2);
    nodes_per_element = size(grid.loc2glb,1);
    nfaces = size(grid.face2glb, 3);
    faces_per_element = size(grid.element2face, 1);
    nodes_per_face = size(grid.face2glb, 1);
    nbids = length(grid.bids);
    
    # Flag each node as inside or outside
    node_inside = zeros(Bool, nnodes);
    if dim == 1
        for i=1:nnodes
            node_inside[i] = subdomain.inside(grid.allnodes[1,i]);
        end
    elseif dim == 2
        for i=1:nnodes
            node_inside[i] = subdomain.inside(grid.allnodes[1,i], grid.allnodes[2,i]);
        end
    elseif dim == 3
        for i=1:nnodes
            node_inside[i] = subdomain.inside(grid.allnodes[1,i], grid.allnodes[2,i], grid.allnodes[3,i]);
        end
    end
    
    # Flag each element, face, and node as kept or removed
    keep_element = zeros(Bool, nel);
    keep_node = zeros(Bool, nnodes);
    keep_face = zeros(Bool, nfaces);
    shift_element = zeros(Int, nel);
    shift_node = zeros(Int, nnodes); # old->new maps
    shift_face = zeros(Int, nfaces);
    removed_elements = 0; # number of elements to remove
    removed_faces = 0; # number of faces to remove
    removed_nodes = 0; # number of nodes to remove
    next_index = 1;
    for ei=1:nel
        in_nodes = 0;
        for ni=1:nodes_per_element
            if node_inside[grid.loc2glb[ni,ei]]
                in_nodes += 1;
            end
        end
        if (in_nodes == nodes_per_element) || (in_nodes > 0 && keep_cuts)
            keep_element[ei] = true;
            shift_element[ei] = next_index;
            next_index += 1;
            for ni=1:nodes_per_element
                keep_node[grid.loc2glb[ni,ei]] = true;
            end
            for fi=1:faces_per_element
                keep_face[grid.element2face[fi,ei]] = true;
            end
        else
            removed_elements += 1;
        end
    end
    next_index = 1;
    for i=1:nnodes
        if !keep_node[i]
            removed_nodes += 1;
        else
            shift_node[i] = next_index;
            next_index += 1;
        end
    end
    next_index = 1;
    for i=1:nfaces
        if !keep_face[i]
            removed_faces += 1;
        else
            shift_face[i] = next_index;
            next_index += 1;
        end
    end
    
    # Make copies of the Grid data with parts removed
    ftype = typeof(grid.allnodes[1,1]);
    new_nnodes = nnodes - removed_nodes;
    new_nel = nel - removed_elements;
    new_nfaces = nfaces - removed_faces;
    bid_node_counts = zeros(Int, nbids);
    bid_face_counts = zeros(Int, nbids);
    
    new_allnodes = zeros(ftype, dim, new_nnodes); #
    
    new_bdry = fill(zeros(Int,0), nbids);
    new_bdryface = fill(zeros(Int,0), nbids);
    new_bdrynorm = fill(zeros(ftype, dim,0), nbids);
    new_bids = grid.bids; #
    new_nodebid = zeros(Int, new_nnodes); #
    
    new_loc2glb = zeros(Int, nodes_per_element, new_nel); #
    new_glbvertex = zeros(Int, size(grid.glbvertex,1), new_nel); #
    
    new_face2glb = zeros(Int, nodes_per_face, size(grid.face2glb,2), new_nfaces); #
    new_element2face = zeros(Int, faces_per_element, new_nel); #
    new_face2element = zeros(Int, 2, new_nfaces); #
    new_facenormals = zeros(ftype, dim, new_nfaces); #
    new_faceRefelInd = zeros(Int, 2, new_nfaces); #
    new_facebid = zeros(Int, new_nfaces); #
    
    # copy node info
    next_index = 1;
    for i=1:nnodes
        if keep_node[i]
            new_allnodes[:, next_index] .= grid.allnodes[:, i];
            new_nodebid[next_index] = grid.nodebid[i];
            if new_nodebid[next_index] > 0
                bid_node_counts[new_nodebid[next_index]] += 1;
            end
            
            next_index += 1;
        end
    end
    
    # copy element info
    next_index = 1;
    for i=1:nel
        if keep_element[i]
            for ni=1:nodes_per_element
                new_loc2glb[ni,next_index] = shift_node[grid.loc2glb[ni,i]];
            end
            for ni=1:size(grid.glbvertex,1)
                new_glbvertex[ni,next_index] = shift_node[grid.glbvertex[ni,i]];
            end
            for fi=1:faces_per_element
                new_element2face[fi,next_index] = shift_face[grid.element2face[fi,i]];
            end
            
            next_index += 1;
        end
    end
    
    # copy face info
    next_index = 1;
    for i=1:nfaces
        if keep_face[i]
            for ni=1:nodes_per_face
                new_face2glb[ni,1,next_index] = shift_node[grid.face2glb[ni,1,i]];
                if size(new_face2glb,2) > 1
                    new_face2glb[ni,2,next_index] = shift_node[grid.face2glb[ni,2,i]];
                end
            end
            # If element on one side has been removed, adjust and make boundary face
            if grid.facebid[i] > 0
                new_face2element[1, next_index] = shift_element[grid.face2element[1, i]];
                new_facebid[next_index] = grid.facebid[i];
                new_facenormals[:,next_index] .= grid.facenormals[:,i];
                new_faceRefelInd[:,next_index] .= grid.faceRefelInd[:,i];
                bid_face_counts[new_facebid[next_index]] += 1;
            elseif !keep_element[grid.face2element[1, i]]
                # This is a new boundary face
                new_face2element[1, next_index] = shift_element[grid.face2element[2,i]];
                new_facebid[next_index] = 1;
                bid_face_counts[1] += 1;
                for ni=1:nodes_per_face
                    if new_nodebid[new_face2glb[ni,1,next_index]] == 0 # new bdry node
                        new_nodebid[new_face2glb[ni,1,next_index]] = 1;
                        bid_node_counts[1] += 1;
                    end
                end
                new_facenormals[:,next_index] .= -grid.facenormals[:,i];
                new_faceRefelInd[1,next_index] = grid.faceRefelInd[2,i];
                new_faceRefelInd[2,next_index] = 0; #grid.faceRefelInd[2,i];
            elseif !keep_element[grid.face2element[2, i]]
                # This is a new boundary face
                new_face2element[1, next_index] = shift_element[grid.face2element[1,i]];
                new_facebid[next_index] = 1;
                bid_face_counts[1] += 1;
                for ni=1:nodes_per_face
                    if new_nodebid[new_face2glb[ni,1,next_index]] == 0 # new bdry node
                        new_nodebid[new_face2glb[ni,1,next_index]] = 1;
                        bid_node_counts[1] += 1;
                    end
                end
                new_facenormals[:,next_index] .= grid.facenormals[:,i];
                new_faceRefelInd[1,next_index] = grid.faceRefelInd[1,i];
                new_faceRefelInd[2,next_index] = 0; #grid.faceRefelInd[1,i];
            else
                # not a boundary face
                new_face2element[1, next_index] = shift_element[grid.face2element[1,i]];
                new_face2element[2, next_index] = shift_element[grid.face2element[2,i]];
                new_facenormals[:,next_index] .= grid.facenormals[:,i];
                new_faceRefelInd[:,next_index] .= grid.faceRefelInd[:,i];
            end
            
            
            next_index += 1;
        end
    end
    
    # Construct boundary info
    # bdry::Vector{Vector{Int}}        # Indices of boundary nodes for each BID (bdry[bid][nodes])*note:array of arrays
    # bdryface::Vector{Vector{Int}}    # Indices of faces touching each BID (bdryface[bid][faces])*note:array of arrays
    # bdrynorm::Vector{Matrix{T}}    # Normal vector for boundary nodes for each BID (bdrynorm[bid][dim, nodes])*note:array of arrays
    for i=1:nbids
        new_bdry[i] = zeros(Int, bid_node_counts[i]);
        new_bdryface[i] = zeros(Int, bid_face_counts[i]);
        new_bdrynorm[i] = zeros(ftype, dim, bid_node_counts[i]);
    end
    next_indices = ones(Int, nbids);
    for i=1:new_nnodes
        bid = new_nodebid[i];
        if bid > 0
            new_bdry[bid][next_indices[bid]] = i;
            next_indices[bid] += 1;
        end
    end
    next_indices = ones(Int, nbids);
    for i=1:new_nfaces
        bid = new_facebid[i];
        if bid > 0
            new_bdryface[bid][next_indices[bid]] = i;
            # also use this face's normal for nodes that share BID
            for ni=1:nodes_per_face
                nodeid = new_face2glb[ni,1,i];
                if new_nodebid[nodeid] == bid # node and face bid match
                    # find its index in bdry
                    for nj = 1:bid_node_counts[bid]
                        if new_bdry[bid][nj] == nodeid
                            new_bdrynorm[bid][:, nj] .= new_facenormals[:, i];
                            break;
                        end
                    end
                end
            end
            next_indices[bid] += 1;
        end
    end
    
    # Partitioned grids need to be reconfigured
    if grid.is_subgrid
        # TODO
        
        # # When partitioning the grid, this stores the ghost info.
        # # Items specifying (for solver type) will be empty/0 for other types.
        # is_subgrid::Bool            # Is this a partition of a greater grid?
        # elemental_order::Vector{Int}     # Order used in elemental loops
        # nel_global::Int             # Number of global elements
        # nel_owned::Int              # Number of elements owned by this partition
        # nel_ghost::Int              # Number of ghost elements (for FV)
        # nface_owned::Int            # Number of faces owned by this partition
        # nface_ghost::Int            # Number of ghost faces that are not owned (for FV)
        # nnodes_global::Int          # Number of global nodes
        # nnodes_borrowed::Int        # Number of nodes borrowed from another partition (for CG)
        # element_owner::Vector{Int}       # The rank of each element's owner or -1 if locally owned (for FV)
        # node_owner::Vector{Int}          # The rank of each node's owner (for FE)
        # partition2global_element::Vector{Int}           # Map from partition elements to global mesh element index
        # partition2global::Vector{Int}    # Global index of nodes (for CG,DG)
        # global_bdry_index::Vector{Int8}   # Index in bids for every global node, or 0 for interior (Only proc 0 holds, only for FE)
        
        # num_neighbor_partitions::Int   # number of partitions that share ghosts with this.
        # neighboring_partitions::Vector{Int} # IDs of neighboring partitions
        # ghost_counts::Vector{Int}           # How many ghost elements for each neighbor (for FV)
        # ghost_index::Vector{Matrix{Int}}     # Lists of ghost elements to send/recv for each neighbor (for FV)
        
        new_grid = grid;
        
    else # not partitioned
        new_grid = Grid(ftype, new_allnodes, new_bdry, new_bdryface, new_bdrynorm, new_bids, new_nodebid, 
                        new_loc2glb, new_glbvertex, new_face2glb, new_element2face, new_face2element, 
                        new_facenormals, new_faceRefelInd, new_facebid);
    end
    
    return new_grid;
end