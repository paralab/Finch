#=
# Contains info about all nodes on the domain
# Unlike MeshData struct, this accounts for interior nodes and corresponds to nodal DOFs.
=#

# struct Grid
#     # nodes
#     allnodes::Array                 # All node coordinates size = (dim, nnodes)
    
#     # boundaries
#     bdry::Vector{Vector}        # Indices of boundary nodes for each BID (bdry[bid][nodes])*note:array of arrays
#     bdryface::Vector{Vector}    # Indices of faces touching each BID (bdryface[bid][faces])*note:array of arrays
#     bdrynorm::Vector{Matrix}    # Normal vector for boundary nodes for each BID (bdrynorm[bid][dim, nodes])*note:array of arrays
#     bids::Vector                # BID corresponding to rows of bdrynodes
#     nodebid::Vector             # BID for every node in allnodes order(interior=0)
    
#     # elements
#     loc2glb::Matrix             # local to global map for each element's nodes (size is (Np, nel))
#     glbvertex::Matrix           # global indices of each elements' vertices (size is (Nvertex, nel))
    
#     # faces (For CG, G=1. For DG, G=2)
#     face2glb::Array             # local to global map for faces (size is (Nfp, G, Nfaces))
#     element2face::Matrix        # face indices for each element (size is (Nfaces, nel))
#     face2element::Matrix        # elements on both sides of a face, 0=boundary (size is (2, Nfaces))
#     facenormals::Matrix         # normal vector for each face
#     faceRefelInd::Matrix        # Index for face within the refel for each side
#     facebid::Vector             # BID of each face (0=interior face)
    
#     # When partitioning the grid, this stores the ghost info.
#     # Items specifying (for solver type) will be empty/0 for other types.
#     is_subgrid::Bool            # Is this a partition of a greater grid?
#     elemental_order::Vector     # Order used in elemental loops
#     nel_global::Int             # Number of global elements
#     nel_owned::Int              # Number of elements owned by this partition
#     nel_ghost::Int              # Number of ghost elements (for FV)
#     nface_owned::Int            # Number of faces owned by this partition
#     nface_ghost::Int            # Number of ghost faces that are not owned (for FV)
#     nnodes_global::Int          # Number of global nodes
#     nnodes_borrowed::Int        # Number of nodes borrowed from another partition (for CG)
#     element_owner::Vector       # The rank of each element's owner or -1 if locally owned (for FV)
#     node_owner::Vector          # The rank of each node's owner (for FE)
#     partition2global_element::Vector           # Map from partition elements to global mesh element index
#     partition2global::Vector    # Global index of nodes (for CG,DG)
#     global_bdry_index::Vector   # Index in bids for every global node, or 0 for interior (Only proc 0 holds, only for FE)
    
#     num_neighbor_partitions::Int   # number of partitions that share ghosts with this.
#     neighboring_partitions::Vector # IDs of neighboring partitions
#     ghost_counts::Vector           # How many ghost elements for each neighbor (for FV)
#     ghost_index::Vector{Array}     # Lists of ghost elements to send/recv for each neighbor (for FV)
    
#     # constructors
#     Grid(allnodes, bdry, bdryfc, bdrynorm, bids, nodebid, loc2glb, glbvertex, f2glb, element2face, 
#          face2element, facenormals, faceRefelInd, facebid) = 
#      new(allnodes, bdry, bdryfc, bdrynorm, bids, nodebid, loc2glb, glbvertex, f2glb, element2face, 
#          face2element, facenormals, faceRefelInd, facebid, 
#          false, Array(1:size(loc2glb,2)), size(loc2glb,2), size(loc2glb,2), 0,size(face2element,2), 0, 0, 0, zeros(Int,0), zeros(Int,0), 
#          zeros(Int,0), zeros(Int,0), zeros(Int8,0), 0, zeros(Int,0), zeros(Int,0), [zeros(Int,2,0)]); # up to facebid only
     
#     Grid(allnodes, bdry, bdryfc, bdrynorm, bids, nodebid, loc2glb, glbvertex, f2glb, element2face, 
#          face2element, facenormals, faceRefelInd, facebid, 
#          ispartitioned, el_order, nel_global, nel_owned, nel_ghost, nface_owned, nface_ghost, nnodes_global, nnodes_borrowed, element_owners, 
#          node_owner, partition2global_element, partition2global, glb_bid, num_neighbors, neighbor_ids, ghost_counts, ghost_ind) = 
#      new(allnodes, bdry, bdryfc, bdrynorm, bids, nodebid, loc2glb, glbvertex, f2glb, element2face, 
#          face2element, facenormals, faceRefelInd, facebid, 
#          ispartitioned, el_order, nel_global, nel_owned, nel_ghost, nface_owned, nface_ghost, nnodes_global, nnodes_borrowed, element_owners, 
#          node_owner, partition2global_element, partition2global, glb_bid, num_neighbors, neighbor_ids, ghost_counts, ghost_ind); # subgrid parts included
         
#     # An empty Grid
#     Grid() = new(
#         [],[],[],[],[],[],[],[],[],[],[],[],[],[],
#         false,[],0,0,0,0,0,0,0,[],[],[],[],[],0,[],[],[]
#     )
# end

# Build a grid from a mesh
# This is for full grids. For partitioned grids see partitioned_grid_from_mesh()
function grid_from_mesh(mesh::MeshData; grid_type=CG, order=1)
    etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5]; # number of vertices for each type
    etypetodim= [1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3]; # dimension of each type
    etypetonf = [2, 3, 4, 4, 6, 5, 5, 2, 3, 4, 4, 6, 5, 5, 1, 4, 6, 5, 5]; # number of faces for element types
    etypetoftype=[0,1, 1, 2, 3, 3, 3, 0, 1, 1, 2, 3, 3, 3, 0, 0, 0, 0, 0]; # type of faces for this element type
    
    config = finch_state.config;
    log_entry("Building full grid from mesh data. Types = "*string(grid_type), 2);
    t_grid_from_mesh = Base.Libc.time();
    int_type = config.index_type;
    float_type = config.float_type;
    if float_type==Float64
        tol = 1e-8;
    elseif float_type==Float32
        tol = 1e-4;
    else
        tol = 1e-4;
    end
    
    dim = config.dimension;
    ord = order;
    nfaces = etypetonf[mesh.etypes[1]];
    #totalfaces = nfaces*mesh.nel;
    totalfaces = size(mesh.normals,2);
    nel = mesh.nel;
    if dim == 1
        facenvtx = 1
    else
        facenvtx = etypetonv[etypetoftype[mesh.etypes[1]]]; # Assumes one element type
    end
    nvtx = etypetonv[mesh.etypes[1]]; # Assumes one element type
    
    log_entry("Building reference element: "*string(dim)*"D, order="*string(ord)*", nfaces="*string(nfaces), 3);
    refel::Refel = build_refel(dim, ord, nfaces, config.elemental_nodes);
    
    if grid_type == DG
        Gness = 2;
    else
        Gness = 1;
    end
    
    vertex_dist_scale = zeros(float_type, nel); # a scale for relative tolerance
    
    Np = refel.Np;                      # number of nodes per element
    bdry = Vector{Vector{int_type}}(undef,0);       # index(in x) of boundary nodes for each BID
    bdryfc = Vector{Vector{int_type}}(undef,0);     # index of faces touching each BID
    bdrynorm = Vector{Matrix{float_type}}(undef,0); # normal at boundary nodes
    bids = collectBIDs(mesh);           # BID list
    nbids = length(bids);
    for i=1:nbids
        push!(bdry, zeros(int_type, 0));
        push!(bdryfc, zeros(int_type,0));
        push!(bdrynorm, zeros(float_type, config.dimension,0));
    end
    loc2glb = zeros(int_type, Np, nel)       # local to global index map for each element's nodes
    glbvertex = zeros(int_type, nvtx, nel);     # local to global for vertices
    f2glb = zeros(int_type, refel.Nfp[1], Gness, totalfaces);  # face node local to global
    element2face = zeros(int_type, nfaces, nel);  # element to face map
    face2element = zeros(int_type, 2, size(mesh.face2element,2));  # face to element map
    facenormals = zeros(float_type, dim, totalfaces); # normal vectors for every face
    faceRefelInd = zeros(int_type, 2, totalfaces); # Index in refel for this face for elements on both sides
    facebid = zeros(int_type, totalfaces); # BID of each face
    
    t_nodes1 = Base.Libc.time();
    
    tmpallnodes = zeros(float_type, dim, mesh.nel*refel.Np);
    n_vert = etypetonv[mesh.etypes[1]]; # Assumes one element type
    e_vert = zeros(dim, n_vert); # Assumes one element type
    e_x = zeros(float_type, refel.Np);
    e_y = zeros(float_type, refel.Np);
    e_z = zeros(float_type, refel.Np);
    
    for ei=1:nel
        # Find this element's nodes
        # n_vert = etypetonv[mesh.etypes[ei]];
        # e_vert = mesh.nodes[1:dim, mesh.elements[1:n_vert, ei]];
        for ni=1:n_vert
            for di=1:dim
                e_vert[di, ni] = mesh.nodes[di, mesh.elements[ni, ei]];
            end
        end
        
        # sample distance betweeen vertices
        for i=2:n_vert
            tmp = 0.0;
            for di=1:dim
                tmp += abs(e_vert[di,i] .- e_vert[di,i-1]);
            end
            vertex_dist_scale[ei] += tmp;
        end
        vertex_dist_scale[ei] /= (n_vert-1);
        
        if dim == 1
            e_x = line_refel_to_x(refel.r[:,1], e_vert);
        elseif dim == 2
            if nfaces == 3 # triangles
                (e_x, e_y) = triangle_refel_to_xy_(refel.r[:,1], refel.r[:,2], e_vert);
            else # quads
                (e_x, e_y) = quad_refel_to_xy(refel.r[:,1], refel.r[:,2], e_vert);
            end
        elseif dim == 3
            if nvtx == 8 # hexes
                (e_x, e_y, e_z) = hex_refel_to_xyz(refel.r[:,1], refel.r[:,2], refel.r[:,3], e_vert);
            else # tets
                # (e_x, e_y, e_z) = tetrahedron_refel_to_xyz(refel.r[:,1], refel.r[:,2], refel.r[:,3], e_vert);
                tetrahedron_refel_to_xyz!(refel.r, e_vert, e_x, e_y, e_z);
            end
        end
        
        # Add them to the tmp global nodes
        tmpallnodes[1, ((ei-1)*Np+1):(ei*Np)] .= e_x;
        if dim > 1
            tmpallnodes[2, ((ei-1)*Np+1):(ei*Np)] .= e_y;
            if dim > 2
                tmpallnodes[3, ((ei-1)*Np+1):(ei*Np)] .= e_z;
            end
        end
        
        # temporary mapping
        # loc2glb[:,ei] = ((ei-1)*Np+1):(ei*Np);
        for ni=1:Np
            loc2glb[ni,ei] = (ei-1)*Np+ni;
        end
    end
    t_nodes1 = Base.Libc.time() - t_nodes1;
    log_entry("Initial node setup: "*string(t_nodes1), 3);
    
    if Gness == 1 # CG
        # Go back and remove duplicate nodes. Adjust loc2glb.
        t_nodes1 = Base.Libc.time();
        (allnodes, loc2glb) = remove_duplicate_nodes(tmpallnodes, loc2glb, tol=tol);
        t_nodes1 = Base.Libc.time() - t_nodes1;
        log_entry("Remove duplicate nodes: "*string(t_nodes1), 3);
        
    else
        allnodes = tmpallnodes; # DG grid is already made
    end
    
    # bid for every node start with 0 for interior
    node_bids = zeros(int_type, size(allnodes,2));
    
    # vertices, faces and boundary
    t_faces1 = Base.Libc.time();
    
    mfids = zeros(Int, nfaces);
    normals = zeros(dim, nfaces);
    el_center = zeros(float_type, dim);
    Nfp = refel.Nfp[1];
    tmpf2glb = zeros(Int, Nfp);
    faceNodesA = zeros(dim, Nfp);
    faceNodesB = zeros(dim, Nfp);
    
    for ei=1:nel
        n_vert = etypetonv[mesh.etypes[ei]];
        # mfids = mesh.element2face[:,ei];
        # normals = mesh.normals[:,mfids];
        for fi=1:nfaces
            mfids[fi] = mesh.element2face[fi, ei];
            for di=1:dim
                normals[di,fi] = mesh.normals[di, mfids[fi]];
            end
        end
        
        # vertices and center
        el_center .= 0.0;
        for ni=1:Np
            for di=1:dim
                el_center[di] += allnodes[di, loc2glb[ni,ei]];
            end
            
            for vi=1:n_vert
                # if is_same_node(mesh.nodes[:, mesh.elements[vi,ei]], allnodes[:,loc2glb[ni,ei]], tol, vertex_dist_scale[ei])
                if (dim == 3 && is_same_node_3d(tol, vertex_dist_scale[ei], mesh.nodes[1, mesh.elements[vi,ei]], allnodes[1,loc2glb[ni,ei]], mesh.nodes[2, mesh.elements[vi,ei]], allnodes[2,loc2glb[ni,ei]], mesh.nodes[3, mesh.elements[vi,ei]], allnodes[3,loc2glb[ni,ei]])) ||
                    (dim == 2 && is_same_node_2d(tol, vertex_dist_scale[ei], mesh.nodes[1, mesh.elements[vi,ei]], allnodes[1,loc2glb[ni,ei]], mesh.nodes[2, mesh.elements[vi,ei]], allnodes[2,loc2glb[ni,ei]])) ||
                    (dim == 1 && is_same_node_1d(tol, vertex_dist_scale[ei], mesh.nodes[1, mesh.elements[vi,ei]], allnodes[1,loc2glb[ni,ei]]))
                    glbvertex[vi, ei] = loc2glb[ni,ei];
                end
            end
        end
        el_center ./= Np;
        
        # f2glb has duplicates. Compare to mesh faces and keep same ordering as mesh.
        # Copy normals and bdry info.
        # Set element2face map.
        if dim == 1
            test_same_face = is_same_node;
        elseif dim == 2
            test_same_face = is_same_line;
            # test_same_face = is_same_face_center; # doesn't always work?
        elseif dim == 3
            test_same_face = is_same_plane;
            # test_same_face = is_same_face_center; # doesn't always work?
        end
        
        for gfi=1:nfaces
            for fpi=1:Nfp
                tmpf2glb[fpi] = loc2glb[refel.face2local[gfi][fpi], ei];
                for di=1:dim
                    faceNodesA[di,fpi] = allnodes[di, tmpf2glb[fpi]];
                end
            end
            
            for mfi=1:nfaces
                thisfaceind = mesh.element2face[mfi, ei];
                for fpi=1:Nfp
                    tmp_nodeid = mesh.face2vertex[fpi,thisfaceind];
                    for di=1:dim
                        faceNodesB[di,fpi] = mesh.nodes[di, tmp_nodeid];
                    end
                end
                if test_same_face(faceNodesB, faceNodesA, tol, vertex_dist_scale[ei])
                    # This mesh face corresponds to this tmpf2glb face
                    # Put the tmpf2glb map into f2glb at the mesh index(thisfaceind).
                    # Move the f2glb[:,1,ind] to f2glb[:,2,ind] first if DG (Gness==2)
                    # Set element2face according to gfi(not mfi)
                    # Move face2element[1] to face2element[2] and put this one in [1]
                    for fpi=1:Nfp
                        if (Gness == 2) f2glb[fpi, 2, thisfaceind] = f2glb[fpi, 1, thisfaceind]; end
                        f2glb[fpi, 1, thisfaceind] = tmpf2glb[fpi];
                    end
                    element2face[gfi, ei] = thisfaceind;
                    face2element[2, thisfaceind] = face2element[1, thisfaceind];
                    face2element[1, thisfaceind] = ei;
                    
                    # Find the normal for every face. The normal points from e1 to e2 or outward for boundary.
                    # Note that the normal stored in mesh_data could be pointing either way.
                    thisnormal = normals[:, mfi];
                    f_center = zeros(float_type, dim);
                    for ni=1:Nfp
                        f_center += faceNodesA[:, ni];
                    end
                    f_center ./= Nfp
                    d1 = norm(f_center + thisnormal - el_center);
                    d2 = norm(f_center - thisnormal - el_center);
                    fdotn = sum((f_center-el_center) .* thisnormal);
                    # if d1 < d2 # normal is pointing in the wrong direction
                    #     thisnormal = -thisnormal;
                    # end
                    if fdotn < 0 # normal is pointing in the wrong direction
                        thisnormal = -thisnormal;
                    end
                    facenormals[:, thisfaceind] = thisnormal;
                    
                    # Copy boundary info: bdry, bdryface, bdrynorm
                    mbid = mesh.bdryID[thisfaceind];
                    gbid = indexin([mbid], bids)[1];
                    if !(gbid === nothing) # This is a boundary face
                        append!(bdry[gbid], tmpf2glb);
                        push!(bdryfc[gbid], thisfaceind);
                        facebid[thisfaceind] = gbid;
                        thisnormal = normals[:, mfi];
                        normchunk = zeros(float_type, config.dimension, Nfp);
                        for ni=1:Nfp
                            normchunk[:,ni] = thisnormal;
                        end
                        bdrynorm[gbid] = hcat(bdrynorm[gbid], normchunk);
                    end
                end
            end
        end
        
    end # element loop
    t_faces1 = Base.Libc.time() - t_faces1;
    log_entry("Vertex, face, bdry setup: "*string(t_faces1), 3);
    
    # There are duplicates in the bdry info. Remove them
    t_faces1 = Base.Libc.time();
    newbdry = Vector{Vector{Int}}(undef,nbids);
    newbdrynorm = Vector{Matrix{float_type}}(undef,nbids);
    for i=1:nbids
        newbdry[i] = zeros(int_type, length(bdry[i]));
        newbdrynorm[i] = zeros(float_type, size(bdrynorm[i]));
    end
    
    for bidi=1:nbids
        nextbdryind = 1;
        for bi=1:length(bdry[bidi])
            found = false;
            for bidj=1:bidi
                for bj=1:(nextbdryind-1)
                    if bdry[bidi][bi] == newbdry[bidj][bj]
                        # bi already exists in newbdry
                        found = true;
                        break;
                    end
                end
                if found
                    break;
                end
            end
            if !found # it's a new bdry node
                newbdry[bidi][nextbdryind] = bdry[bidi][bi];
                for di=1:dim
                    newbdrynorm[bidi][di,nextbdryind] = bdrynorm[bidi][di,bi];
                end
                nextbdryind += 1;
                
                node_bids[bdry[bidi][bi]] = bids[bidi];
            end
        end
        newbdry[bidi] = newbdry[bidi][1:(nextbdryind-1)];
        newbdrynorm[bidi] = newbdrynorm[bidi][:, 1:(nextbdryind-1)];
    end
    bdry = newbdry;
    bdrynorm = newbdrynorm;
    
    t_faces1 = Base.Libc.time() - t_faces1;
    log_entry("Remove duplicate bdry nodes: "*string(t_faces1), 3);
    
    # Refel index for each face
    for fi=1:totalfaces
        eL = face2element[1,fi];
        eR = face2element[2,fi];
        
        # Check f2glb against the face2local in refel
        for fpi=1:Nfp
            for di=1:dim
                faceNodesA[di,fpi] = allnodes[di, f2glb[fpi,1,fi]];
            end
        end
        for fi=1:refel.Nfaces
            for fpi=1:Nfp
                for di=1:dim
                    faceNodesB[di,fpi] = allnodes[di, loc2glb[refel.face2local[fi][fpi], eL]];
                end
            end
            
            if is_same_face(faceNodesA, faceNodesB, dim, tol, vertex_dist_scale[eL])
                faceRefelInd[1,fi] = fi;
                break;
            end
        end
        
        # Then do the same things for the other side of the face
        if eR > 0
            # Check f2glb against the face2local in refel
            for fi=1:refel.Nfaces
                for fpi=1:Nfp
                    for di=1:dim
                        faceNodesB[di,fpi] = allnodes[di, loc2glb[refel.face2local[fi][fpi], eR]];
                    end
                end
                
                if is_same_face(faceNodesA, faceNodesB, dim, tol, vertex_dist_scale[eL])
                    faceRefelInd[2,fi] = fi;
                    break;
                end
            end
        end
    end
    
    t_grid_from_mesh = Base.Libc.time() - t_grid_from_mesh;
    log_entry("Total grid building time: "*string(t_grid_from_mesh), 2);
    
    return (refel, Grid(finch_state.config.float_type, allnodes, bdry, bdryfc, bdrynorm, bids, node_bids, loc2glb, glbvertex, f2glb, 
                        element2face, face2element, facenormals, faceRefelInd, facebid));
end

#######################################################################################################

######       ###      ######   ######  ######  ######  ######    #####    ##    ##  #######  ######    
##   ##     ## ##     ##   ##    ##      ##      ##      ##    ###   ###  ###   ##  ##       ##    ##  
######     ##   ##    ######     ##      ##      ##      ##    ##     ##  ## ## ##  ######   ##    ##  
##        #########   ##  ##     ##      ##      ##      ##    ###   ###  ##   ###  ##       ##    ##  
##       ##       ##  ##   ##    ##    ######    ##    ######    #####    ##    ##  #######  ######    

#######################################################################################################
# This takes a full mesh and partition information in a format supplied by METIS.
# It constructs a grid that only holds this partition and ghosted neighbor elements.
#   
# Element status: 0 = owned, 1 = ghost sharing face, 2 = ghost sharing vertex, but not face
#  -------------
#  | 2 | 1 | 2 |
#  |---+---+---|
#  | 1 | 0 | 1 |
#  |---+---+---|
#  | 2 | 1 | 2 |
#  -------------
#
# For FV: only status 1 ghosts are needed.
# For FE: all ghosts are needed.
function partitioned_grid_from_mesh(mesh, epart; grid_type=CG, order=1)
    etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5]; # number of vertices for each type
    etypetodim= [1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3]; # dimension of each type
    etypetonf = [2, 3, 4, 4, 6, 5, 5, 2, 3, 4, 4, 6, 5, 5, 1, 4, 6, 5, 5]; # number of faces for element types
    etypetoftype=[0,1, 1, 2, 3, 3, 3, 0, 1, 1, 2, 3, 3, 3, 0, 0, 0, 0, 0]; # type of faces for this element type
    
    config = finch_state.config;
    log_entry("Building partitioned grid from mesh data", 2);
    t_grid_from_mesh = Base.Libc.time();
    int_type = config.index_type;
    float_type = config.float_type;
    if float_type==Float64
        tol = 1e-8;
    elseif float_type==Float32
        tol = 1e-4;
    else
        tol = 1e-4;
    end
    
    dim = config.dimension;
    ord = order;
    nfaces = etypetonf[mesh.etypes[1]]; # faces per elements
    nnodes_global = mesh.nx; # This will be updated for FE, but correct for FV
    
    # Count the owned faces and owned/ghost elements
    owned_faces = 0;
    ghost_faces = 0;
    nel_global = mesh.nel;
    nel_owned = 0;
    nel_face_ghost = 0;
    nel_node_ghost = 0; # ghosts that share nodes, but not faces
    
    element_status = fill(-1, mesh.nel); # 2 for node ghosts, 1 for face ghosts, 0 for owned, -1 for other
    face_status = fill(-1, size(mesh.normals,2)); # 1 for face to ghosts, 0 for owned, -1 for other
    mesh2grid_face = fill(-1, size(mesh.normals,2)); # -1 if this face is not in the partition, otherwise the index of the grid face
    partition2global = zeros(int_type,0); # global index of nodes
    num_neighbors = 0; # number of partitions neighboring this
    neighbor_ids = zeros(int_type,0); # ID of neighbors
    ghost_counts = zeros(int_type,0); # number of ghosts per neighbor
    # ghost_inds = [zeros(int_type,2,0)]; # local index of ghost elements
    
    # - Label each vertex with the lowest partition index that touches it.
    # - Label each element with its lowest and highest labeled vertex.
    # - Elements with lowest labels equal to this partition are owned by this partition(all nodes therein are owned).
    #    |__[1. Elements with both labels equal to this partition are fully interior to this partition.
    #       [2. Elements with highest labels > this partition have shared nodes, but are owned for indexing.
    # - 3. Elements with lowest labels < this partition have nodes that are borrowed from another partition.
    #
    # Global indices of owned nodes are set by this partition.
    # Global indices of borrowed nodes need to be determined.
    
    vertex_labels = fill(-1, mesh.nx);
    vertex_touching = fill(false, mesh.nx); # set to true if vertex touches this partition
    element_high_labels = fill(-1, mesh.nel);
    element_low_labels = fill(-1, mesh.nel);
    element_touching = fill(false, mesh.nel); # set to true if element touches this partition
    
    for ei=1:mesh.nel
        n_vert = etypetonv[mesh.etypes[ei]];
        for ni=1:n_vert
            vtx = mesh.elements[ni,ei];
            if vertex_labels[vtx] < 0 || vertex_labels[vtx] > epart[ei]
                vertex_labels[vtx] = epart[ei];
            end
            if epart[ei] == config.partition_index
                vertex_touching[vtx] = true;
            end
        end
    end
    for ei=1:mesh.nel
        n_vert = etypetonv[mesh.etypes[ei]];
        for ni=1:n_vert
            vtx = mesh.elements[ni,ei];
            if element_high_labels[ei] < 0
                element_high_labels[ei] = vertex_labels[vtx];
                element_low_labels[ei] = vertex_labels[vtx];
            else
                if vertex_labels[vtx] < element_low_labels[ei]
                    element_low_labels[ei] = vertex_labels[vtx];
                end
                if vertex_labels[vtx] > element_high_labels[ei]
                    element_high_labels[ei] = vertex_labels[vtx];
                end
            end
            if vertex_touching[vtx]
                element_touching[ei] = true;
            end
        end
    end
    
    # Set element status: owned, ghost, or other
    for ei=1:mesh.nel
        if epart[ei] == config.partition_index
            nel_owned += 1;
            element_status[ei] = 0;
        elseif element_touching[ei]
            element_status[ei] = 2; # may be set to 1 later
            nel_node_ghost += 1;
        end
    end
    
    # Find owned faces and face ghost elements
    for fi=1:size(mesh.normals,2)
        e1 = mesh.face2element[1,fi];
        e2 = mesh.face2element[2,fi];
        if element_status[e1] == 0 || (e2 > 0 && element_status[e2] == 0)
            owned_faces += 1;
            if !(element_status[e1] == 0)
                if element_status[e1] == -1
                    printerr("mixed up ghost? eid="*string(ei)*" ghost info may have errors");
                elseif element_status[e1] == 2
                    nel_face_ghost += 1;
                    nel_node_ghost -= 1;
                end
                element_status[e1] = 1;
                face_status[fi] = 1;
            elseif e2 > 0 && !(element_status[e2] == 0)
                if element_status[e2] == -1
                    printerr("mixed up ghost? eid="*string(ei)*" ghost info may have errors");
                elseif element_status[e2] == 2
                    nel_face_ghost += 1;
                    nel_node_ghost -= 1;
                end
                element_status[e2] = 1;
                face_status[fi] = 1;
            else
                face_status[fi] = 0;
            end
            mesh2grid_face[fi] = owned_faces;
        end
    end
    if grid_type == FV
        nel = nel_owned + nel_face_ghost;
    else
        nel = nel_owned;
    end
    
    # Now that ghosts are known, add their faces that are not owned
    # Only faces of face ghosts for FV are added.
    for ei=1:mesh.nel
        if element_status[ei] == 1
            for i=1:size(mesh.element2face,1)
                fi = mesh.element2face[i,ei];
                if fi > 0 && mesh2grid_face[fi] < 0 # a valid ghost face that is not owned
                    ghost_faces += 1;
                    mesh2grid_face[fi] = owned_faces + ghost_faces;
                end
            end
        end
    end
    if grid_type == FV
        totalfaces = owned_faces + ghost_faces;
    else
        totalfaces = owned_faces;
    end
    
    if dim == 1
        facenvtx = 1
    else
        facenvtx = etypetonv[etypetoftype[mesh.etypes[1]]]; # Assumes one element type
    end
    nvtx = etypetonv[mesh.etypes[1]]; # Assumes one element type
    
    log_entry("Building reference element: "*string(dim)*"D, order="*string(ord)*", nfaces="*string(nfaces), 3);
    refel = build_refel(dim, ord, nfaces, config.elemental_nodes);
    
    if grid_type == DG
        Gness = 2;
    else
        Gness = 1;
    end
    
    vertex_dist_scale = zeros(float_type, mesh.nel); # a scale for relative tolerance
    
    Np = refel.Np;                      # number of nodes per element
    bdry = [];                          # index(in x) of boundary nodes for each BID
    bdryfc = [];                        # index of faces touching each BID
    bdrynorm = [];                      # normal at boundary nodes
    bids = collectBIDs(mesh);           # BID list
    nbids = length(bids);
    for i=1:nbids
        push!(bdry, zeros(int_type, 0));
        push!(bdryfc, zeros(int_type,0));
        push!(bdrynorm, zeros(float_type, config.dimension,0));
    end
    loc2glb = zeros(int_type, Np, nel)           # local to global index map for each element's nodes
    glbvertex = zeros(int_type, nvtx, nel);      # local to global for vertices
    f2glb = zeros(int_type, refel.Nfp[1], Gness, totalfaces);  # face node local to global
    element2face = zeros(int_type, nfaces, nel); # element to face map
    face2element = zeros(int_type, 2, totalfaces); # face to element map
    facenormals = zeros(float_type, dim, totalfaces);   # normal vectors for every face
    faceRefelInd = zeros(int_type, 2, totalfaces); # Index in refel for this face for elements on both sides
    facebid = zeros(int_type, totalfaces);       # BID of each face
    
    partition2global_element = zeros(int_type, nel); # maps partition elements to global mesh elements
    
    if grid_type == FV
        tmpallnodes = zeros(float_type, dim, nel*refel.Np);
        element_owners = fill(-1, nel) # partition number of each ghost, or -1 for owned elements
    else
        tmpallnodes = zeros(float_type, dim, (nel + nel_node_ghost + nel_face_ghost) * refel.Np);
        loc2glb = zeros(int_type, Np, nel + nel_node_ghost + nel_face_ghost); # This will be trimmed later
        element_owners = fill(-1, nel + nel_node_ghost + nel_face_ghost) # partition number of each element
    end
    
    # compute node coordinates
    next_e_index = 1;
    next_g_index = nel_owned + 1; #index for next ghost
    t_nodes1 = Base.Libc.time();
    for ei=1:mesh.nel
        if (element_status[ei] == 0 || # owned
            (element_status[ei] == 1 && grid_type == FV) || # face_ghost && FV
            (element_status[ei] > 0 && !(grid_type == FV))) # any ghost && FE
            
            if element_status[ei] == 0 # owned
                next_index = next_e_index;
                next_e_index += 1;
            else # ghost
                next_index = next_g_index;
                next_g_index += 1;
            end
            
            # Find this element's nodes
            n_vert = etypetonv[mesh.etypes[ei]];
            e_vert = mesh.nodes[1:dim, mesh.elements[1:n_vert, ei]];
            
            # sample distance betweeen vertices
            for i=2:n_vert
                vertex_dist_scale[next_index] += sum(abs.(e_vert[:,i] .- e_vert[:,i-1]));
            end
            vertex_dist_scale[next_index] /= (n_vert-1);
            
            if dim == 1
                e_x = line_refel_to_x(refel.r[:,1], e_vert);
            elseif dim == 2
                if nfaces == 3 # triangles
                    (e_x, e_y) = triangle_refel_to_xy_(refel.r[:,1], refel.r[:,2], e_vert);
                else # quads
                    (e_x, e_y) = quad_refel_to_xy(refel.r[:,1], refel.r[:,2], e_vert);
                end
            elseif dim == 3
                if nvtx == 8 # hexes
                    (e_x, e_y, e_z) = hex_refel_to_xyz(refel.r[:,1], refel.r[:,2], refel.r[:,3], e_vert);
                else # tets
                    (e_x, e_y, e_z) = tetrahedron_refel_to_xyz(refel.r[:,1], refel.r[:,2], refel.r[:,3], e_vert);
                end
            end
            
            
            # Add them to the tmp global nodes
            tmpallnodes[1, ((next_index-1)*Np+1):(next_index*Np)] = e_x;
            if dim > 1
                tmpallnodes[2, ((next_index-1)*Np+1):(next_index*Np)] = e_y;
                if dim > 2
                    tmpallnodes[3, ((next_index-1)*Np+1):(next_index*Np)] = e_z;
                end
            end
            
            # temporary mapping
            loc2glb[:,next_index] = ((next_index-1)*Np+1):(next_index*Np);
            
            if grid_type == FV
                partition2global_element[next_index] = ei;
                
                # ghost info
                if element_status[ei] == 1 # ghost
                    element_owners[next_index] = epart[ei];
                    list_ind = 0;
                    for i=1:num_neighbors
                        if neighbor_ids[i] == epart[ei]
                            list_ind += i; # += used for scope of list_ind
                        end
                    end
                    if list_ind == 0
                        # A new neighbor is found
                        num_neighbors += 1;
                        push!(neighbor_ids, epart[ei]);
                    end
                end
                
            else ### FE ###
                element_owners[next_index] = epart[ei];
                if element_status[ei] == 0 # owned
                    partition2global_element[next_index] = ei;
                end
            end
        end
    end
    t_nodes1 = Base.Libc.time() - t_nodes1;
    log_entry("Initial node setup: "*string(t_nodes1), 3);
    
    if Gness == 1 # CG
        # Go back and remove duplicate nodes. Adjust loc2glb.
        t_nodes1 = Base.Libc.time();
        (allnodes, loc2glb) = remove_duplicate_nodes(tmpallnodes, loc2glb, tol=tol);
        t_nodes1 = Base.Libc.time() - t_nodes1;
        log_entry("Remove duplicate nodes: "*string(t_nodes1), 3);
    else
        allnodes = tmpallnodes; # DG grid is already made
    end
    
    # bid for every node start with 0 for interior
    node_bids = zeros(int_type, size(allnodes,2));
    
    # vertices, faces and boundary
    # first make a map from mesh faces to grid faces
    next_e_index = 1;
    next_g_index = nel_owned+1; #index for next ghost
    t_faces1 = Base.Libc.time();
    for ei=1:mesh.nel
        if element_status[ei] == 0 || (element_status[ei] == 1 && grid_type == FV)# owned or (ghost && FV) element
            n_vert = etypetonv[mesh.etypes[ei]];
            mfids = mesh.element2face[:,ei];
            normals = mesh.normals[:,mfids];
            el_center = zeros(float_type, dim);
            
            if element_status[ei] == 0 # owned
                next_index = next_e_index;
                next_e_index += 1;
            else # ghost
                next_index = next_g_index;
                next_g_index += 1;
            end
            
            # vertices and center
            for ni=1:Np
                el_center += allnodes[:,loc2glb[ni,next_index]];
                
                for vi=1:n_vert
                    if is_same_node(mesh.nodes[:, mesh.elements[vi,ei]], allnodes[:,loc2glb[ni,next_index]], tol, vertex_dist_scale[next_index])
                        glbvertex[vi, next_index] = loc2glb[ni,next_index];
                    end
                end
            end
            el_center ./= Np;
            
            #  Compare to mesh faces and keep same ordering as mesh.
            # Copy normals and bdry info.
            # Set element2face map.
            meshfaces = mesh.element2face[:,ei];    # index from mesh
            if dim == 1
                test_same_face = is_same_node;
            elseif dim == 2
                # test_same_face = is_same_line;
                test_same_face = is_same_face_center;
            elseif dim == 3
                # test_same_face = is_same_plane;
                test_same_face = is_same_face_center;
            end
            
            for gfi=1:nfaces
                tmpf2glb = loc2glb[refel.face2local[gfi], next_index];
                found_face = false;
                for mfi=1:nfaces
                    if !found_face
                        meshfaceind = meshfaces[mfi];
                        gridfaceind = mesh2grid_face[meshfaceind]; # the face index in the grid
                        if gridfaceind > 0 && test_same_face(mesh.nodes[:,mesh.face2vertex[:,meshfaceind]], allnodes[:, tmpf2glb], tol, vertex_dist_scale[next_index])
                            found_face = true;
                            # This mesh face corresponds to this tmpf2glb face
                            # Put the tmpf2glb map into f2glb at the mesh index(meshfaceind).
                            # Move the f2glb[:,1,ind] to f2glb[:,2,ind] first if DG (Gness==2)
                            # Set element2face according to gfi(not mfi)
                            # Move face2element[1] to face2element[2] and put this one in [1]
                            if (Gness == 2) f2glb[:, 2, gridfaceind] = f2glb[:, 1, gridfaceind]; end
                            f2glb[:, 1, gridfaceind] = tmpf2glb;
                            element2face[gfi, next_index] = gridfaceind;
                            face2element[2, gridfaceind] = face2element[1, gridfaceind];
                            face2element[1, gridfaceind] = next_index;
                            
                            # Find the normal for every face. The normal points from e1 to e2 or outward for boundary.
                            # Note that the normal stored in mesh_data could be pointing either way.
                            thisnormal = normals[:, mfi];
                            f_center = zeros(float_type, dim);
                            for ni=1:length(tmpf2glb)
                                f_center += allnodes[:, tmpf2glb[ni]];
                            end
                            f_center ./= length(tmpf2glb)
                            fdotn = sum((f_center-el_center) .* thisnormal);
                            if fdotn < 0 # normal is pointing in the wrong direction
                                thisnormal = -thisnormal;
                            end
                            facenormals[:, gridfaceind] = thisnormal;
                            
                            # Copy boundary info: bdry, bdryface, bdrynorm
                            mbid = mesh.bdryID[meshfaceind];
                            gbid = indexin([mbid], bids)[1];
                            nfacenodes = length(tmpf2glb);
                            if !(gbid === nothing) && gridfaceind <= owned_faces # This is a boundary face of an owned element
                                append!(bdry[gbid], tmpf2glb);
                                push!(bdryfc[gbid], gridfaceind);
                                facebid[gridfaceind] = gbid;
                                thisnormal = normals[:, mfi];
                                normchunk = zeros(float_type, config.dimension, nfacenodes);
                                for ni=1:nfacenodes
                                    normchunk[:,ni] = thisnormal;
                                end
                                bdrynorm[gbid] = hcat(bdrynorm[gbid], normchunk);
                            end
                        end
                    end
                end
                if !found_face
                    printerr("couldn't find a match for a face");
                end
            end
        end
    end # element loop
    t_faces1 = Base.Libc.time() - t_faces1;
    log_entry("Vertex, face, bdry setup: "*string(t_faces1), 3);
    
    # There are duplicates in the bdry info. Remove them
    t_faces1 = Base.Libc.time();
    newbdry = similar(bdry);
    newbdrynorm = similar(bdrynorm);
    for i=1:length(bdry)
        newbdry[i] = zeros(int_type, length(bdry[i]));
        newbdrynorm[i] = zeros(float_type, size(bdrynorm[i]));
    end
    
    for bidi=1:length(bids)
        nextbdryind = 1;
        for bi=1:length(bdry[bidi])
            found = false;
            for bidj=1:bidi
                for bj=1:length(newbdry[bidj])
                    if bdry[bidi][bi] == newbdry[bidj][bj]
                        # bi already exists in newbdry
                        found = true;
                        break;
                    end
                end
                if found
                    break;
                end
            end
            if !found # it's a new bdry node
                newbdry[bidi][nextbdryind] = bdry[bidi][bi];
                newbdrynorm[bidi][:,nextbdryind] = bdrynorm[bidi][:,bi];
                nextbdryind += 1;
                
                node_bids[bdry[bidi][bi]] = bids[bidi];
            end
        end
        newbdry[bidi] = newbdry[bidi][1:nextbdryind-1];
        newbdrynorm[bidi] = newbdrynorm[bidi][:, 1:nextbdryind-1];
    end
    bdry = newbdry;
    bdrynorm = newbdrynorm;
    
    t_faces1 = Base.Libc.time() - t_faces1;
    log_entry("Remove duplicate bdry nodes: "*string(t_faces1), 3);
    
    # Refel index for each face
    for fi=1:totalfaces
        eL = face2element[1,fi];
        eR = face2element[2,fi];
        faceRefelInd[1,fi] = which_refel_face(f2glb[:,1,fi], allnodes, refel, loc2glb[:,eL], tol, vertex_dist_scale[eL]);
        if eR > 0
            faceRefelInd[2,fi] = which_refel_face(f2glb[:,Gness,fi], allnodes, refel, loc2glb[:,eR], tol, vertex_dist_scale[eR]);
        end
    end
    
    #### GHOST INFO ######################################################################
    # if more than one proc uses the same mesh partition, neighbor ID adjustment is needed.
    # 0, 1, 2, ...  ->   0, 0, 0, 1, 1, 1, 2, 2, 2, ...
    # |__^               |________^     
    procs_per_partition = 1;
    proc_offset = 0;
    if config.num_procs >= (config.num_partitions * 2)
        procs_per_partition = Int(floor(config.num_procs / config.num_partitions));
        proc_offset = config.proc_rank - config.partition_index * procs_per_partition;
    end
    
    #### FV ONLY ####
    # Form ghost pairs for send/recv
    # First count how many pairs are needed for each neighbor
    if grid_type == FV
        ghost_counts = zeros(int_type, num_neighbors)
        for fi=1:owned_faces
            e1 = face2element[1,fi];
            e2 = face2element[2,fi];
            eg = e2;
            eo = e1;
            neigh = -1;
            if element_owners[e1] < 0
                if e2 > 0 && element_owners[e2] >= 0
                    # e1 owned, e2 ghost
                    neigh = element_owners[e2];
                end
            elseif e2 > 0 && element_owners[e2] < 0
                # e2 owned, e1 ghost
                eg = e1;
                eo = e2;
                neigh = element_owners[e1];
            end
            
            if neigh >= 0 # It is a face between neighbors
                list_ind = 0;
                for ni=1:num_neighbors
                    if neighbor_ids[ni] == neigh;
                        list_ind += ni; # += used for scope of list_ind
                            break;
                    end
                end
                if list_ind > 0
                    ghost_counts[list_ind] += 1;
                end
            end
        end
        
        # Then build the lists of ghost pairs for each neighbor
        ghost_inds = Vector{Array{Int,2}}(undef, num_neighbors);
        for ni=1:num_neighbors
            ghost_inds[ni] = zeros(int_type, 2, ghost_counts[ni]);
        end
        # Need to loop over mesh faces to be sure they are built in the same order on each partition
        for i=1:length(mesh2grid_face)
            fi = mesh2grid_face[i];
            if fi > 0  && fi <= owned_faces # exists in this partition and at least one of e1,e2 are owned
                e1 = face2element[1,fi];
                e2 = face2element[2,fi];
                if element_owners[e1] >= 0 || (e2>0 && element_owners[e2] >= 0) # It is a face between owned and ghost
                    if element_owners[e1] >= 0 # e1 is ghost
                        eg = e1; # ghost
                        eo = e2; # owned
                    elseif element_owners[e2] >= 0 # e2 is ghost
                        eg = e2;
                        eo = e1;
                    end
                    list_ind = 0;
                    for ni=1:num_neighbors
                        if neighbor_ids[ni] == element_owners[eg];
                            list_ind += ni; # += used for scope of list_ind
                        end
                    end
                    if list_ind > 0
                        for j=1:ghost_counts[list_ind]
                            if ghost_inds[list_ind][1,j] == 0
                                ghost_inds[list_ind][:,j] = [eg, eo];
                                break;
                            end
                        end
                    else
                        # There was a problem...
                        printerr("Problem with building ghost lists. Expect an error.")
                    end
                end
            end
        end
        # Adjust neighbor IDs from partition index to proc rank
        for ni=1:num_neighbors
            neighbor_ids[ni] = neighbor_ids[ni] * procs_per_partition + proc_offset;
        end
    end#### FV ONLY ####
    
    #### FE ONLY ####
    if !(grid_type == FV)
        # Build the nodal partition to global map.
        # This requires finding the number of nodes for all lower number partitions, which is a significant task.
        #
        # Global indices of owned nodes are set by this partition.
        # Global indices of borrowed nodes need to be determined.
        
        nnodes_with_ghosts = size(allnodes, 2);
        node_owner_index = fill(config.num_partitions+1, nnodes_with_ghosts); # will be set to lowest partition index that has each node
        node_shared_flag = fill(false, nnodes_with_ghosts); # When a node is shared by another partition, flag it (1)
        for ei=1:size(loc2glb,2)
            for ni=1:size(loc2glb,1)
                nid = loc2glb[ni,ei];
                node_owner_index[nid] = min(node_owner_index[nid], element_owners[ei]);
                if !(element_owners[ei] == config.partition_index)
                    node_shared_flag[nid] = true;
                end
            end
        end
        
        # Need to trim the extra nodes off of allnodes and discard ghosts in loc2glb
        node_keep_flag = fill(false, nnodes_with_ghosts);
        highest_kept = 0;
        for ei=1:nel # only owned elements
            if !(element_owners[ei] == config.partition_index)
                printerr("Ghost element mixed in with owned elements? eid = "*string(ei)*" part = "*string(config.partition_index)*", owner = "*string(element_owners[ei]))
            end
            for ni=1:size(loc2glb,1)
                nid = loc2glb[ni,ei];
                node_keep_flag[nid] = true;
                highest_kept = max(highest_kept, nid);
            end
        end
        # Make sure the kept nodes are all before the discards
        for ni=1:highest_kept
            if !node_keep_flag[ni]
                printerr("Nodes are out of order while trimming ghosts. discard = "*string(ni)*", highest kept = "*string(highest_kept), fatal=true);
            end
            if node_owner_index[ni] > config.partition_index
                printerr("Kept node had owner > part. node = "*string(ni)*" part = "*string(config.partition_index)*", owner = "*string(node_owner_index[ni]));
            end
        end
        # If that worked, we can safely trim allnodes and loc2glb
        allnodes = allnodes[:,1:highest_kept];
        loc2glb = loc2glb[:,1:nel];
        element_owners = element_owners[1:nel];
        node_owner_index = node_owner_index[1:highest_kept];
        node_bids = node_bids[1:highest_kept];
        
        # Count nodes that are owned, shared or borrowed
        nnodes_owned = highest_kept;
        nnodes_shared = 0;
        nnodes_borrowed = 0;
        for ni=1:highest_kept
            if node_shared_flag[ni]
                nnodes_owned -= 1;
                if node_owner_index[ni] == config.partition_index
                    nnodes_shared += 1;
                else
                    nnodes_borrowed += 1;
                end
            end
        end
        
        # The ultimate goal
        partition2global = zeros(int_type, highest_kept); # zero means it hasn't been determined yet
        
        # Have all processes gather their [partition index, nodes_owned, nodes_shared, nodes_borrowed].
        p_data = zeros(Int, config.num_procs * 4);
        p_data_in = MPI.UBuffer(p_data, 4, config.num_procs, MPI.Datatype(Int));
        p_data_out = [config.partition_index, nnodes_owned, nnodes_shared, nnodes_borrowed];
        MPI.Allgather!(p_data_out, p_data_in, MPI.COMM_WORLD);
        
        start_offset = 0;
        nnodes_global = 0;
        owned_nodes_per_partition = zeros(int_type, config.num_partitions);
        shared_nodes_per_partition = zeros(int_type, config.num_partitions);
        borrowed_nodes_per_partition = zeros(int_type, config.num_partitions);
        nodes_to_share_per_proc = zeros(int_type, config.num_procs);
        proc2partition = zeros(int_type, config.num_procs);
        total_shared_node_size = 0;
        i_need_to_send_shared = false;
        for proc_i = 1:config.num_procs
            part_id = p_data[(proc_i-1)*4+1];
            proc2partition[proc_i] = part_id;
            
            if shared_nodes_per_partition[part_id + 1] == 0
                total_shared_node_size += p_data[(proc_i-1)*4+3]; # The number of all shared nodes to communicate
                if proc_i == config.proc_rank+1
                    scoper = i_need_to_send_shared;
                    i_need_to_send_shared = true; # This proc will be responsible for sharing for this partition
                end
                
                nodes_to_share_per_proc[proc_i] = p_data[(proc_i-1)*4+3];
            else
                nodes_to_share_per_proc[proc_i] = 1; # a dummy to make things work
                total_shared_node_size += 1;
            end
            
            owned_nodes_per_partition[part_id + 1] = p_data[(proc_i-1)*4+2];
            shared_nodes_per_partition[part_id + 1] = p_data[(proc_i-1)*4+3];
            borrowed_nodes_per_partition[part_id + 1] = p_data[(proc_i-1)*4+4];
        end
        # Determine the starting index of this partition's nodes
        for part_i = 1:config.num_partitions
            if part_i < config.partition_index + 1
                start_offset += owned_nodes_per_partition[part_i] + shared_nodes_per_partition[part_i];
            end
            nnodes_global += owned_nodes_per_partition[part_i] + shared_nodes_per_partition[part_i];
        end
        
        # The global index for owned nodes will simply add the start_offset
        next_index = start_offset+1;
        for ni=1:size(allnodes, 2)
            if node_owner_index[ni] == config.partition_index
                partition2global[ni] = next_index;
                next_index += 1;
            end
        end
        
        # Now only the shared nodes' global indices have not been set.
        # ...
        # Collect all of the shared nodes coords and indices.
        # Then each process will search for their borrowed nodes and set the global index.
        # The buffer must be total_shared_node_size * (dimension+1)
        chunk_sizes = nodes_to_share_per_proc .* (config.dimension + 1);
        displacements = zeros(int_type, config.num_procs); # for the irregular gatherv
        d = 0;
        for proc_i=1:config.num_procs
            displacements[proc_i] = d;
            d += nodes_to_share_per_proc[proc_i] * (config.dimension + 1);
        end
        p_data = zeros(total_shared_node_size * (config.dimension + 1));
        p_data_in = MPI.VBuffer(p_data, chunk_sizes, displacements, MPI.Datatype(Float64));
        if i_need_to_send_shared
            p_data_out = zeros(float_type, (config.dimension + 1), nnodes_shared);
            next_index = 1;
            for ni=1:size(allnodes, 2)
                if node_shared_flag[ni] && node_owner_index[ni] == config.partition_index
                    p_data_out[1:config.dimension, next_index] = allnodes[:,ni];
                    p_data_out[config.dimension+1, next_index] = partition2global[ni];
                    next_index += 1;
                end
            end
            
        else
            p_data_out = fill(-1, (config.dimension + 1), 1);
        end
        
        MPI.Allgatherv!(p_data_out, p_data_in, MPI.COMM_WORLD);
        
        # Now search for shared coordinates of borrowed nodes and set the global index
        for ni=1:size(allnodes, 2)
            if node_owner_index[ni] < config.partition_index
                these_coords = allnodes[:,ni];
                found_index = -1; # will be set to the index when found
                for proc_i=1:config.num_procs
                    if node_owner_index[ni] == proc2partition[proc_i]
                        p_offset = displacements[proc_i] + 1;
                        for nj=1:nodes_to_share_per_proc[proc_i]
                            n_offset = p_offset + (nj-1) * (config.dimension + 1);
                            those_coords = p_data[n_offset:(n_offset+config.dimension-1)];
                            # Are they the same node?
                            if p_data[n_offset+config.dimension] >= 0 && sum(abs.(these_coords-those_coords)) < 1e-10
                                found_index = p_data[n_offset+config.dimension];
                                break;
                            end
                        end
                        if found_index >= 0
                            break;
                        end
                    end
                end
                if found_index < 0
                    printerr("proc "*string(config.proc_rank)*" Didn't find a global index for "*string(ni)*" owned by "*string(node_owner_index[ni])*
                                " coords= "*string(these_coords), fatal=true);
                end
                partition2global[ni] = Int(found_index);
            end
        end
        
        # Finally, armed with a complete partition2global map, flag each boundary node 
        # so that boundary conditions don't get messed up when gathering the system.
        nbdrynodes = 0;
        for bi=1:length(bids)
            nbdrynodes += length(bdry[bi]);
        end
        bidmap = zeros(int_type, 2, nbdrynodes);
        next_ind = 1;
        for bi=1:length(bids)
            for ni=1:length(bdry[bi])
                bidmap[1,next_ind] = partition2global[bdry[bi][ni]];
                bidmap[2,next_ind] = bi;
                next_ind += 1;
            end
        end
        
        # gather the numbers of bdry nodes
        p_data = zeros(int_type, config.num_procs);
        p_data_in = MPI.UBuffer(p_data, 1, config.num_procs, MPI.Datatype(Int));
        p_data_out = [nbdrynodes * 2];
        MPI.Allgather!(p_data_out, p_data_in, MPI.COMM_WORLD);
        
        # gather the bidmaps
        chunk_sizes = p_data;
        displacements = zeros(int_type, config.num_procs); # for the irregular gatherv
        d = 0;
        for proc_i=1:config.num_procs
            displacements[proc_i] = d;
            d += p_data[proc_i];
        end
        p_data = zeros(int_type, d);
        p_data_in = MPI.VBuffer(p_data, chunk_sizes, displacements, MPI.Datatype(Int));
        p_data_out = bidmap[:];
        
        MPI.Allgatherv!(p_data_out, p_data_in, MPI.COMM_WORLD);
        
        global_bdry_flag = zeros(Int8, nnodes_global);
        for ni=1:2:d
            global_bdry_flag[p_data[ni]] = max(global_bdry_flag[p_data[ni]], p_data[ni+1]); # keep highest values
        end
        
        # There is a chance that a partition owns a boundary node without owning an 
        # adjoining boundary face.
        bdry_adjustment = zeros(int_type,0);
        for ni=1:length(partition2global)
            if global_bdry_flag[partition2global[ni]] > 0
                # make sure it is in bdry
                foundit = false;
                for bi=1:length(bdry)
                    for bni=1:length(bdry[bi])
                        if bdry[bi][bni] == ni
                            scoper = foundit;
                            foundit = true;
                            break;
                        end
                    end
                    if foundit
                        break;
                    end
                end
                
                if !foundit
                    log_entry("Part "*string(config.partition_index)*" has a boundary node without adjoining boundary face!");
                    # Figure out who owns a boundary face with this node and grant it unto them.
                    for proc_i=1:config.num_procs
                        for pni=(displacements[proc_i]+1):(displacements[proc_i]+chunk_sizes[proc_i]-1)
                            if p_data[pni] == partition2global[ni]
                                # Found a proc that knows this boundary
                                push!(bdry_adjustment, proc2partition[proc_i]); # The partition to give it to
                                push!(bdry_adjustment, partition2global[ni]);   # The global index to give
                                scoper = foundit;
                                foundit = true;
                                break;
                            end
                        end
                        if foundit
                            break;
                        end
                    end
                    if foundit
                        # change node owner
                        node_owner_index[ni] = bdry_adjustment[end-1];
                        nnodes_shared -= 1;
                        nnodes_borrowed += 1;
                    end
                end
            end
        end
        
        # Share bdry_adjustment and change owners as needed
        # gather the numbers of adjusted nodes
        p_data = zeros(int_type, config.num_procs);
        p_data_in = MPI.UBuffer(p_data, 1, config.num_procs, MPI.Datatype(Int));
        p_data_out = [length(bdry_adjustment)];
        MPI.Allgather!(p_data_out, p_data_in, MPI.COMM_WORLD);
        
        # gather the adjustments
        chunk_sizes = p_data;
        displacements = zeros(int_type, config.num_procs); # for the irregular gatherv
        d = 0;
        for proc_i=1:config.num_procs
            displacements[proc_i] = d;
            d += p_data[proc_i];
        end
        p_data = zeros(int_type, d);
        p_data_in = MPI.VBuffer(p_data, chunk_sizes, displacements, MPI.Datatype(Int));
        p_data_out = bdry_adjustment;
        
        MPI.Allgatherv!(p_data_out, p_data_in, MPI.COMM_WORLD);
        
        for i=1:2:d
            if p_data[i] == config.partition_index
                # This node was given unto me
                # Find the local index and change owner
                setit = false;
                for ni=1:length(partition2global)
                    if partition2global[ni] == p_data[ni+1] && node_owner_index[ni] != config.partition_index
                        node_owner_index[ni] = config.partition_index;
                        nnodes_shared += 1;
                        nnodes_borrowed -= 1;
                        scoper = setit;
                        setit = true;
                        break;
                    end
                end
                if !setit
                    #printerr("Problem with adjusting boundary nodes. Couldn't find one that was given to me.")
                end
            end
        end
        
    end#### FE ONLY ####
    
    t_grid_from_mesh = Base.Libc.time() - t_grid_from_mesh;
    log_entry("Partitioned grid building time: "*string(t_grid_from_mesh), 2);
    
    if grid_type == FV
        return (refel, Grid(finch_state.config.float_type, allnodes, bdry, bdryfc, bdrynorm, bids, node_bids, loc2glb, glbvertex, f2glb, element2face, 
            face2element, facenormals, faceRefelInd, facebid, 
            true, Array(1:nel_owned), nel_global, nel_owned, nel_face_ghost, owned_faces, ghost_faces, nnodes_global, 0, element_owners, zeros(int_type,0), partition2global_element, zeros(int_type,0), 
            zeros(Int8, 0), num_neighbors, neighbor_ids, ghost_counts, ghost_inds));
    else
        return (refel, Grid(finch_state.config.float_type, allnodes, bdry, bdryfc, bdrynorm, bids, node_bids, loc2glb, glbvertex, f2glb, element2face, 
            face2element, facenormals, faceRefelInd, facebid, 
            true, Array(1:nel_owned), nel_global, nel_owned, 0, owned_faces, 0, nnodes_global, nnodes_borrowed, zeros(int_type,0), node_owner_index, partition2global_element, partition2global, 
            global_bdry_flag, num_neighbors, neighbor_ids, zeros(int_type,0), [zeros(int_type,0,0)]));
    end
end

#########################################################################

##    ##  ######  ######  ##       #####
##    ##    ##      ##    ##      ###   #
##    ##    ##      ##    ##        ###
##    ##    ##      ##    ##      #   ###
  ####      ##    ######  ######   #####
  
#########################################################################

function collectBIDs(mesh::MeshData)
    bids = zeros(Int,0);
    N = length(mesh.bdryID);
    nbids = 0;
    for i=1:N
        if mesh.bdryID[i] > 0
            already = false;
            for j=1:nbids
                if mesh.bdryID[i] == bids[j]
                    already = true;
                end
            end
            if !already
                push!(bids, mesh.bdryID[i]);
                nbids += 1;
            end
        end
    end
    return bids;
end

# Removes duplicate nodes and updates local to global maps
function remove_duplicate_nodes(nodes::Matrix, loc2glb::Matrix; tol=1e-12, scale=1, depth=5, mincount=50, other2glb=[])
    # defaults
    # depth = 5; # 32768 for 3D
    # mincount = 50; # don't subdivide if less than this
    int_type = finch_state.config.index_type;
    float_type = finch_state.config.float_type;
    
    # Strategy: divide nodes into bins compare against nodes in bin
    tmpNnodes = size(nodes,2);
    dim = size(nodes,1);
    replace_with = zeros(Int, tmpNnodes); # holds index of replacement node being kept. 0 if this is kept.
    abins = zeros(Int, tmpNnodes);
    bin_ends = [tmpNnodes]; # The last index of each bin
    
    # init bins and find extrema
    xlim = [1e10, -1e10];
    ylim = [1e10, -1e10];
    zlim = [1e10, -1e10];
    for i=1:tmpNnodes
        abins[i] = i;
        x = nodes[1,i];
        y = nodes[2,i];
        z = nodes[3,i];
        xlim[1] = min(xlim[1], x);
        xlim[2] = max(xlim[2], x);
        if dim>1
            ylim[1] = min(ylim[1], y);
            ylim[2] = max(ylim[2], y);
            if dim > 2
                zlim[1] = min(zlim[1], z);
                zlim[2] = max(zlim[2], z);
            end
        end
    end
    
    if dim == 1
        (abins, bin_ends, cbin) = partition_nodes_in_bins_1d(nodes, abins, xlim, depth, mincount);
    elseif dim == 2
        (abins, bin_ends, cbin) = partition_nodes_in_bins_2d(nodes, abins, xlim, ylim, depth, mincount);
    else # dim==3
        (abins, bin_ends, cbin) = partition_nodes_in_bins_3d(nodes, abins, xlim, ylim, zlim, depth, mincount);
    end
    
    nbins::Int = length(bin_ends);
    remove_count = 0;
    startind = 1;
    lastind = 1;
    # Loop over nodes in bins and put numbers in replace_with where duplicates will be removed.
    for bini = 1:nbins
        lastind = bin_ends[bini];
        for ni=startind:lastind
            if replace_with[abins[ni]] == 0 # It may have already been handled
                for nj=startind:(ni-1)
                    # If node[nj] == node[ni], keep nj, remove ni
                    if (dim == 3 && is_same_node_3d(tol, scale, nodes[1,abins[ni]], nodes[1,abins[nj]], nodes[2,abins[ni]], nodes[2,abins[nj]], nodes[3,abins[ni]], nodes[3,abins[nj]])) ||
                        (dim == 2 && is_same_node_2d(tol, scale, nodes[1,abins[ni]], nodes[1,abins[nj]], nodes[2,abins[ni]], nodes[2,abins[nj]])) ||
                        (dim == 1 && is_same_node_1d(tol, scale, nodes[1,abins[ni]], nodes[1,abins[nj]]))
                        
                        remove_count += 1;
                        replace_with[abins[ni]] = abins[nj];
                        break;
                    end
                end
            end
        end
        startind = bin_ends[bini] + 1;
    end
    # Then check cbin against itself for the rare case of split nodes
    # println("part "*string(finch_state.config.partition_index)*" cbin "*string(length(cbin)));
    ncbin = length(cbin);
    for ni=2:ncbin
        if replace_with[cbin[ni]] == 0 # It may have already been handled
            for nj=1:ni-1
                if (dim == 3 && is_same_node_3d(tol, scale, nodes[1,cbin[ni]], nodes[1,cbin[nj]], nodes[2,cbin[ni]], nodes[2,cbin[nj]], nodes[3,cbin[ni]], nodes[3,cbin[nj]])) ||
                    (dim == 2 && is_same_node_2d(tol, scale, nodes[1,cbin[ni]], nodes[1,cbin[nj]], nodes[2,cbin[ni]], nodes[2,cbin[nj]])) ||
                    (dim == 1 && is_same_node_1d(tol, scale, nodes[1,cbin[ni]], nodes[1,cbin[nj]]))
                    # make sure I'm not replacing with something replaced by this
                    if replace_with[cbin[nj]] > 0
                        tmp = replace_with[cbin[nj]];
                        cycle_counter = 0;
                        while replace_with[tmp] > 0 && cycle_counter < 100
                            tmp = replace_with[tmp];
                            cycle_counter += 1;
                        end
                        if tmp == cbin[ni]
                            continue; # don't replace i if j is replaced by i...
                        else
                            remove_count += 1;
                            replace_with[cbin[ni]] = tmp;
                            break;
                        end
                    else
                        remove_count += 1;
                        replace_with[cbin[ni]] = cbin[nj];
                        break;
                    end
                end
            end
        end
    end
    # Since cbin could have caused a chain of replacements, check each one
    for i=1:tmpNnodes
        if replace_with[i] > 0
            tmp = replace_with[i];
            cycle_counter = 0;
            while replace_with[tmp] > 0 && cycle_counter < 100
                tmp = replace_with[tmp];
                cycle_counter += 1;
            end
            if cycle_counter == 100 # a cycle was detected, let this one be kept
                replace_with[i] = 0;
                remove_count -= 1;
            else
                replace_with[i] = tmp;
            end
        end
    end
    
    # println(string(length(bin_ends))*" bins: " )
    # println(bin_ends);
    # println("Nnodes before: "*string(tmpNnodes)*" after: "*string(tmpNnodes-remove_count));
    
    # Loop over nodes and place in new allnodes array while adjusting replace_with
    Nnodes = tmpNnodes-remove_count;
    next_ind = 1;
    newnodes = zeros(float_type, dim, Nnodes);
    new_homes = zeros(Int, tmpNnodes);
    for i=1:tmpNnodes
        if replace_with[i] == 0
            # A node to keep
            for j=1:dim
                newnodes[j,next_ind] = nodes[j,i];
            end
            new_homes[i] = next_ind;
            next_ind += 1;
        end
    end
    for i=1:tmpNnodes
        if replace_with[i] > 0
            # A node that was removed. update replace_with with new_homes
            replace_with[i] = new_homes[replace_with[i]];
        end
    end
    
    # Loop over each element's nodes and adjust loc2glb
    nel = size(loc2glb,2);
    Np = size(loc2glb,1);
    for ei=1:nel
        for ni=1:Np
            ind = loc2glb[ni,ei];
            if replace_with[ind] > 0
                loc2glb[ni,ei] = replace_with[ind];
            else
                loc2glb[ni,ei] = new_homes[ind];
            end
        end
    end
    
    # If another 2glb map is present, do the same to it
    if length(other2glb) > 0
        for i=1:length(other2glb) # 
            ind = other2glb[i];
            if replace_with[ind] > 0
                other2glb[i] = replace_with[ind];
            else
                other2glb[i] = new_homes[ind];
            end
        end
        
        return (newnodes, loc2glb, other2glb);
    end
    
    return (newnodes, loc2glb);
end

# recursively subdivide bins of nodes
function partition_nodes_in_bins_1d(nodes, abin, xlim, depth, mincount)
    halfx = (xlim[1] + xlim[2])/2 + 1e-11;
    N = length(abin);
    
    bbin = similar(abin); # Temporary storage
    cbin = zeros(Int, 0); # These are ones too close to the cutoff that will be added to an extra bin
    lcount = 0;
    rcount = 0;
    for i=1:N
        if nodes[1,abin[i]] <= halfx
            lcount += 1;
            bbin[lcount] = abin[i];
        else
            rcount += 1;
            bbin[N-rcount+1] = abin[i];
        end
        if abs(nodes[1,abin[i]] - halfx) < 1e-12
            push!(cbin, abin[i]);
        end
    end
    
    # Only recurse if depth and mincount allow
    if depth > 0
        if lcount >= mincount
            (abin[1:lcount], lbinends, cbini) = partition_nodes_in_bins_1d(nodes, bbin[1:lcount], [xlim[1], halfx], depth-1, mincount);
            Ncbin = length(cbin);
            next_ind = 1;
            for ci=1:length(cbini)
                for cj=1:Ncbin
                    if cbini[ci] == cbin[cj]
                        break;
                    elseif cj == Ncbin
                        cbini[next_ind] = cbini[ci];
                        next_ind += 1;
                    end
                end
            end
            append!(cbin, cbini[1:(next_ind-1)]);
        else
            abin[1:lcount] = bbin[1:lcount];
            lbinends = [lcount];
        end
        if rcount >= mincount
            (abin[(lcount+1):N], rbinends, cbini) = partition_nodes_in_bins_1d(nodes, bbin[(lcount+1):N], [halfx, xlim[2]], depth-1, mincount);
            Ncbin = length(cbin);
            next_ind = 1;
            for ci=1:length(cbini)
                for cj=1:Ncbin
                    if cbini[ci] == cbin[cj]
                        break;
                    elseif cj == Ncbin
                        cbini[next_ind] = cbini[ci];
                        next_ind += 1;
                    end
                end
            end
            append!(cbin, cbini[1:(next_ind-1)]);
        else
            abin[(lcount+1):N] = bbin[(lcount+1):N];
            rbinends = [rcount];
        end
        binends = append!(lbinends, rbinends .+ lcount);
        
    else
        abin = bbin;
        binends = [lcount, N];
    end
    
    return (abin, binends, cbin);
end

# recursively subdivide bins of nodes
function partition_nodes_in_bins_2d(nodes, abin, xlim, ylim, depth, mincount)
    halfx = (xlim[1] + xlim[2])/2 + 1e-11;
    halfy = (ylim[1] + ylim[2])/2 + 1e-11;
    N = length(abin);
    
    which_bin = zeros(Int, N) # Which bin will the node go in
    bin_counts = zeros(Int, 4) # How many go in each
    cbin = zeros(Int, 0); # These are ones too close to the cutoff that will be added to an extra bin
    for i=1:N
        if nodes[1,abin[i]] <= halfx
            if nodes[2,abin[i]] <= halfy
                bin_counts[1] += 1;
                which_bin[i] = 1;
            else
                bin_counts[2] += 1;
                which_bin[i] = 2;
            end
            
        else
            if nodes[2,abin[i]] <= halfy
                bin_counts[3] += 1;
                which_bin[i] = 3;
            else
                bin_counts[4] += 1;
                which_bin[i] = 4;
            end
        end
        if abs(nodes[1,abin[i]] - halfx) < 1e-12 || abs(nodes[2,abin[i]] - halfy) < 1e-12
            push!(cbin, abin[i]);
        end
    end
    
    # Make bbin arrays and fill them.
    bbins = [zeros(Int, bin_counts[1]), zeros(Int, bin_counts[2]), zeros(Int, bin_counts[3]), zeros(Int, bin_counts[4])];
    next_inds = [1,1,1,1];
    for i=1:N
        ind = which_bin[i];
        bbins[ind][next_inds[ind]] = abin[i];
        next_inds[ind] += 1;
    end
    
    # convenient offsets and limits
    starts = [1, 0,0,0];
    ends   = [bin_counts[1], 0,0,0];
    for i=2:4
        starts[i] = ends[i-1] + 1;
        ends[i] = bin_counts[i] + ends[i-1];
    end
    xlims = [[xlim[1], halfx], [xlim[1], halfx], [halfx, xlim[2]], [halfx, xlim[2]]];
    ylims = [[ylim[1], halfy], [halfy, ylim[2]], [ylim[1], halfy], [halfy, ylim[2]]];
    
    # Only recurse if depth and mincount allow
    if depth > 0
        binends = [];
        for i=1:4
            if bin_counts[i] >= mincount
                (abin[starts[i]:ends[i]], tmpbinends, cbini) = partition_nodes_in_bins_2d(nodes, bbins[i], xlims[i], ylims[i], depth-1, mincount);
                if i>1
                    append!(binends, tmpbinends .+ ends[i-1]);
                else
                    binends = tmpbinends;
                end
                Ncbin = length(cbin);
                next_ind = 1;
                for ci=1:length(cbini)
                    for cj=1:Ncbin
                        if cbini[ci] == cbin[cj]
                            break;
                        elseif cj == Ncbin
                            cbini[next_ind] = cbini[ci];
                            next_ind += 1;
                        end
                    end
                end
                append!(cbin, cbini[1:(next_ind-1)]);
                
            else
                abin[starts[i]:ends[i]] = bbins[i];
                append!(binends, [ends[i]]);
            end
        end
        
    else
        abin = [bbins[1]; bbins[2]; bbins[3]; bbins[4]];
        binends = ends;
    end
    
    return (abin, binends, cbin);
end

# recursively subdivide bins of nodes
function partition_nodes_in_bins_3d(nodes, abin, xlim, ylim, zlim, depth, mincount)
    halfx = (xlim[1] + xlim[2])/2 + 1e-11;
    halfy = (ylim[1] + ylim[2])/2 + 1e-11;
    halfz = (zlim[1] + zlim[2])/2 + 1e-11;
    N = length(abin);
    
    which_bin = zeros(Int, N) # Which bin will the node go in
    bin_counts = zeros(Int, 8) # How many go in each
    cbin = zeros(Int, 0); # These are ones too close to the cutoff that will be added to an extra bin
    for i=1:N
        if nodes[1,abin[i]] <= halfx
            if nodes[2,abin[i]] <= halfy
                if nodes[2,abin[i]] <= halfz
                    bin_counts[1] += 1;
                    which_bin[i] = 1;
                else
                    bin_counts[2] += 1;
                    which_bin[i] = 2;
                end
            else
                if nodes[2,abin[i]] <= halfy
                    bin_counts[3] += 1;
                    which_bin[i] = 3;
                else
                    bin_counts[4] += 1;
                    which_bin[i] = 4;
                end
            end
            
        else
            if nodes[2,abin[i]] <= halfy
                if nodes[2,abin[i]] <= halfz
                    bin_counts[5] += 1;
                    which_bin[i] = 5;
                else
                    bin_counts[6] += 1;
                    which_bin[i] = 6;
                end
            else
                if nodes[2,abin[i]] <= halfy
                    bin_counts[7] += 1;
                    which_bin[i] = 7;
                else
                    bin_counts[8] += 1;
                    which_bin[i] = 8;
                end
            end
        end
        if abs(nodes[1,abin[i]] - halfx) < 1e-12 || abs(nodes[2,abin[i]] - halfy) < 1e-12 || abs(nodes[3,abin[i]] - halfz) < 1e-12
            push!(cbin, abin[i]);
        end
    end
    
    # Make bbin arrays and fill them.
    bbins = [zeros(Int, bin_counts[1]), zeros(Int, bin_counts[2]), zeros(Int, bin_counts[3]), zeros(Int, bin_counts[4]), 
             zeros(Int, bin_counts[5]), zeros(Int, bin_counts[6]), zeros(Int, bin_counts[7]), zeros(Int, bin_counts[8])];
    next_inds = [1,1,1,1,1,1,1,1];
    for i=1:N
        ind = which_bin[i];
        bbins[ind][next_inds[ind]] = abin[i];
        next_inds[ind] += 1;
    end
    
    # convenient offsets and limits
    starts = [1, 0,0,0,0,0,0,0];
    ends   = [bin_counts[1], 0,0,0,0,0,0,0];
    for i=2:8
        starts[i] = ends[i-1] + 1;
        ends[i] = bin_counts[i] + ends[i-1];
    end
    xlims = [[xlim[1], halfx], [xlim[1], halfx], [xlim[1], halfx], [xlim[1], halfx], 
             [halfx, xlim[2]], [halfx, xlim[2]], [halfx, xlim[2]], [halfx, xlim[2]]];
    ylims = [[ylim[1], halfy], [ylim[1], halfy], [halfy, ylim[2]], [halfy, ylim[2]], 
             [ylim[1], halfy], [ylim[1], halfy], [halfy, ylim[2]], [halfy, ylim[2]]];
    zlims = [[zlim[1], halfz], [halfz, zlim[2]], [zlim[1], halfz], [halfz, zlim[2]], 
             [zlim[1], halfz], [halfz, zlim[2]], [zlim[1], halfz], [halfz, zlim[2]]];
    
    # Only recurse if depth and mincount allow
    if depth > 0
        binends = [];
        for i=1:8
            if bin_counts[i] >= mincount
                (abin[starts[i]:ends[i]], tmpbinends, cbini) = partition_nodes_in_bins_3d(nodes, bbins[i], xlims[i], ylims[i], zlims[i], depth-1, mincount);
                if i>1
                    append!(binends, tmpbinends .+ ends[i-1]);
                else
                    binends = tmpbinends;
                end
                Ncbin = length(cbin);
                next_ind = 1;
                for ci=1:length(cbini)
                    for cj=1:Ncbin
                        if cbini[ci] == cbin[cj]
                            break;
                        elseif cj == Ncbin
                            cbini[next_ind] = cbini[ci];
                            next_ind += 1;
                        end
                    end
                end
                append!(cbin, cbini[1:(next_ind-1)]);
            else
                abin[starts[i]:ends[i]] = bbins[i];
                append!(binends, [ends[i]]);
            end
        end
        
    else
        abin = [bbins[1]; bbins[2]; bbins[3]; bbins[4]; bbins[5]; bbins[6]; bbins[7]; bbins[8]];
        binends = ends;
    end
    
    return (abin, binends, cbin);
end

#Extra remove later 

function triangle_element_nodes_(refel, v)
    return  triangle_refel_to_xy_(refel.r[:,1], refel.r[:,2], v);
end

function line_refel_to_x(r, v)
    x = v[1] .+ 0.5 .* (v[2]-v[1]) .* (1 .+ r);
    
    return x;
end

function triangle_refel_to_xy_(r, s, v)
    x = 0.5 * (-(r .+ s) * v[1,1] .+ (1 .+ r) * v[1,2] .+ (1 .+ s) * v[1,3]);
    y = 0.5 * (-(r .+ s) * v[2,1] .+ (1 .+ r) * v[2,2] .+ (1 .+ s) * v[2,3]);
    
    return (x, y);
end

function quad_refel_to_xy(r, s, v)
    vx = v[1,:];
    vy = v[2,:];
    
    rp = 1 .+ r;
    rm = 1 .- r;
    sp = 1 .+ s;
    sm = 1 .- s;
    x = 0.25 .* (rm .* sm .* vx[1] .+ rp .* sm .* vx[2] .+ rp .* sp .* vx[3] .+ rm .* sp .* vx[4]); 
    y = 0.25 .* (rm .* sm .* vy[1] .+ rp .* sm .* vy[2] .+ rp .* sp .* vy[3] .+ rm .* sp .* vy[4]); 
    
    return (x, y);
end

function hex_refel_to_xyz(r, s, t, v)
    vx = v[1,:];
    vy = v[2,:];
    vz = v[3,:];
    
    rp = 1 .+ r;
    rm = 1 .- r;
    sp = 1 .+ s;
    sm = 1 .- s;
    tp = 1 .+ t;
    tm = 1 .- t;
    x = 0.125 .* (rm .* sm .* tm .* vx[1] .+ rp .* sm .* tm .* vx[2] .+ rp .* sp .* tm .* vx[3] .+ rm .* sp .* tm .* vx[4]
                    .+ rm .* sm .* tp .* vx[5] .+ rp .* sm .* tp .* vx[6] .+ rp .* sp .* tp .* vx[7] .+ rm .* sp .* tp .* vx[8]);
    y = 0.125 .* (rm .* sm .* tm .* vy[1] .+ rp .* sm .* tm .* vy[2] .+ rp .* sp .* tm .* vy[3] .+ rm .* sp .* tm .* vy[4]
                    .+ rm .* sm .* tp .* vy[5] .+ rp .* sm .* tp .* vy[6] .+ rp .* sp .* tp .* vy[7] .+ rm .* sp .* tp .* vy[8]);
    z = 0.125 .* (rm .* sm .* tm .* vz[1] .+ rp .* sm .* tm .* vz[2] .+ rp .* sp .* tm .* vz[3] .+ rm .* sp .* tm .* vz[4]
                    .+ rm .* sm .* tp .* vz[5] .+ rp .* sm .* tp .* vz[6] .+ rp .* sp .* tp .* vz[7] .+ rm .* sp .* tp .* vz[8]);
    
    return (x, y, z);
end

function tetrahedron_refel_to_xyz(r, s, t, v)
    # Check the orientation of the vertices.
    # If they don't match refel, the jacobian will be bad.
    # (p2-p1) X (p3-p1) points toward p4
    # If not, swap p3 and p4
    p1 = v[:,1];
    p2 = v[:,2];
    p3 = v[:,3];
    p4 = v[:,4];
    e1 = p2 - p1;
    e2 = p3 - p1;
    e1xe2 = [e1[2]*e2[3] - e1[3]*e2[2], e1[3]*e2[1] - e1[1]*e2[3], e1[1]*e2[2] - e1[2]*e2[1]];
    # If dist(p1+e1xe2, p4) < dist(p1-e1xe2, p4), it is pointing toward p4
    d1 = (p4[1] - (p1[1]+e1xe2[1]))^2 + (p4[2] - (p1[2]+e1xe2[2]))^2 + (p4[3] - (p1[3]+e1xe2[3]))^2;
    d2 = (p4[1] - (p1[1]-e1xe2[1]))^2 + (p4[2] - (p1[2]-e1xe2[2]))^2 + (p4[3] - (p1[3]-e1xe2[3]))^2;
    if d1 > d2
        #3 and 4 need to be swapped
        p3 = v[:,4];
        p4 = v[:,3];
    end
    
    A = [p2-p1   p3-p1   p4-p1];
    
    np = length(r);
    mv = zeros(finch_state.config.float_type, 3,np);
    for i=1:np
        tmp = [(r[i]+1)/2, (s[i]+1)/2, (t[i]+1)/2];
        mv[:,i] = A*tmp + p1;
    end
    
    x = mv[1,:]
    y = mv[2,:]
    z = mv[3,:]
    
    return (x, y, z);
end

function tetrahedron_refel_to_xyz!(rst::Matrix, v::Matrix, x::Vector, y::Vector, z::Vector)
    # Check the orientation of the vertices.
    # If they don't match refel, the jacobian will be bad.
    # (p2-p1) X (p3-p1) points toward p4
    # If not, swap p3 and p4
    p1 = v[1:3,1];
    p2 = v[1:3,2];
    p3 = v[1:3,3];
    p4 = v[1:3,4];
    e1 = p2 - p1;
    e2 = p3 - p1;
    e1xe2 = [e1[2]*e2[3] - e1[3]*e2[2], e1[3]*e2[1] - e1[1]*e2[3], e1[1]*e2[2] - e1[2]*e2[1]];
    # If dist(p1+e1xe2, p4) < dist(p1-e1xe2, p4), it is pointing toward p4
    d1 = (p4[1] - (p1[1]+e1xe2[1]))^2 + (p4[2] - (p1[2]+e1xe2[2]))^2 + (p4[3] - (p1[3]+e1xe2[3]))^2;
    d2 = (p4[1] - (p1[1]-e1xe2[1]))^2 + (p4[2] - (p1[2]-e1xe2[2]))^2 + (p4[3] - (p1[3]-e1xe2[3]))^2;
    if d1 > d2
        #3 and 4 need to be swapped
        p3 = v[1:3,4];
        p4 = v[1:3,3];
    end
    
    np = size(rst,1);
    # mv = zeros(finch_state.config.float_type, 3, np);
    # A = [p2-p1   p3-p1   p4-p1];
    for i=1:np
        # Leave this to remember what this does
        # tmp = [(rst[i,1]+1)/2, (rst[i,2]+1)/2, (rst[i,3]+1)/2];
        # mv[:,i] = A*tmp + p1;
        # x[i] = mv[1,i];
        # y[i] = mv[2,i];
        # z[i] = mv[3,i];
        
        x[i] = p1[1] + ( (p2[1]-p1[1])*(rst[i,1]+1) + 
                         (p3[1]-p1[1])*(rst[i,2]+1) + 
                         (p4[1]-p1[1])*(rst[i,3]+1) )*0.5;
        y[i] = p1[2] + ( (p2[2]-p1[2])*(rst[i,1]+1) + 
                         (p3[2]-p1[2])*(rst[i,2]+1) + 
                         (p4[2]-p1[2])*(rst[i,3]+1) )*0.5;
        z[i] = p1[3] + ( (p2[3]-p1[3])*(rst[i,1]+1) + 
                         (p3[3]-p1[3])*(rst[i,2]+1) + 
                         (p4[3]-p1[3])*(rst[i,3]+1) )*0.5;
    end
end

# Returns true if the nodes are within tol of each other.
# If scale is provided, it is a relative tolerance.
function is_same_node(x1, x2, tol, scale=1)
    return (sum(abs.(x1 - x2)) / scale) < tol
end
function is_same_node_1d(tol, scale, x1, x2)
    return (abs(x1-x2) / scale) < tol
end
function is_same_node_2d(tol, scale, x1, x2, y1, y2)
    return ((abs(x1-x2) + abs(y1-y2)) / scale) < tol
end
function is_same_node_3d(tol, scale, x1, x2, y1, y2, z1, z2)
    return ((abs(x1-x2) + abs(y1-y2) + abs(z1-z2)) / scale) < tol
end

# Returns true if the centroid of each line is close enough.
function is_same_face_center(l1, l2, tol, scale=1)
    n1 = size(l1,2);
    n2 = size(l2,2);
    center1 = zeros(finch_state.config.float_type, size(l1,1));
    center2 = zeros(finch_state.config.float_type, size(l2,1));
    for i=1:n1
        center1 = center1 .+ l1[:,i];
    end
    center1 = center1 ./ n1;
    
    for i=1:n2
        center2 = center2 .+ l2[:,i];
    end
    center2 = center2 ./ n2;
    
    return is_same_node(center1, center2, tol, scale);
end

# Returns the distance between face centers.
function face_center_distance(l1, l2)
    n1 = size(l1,2);
    n2 = size(l2,2);
    center1 = zeros(finch_state.config.float_type, size(l1,1));
    center2 = zeros(finch_state.config.float_type, size(l2,1));
    for i=1:n1
        center1 = center1 .+ l1[:,i];
    end
    center1 = center1 ./ n1;
    
    for i=1:n2
        center2 = center2 .+ l2[:,i];
    end
    center2 = center2 ./ n2;
    
    return sum(abs.(center1 - center2));
end

# Returns true if the two node lists have at least two of the same nodes.
function is_same_line(l1, l2, tol, scale=1)
    found = 0;
    n1 = size(l1,2);
    n2 = size(l2,2);
    dim = size(l1,1);
    
    if dim == 2
        for i=1:n1
            for j=1:n2
                if is_same_node_2d(tol, scale, l1[1,i], l2[1,j], l1[2,i], l2[2,j]) # is_same_node(l1[:,i], l2[:,j], tol, scale)
                    found += 1;
                end
                if found >= 2
                    return true;
                end
            end
        end
    else
        for i=1:n1
            for j=1:n2
                if is_same_node_3d(tol, scale, l1[1,i], l2[1,j], l1[2,i], l2[2,j], l1[3,i], l2[3,j]) # is_same_node(l1[:,i], l2[:,j], tol, scale)
                    found += 1;
                end
                if found >= 2
                    return true;
                end
            end
        end
    end
    
    return false;
end

# Returns true if the two node lists have at least three of the same nodes.
function is_same_plane(p1, p2, tol, scale=1)
    found = 0;
    n1 = size(p1,2);
    n2 = size(p2,2);
    for i=1:n1
        for j=1:n2
            if is_same_node_3d(tol, scale, p1[1,i], p2[1,j], p1[2,i], p2[2,j], p1[3,i], p2[3,j])# is_same_node(p1[:,i], p2[:,j], tol, scale)
                found += 1;
            end
            if found >= 3
                return true;
            end
        end
    end
    
    return false;
end

function which_refel_face(f2glb, nodes, refel, elem, tol, scale)
    # Check f2glb against the face2local in refel
    fnodes = nodes[:,f2glb];
    for fi=1:refel.Nfaces
        refnodes = nodes[:, elem[refel.face2local[fi]]];
        
        if is_same_face(fnodes, refnodes, refel.dim, tol, scale)
            return fi;
        end
    end
    
    printerr("Couldn't match face when building grid (see which_refel_face() in grid.jl)");
end

function is_same_face(f1, f2, dim, tol, scale=1)
    if dim == 1 # one point
        return is_same_node_1d(tol, scale, f1[1], f2[1]);
    elseif dim == 2 # same line(two same points)
        return is_same_line(f1, f2, tol, scale);
    elseif dim == 3 # same plane(three same points)
        return is_same_plane(f1, f2, tol, scale);
    end
end

# Adds a boundary ID to some region. Find boundary points satifying on_bdry and moves them to a new set for this bid.
function add_boundary_ID_to_grid(bid, on_bdry, grid)
    int_type = finch_state.config.index_type;
    float_type = finch_state.config.float_type;
    dim = finch_state.config.dimension;
    # Find if this bid exists. If so, just add points to it, removing from others.
    ind = indexin([bid], grid.bids)[1];
    nbids = length(grid.bids);
    if ind === nothing
        # This is a new bid, add to bids, bdry, bdryface, bdrynorm
        ind = nbids + 1;
        nbids += 1;
        push!(grid.bids, bid);
        push!(grid.bdry, zeros(int_type, 0));
        push!(grid.bdryface, zeros(int_type, 0));
        push!(grid.bdrynorm, zeros(float_type, dim, 0));
    end
    
    # If it is an empty Grid, do nothing else
    if size(grid.allnodes,2) == 0
        log_entry("Added boundary ID: "*string(bid));
        return;
    end
    
    # Search all other bids for nodes and faces on this segment. Remove them there and add them here.
    # First find indices and count them. Then move.
    move_nodes = Array{Array{Int,1},1}(undef,nbids);
    node_count = zeros(int_type, nbids);
    move_faces = Array{Array{Int,1},1}(undef,nbids);
    face_count = zeros(int_type, nbids);
    for i=1:nbids
        bi = grid.bids[i];
        move_nodes[i] = [];
        move_faces[i] = [];
        if bi != bid
            # First the nodes
            for j=1:length(grid.bdry[i])
                nj = grid.bdry[i][j];
                if dim == 1
                    if Base.invokelatest(on_bdry, grid.allnodes[1, nj])
                        push!(move_nodes[i], nj);
                        node_count[i] += 1;
                    end
                elseif dim == 2
                    if Base.invokelatest(on_bdry, grid.allnodes[1, nj], grid.allnodes[2, nj])
                        push!(move_nodes[i], nj);
                        node_count[i] += 1;
                    end
                elseif dim == 3
                    if Base.invokelatest(on_bdry, grid.allnodes[1, nj], grid.allnodes[2, nj], grid.allnodes[3, nj])
                        push!(move_nodes[i], nj);
                        node_count[i] += 1;
                    end
                end
            end
            # Then the faces
            # Base this decision on the center of the face
            for j=1:length(grid.bdryface[i])
                fj = grid.bdryface[i][j];
                nfp = size(grid.face2glb,1)
                isbdryface = true
                # find the center
                fcenter = zeros(float_type, size(grid.allnodes,1));
                for ni=1:nfp
                    fcenter = fcenter + grid.allnodes[:,grid.face2glb[ni,1,fj]];
                end
                fcenter = fcenter./nfp;
                
                if dim == 1
                    isbdryface = Base.invokelatest(on_bdry, fcenter[1]);
                elseif dim == 2
                    isbdryface = Base.invokelatest(on_bdry, fcenter[1], fcenter[2]);
                elseif dim == 3
                    isbdryface = Base.invokelatest(on_bdry, fcenter[1], fcenter[2], fcenter[3]);
                end
                
                if isbdryface
                    push!(move_faces[i], fj);
                    face_count[i] += 1;
                end
            end
        end
    end # find indices
    
    # Move things from other bids to this one
    for i=1:nbids
        if i != ind
            # Add to this bid
            append!(grid.bdry[ind], move_nodes[i]);
            append!(grid.bdryface[ind], move_faces[i]);
            grid.bdrynorm[ind] = hcat(grid.bdrynorm[ind], grid.bdrynorm[i][:,indexin(move_nodes[i], grid.bdry[i])]);
            
            # Make sure all of the norms correspond to the face on this bdry
            startnodeind = length(grid.bdry[ind]) - length(move_nodes[i]);
            for ni=1:length(move_nodes[i])
                for fi=1:length(move_faces[i])
                    # Does this node lie on this face?
                    facenodeindex = indexin(move_nodes[i][ni], grid.face2glb[:,1,move_faces[i][fi]]);
                    if length(facenodeindex) > 0
                        grid.bdrynorm[ind][:,startnodeind + ni] = grid.facenormals[:,move_faces[i][fi]];
                        break;
                    end
                end
            end
            
            # Remove things from other bids
            # Remove bdrynorm and bdryfacenorm
            numremove = length(move_nodes[i]);
            if numremove > 0
                newbdrynorm = zeros(float_type, dim, size(grid.bdrynorm[i],2) - numremove);
                nextind = 1;
                for j=1:length(grid.bdry[i])
                    keepit = true;
                    for k=1:numremove
                        if grid.bdry[i][j] == move_nodes[i][k]
                            keepit = false;
                            break;
                        end
                    end
                    if keepit
                        newbdrynorm[:,nextind] = grid.bdrynorm[i][:,j];
                        nextind += 1;
                    end
                end
                grid.bdrynorm[i] = newbdrynorm;
            end
            
            # Remove nodes
            deleteat!(grid.bdry[i], indexin(move_nodes[i], grid.bdry[i]));
            
            # Remove bdryface
            deleteat!(grid.bdryface[i], indexin(move_faces[i], grid.bdryface[i]));
            
            # Change the bid of moved faces in facebid
            for fi=1:length(move_faces[i])
                grid.facebid[move_faces[i][fi]] = bid;
            end
        end
    end
    
    # update nodebid
    for i=1:length(grid.bdry[ind])
        grid.nodebid[grid.bdry[ind][i]] = bid;
    end
    
    log_entry("Added boundary ID: "*string(bid)*" including "*string(node_count)*" nodes, "*string(face_count)*" faces.");
    
end

# Reorder the nodes and update the maps accordingly
# new_map maps old to new index 
# newnodes[:,new_map[i]] = oldnodes[:,i]
# newloc2glb[ni,ei] = new_map[oldloc2glb[ni,ei]]
function reorder_grid_nodes!(grid, new_map)
    # Need to change:
    # - allnodes
    # - bdry
    # - loc2glb
    # - glbvertex
    # - face2glb
    # These are temporary copies
    copy_nodes = copy(grid.allnodes);
    copy_nodebid = copy(grid.nodebid);
    copy_bdry = copy(grid.bdry);
    copy_loc2glb = copy(grid.loc2glb);
    copy_glbvertex = copy(grid.glbvertex);
    copy_face2glb = copy(grid.face2glb);
    
    N = size(grid.allnodes,2);
    nel = size(grid.loc2glb,2);
    nface = size(grid.face2glb,3);
    
    # allnodes
    for ni=1:N
        grid.allnodes[:,new_map[ni]] = copy_nodes[:,ni];
        grid.nodebid[new_map[ni]] = copy_nodebid[ni];
    end
    
    # bdry
    for bi=1:length(grid.bids)
        for ni=1:length(grid.bdry[bi])
            grid.bdry[bi][ni] = new_map[copy_bdry[bi][ni]];
        end
    end
    
    # loc2glb, glbvertex
    for ei=1:nel
        for ni=1:size(grid.loc2glb,1)
            grid.loc2glb[ni,ei] = new_map[copy_loc2glb[ni,ei]];
        end
        for vi=1:size(grid.glbvertex,1)
            grid.glbvertex[vi,ei] = new_map[copy_glbvertex[vi,ei]];
        end
    end
    
    # face2glb
    for fi=1:nface
        for gi=1:size(grid.face2glb,2)
            for ni=1:size(grid.face2glb,1)
                grid.face2glb[ni,gi,fi] = new_map[copy_face2glb[ni,gi,fi]];
            end
        end
    end
end
