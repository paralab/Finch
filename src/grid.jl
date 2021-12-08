#=
# Contains info about all nodes on the domain
# Unlike MeshData struct, this accounts for interior nodes and corresponds to nodal DOFs.
=#

struct Grid
    allnodes::Array{Float64}        # All node coordinates size = (dim, nnodes)
    # boundaries
    bdry::Array{Array{Int,1},1}     # Indices of boundary nodes for each BID (bdry[bid][nodes])*note:array of arrays
    bdryface::Array{Array{Int,1},1} # Indices of faces touching each BID (bdryface[bid][faces])*note:array of arrays
    bdrynorm::Array{Array{Float64,2},1} # Normal vector for boundary nodes for each BID (bdrynorm[bid][dim, nodes])*note:array of arrays
    bids::Array{Int,1}              # BID corresponding to rows of bdrynodes
    # elements
    loc2glb::Array{Int,2}           # local to global map for each element's nodes (size is (Np, nel))
    glbvertex::Array{Int,2}         # global indices of each elements' vertices (size is (Nvertex, nel))
    # faces (For CG, G=1. For DG, G=2)
    face2glb::Array{Int,3}          # local to global map for faces (size is (Nfp, G, Nfaces))
    element2face::Array{Int,2}      # face indices for each element (size is (Nfaces, nel))
    face2element::Array{Int,2}      # elements on both sides of a face, 0=boundary (size is (2, Nfaces))
    facenormals::Array{Float64,2}   # normal vector for each face
    faceRefelInd::Array{Int,2}      # Index for face within the refel for each side
    
    facebid::Array{Int,1}           # BID of each face (0=interior face)
    
    # When partitioning the grid, this stores the ghost info.
    is_subgrid::Bool                # Is this a partition of a greater grid?
    nel_owned::Int                  # Number of elements owned by this partition
    nel_ghost::Int                  # Number of ghost elements
    element_owner::Vector{Int}      # The rank of each ghost element's owner or -1 if locally owned
    grid2mesh::Vector{Int}          # Map from partition elements to global mesh element index
    
    num_neighbor_partitions::Int        # number of partitions that share ghosts with this.
    neighboring_partitions::Vector{Int} # IDs of neighboring partitions
    ghost_counts::Vector{Int}           # How many ghosts for each neighbor
    ghost_index::Vector{Array{Int}}   # Lists of ghost elements to send/recv for each neighbor
    
    # constructors
    Grid(allnodes, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, element2face, 
         face2element, facenormals, faceRefelInd, facebid) = 
     new(allnodes, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, element2face, 
         face2element, facenormals, faceRefelInd, facebid, 
         false, size(loc2glb,2), 0, zeros(Int,0), zeros(Int,0), 0, zeros(Int,0), zeros(Int,0), [zeros(Int,2,0)]); # up to facebid only
     
    Grid(allnodes, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, element2face, 
         face2element, facenormals, faceRefelInd, facebid, 
         ispartitioned, nel_owned, nel_ghost, element_owners, grid2mesh, num_neighbors, neighbor_ids, ghost_counts, ghost_ind) = 
     new(allnodes, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, element2face, 
         face2element, facenormals, faceRefelInd, facebid, 
         ispartitioned, nel_owned, nel_ghost, element_owners, grid2mesh, num_neighbors, neighbor_ids, ghost_counts, ghost_ind); # subgrid parts included
end

etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5]; # number of vertices for each type
etypetodim= [1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3]; # dimension of each type
etypetonf = [2, 3, 4, 4, 6, 5, 5, 2, 3, 4, 4, 6, 5, 5, 1, 4, 6, 5, 5]; # number of faces for element types
etypetoftype=[0,1, 1, 2, 3, 3, 3, 0, 1, 1, 2, 3, 3, 3, 0, 0, 0, 0, 0]; # type of faces for this element type

# Build a grid from a mesh
# This is for full grids. For partitioned grids see partitioned_grid_from_mesh()
function grid_from_mesh(mesh)
    log_entry("Building full grid from mesh data", 2);
    t_grid_from_mesh = Base.Libc.time();
    dim = config.dimension;
    ord = config.basis_order_min;
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
    refel = build_refel(dim, ord, nfaces, config.elemental_nodes);
    
    if config.solver_type == DG # || config.solver_type == FV ????
        Gness = 2;
    else
        Gness = 1;
    end
    
    tol = 1e-12; # tolerance for is_same_node
    
    Np = refel.Np;                      # number of nodes per element
    bdry = [];                          # index(in x) of boundary nodes for each BID
    bdryfc = [];                        # index of faces touching each BID
    bdrynorm = [];                      # normal at boundary nodes
    bids = collectBIDs(mesh);           # BID list
    nbids = length(bids);
    for i=1:nbids
        push!(bdry, zeros(Int, 0));
        push!(bdryfc, zeros(Int,0));
        push!(bdrynorm, zeros(config.dimension,0));
    end
    loc2glb = zeros(Int, Np, nel)       # local to global index map for each element's nodes
    glbvertex = zeros(Int, nvtx, nel);     # local to global for vertices
    f2glb = zeros(Int, refel.Nfp[1], Gness, totalfaces);  # face node local to global
    element2face = zeros(Int, nfaces, nel);  # element to face map
    face2element = zeros(Int, 2, size(mesh.face2element,2));  # face to element map
    facenormals = zeros(dim, totalfaces); # normal vectors for every face
    faceRefelInd = zeros(Int, 2, totalfaces); # Index in refel for this face for elements on both sides
    facebid = zeros(Int, totalfaces); # BID of each face
    
    tmpallnodes = zeros(dim, mesh.nel*refel.Np);
    t_nodes1 = Base.Libc.time();
    for ei=1:nel
        # Find this element's nodes
        n_vert = etypetonv[mesh.etypes[ei]];
        e_vert = mesh.nodes[1:dim, mesh.elements[1:n_vert, ei]];
        
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
        tmpallnodes[1, ((ei-1)*Np+1):(ei*Np)] = e_x;
        if dim > 1
            tmpallnodes[2, ((ei-1)*Np+1):(ei*Np)] = e_y;
            if dim > 2
                tmpallnodes[3, ((ei-1)*Np+1):(ei*Np)] = e_z;
            end
        end
        
        # temporary mapping
        loc2glb[:,ei] = ((ei-1)*Np+1):(ei*Np);
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
    
    # vertices, faces and boundary
    t_faces1 = Base.Libc.time();
    for ei=1:nel
        n_vert = etypetonv[mesh.etypes[ei]];
        mfids = mesh.element2face[:,ei];
        normals = mesh.normals[:,mfids];
        el_center = zeros(dim);
        
        # vertices and center
        for ni=1:Np
            el_center += allnodes[:,loc2glb[ni,ei]];
            
            for vi=1:n_vert
                if is_same_node(mesh.nodes[:, mesh.elements[vi,ei]], allnodes[:,loc2glb[ni,ei]], tol)
                    glbvertex[vi, ei] = loc2glb[ni,ei];
                end
            end
        end
        el_center ./= Np;
        
        # f2glb has duplicates. Compare to mesh faces and keep same ordering as mesh.
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
            tmpf2glb = loc2glb[refel.face2local[gfi], ei];
            
            for mfi=1:nfaces
                thisfaceind = meshfaces[mfi];
                if test_same_face(mesh.nodes[:,mesh.face2vertex[:,thisfaceind]], allnodes[:, tmpf2glb], tol)
                    # This mesh face corresponds to this tmpf2glb face
                    # Put the tmpf2glb map into f2glb at the mesh index(thisfaceind).
                    # Move the f2glb[:,1,ind] to f2glb[:,2,ind] first if DG (Gness==2)
                    # Set element2face according to gfi(not mfi)
                    # Move face2element[1] to face2element[2] and put this one in [1]
                    if (Gness == 2) f2glb[:, 2, thisfaceind] = f2glb[:, 1, thisfaceind]; end
                    f2glb[:, 1, thisfaceind] = tmpf2glb;
                    element2face[gfi, ei] = thisfaceind;
                    face2element[2, thisfaceind] = face2element[1, thisfaceind];
                    face2element[1, thisfaceind] = ei;
                    
                    # Find the normal for every face. The normal points from e1 to e2 or outward for boundary.
                    # Note that the normal stored in mesh_data could be pointing either way.
                    thisnormal = normals[:, mfi];
                    f_center = zeros(dim);
                    for ni=1:length(tmpf2glb)
                        f_center += allnodes[:, tmpf2glb[ni]];
                    end
                    f_center ./= length(tmpf2glb)
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
                    nfacenodes = length(tmpf2glb);
                    if !(gbid === nothing) # This is a boundary face
                        append!(bdry[gbid], tmpf2glb);
                        push!(bdryfc[gbid], thisfaceind);
                        facebid[thisfaceind] = gbid;
                        thisnormal = normals[:, mfi];
                        normchunk = zeros(config.dimension, nfacenodes);
                        for ni=1:nfacenodes
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
    newbdry = similar(bdry);
    newbdrynorm = similar(bdrynorm);
    for i=1:length(bdry)
        newbdry[i] = zeros(length(bdry[i]));
        newbdrynorm[i] = zeros(size(bdrynorm[i]));
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
        faceRefelInd[1,fi] = which_refel_face(f2glb[:,1,fi], allnodes, refel, loc2glb[:,face2element[1,fi]]);
        if face2element[2,fi] > 0
            faceRefelInd[2,fi] = which_refel_face(f2glb[:,Gness,fi], allnodes, refel, loc2glb[:,face2element[2,fi]]);
        end
    end
    
    t_grid_from_mesh = Base.Libc.time() - t_grid_from_mesh;
    log_entry("Total grid building time: "*string(t_grid_from_mesh), 2);
    
    return (refel, Grid(allnodes, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, element2face, face2element, facenormals, faceRefelInd, facebid));
end

# This takes a full mesh and partition information in a format supplied by METIS.
# It constructs a grid that only holds this partition and ghosted neighbor elements.
function partitioned_grid_from_mesh(mesh, epart)
    log_entry("Building partitioned grid from mesh data", 2);
    t_grid_from_mesh = Base.Libc.time();
    dim = config.dimension;
    ord = config.basis_order_min;
    nfaces = etypetonf[mesh.etypes[1]]; # faces per elements
    
    # Count the owned faces and owned/ghost elements
    totalfaces = 0;
    nel_owned = 0;
    nel_ghost = 0;
    element_status = fill(-1, mesh.nel); # 1 for ghosts, 0 for owned, -1 for other
    face_status = fill(-1, size(mesh.normals,2)); # 1 for face to ghosts, 0 for owned, -1 for other
    mesh2grid_face = fill(-1, size(mesh.normals,2)); # -1 if this face is not in the partition, otherwise the index of the grid face
    num_neighbors = 0; # number of partitions neighboring this
    neighbor_ids = zeros(Int,0); # ID of neighbors
    ghost_counts = zeros(Int,0); # number of ghosts per neighbor
    # ghost_inds = [zeros(Int,2,0)]; # local index of ghost elements
    
    for ei=1:mesh.nel
        if epart[ei] == config.partition_index
            nel_owned += 1;
            element_status[ei] = 0;
        end
    end
    for fi=1:size(mesh.normals,2)
        e1 = mesh.face2element[1,fi];
        e2 = mesh.face2element[2,fi];
        if element_status[e1] == 0 || (e2 > 0 && element_status[e2] == 0)
            totalfaces += 1;
            if !(element_status[e1] == 0)
                if element_status[e1] == -1
                    nel_ghost += 1;
                end
                element_status[e1] = 1;
                face_status[fi] = 1;
            elseif e2 > 0 && !(element_status[e2] == 0)
                if element_status[e2] == -1
                    nel_ghost += 1;
                end
                element_status[e2] = 1;
                face_status[fi] = 1;
            else
                face_status[fi] = 0;
            end
            mesh2grid_face[fi] = totalfaces;
        end
    end
    nel = nel_owned + nel_ghost;
    
    if dim == 1
        facenvtx = 1
    else
        facenvtx = etypetonv[etypetoftype[mesh.etypes[1]]]; # Assumes one element type
    end
    nvtx = etypetonv[mesh.etypes[1]]; # Assumes one element type
    
    log_entry("Building reference element: "*string(dim)*"D, order="*string(ord)*", nfaces="*string(nfaces), 3);
    refel = build_refel(dim, ord, nfaces, config.elemental_nodes);
    
    if config.solver_type == DG # || config.solver_type == FV ????
        Gness = 2;
    else
        Gness = 1;
    end
    
    tol = 1e-12; # tolerance for is_same_node
    
    Np = refel.Np;                      # number of nodes per element
    bdry = [];                          # index(in x) of boundary nodes for each BID
    bdryfc = [];                        # index of faces touching each BID
    bdrynorm = [];                      # normal at boundary nodes
    bids = collectBIDs(mesh);           # BID list
    nbids = length(bids);
    for i=1:nbids
        push!(bdry, zeros(Int, 0));
        push!(bdryfc, zeros(Int,0));
        push!(bdrynorm, zeros(config.dimension,0));
    end
    loc2glb = zeros(Int, Np, nel)       # local to global index map for each element's nodes
    glbvertex = zeros(Int, nvtx, nel);     # local to global for vertices
    f2glb = zeros(Int, refel.Nfp[1], Gness, totalfaces);  # face node local to global
    element2face = zeros(Int, nfaces, nel);  # element to face map
    face2element = zeros(Int, 2, totalfaces);  # face to element map
    facenormals = zeros(dim, totalfaces); # normal vectors for every face
    faceRefelInd = zeros(Int, 2, totalfaces); # Index in refel for this face for elements on both sides
    facebid = zeros(Int, totalfaces); # BID of each face
    
    element_owners = fill(-1, nel) # partition number of each ghost, or -1 for owned elements
    grid2mesh = zeros(Int, nel); # maps partition elements to global mesh elements
    
    # compute node coordinates
    tmpallnodes = zeros(dim, nel*refel.Np);
    next_e_index = 1;
    next_g_index = nel_owned+1; #index for next ghost
    t_nodes1 = Base.Libc.time();
    for ei=1:mesh.nel
        if element_status[ei] >= 0 # owned or ghost element
            # Find this element's nodes
            n_vert = etypetonv[mesh.etypes[ei]];
            e_vert = mesh.nodes[1:dim, mesh.elements[1:n_vert, ei]];
            
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
            if element_status[ei] == 0 # owned
                next_index = next_e_index;
                next_e_index += 1;
            else # ghost
                next_index = next_g_index;
                next_g_index += 1;
            end
            
            tmpallnodes[1, ((next_index-1)*Np+1):(next_index*Np)] = e_x;
            if dim > 1
                tmpallnodes[2, ((next_index-1)*Np+1):(next_index*Np)] = e_y;
                if dim > 2
                    tmpallnodes[3, ((next_index-1)*Np+1):(next_index*Np)] = e_z;
                end
            end
            
            # temporary mapping
            loc2glb[:,next_index] = ((next_index-1)*Np+1):(next_index*Np);
            
            grid2mesh[next_index] = ei;
            
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
    
    # vertices, faces and boundary
    # first make a map from mesh faces to grid faces
    next_e_index = 1;
    next_g_index = nel_owned+1; #index for next ghost
    t_faces1 = Base.Libc.time();
    for ei=1:mesh.nel
        if element_status[ei] >= 0 # owned or ghost element
            n_vert = etypetonv[mesh.etypes[ei]];
            mfids = mesh.element2face[:,ei];
            normals = mesh.normals[:,mfids];
            el_center = zeros(dim);
            
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
                    if is_same_node(mesh.nodes[:, mesh.elements[vi,ei]], allnodes[:,loc2glb[ni,next_index]], tol)
                        glbvertex[vi, next_index] = loc2glb[ni,next_index];
                    end
                end
            end
            el_center ./= Np;
            
            # f2glb has duplicates. Compare to mesh faces and keep same ordering as mesh.
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
                
                for mfi=1:nfaces
                    thisfaceind = meshfaces[mfi];
                    gridfaceind = mesh2grid_face[thisfaceind]; # the face index in the grid
                    if gridfaceind > 0 && test_same_face(mesh.nodes[:,mesh.face2vertex[:,thisfaceind]], allnodes[:, tmpf2glb], tol)
                        # This mesh face corresponds to this tmpf2glb face
                        # Put the tmpf2glb map into f2glb at the mesh index(thisfaceind).
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
                        f_center = zeros(dim);
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
                        mbid = mesh.bdryID[thisfaceind];
                        gbid = indexin([mbid], bids)[1];
                        nfacenodes = length(tmpf2glb);
                        if !(gbid === nothing) # This is a boundary face
                            append!(bdry[gbid], tmpf2glb);
                            push!(bdryfc[gbid], gridfaceind);
                            facebid[gridfaceind] = gbid;
                            thisnormal = normals[:, mfi];
                            normchunk = zeros(config.dimension, nfacenodes);
                            for ni=1:nfacenodes
                                normchunk[:,ni] = thisnormal;
                            end
                            bdrynorm[gbid] = hcat(bdrynorm[gbid], normchunk);
                        end
                    end
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
        newbdry[i] = zeros(length(bdry[i]));
        newbdrynorm[i] = zeros(size(bdrynorm[i]));
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
        faceRefelInd[1,fi] = which_refel_face(f2glb[:,1,fi], allnodes, refel, loc2glb[:,face2element[1,fi]]);
        if face2element[2,fi] > 0
            faceRefelInd[2,fi] = which_refel_face(f2glb[:,Gness,fi], allnodes, refel, loc2glb[:,face2element[2,fi]]);
        end
    end
    
    # Form ghost pairs for send/recv
    # First count how many pairs are needed for each neighbor
    ghost_counts = zeros(Int, num_neighbors)
    for fi=1:totalfaces
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
        ghost_inds[ni] = zeros(Int, 2, ghost_counts[ni]);
    end
    # Need to loop over mesh faces to be sure they are built in the same order on each partition
    for i=1:length(mesh2grid_face)
        fi = mesh2grid_face[i];
        if fi > 0 # exists in this partition (at least one of e1,e2 are owned)
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
    
    t_grid_from_mesh = Base.Libc.time() - t_grid_from_mesh;
    log_entry("Total grid building time: "*string(t_grid_from_mesh), 2);
    
    return (refel, Grid(allnodes, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, element2face, 
            face2element, facenormals, faceRefelInd, facebid, 
            true, nel_owned, nel_ghost, element_owners, grid2mesh, num_neighbors, neighbor_ids, ghost_counts, ghost_inds));
end

### Utilities ###

function collectBIDs(mesh)
    bids = [];
    for i=1:length(mesh.bdryID)
        if mesh.bdryID[i] > 0
            already = false;
            for j=1:length(bids)
                if mesh.bdryID[i] == bids[j]
                    already = true;
                end
            end
            if !already
                push!(bids, mesh.bdryID[i]);
            end
        end
    end
    return bids;
end

# Removes duplicate nodes and updates local to global maps
function remove_duplicate_nodes(nodes, loc2glb; tol=1e-12, depth=5, mincount=50, other2glb=[])
    # defaults
    # depth = 5; # 32768 for 3D
    # mincount = 50; # don't subdivide if less than this
    
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
        xyz=nodes[:,i];
        xlim[1] = min(xlim[1], xyz[1]);
        xlim[2] = max(xlim[2], xyz[1]);
        if dim>1
            ylim[1] = min(ylim[1], xyz[2]);
            ylim[2] = max(ylim[2], xyz[2]);
            if dim > 2
                zlim[1] = min(zlim[1], xyz[3]);
                zlim[2] = max(zlim[2], xyz[3]);
            end
        end
    end
    
    if dim == 1
        (abins, bin_ends) = partition_nodes_in_bins_1d(nodes, abins, xlim, depth, mincount);
    elseif dim == 2
        (abins, bin_ends) = partition_nodes_in_bins_2d(nodes, abins, xlim, ylim, depth, mincount);
    else # dim==3
        (abins, bin_ends) = partition_nodes_in_bins_3d(nodes, abins, xlim, ylim, zlim, depth, mincount);
    end
    
    remove_count = 0;
    startind = 1;
    # Loop over nodes in bins and put numbers in replace_with where duplicates will be removed.
    for bini = 1:length(bin_ends)
        for ni=startind:bin_ends[bini]
            for nj=startind:(ni-1)
                # If node[nj] == node[ni], keep nj, remove ni
                if is_same_node(nodes[:,abins[ni]], nodes[:,abins[nj]], tol)
                    remove_count += 1;
                    replace_with[abins[ni]] = abins[nj];
                    break;
                end
            end
        end
        startind = bin_ends[bini] + 1;
    end
    
    # println(string(length(bin_ends))*" bins: " )
    # println(bin_ends);
    # println("Nnodes before: "*string(tmpNnodes)*" after: "*string(tmpNnodes-remove_count));
    
    # Loop over nodes and place in new allnodes array while adjusting replace_with
    Nnodes = tmpNnodes-remove_count;
    next_ind = 1;
    newnodes = zeros(size(nodes,1), Nnodes);
    new_homes = zeros(Int, tmpNnodes);
    for i=1:tmpNnodes
        if replace_with[i] == 0
            # A node to keep
            newnodes[:,next_ind] = nodes[:,i];
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
    for ei=1:size(loc2glb,2) # nel
        for ni=1:size(loc2glb,1) # Np
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
    halfx = (xlim[1] + xlim[2])/2;
    N = length(abin);
    
    bbin = similar(abin); # Temporary storage
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
    end
    
    # Only recurse if depth and mincount allow
    if depth > 0
        if lcount >= mincount
            (abin[1:lcount], lbinends) = partition_nodes_in_bins_1d(nodes, bbin[1:lcount], [xlim[1], halfx], depth-1, mincount);
        else
            abin[1:lcount] = bbin[1:lcount];
            lbinends = [lcount];
        end
        if rcount >= mincount
            (abin[(lcount+1):N], rbinends) = partition_nodes_in_bins_1d(nodes, bbin[(lcount+1):N], [halfx, xlim[2]], depth-1, mincount);
        else
            abin[(lcount+1):N] = bbin[(lcount+1):N];
            rbinends = [rcount];
        end
        binends = append!(lbinends, rbinends .+ lcount);
        
    else
        abin = bbin;
        binends = [lcount, N];
    end
    
    return (abin, binends);
end

# recursively subdivide bins of nodes
function partition_nodes_in_bins_2d(nodes, abin, xlim, ylim, depth, mincount)
    halfx = (xlim[1] + xlim[2])/2;
    halfy = (ylim[1] + ylim[2])/2;
    N = length(abin);
    
    which_bin = zeros(Int, N) # Which bin will the node go in
    bin_counts = zeros(Int, 4) # How many go in each
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
                (abin[starts[i]:ends[i]], tmpbinends) = partition_nodes_in_bins_2d(nodes, bbins[i], xlims[i], ylims[i], depth-1, mincount);
                if i>1
                    append!(binends, tmpbinends .+ ends[i-1]);
                else
                    binends = tmpbinends;
                end
                
            else
                abin[starts[i]:ends[i]] = bbins[i];
                append!(binends, [ends[i]]);
            end
        end
        
    else
        abin = [bbins[1]; bbins[2]; bbins[3]; bbins[4]];
        binends = ends;
    end
    
    return (abin, binends);
end

# recursively subdivide bins of nodes
function partition_nodes_in_bins_3d(nodes, abin, xlim, ylim, zlim, depth, mincount)
    halfx = (xlim[1] + xlim[2])/2;
    halfy = (ylim[1] + ylim[2])/2;
    halfz = (zlim[1] + zlim[2])/2;
    N = length(abin);
    
    which_bin = zeros(Int, N) # Which bin will the node go in
    bin_counts = zeros(Int, 8) # How many go in each
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
                (abin[starts[i]:ends[i]], tmpbinends) = partition_nodes_in_bins_3d(nodes, bbins[i], xlims[i], ylims[i], zlims[i], depth-1, mincount);
                if i>1
                    append!(binends, tmpbinends .+ ends[i-1]);
                else
                    binends = tmpbinends;
                end
            else
                abin[starts[i]:ends[i]] = bbins[i];
                append!(binends, [ends[i]]);
            end
        end
        
    else
        abin = [bbins[1]; bbins[2]; bbins[3]; bbins[4]; bbins[5]; bbins[6]; bbins[7]; bbins[8]];
        binends = ends;
    end
    
    return (abin, binends);
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
    d = v[:,1];
    A = [v[:,2].-d v[:,3].-d v[:,4].-d];
    
    np = length(r);
    mv = zeros(3,np);
    for i=1:np
        tmp = [(r[i]+1)/2, (s[i]+1)/2, (t[i]+1)/2];
        mv[:,i] = A*tmp + d;
    end
    
    x = mv[1,:]
    y = mv[2,:]
    z = mv[3,:]
    
    return (x, y, z);
end

# Returns true if the nodes are within tol of each other.
function is_same_node(x1, x2, tol)
    return sum(abs.(x1 - x2)) < tol
end

# Returns true if the centroid of each line is close enough.
function is_same_face_center(l1, l2, tol)
    n1 = size(l1,2);
    n2 = size(l2,2);
    center1 = zeros(size(l1,1));
    center2 = zeros(size(l2,1));
    for i=1:n1
        center1 = center1 .+ l1[:,i];
    end
    center1 = center1 ./ n1;
    
    for i=1:n2
        center2 = center2 .+ l2[:,i];
    end
    center2 = center2 ./ n2;
    
    return is_same_node(center1, center2, tol);
end

# Returns true if the two node lists have at least two of the same nodes.
function is_same_line(l1, l2, tol)
    found = 0;
    n1 = size(l1,2);
    n2 = size(l2,2);
    for i=1:n1
        for j=1:n2
            if is_same_node(l1[:,i], l2[:,j], tol)
                found += 1;
            end
            if found >= 2
                return true;
            end
        end
    end
    
    return false;
end

# Returns true if the two node lists have at least three of the same nodes.
function is_same_plane(p1, p2, tol)
    found = 0;
    n1 = size(p1,2);
    n2 = size(p2,2);
    for i=1:n1
        for j=1:n2
            if is_same_node(p1[:,i], p2[:,j], tol)
                found += 1;
            end
            if found >= 3
                return true;
            end
        end
    end
    
    return false;
end

function which_refel_face(f2glb, nodes, ref, elem)
    # Check f2glb against the face2local in refel
    fnodes = nodes[:,f2glb];
    for fi=1:ref.Nfaces
        refnodes = nodes[:, elem[ref.face2local[fi]]];
        
        if is_same_face(fnodes, refnodes, ref.dim)
            return fi;
        end
    end
    
    printerr("Couldn't match face when building grid (see which_refel_face() in grid.jl)");
end

function is_same_face(f1, f2, dim)
    tol = 1e-12;
    if dim == 1 # one point
        return is_same_node(f1, f2, tol);
    elseif dim == 2 # same line(two same points)
        return is_same_line(f1, f2, tol);
    elseif dim == 3 # same plane(three same points)
        return is_same_plane(f1, f2, tol);
    end
end

# Adds a boundary ID to some region. Find boundary points satifying on_bdry and moves them to a new set for this bid.
function add_boundary_ID_to_grid(bid, on_bdry, grid)
    # Find if this bid exists. If so, just add points to it, removing from others.
    ind = indexin([bid], grid.bids)[1];
    nbids = length(grid.bids);
    if ind === nothing
        # This is a new bid, add to bids, bdry, bdryface, bdrynorm
        ind = nbids + 1;
        nbids += 1;
        push!(grid.bids, bid);
        push!(grid.bdry, zeros(Int, 0));
        push!(grid.bdryface, zeros(Int, 0));
        push!(grid.bdrynorm, zeros(config.dimension, 0));
    end
    
    # Search all other bids for nodes and faces on this segment. Remove them there and add them here.
    # First find indices and count them. Then move.
    move_nodes = Array{Array{Int,1},1}(undef,nbids);
    node_count = zeros(Int, nbids);
    move_faces = Array{Array{Int,1},1}(undef,nbids);
    face_count = zeros(Int, nbids);
    for i=1:nbids
        bi = grid.bids[i];
        move_nodes[i] = [];
        move_faces[i] = [];
        if bi != bid
            # First the nodes
            for j=1:length(grid.bdry[i])
                nj = grid.bdry[i][j];
                if config.dimension == 1
                    if on_bdry(grid.allnodes[1, nj])
                        push!(move_nodes[i], nj);
                        node_count[i] += 1;
                    end
                elseif config.dimension == 2
                    if on_bdry(grid.allnodes[1, nj], grid.allnodes[2, nj])
                        push!(move_nodes[i], nj);
                        node_count[i] += 1;
                    end
                elseif config.dimension == 3
                    if on_bdry(grid.allnodes[1, nj], grid.allnodes[2, nj], grid.allnodes[3, nj])
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
                fcenter = zeros(size(grid.allnodes,1));
                for ni=1:nfp
                    fcenter = fcenter + grid.allnodes[:,grid.face2glb[ni,1,fj]];
                end
                fcenter = fcenter./nfp;
                
                if config.dimension == 1
                    isbdryface = on_bdry(fcenter[1]);
                elseif config.dimension == 2
                    isbdryface = on_bdry(fcenter[1], fcenter[2]);
                elseif config.dimension == 3
                    isbdryface = on_bdry(fcenter[1], fcenter[2], fcenter[3]);
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
                newbdrynorm = zeros(config.dimension, size(grid.bdrynorm[i],2) - numremove);
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