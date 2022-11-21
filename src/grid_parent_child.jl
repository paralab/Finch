#=
Subdivide each element in the grid into several children, making a structured parent/child relation.
This is for higher order FV use.
Important: Since this is for FV, only first order elements should be used. Only vertex nodes will be 
considered and the resulting grid will only have vertex nodes.
=#

# # Struct for keeping the parent info.
# struct ParentMaps
#     num_parents::Int                # Number of parents including ghosts
#     children_per_parent::Int        # Number of children per parent assuming similar element types
#     faces_per_parent::Int           # Number of child faces in a parent
#     patch_size::Int                 # Number of children in a patch (parent plus neighboring parents)
    
#     child2parent::Array{Int,2}      # size = (2, allChildren) 1 = index of parent, 2 = location within parent
#     parent2child::Array{Int,2}      # size = (myChildren, allParents) global index of each child in each parent
#     parent2face::Array{Int,2}       # size = (myFaces, allParents) global index of each face in each parent
#     cface2pface::Array{Int,3}       # size = (elementFaces, myChildren, allParents) relative face index in parent
#     parent2neighbor::Array{Int,2}   # size = (outerFaces, allParents) index of neighboring parents
    
#     patches::Array{Int, 2}          # size = (outerfaces*neighborChildren) local patch around each parent
#     leftCells::Vector               # Patch indices for left and right cell groups for each face
#     rightCells::Vector              #
    
#     face_neighborhoods::Matrix      # indices of a group of elements to the left and right of each face
# end

# Divides all elements in a grid
function divide_parent_grid(grid, order)
    # Types used in grid
    int_type = finch_state.config.index_type;
    float_type = finch_state.config.float_type;
    
    # The convention for n numbers is: N->global, n->local
    dim = size(grid.allnodes,1);
    Nparent = size(grid.glbvertex,2);
    nvertex = size(grid.glbvertex,1); # number of vertices for each element assuming one element type.
    nneighbor = size(grid.element2face,1);
    Npnodes = size(grid.allnodes,2);
    Npfaces = size(grid.face2element,2);
    nbids = length(grid.bdry);
    nchildren = 0; # will be set below
    
    if order < 2 && dim==1
        #nothing needs to happen to the grid, but we should make the maps anyway
        nchildren = 1;
        nfaces = 2;
        Ncfaces = Npfaces;
        ncells_in_patch = (nneighbor+1)*nchildren;
        
        c2p = ones(Int, 2, Nparent);
        c2p[1,:] = 1:Nparent;
        p2c = Array(1:Nparent)';
        p2f = grid.element2face;
        cf2pf = zeros(Int, nneighbor, 1, Nparent);
        p2n = zeros(Int, nneighbor, Nparent);
        for i=1:Nparent
            cf2pf[:,1,i] .= 1:Ncfaces;
            for j=1:nneighbor
                fid = grid.element2face[j,i];
                if grid.face2element[1,fid] == i
                    p2n[j,i] = grid.face2element[2,fid];
                else
                    p2n[j,i] = grid.face2element[1,fid];
                end
            end
        end
        
        tmp_parent_maps = ParentMaps(Nparent, nchildren, nfaces, ncells_in_patch, c2p, p2c, p2f, cf2pf, p2n, 
            zeros(Int,0,0),Vector{Vector{Int}}(undef,0),Vector{Vector{Int}}(undef,0),Matrix{Vector{Int}}(undef,0,0)); # no patches built yet
        
        # build patches
        patches = zeros(Int, ncells_in_patch, Nparent);
        for ei=1:Nparent
            patches[:, ei] = build_local_patch(tmp_parent_maps, grid, ei);
        end
        (left_cells, right_cells) = get_left_right_cells(dim, nfaces, nneighbor, order);
        
        # Build left and right neighborhoods for each face
        f_neighborhoods = fill(zeros(Int,0), 2, Ncfaces);
        for parentid=1:Nparent
            for fi=1:nfaces
                fid = p2f[fi,parentid];
                f_neighborhoods[1,fid] = patches[left_cells[fi], parentid];
                f_neighborhoods[2,fid] = patches[right_cells[fi], parentid];
                
                # remove zeros
                first_zero = 1;
                for i=1:length(f_neighborhoods[1,fid])
                    if f_neighborhoods[1,fid][i] == 0
                        break;
                    end
                    first_zero += 1;
                end
                f_neighborhoods[1,fid] = f_neighborhoods[1,fid][1:(first_zero-1)];
                
                first_zero = 1;
                for i=1:length(f_neighborhoods[2,fid])
                    if f_neighborhoods[2,fid][i] == 0
                        break;
                    end
                    first_zero += 1;
                end
                f_neighborhoods[2,fid] = f_neighborhoods[2,fid][1:(first_zero-1)];
            end
        end
        
        parent_maps = ParentMaps(Nparent, nchildren, Ncfaces, ncells_in_patch, c2p, p2c, p2f, cf2pf, p2n, 
                                patches, left_cells, right_cells, f_neighborhoods);
        return (parent_maps, grid);
    end
    
    if dim==1
        level = order-1;
        nchildren = 1 + level;
        nfaces = 2 + level;
        Ncnodes = Nparent * nchildren + 1;
        Ncfaces = Npfaces + level*Nparent;
        ncvertex = 2;
        nfacevertex = 1;
        nfaceperchild = 2;
        cfaceperpface = 1;
    elseif dim==2
        if order > 6
            printerr("Orders greater then 6 are not ready for 2D FV. Changing to 6.")
            order = 6;
        end
        nchildren = 4; # same for triangles and quads
        if nvertex==3 # triangle
            nfaces = 9;
            Ncnodes = Npnodes + Npfaces; # add one node per face
            Ncfaces = Npfaces*2 + Nparent*3; # each parent face becomes 2, plus 3 internal
            ncvertex = 3;
            nfacevertex = 2;
            nfaceperchild = 3;
            cfaceperpface = 2;
        elseif nvertex==4 # quad
            nfaces = 12;
            Ncnodes = Npnodes + Npfaces + Nparent; # add one per face plus one per parent(center)
            Ncfaces = Npfaces*2 + Nparent*4; # each parent face becomes 2, plus 4 internal
            ncvertex = 4;
            nfacevertex = 2;
            nfaceperchild = 4;
            cfaceperpface = 2;
        end
    else #dim==3
        ### NOT READY
        printerr("Not ready to divide 3D grid for higher order FV. TODO")
        return (nothing, grid);
        if nvertex==4 # tet
            nchildren = 4; # divide into 4 hexes
            nfaces = 18;
            Ncnodes = Npnodes; # TODO
            Ncfaces = Npfaces*3 + Nparent*6; # each parent face becomes 3, plus 6 internal
            ncvertex = 8;
            nfacevertex = 4;
            nfaceperchild = 6;
            cfaceperpface = 3;
        elseif nvertex==8 # hex
            nchildren = 8;
            nfaces = 36;
            Ncnodes = Npnodes; # TODO
            Ncfaces = Npfaces*4 + Nparent*12; # each parent face becomes 4, plus 12 internal
            ncvertex = 8;
            nfacevertex = 4;
            nfaceperchild = 6;
            cfaceperpface = 4;
        end
    end
    Nchildren = nchildren * Nparent;
    
    # Pieces needed for the ParentGrid struct
    c2p = zeros(Int, 2, Nchildren);         # child to parent
    p2c = zeros(Int, nchildren, Nparent);   # parent to child
    p2f = zeros(Int, nfaces, Nparent);      # parent to face
    cf2pf = zeros(Int, nneighbor, nchildren, Nparent); # cface to pface
    p2n = zeros(Int, nneighbor, Nparent);            # parent to neighbor
    
    # Pieces needed for the new child grid
    # Refer to grid.jl for details
    # allnodes = zeros(dim, Ncnodes); # First build DG grid, then remove duplicates
    bdry = Vector{Vector{int_type}}(undef,nbids);
    bdryface = Vector{Vector{int_type}}(undef,nbids);
    bdrynorm = Vector{Matrix{float_type}}(undef,nbids);
    bids = copy(grid.bids);
    loc2glb = zeros(int_type, ncvertex, Nchildren);
    # glbvertex = zeros(Int, ncvertex, Nchildren); # will be identical to loc2glb
    face2glb = zeros(int_type, nfacevertex, 1, Ncfaces);
    element2face = zeros(int_type, nfaceperchild, Nchildren);
    face2element = zeros(int_type, 2, Ncfaces);
    facenormals = zeros(dim, Ncfaces);
    faceRefelInd = zeros(int_type, 2, Ncfaces);
    facebid = zeros(int_type, Ncfaces);
    nodebid = zeros(int_type, Ncnodes);
    
    if grid.is_subgrid
        nel_owned = grid.nel_owned * nchildren;
        nel_ghost = grid.nel_ghost * nchildren;
        nface_owned = 0; # TBD
        nface_ghost = 0; # TBD
        element_owner = fill(-1, Nchildren); # same owner as parent
        grid2mesh = fill(-1, Nchildren); # This will now refer to the parent element in the mesh.
        
        num_neighbor_partitions = grid.num_neighbor_partitions;
        neighboring_partitions = grid.neighboring_partitions;
        ghost_counts = grid.ghost_counts .* nchildren;
        ghost_index = fill(zeros(int_type,2,0), num_neighbor_partitions); # TBD
    end
    
    for i=1:nbids
        if dim==1
            bdry[i] = zeros(int_type, length(grid.bdry[i])); # same number
        elseif dim==2
            bdry[i] = zeros(int_type, length(grid.bdry[i]) + length(grid.bdryface[i])); # add one for each bdry face
        else
            bdry[i] = zeros(int_type, length(grid.bdry[i])); # TODO
        end
        bdryface[i] = zeros(int_type, length(grid.bdryface[i]) * cfaceperpface);
        bdrynorm[i] = zeros(dim, length(bdry[i]));
    end
    
    # All of the arrays are set up. Now fill them.
    # 1D is simple
    if dim==1
        # Just directly make the correct allnodes
        allnodes = zeros(1, Ncnodes);
        
        face_done = zeros(Int, Npfaces); # set to fid after adding it
        next_node = 1;
        next_child = 1;
        next_face = 1;
        next_bdryface = ones(Int, length(grid.bids));
        for ei=1:Nparent
            leftnode = grid.glbvertex[1,ei];
            rightnode = grid.glbvertex[2,ei];
            leftx = grid.allnodes[1,leftnode];
            rightx = grid.allnodes[1,rightnode];
            leftface = grid.element2face[1,ei];
            rightface = grid.element2face[2,ei];
            
            # child/parent maps
            for ci=1:nchildren
                p2c[ci,ei] = next_child;
                c2p[:,next_child] = [ei, ci];
                next_child += 1;
            end
            
            # neighbors
            if grid.face2element[1,leftface] == ei
                p2n[1,ei] = grid.face2element[2,leftface];
            else
                p2n[1,ei] = grid.face2element[1,leftface];
            end
            if grid.face2element[1,rightface] == ei
                p2n[2,ei] = grid.face2element[2,rightface];
            else
                p2n[2,ei] = grid.face2element[1,rightface];
            end
            
            # nodes and faces
            dx = (rightx-leftx)/nchildren;
            # The left face
            if face_done[leftface] == 0
                face_done[leftface] = next_face;
                p2f[1,ei] = next_face;
                
                allnodes[1, next_node] = leftx;
                loc2glb[1,p2c[1,ei]] = next_node;
                face2glb[1,1,next_face] = next_node;
                element2face[1, p2c[1,ei]] = next_face;
                face2element[1, next_face] = p2c[1,ei];
                facenormals[1,next_face] = leftx < rightx ? -1 : 1 #grid.facenormals[1,leftface];
                faceRefelInd[1, next_face] = 1; # left of element
                facebid[next_face] = grid.facebid[leftface];
                nodebid[next_node] = grid.nodebid[leftnode];
                
                fbid = facebid[next_face];
                if fbid > 0
                    bdryface[fbid][next_bdryface[fbid]] = next_face;
                    # since a 1d bdry face is only one node
                    bdry[fbid][next_bdryface[fbid]] = next_node;
                    bdrynorm[fbid][next_bdryface[fbid]] = facenormals[1,next_face];
                    
                    next_bdryface[fbid] += 1;
                end
                
                next_node += 1;
                next_face += 1;
                
            else # parent face already set
                fid = face_done[leftface];
                p2f[1,ei] = fid;
                loc2glb[1,p2c[1,ei]] = face2glb[1,1,fid];
                element2face[1, p2c[1,ei]] = fid;
                face2element[2, fid] = p2c[1,ei];
                faceRefelInd[2, fid] = 1; # left of element
            end
            
            # The interior
            for ni=1:nchildren-1
                p2f[ni+1,ei] = next_face;
                
                allnodes[1, next_node] = leftx + ni*dx; 
                loc2glb[2,p2c[ni,ei]] = next_node;
                loc2glb[1,p2c[ni+1,ei]] = next_node;
                face2glb[1,1,next_face] = next_node;
                element2face[2, p2c[ni,ei]] = next_face;
                element2face[1, p2c[ni+1,ei]] = next_face;
                face2element[1, next_face] = p2c[ni,ei];
                face2element[2, next_face] = p2c[ni+1,ei];
                facenormals[1,next_face] = leftx < rightx ? 1 : -1  # will point right because building from left
                faceRefelInd[1, next_face] = 2; # right of element
                faceRefelInd[2, next_face] = 1; # left of element
                facebid[next_face] = 0;
                nodebid[next_node] = 0;
                
                next_node += 1;
                next_face += 1;
            end
            
            # The right face
            if face_done[rightface] == 0
                face_done[rightface] = next_face;
                p2f[nchildren+1,ei] = next_face;
                
                allnodes[1, next_node] = rightx;
                loc2glb[2,p2c[nchildren,ei]] = next_node;
                face2glb[1,1,next_face] = next_node;
                element2face[2, p2c[nchildren,ei]] = next_face;
                face2element[1, next_face] = p2c[nchildren,ei];
                facenormals[1,next_face] = leftx < rightx ? 1 : -1 
                faceRefelInd[1, next_face] = 2; # right of element
                facebid[next_face] = grid.facebid[rightface];
                nodebid[next_node] = grid.nodebid[rightnode];
                
                fbid = facebid[next_face];
                if fbid > 0
                    bdryface[fbid][next_bdryface[fbid]] = next_face;
                    # since a 1d bdry face is only one node
                    bdry[fbid][next_bdryface[fbid]] = next_node;
                    bdrynorm[fbid][next_bdryface[fbid]] = facenormals[1,next_face];
                    
                    next_bdryface[fbid] += 1;
                end
                
                next_node += 1;
                next_face += 1;
                
            else # parent face already set
                fid = face_done[rightface];
                p2f[nchildren,ei] = fid;
                loc2glb[2,p2c[nchildren,ei]] = face2glb[1,1,fid];
                element2face[2, p2c[nchildren,ei]] = fid;
                face2element[2, fid] = p2c[nchildren,ei];
                faceRefelInd[2, fid] = 2; # right of element
            end
            
        end # parent loop
        
    elseif dim==2
        # Make the DG nodes, then remove duplicates
        if length(grid.element2face[:,1]) == 3 # triangle
            tmpallnodes = zeros(2, Nparent * 6);
        else
            tmpallnodes = zeros(2, Nparent * 9);
        end
        
        face_done = zeros(Int, 2, Npfaces); # set to fid after adding it (each parent face -> 2 child faces)
        vertex_done = zeros(Int, Npnodes); # set to node index after adding it
        next_node = 1;
        next_child = 1;
        next_face = 1;
        next_bdryface = ones(Int, length(grid.bids));
        for ei=1:Nparent
            pnodes = grid.glbvertex[:,ei];
            pnodex = grid.allnodes[:, pnodes];
            pfaces = grid.element2face[:,ei];
            
            # child/parent maps
            #          1
            #         / \
            #    n1  /c1 \  n3
            #       /_____\
            #      / \c4 / \
            #     /c2 \ /c3 \ 
            #    2_____V_____3 
            #         n2
            #       
            #        n3
            #    4---------3
            #    | c4 | c3 |
            # n4 |---------| n2
            #    | c1 | c2 | 
            #    1---------2
            #        n1
            for ci=1:nchildren
                p2c[ci,ei] = next_child;
                c2p[:,next_child] = [ei, ci];
                next_child += 1;
            end
            
            # neighbors    1
            local_face_ind = zeros(Int, length(pfaces));
            for fi=1:length(pfaces)
                # figure out the local index of the corresponding face in grid
                lfi = 0;
                for fj=1:length(pfaces)
                    if pnodes[fj] in grid.face2glb[:,1,pfaces[fi]]
                        if fj < length(pfaces) && pnodes[fj+1] in grid.face2glb[:,1,pfaces[fi]]
                            lfi = fj;
                        elseif fj == length(pfaces) && pnodes[1] in grid.face2glb[:,1,pfaces[fi]]
                            lfi = fj;
                        else
                            lfi = fj>1 ? fj-1 : length(pfaces);
                        end
                    end
                end
                if lfi == 0
                    println("error: local face index was 0");
                    lfi = fi;
                end
                local_face_ind[fi] = lfi;
                
                if grid.face2element[1,pfaces[fi]] == ei
                    p2n[lfi,ei] = grid.face2element[2,pfaces[fi]];
                else
                    p2n[lfi,ei] = grid.face2element[1,pfaces[fi]];
                end
            end
            # reorder pfaces to match local
            tmp = zeros(Int, length(pfaces))
            for fi=1:length(pfaces)
                tmp[local_face_ind[fi]] = pfaces[fi];
            end
            pfaces = tmp;
            
            # println("lfi: "*string(local_face_ind));
            # println("pfaces: "*string(pfaces));
            # println("pnodes: "*string(pnodes));
            
            # nodes and faces
            # Here things have to be done separately for triangles and quads
            if length(pfaces) == 3 # triangle
                child_nodes = zeros(2,6);
                child2loc = zeros(Int, 3, 4);
                child2face = zeros(Int, 3, 4);
                
                child_nodes[:,1:3] = pnodex[:,1:3];
                child_nodes[:,4] = (pnodex[:,1] + pnodex[:,2])/2;
                child_nodes[:,5] = (pnodex[:,2] + pnodex[:,3])/2;
                child_nodes[:,6] = (pnodex[:,1] + pnodex[:,3])/2;
                tmpallnodes[:,((ei-1)*6 + 1):(ei*6)] = child_nodes;
                
                child2loc[:,1] = [1, 4, 6];
                child2loc[:,2] = [2, 5, 4];
                child2loc[:,3] = [3, 6, 5];
                child2loc[:,4] = [4, 5, 6];
                child2face[:,1] = [1, 7, 6];
                child2face[:,2] = [3, 8, 2];
                child2face[:,3] = [5, 9, 4];
                child2face[:,4] = [8, 9, 7];
                face2loc = [1 4; 4 2; 2 5; 5 3; 3 6; 6 1; 4 6; 4 5; 5 6]; # local index of face nodes
                face2child = [1 0; 2 0; 2 0; 3 0; 3 0; 1 0; 1 4; 2 4; 3 4]; # elements having this face (local index)
                
                parent2glb = ((ei-1)*6+1):(ei*6); # global index of parent nodes
                parentface2glb = zeros(2,9); # global index of all face nodes in parent
                for fi=1:9
                    parentface2glb[:,fi] = parent2glb[face2loc[fi,:]];
                end
                
                # The external faces
                for fi=1:3
                    if face_done[1,pfaces[fi]] == 0
                        face_done[1,pfaces[fi]] = next_face;
                        face_done[2,pfaces[fi]] = next_face+1;
                        
                        p2f[fi*2-1,ei] = next_face;
                        p2f[fi*2,ei] = next_face+1;
                        
                        face2element[1, next_face] = p2c[face2child[2*fi-1,1],ei];
                        face2element[1, next_face+1] = p2c[face2child[2*fi,1],ei];
                        # element2face[1, p2c[1,ei]] = next_face; # do later
                        facenormals[:,next_face] = grid.facenormals[:,pfaces[fi]];
                        facenormals[:,next_face+1] = grid.facenormals[:,pfaces[fi]];
                        faceRefelInd[1, next_face] = 1;
                        faceRefelInd[1, next_face+1] = 3;
                        facebid[next_face] = grid.facebid[pfaces[fi]];
                        facebid[next_face+1] = grid.facebid[pfaces[fi]];
                        
                        fbid = facebid[next_face];
                        if fbid > 0
                            bdryface[fbid][next_bdryface[fbid]] = next_face;
                            bdryface[fbid][next_bdryface[fbid]+1] = next_face+1;
                            
                            next_bdryface[fbid] += 2;
                        end
                        
                        next_face += 2;
                        
                    else # parent face already set
                        same_orientation = faceRefelInd[face_done[1,pfaces[fi]]] == 1;
                        if same_orientation
                            fid1 = face_done[1,pfaces[fi]];
                            fid2 = face_done[2,pfaces[fi]];
                        else
                            fid1 = face_done[2,pfaces[fi]];
                            fid2 = face_done[1,pfaces[fi]];
                        end
                        
                        p2f[fi*2-1,ei] = fid1;
                        p2f[fi*2,ei] = fid2;
                        
                        face2element[2, fid1] = p2c[face2child[2*fi-1,1],ei];
                        face2element[2, fid2] = p2c[face2child[2*fi,1],ei];
                        faceRefelInd[2, fid1] = 1;
                        faceRefelInd[2, fid2] = 3;
                    end
                end # exterior faces
                
                # Interior faces
                p2f[7:9,ei] = next_face:(next_face+2);
                face2element[1, next_face:(next_face+2)] = p2c[face2child[7:9,1],ei];
                face2element[2, next_face:(next_face+2)] = p2c[face2child[7:9,2],ei];
                facenormals[:,next_face] = grid.facenormals[:,pfaces[2]];
                facenormals[:,next_face+1] = grid.facenormals[:,pfaces[3]];
                facenormals[:,next_face+2] = grid.facenormals[:,pfaces[1]];
                faceRefelInd[:, next_face] = [2,3];
                faceRefelInd[:, next_face+1] = [2,1];
                faceRefelInd[:, next_face+2] = [2,2];
                facebid[next_face:(next_face+2)] = [0,0,0];
                
                next_face += 3;
                
                # Still need to do: loc2glb, face2glb, element2face, bdry, bdrynorm
                for ci=1:4
                    loc2glb[:, p2c[ci,ei]] = parent2glb[child2loc[:,ci]];
                    element2face[:, p2c[ci,ei]] = p2f[child2face[:,ci],ei];
                end
                for fi=1:9
                    face2glb[:, 1, p2f[fi,ei]] = parent2glb[face2loc[fi,:]];
                end
                
                # println("p2glb for "*string(ei)*": "*string(parent2glb))
                    
            else # quad
                child_nodes = zeros(2,9);
                child2loc = zeros(Int, 4, 4);
                child2face = zeros(Int, 4, 4);
                
                #       6    5 
                #    4----7----3
                # 7  | c4 | c3 | 4
                #    8----9----6
                # 8  | c1 | c2 | 3
                #    1----5----2
                #       1    2
                child_nodes[:,1:4] = pnodex[:,1:4];
                child_nodes[:,5] = (pnodex[:,1] + pnodex[:,2])/2;
                child_nodes[:,6] = (pnodex[:,2] + pnodex[:,3])/2;
                child_nodes[:,7] = (pnodex[:,3] + pnodex[:,4])/2;
                child_nodes[:,8] = (pnodex[:,1] + pnodex[:,4])/2;
                child_nodes[:,9] = (pnodex[:,1] + pnodex[:,3])/2;
                tmpallnodes[:,((ei-1)*9 + 1):(ei*9)] = child_nodes;
                
                child2loc[:,1] = [1, 5, 9, 8];
                child2loc[:,2] = [2, 6, 9, 5];
                child2loc[:,3] = [3, 7, 9, 6];
                child2loc[:,4] = [4, 8, 9, 7];
                child2face[:,1] = [1, 9, 12, 8];
                child2face[:,2] = [3, 10, 9, 2];
                child2face[:,3] = [5, 11, 10, 4];
                child2face[:,4] = [7, 12, 11, 6];
                face2loc = [1 5; 5 2; 2 6; 6 3; 3 7; 7 4; 4 8; 8 1; 5 9; 6 9; 7 9; 8 9]; # local index of face nodes
                face2child = [1 0; 2 0; 2 0; 3 0; 3 0; 4 0; 4 0; 1 0; 1 2; 2 3; 3 4; 4 1]; # elements having this face (local index)
                
                parent2glb = ((ei-1)*9+1):(ei*9); # global index of parent nodes
                parentface2glb = zeros(2,12); # global index of all face nodes in parent
                for fi=1:12
                    parentface2glb[:,fi] = parent2glb[face2loc[fi,:]];
                end
                
                # The external faces
                tmpfri = [2, 3, 4, 1];
                for fi=1:4
                    if face_done[1,pfaces[fi]] == 0
                        face_done[1,pfaces[fi]] = next_face;
                        face_done[2,pfaces[fi]] = next_face+1;
                        
                        p2f[fi*2-1,ei] = next_face;
                        p2f[fi*2,ei] = next_face+1;
                        
                        face2element[1, next_face] = p2c[face2child[2*fi-1,1],ei];
                        face2element[1, next_face+1] = p2c[face2child[2*fi,1],ei];
                        # element2face[1, p2c[1,ei]] = next_face; # do later
                        facenormals[:,next_face] = grid.facenormals[:,pfaces[fi]];
                        facenormals[:,next_face+1] = grid.facenormals[:,pfaces[fi]];
                        faceRefelInd[1, next_face] = tmpfri[fi];
                        faceRefelInd[1, next_face+1] = tmpfri[fi];
                        facebid[next_face] = grid.facebid[pfaces[fi]];
                        facebid[next_face+1] = grid.facebid[pfaces[fi]];
                        
                        fbid = facebid[next_face];
                        if fbid > 0
                            bdryface[fbid][next_bdryface[fbid]] = next_face;
                            bdryface[fbid][next_bdryface[fbid]+1] = next_face+1;
                            
                            next_bdryface[fbid] += 2;
                        end
                        
                        next_face += 2;
                        
                    else # parent face already set
                        same_orientation = false;
                        # Check if the first corner matches the first child's face
                        first_face_nodes = tmpallnodes[:,face2glb[:,1,face_done[1,pfaces[fi]]]];
                        if fi==1
                            corner = child_nodes[:,1];
                        elseif fi==2
                            corner = child_nodes[:,2];
                        elseif fi==3
                            corner = child_nodes[:,3];
                        else
                            corner = child_nodes[:,4];
                        end
                        for ffni=1:size(first_face_nodes,2)
                            if norm(corner - first_face_nodes[:,ffni]) < 1e-16
                                same_orientation = true;
                            end
                        end
                        
                        if same_orientation
                            fid1 = face_done[1,pfaces[fi]];
                            fid2 = face_done[2,pfaces[fi]];
                        else
                            fid1 = face_done[2,pfaces[fi]];
                            fid2 = face_done[1,pfaces[fi]];
                        end
                        
                        p2f[fi*2-1,ei] = fid1;
                        p2f[fi*2,ei] = fid2;
                        
                        face2element[2, fid1] = p2c[face2child[2*fi-1,1],ei];
                        face2element[2, fid2] = p2c[face2child[2*fi,1],ei];
                        faceRefelInd[2, fid1] = tmpfri[fi];
                        faceRefelInd[2, fid2] = tmpfri[fi];
                    end
                end # exterior faces
                
                # Interior faces
                p2f[9:12,ei] = next_face:(next_face+3);
                face2element[1, next_face:(next_face+3)] = p2c[face2child[9:12,1],ei];
                face2element[2, next_face:(next_face+3)] = p2c[face2child[9:12,2],ei];
                
                p2 = child_nodes[:,9];
                for inter=0:3
                    p1 = child_nodes[:,5+inter];
                    p1p2_dist = sqrt((p2[1]-p1[1])^2 + (p2[2]-p1[2])^2);
                    facenormals[:,next_face+inter] = [(p2[2]-p1[2])/p1p2_dist, (p1[1]-p2[1])/p1p2_dist];
                end
                
                faceRefelInd[:, next_face] = [3,1];
                faceRefelInd[:, next_face+1] = [4,2];
                faceRefelInd[:, next_face+2] = [1, 3];
                faceRefelInd[:, next_face+3] = [2, 4];
                facebid[next_face:(next_face+3)] = [0,0,0,0];
                
                next_face += 4;
                
                # Still need to do: loc2glb, face2glb, element2face, bdry, bdrynorm
                for ci=1:4
                    loc2glb[:, p2c[ci,ei]] = parent2glb[child2loc[:,ci]];
                    element2face[:, p2c[ci,ei]] = p2f[child2face[:,ci],ei];
                end
                for fi=1:12
                    face2glb[:, 1, p2f[fi,ei]] = parent2glb[face2loc[fi,:]];
                end
            end
            
        end # parent loop
        
        # println(face2glb[:,1,:])
        
        # Since a DG grid was made, remove duplicate nodes and update loc2glb and face2glb
        (allnodes, loc2glb, face2glb) = remove_duplicate_nodes(tmpallnodes, loc2glb, other2glb=face2glb);
        
    elseif dim==3
        #TODO
        printerr("Sorry, higher order FV is not ready for 3D (TODO: build parent/child grid)");
    end
    
    if grid.is_subgrid
        # element owner is same as parent
        for i=1:Nparent
            for j=1:nchildren
                ei = p2c[j,i];
                element_owner[ei] = grid.element_owner[i];
                grid2mesh[ei] = grid.grid2mesh[i];
            end
        end
        
        # Count owned/ghost faces
        for fi=1:Ncfaces
            e1 = face2element[1,fi];
            e2 = face2element[2,fi];
            if element_owner[e1] == -1
                nface_owned += 1;
            elseif e2 > 0 && element_owner[e2] == -1
                nface_owned += 1;
            else
                nface_ghost += 1;
            end
        end
        
        # Build ghost indices for exchange
        # Loop over parent ghost indices and add their children in the p2c order
        # This should be the same for both partitions, but needs to be checked.
        for i=1:grid.num_neighbor_partitions
            for j=1:grid.ghost_counts[i]
                gij=grid.ghost_index[i][1,j];
                oij=grid.ghost_index[i][2,j];
                block = zeros(Int, 2, nchildren);
                for ei=1:nchildren
                    block[:,ei] = [p2c[ei,gij], p2c[ei,oij]];
                end
                ghost_index[i] = hcat(ghost_index[i], block);
            end
        end
        
        child_grid = Grid(finch_state.config.float_type, allnodes, bdry, bdryface, bdrynorm, bids, nodebid, loc2glb, loc2glb, face2glb, 
                        element2face, face2element, facenormals, faceRefelInd, facebid, 
                        true, Array(1:nel_owned), grid.nel_global*nchildren, nel_owned, nel_ghost, 
                        nface_owned, nface_ghost, 0, 0, element_owner, 
                        zeros(Int,0), grid2mesh, zeros(Int,0), zeros(Int8, 0), 
                        num_neighbor_partitions, neighboring_partitions, ghost_counts, ghost_index);
        
    else # not a subgrid
        child_grid = Grid(finch_state.config.float_type, allnodes, bdry, bdryface, bdrynorm, bids, nodebid, loc2glb, loc2glb, face2glb, element2face, 
                            face2element, facenormals, faceRefelInd, facebid);
    end
    
    # Map child faces to parent faces
    for i=1:Nparent
        for j=1:nchildren
            eid = p2c[j,i];
            for k=1:nneighbor
                fid = child_grid.element2face[k, eid];
                for m=1:Ncfaces
                    if p2f[m,i] == fid
                        cf2pf[k,j,i] = m;
                        break;
                    end
                end
            end
        end
    end
    
    ncells_in_patch = (nneighbor+1)*nchildren;
    tmp_parent_maps = ParentMaps(Nparent, nchildren, Ncfaces, ncells_in_patch, c2p, p2c, p2f, cf2pf, p2n, 
        zeros(Int,0,0),Vector{Vector{Int}}(undef,0),Vector{Vector{Int}}(undef,0),Matrix{Vector{Int}}(undef,0,0)); # no patches built yet
    
    # Normal vectors may be pointing the wrong way. Reorient them.
    check_normal_vectors!(child_grid);
    
    # build patches
    patches = zeros(Int, ncells_in_patch, Nparent);
    for ei=1:Nparent
        patches[:, ei] = build_local_patch(tmp_parent_maps, child_grid, ei);
    end
    
    # patch indices of left and right cells for patch faces
    (left_cells, right_cells) = get_left_right_cells(dim, nfaces, nneighbor, order);
    
    # Build left and right neighborhoods for each face
    f_neighborhoods = fill(zeros(Int,0), 2, Ncfaces);
    for parentid=1:Nparent
        for fi=1:nfaces
            fid = p2f[fi,parentid];
            f_neighborhoods[1,fid] = patches[left_cells[fi], parentid];
            f_neighborhoods[2,fid] = patches[right_cells[fi], parentid];
            
            # remove zeros
            first_zero = 1;
            for i=1:length(f_neighborhoods[1,fid])
                if f_neighborhoods[1,fid][i] == 0
                    break;
                end
                first_zero += 1;
            end
            f_neighborhoods[1,fid] = f_neighborhoods[1,fid][1:(first_zero-1)];
            
            first_zero = 1;
            for i=1:length(f_neighborhoods[2,fid])
                if f_neighborhoods[2,fid][i] == 0
                    break;
                end
                first_zero += 1;
            end
            f_neighborhoods[2,fid] = f_neighborhoods[2,fid][1:(first_zero-1)];
        end
    end
    
    parent_maps = ParentMaps(Nparent, nchildren, Ncfaces, ncells_in_patch, c2p, p2c, p2f, cf2pf, p2n, 
                            patches, left_cells, right_cells, f_neighborhoods);
    
    return (parent_maps, child_grid);
end

# Set normal direction so that it is pointing from 1->2 matching face2element.
function check_normal_vectors!(grid)
    dim = size(grid.allnodes,1);
    nface = size(grid.face2element,2);
    for fi=1:nface
        els=grid.face2element[:,fi];
        nfp = size(grid.face2glb,1);
        np = size(grid.loc2glb,1);
        fcenter = grid.allnodes[:,grid.face2glb[1,1,fi]] ./ nfp; # center of face
        for i=2:nfp
            fcenter += grid.allnodes[:,grid.face2glb[i,1,fi]] ./ nfp;
        end
        c1 = grid.allnodes[:,grid.loc2glb[1,els[1]]] ./ np; # center of element 1
        for i=2:np
            c1 += grid.allnodes[:,grid.loc2glb[i,els[1]]] ./ np;
        end
        if els[2] == 0 # boundary
            c2 = fcenter; # center of face for boundary
        else
            c2 = grid.allnodes[:,grid.loc2glb[1,els[2]]] ./ np; # center of element 2
            for i=2:np
                c2 += grid.allnodes[:,grid.loc2glb[i,els[2]]] ./ np;
            end
        end
        
        # check and reverse if needed
        if dim==1
            grid.facenormals[1,fi] = c2[1]>c1[1] ? 1.0 : -1.0
        else
            # which is further from cell center 1, (face center + normal) or (face center - normal)
            d1 = fcenter + grid.facenormals[:,fi] - c1;
            d2 = fcenter - grid.facenormals[:,fi] - c1;
            if norm(d2) > norm(d1)
                grid.facenormals[:,fi] = -grid.facenormals[:,fi]
            end
        end
    end
end

# Builds the local patch given a central parent.
function build_local_patch(maps, grid, center)
    dim = size(grid.allnodes,1);
    if dim == 1
        # | neighbor 1 | center | neighbor 2| -> | n+1..2n | 1..n | 2n+1..3n |
        nchildren = size(maps.parent2child,1); # (assumes one element type)
        patch = zeros(Int, nchildren*3); # (assumes one element type)
        patch[1:nchildren] = maps.parent2child[:,center];
        for neighborid=1:2
            if maps.parent2neighbor[neighborid, center] == 0
                # This is a boundary face with no neighbor. Set these cells to zero.
                patch[(neighborid*nchildren+1):((neighborid+1)*nchildren)] = zeros(Int, nchildren);
            else
                # Check for orientation
                if maps.parent2face[1, center] == maps.parent2face[end, maps.parent2neighbor[neighborid,center]] || maps.parent2face[end, center] == maps.parent2face[1, maps.parent2neighbor[neighborid,center]]
                    patch[(neighborid*nchildren+1):((neighborid+1)*nchildren)] = maps.parent2child[:,maps.parent2neighbor[neighborid,center]];
                else # reversed orientation
                    patch[((neighborid+1)*nchildren):-1:(neighborid*nchildren+1)] = maps.parent2child[:,maps.parent2neighbor[neighborid,center]];
                end
            end
        end
        
    elseif dim == 2
        # Triangles and quads will be done separately
        if size(maps.parent2neighbor,1) == 3 # triangle
            #          1
            #         / \
            #    n1  /c1 \  n3
            #       /_____\
            #      / \c4 / \
            #     /c2 \ /c3 \ 
            #    2_____V_____3 
            #         n2
            # Neighbor orientation could have six configurations.
            # Determine this by face ID.
            orientation_table = [
                [1, 3, 2, 4],
                [2, 3, 1, 4],
                [2, 1, 3, 4],
                [3, 1, 2, 4],
                [3, 2, 1, 4],
                [1, 2, 3, 4]
            ]
            
            patch = zeros(Int, 4 * 4);
            patch[1:4] = maps.parent2child[:,center];
            
            neighbors = maps.parent2neighbor[:,center];
            cfaces = maps.parent2face[[1,3,5],center]; # first face on each side of parent
            orientation = [0,0,0]; # orientation index for each neighbor
            for ni=1:3
                if neighbors[ni] > 0 # not a boundary
                    for fi=1:6
                        if maps.parent2face[fi,neighbors[ni]] == cfaces[ni]
                            orientation[ni] = fi; # This is the index of the face touching cfaces[ni]
                            break;
                        end
                    end
                    # Use the orientation table to populate the patch
                    patch[(ni*4+1):((ni+1)*4)] = maps.parent2child[orientation_table[orientation[ni]], neighbors[ni]];
                end
            end
            
        else # quad
            #        n3
            #    4---------3
            #    | c4 | c3 |
            # n4 |---------| n2
            #    | c1 | c2 | 
            #    1---------2
            #        n1
            # Neighbor orientation could have eight configurations.
            # Determine this by face ID.
            orientation_table = [
                [1, 4, 3, 2],
                [2, 3, 4, 1],
                [2, 1, 4, 3],
                [3, 4, 1, 2],
                [3, 2, 1, 4],
                [4, 1, 2, 3],
                [4, 3, 2, 1],
                [1, 2, 3, 4]
            ]
            
            patch = zeros(Int, 4 * 5);
            patch[1:4] = maps.parent2child[:,center];
            
            neighbors = maps.parent2neighbor[:,center];
            cfaces = maps.parent2face[[1,3,5,7],center]; # first face on each side of parent
            orientation = [0,0,0,0]; # orientation index for each neighbor
            for ni=1:4
                if neighbors[ni] > 0 # not a boundary
                    for fi=1:8
                        if maps.parent2face[fi,neighbors[ni]] == cfaces[ni]
                            orientation[ni] = fi; # This is the index of the face touching cfaces[ni]
                            break;
                        end
                    end
                    # Use the orientation table to populate the patch
                    patch[(ni*4+1):((ni+1)*4)] = maps.parent2child[orientation_table[orientation[ni]], neighbors[ni]];
                end
            end
            
        end
        
    elseif dim == 3
        # TODO
    end
    
    return patch;
end

# Build maps for left and right cells for each face
function get_left_right_cells(dim::Int, faces::Int, neighbors::Int, order::Int)
    # Provide all of the cells that can be used for this face in a particular order.
    # Make two arrays, one for left one for right, starting from the nearest neighbor.
    left_cells = fill(zeros(Int,0), faces);
    right_cells = fill(zeros(Int,0), faces);
    if dim == 1
        left_cell_table = [
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
        right_cell_table = [
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
        for i=1:faces
            left_cells[i] = left_cell_table[order][i];
            right_cells[i] = right_cell_table[order][i];
        end
        
    elseif dim == 2
        if neighbors == 3 # triangles
            # Triangle parents have 9 faces, patches have 16 cells
            # Here Left means toward the center of the central parent
            left_cells = [
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
            right_cells = [
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
            
        else # quads
            # Quad parents have 12 faces, patches have 20 cells
            # Here Left means toward the center of the central parent
            left_cells = [
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
            right_cells = [
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
        end
        
    elseif dim == 3
        #TODO
    end
    
    # Shorten lists depending on order
    for fi=1:faces
        if length(left_cells[fi]) > order
            left_cells[fi] = left_cells[fi][1:order];
        end
        if length(right_cells[fi]) > order
            right_cells[fi] = right_cells[fi][1:order];
        end
    end
    
    return (left_cells, right_cells);
end