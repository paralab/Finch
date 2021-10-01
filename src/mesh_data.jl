#=
# A struct containing mesh information
# This does not include any information about the nodal placement inside elements.
# It is only the basic elemental mesh with vertices and faces.
=#
export MeshData, build_faces, find_boundaries, find_normals

struct MeshData
    #### Minimal required information ####
    # Nodes
    nx::Int;                    # Number of vertices
    nodes::Array{Float64,2};    # vertex locations (array has size (dim,nx))
    indices::Array{Int,1};      # vertex indices may not be in order
    # Elements
    nel::Int;                   # Number of elements
    elements::Array{Int,2};     # Element vertex mapping (array has size (Np, nel))*assumes only one element type
    etypes::Array{Int,1};       # Element types as defined by GMSH
    nv::Array{Int,1};           # Number of vertices for each element. Only different if they are different types,
    
    #### Optional information that will be built if not provided ####
    invind::Array{Int,1}        # Inverse of indices, maps vertex index to position in nodes array (invind[indices[i]] = i)
    face2vertex::Array{Int,2}   # Vertices defining each face (array has size (Nfp, Nfaces))
    face2element::Array{Int,2}  # Indices of elements on each side of the face. If 0, it is a boundary face. (size is (2,Nfaces))
    element2face::Array{Int,2}  # Indices of faces on each side of the element. (size is (NfacesPerElement, nel))
    normals::Array{Float64,2}   # Normal vectors for each face pointing from first to second in face2element order (size is (dim, Nfaces))
    bdryID::Array{Int,1};       # Boundary ID for each face (0=interior face)
    
    # The minimal constructor needs to build the optional information.
    # Note: Must uncomment to build.
    MeshData(n, x, ind, ne, el, et, v) = (
        # inv = invert_index(ind);
        # face2v = Array{Int,2}(undef,0,0);
        # face2e = Array{Int,2}(undef,0,0);
        # e2face = Array{Int,2}(undef,0,0);
        # norms = Array{Float64,2}(undef,0,0);
        # bdry = Array{Int,1}(undef,0);
        
        # uncomment these to compute. WARNING: can be slow
        inv = invert_index(ind);
        (face2v, face2e, e2face) = build_faces(ne, el, et);
        norms = find_normals(face2v, x);
        bdry = find_boundaries(face2e);
        new(n, x, ind, ne, el, et, v, inv, face2v, face2e, e2face, norms, bdry);
    )
    # The complete constructor
    MeshData(n, x, ind, ne, el, et, v, inv, face2v, face2e, e2face, norms, bdry) = (
        new(n, x, ind, ne, el, et, v, inv, face2v, face2e, e2face, norms, bdry);
    )
end

# numbers of nodes and faces for first and second order elements as defined by GMSH
# line, triangle, quad, tet, hex, prism, 5-pyramid
etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5]; # number of vertices
etypetonf = [2, 3, 4, 4, 6, 5, 5, 2, 3, 4, 4, 6, 5, 5, 1, 4, 6, 5, 5]; # number of faces
etypetonfn= [1, 2, 2, 3, 4, 4, 4, 1, 2, 2, 3, 4, 4, 4, 1, 2, 2, 4, 4]; # number of vertices for each face (except prism and 5-pyramids!)
etypetodim= [1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3]; # dimension of each type

# Builds invind.
function invert_index(ind)
    invind = zeros(Int, size(ind));
    for i=1:length(ind)
        invind[ind[i]] = i;
    end
    return invind;
end

# Builds faces
# For now assumes only one type of element.
function build_faces(nel, elements, etypes)
    NfacesPerElement = etypetonf[etypes[1]]; # For now assumes only one type of element.
    Nfp = etypetonfn[etypes[1]]; # For now assumes only one type of element.
    
    Nfaces = 0; # will be incremented as discovered
    e2face = zeros(Int, NfacesPerElement, nel);
    
    face2v = zeros(Int, Nfp, NfacesPerElement * nel);
    face2e = zeros(Int, 2, NfacesPerElement * nel);
    
    for ei=1:nel
        if etypes[ei] == 1 # line
            face2v[1,Nfaces+1] = elements[1, ei];
            face2v[1,Nfaces+2] = elements[2, ei];
            face2e[1,Nfaces+1] = ei; #
            face2e[1,Nfaces+2] = ei; #
            e2face[1,ei] = Nfaces+1;
            e2face[2,ei] = Nfaces+2;
            Nfaces += 2;
        elseif etypes[ei] == 2 # triangle
            face2v[:,Nfaces+1] = elements[1:2, ei];
            face2v[:,Nfaces+2] = elements[2:3, ei];
            face2v[:,Nfaces+3] = elements[[3,1], ei];
            face2e[1,Nfaces+1] = ei; #
            face2e[1,Nfaces+2] = ei; #
            face2e[1,Nfaces+3] = ei; #
            e2face[1,ei] = Nfaces+1;
            e2face[2,ei] = Nfaces+2;
            e2face[3,ei] = Nfaces+3;
            Nfaces += 3;
        elseif etypes[ei] == 3 # quad
            face2v[:,Nfaces+1] = elements[[4,1], ei];
            face2v[:,Nfaces+2] = elements[1:2, ei];
            face2v[:,Nfaces+3] = elements[2:3, ei];
            face2v[:,Nfaces+4] = elements[[4,3], ei];
            face2e[1,Nfaces+1] = ei; #
            face2e[1,Nfaces+2] = ei; #
            face2e[1,Nfaces+3] = ei; #
            face2e[1,Nfaces+4] = ei; #
            e2face[1,ei] = Nfaces+1;
            e2face[2,ei] = Nfaces+2;
            e2face[3,ei] = Nfaces+3;
            e2face[4,ei] = Nfaces+4;
            Nfaces += 4;
        elseif etypes[ei] == 4 # tet
            face2v[:,Nfaces+1] = elements[[1,3,2], ei];
            face2v[:,Nfaces+2] = elements[[2,3,4], ei];
            face2v[:,Nfaces+3] = elements[[1,2,4], ei];
            face2v[:,Nfaces+4] = elements[[1,4,3], ei];
            face2e[1,Nfaces+1] = ei; #
            face2e[1,Nfaces+2] = ei; #
            face2e[1,Nfaces+3] = ei; #
            face2e[1,Nfaces+4] = ei; #
            e2face[1,ei] = Nfaces+1;
            e2face[2,ei] = Nfaces+2;
            e2face[3,ei] = Nfaces+3;
            e2face[4,ei] = Nfaces+4;
            Nfaces += 4;
        elseif etypes[ei] == 5 # hex
            face2v[:,Nfaces+1] = elements[[1,5,8,4], ei];
            face2v[:,Nfaces+2] = elements[[2,3,7,6], ei];
            face2v[:,Nfaces+3] = elements[[1,2,6,5], ei];
            face2v[:,Nfaces+4] = elements[[3,4,8,7], ei];
            face2v[:,Nfaces+5] = elements[[1,4,3,2], ei];
            face2v[:,Nfaces+6] = elements[[5,6,7,8], ei];
            face2e[1,Nfaces+1] = ei; #
            face2e[1,Nfaces+2] = ei; #
            face2e[1,Nfaces+3] = ei; #
            face2e[1,Nfaces+4] = ei; #
            face2e[1,Nfaces+5] = ei; #
            face2e[1,Nfaces+6] = ei; #
            e2face[1,ei] = Nfaces+1;
            e2face[2,ei] = Nfaces+2;
            e2face[3,ei] = Nfaces+3;
            e2face[4,ei] = Nfaces+4;
            e2face[5,ei] = Nfaces+5;
            e2face[6,ei] = Nfaces+6;
            Nfaces += 6;
        elseif etypes[ei] == 6 # prism
            face2v[1:3,Nfaces+1] = elements[[1,3,2], ei];
            face2v[1:3,Nfaces+2] = elements[[4,5,6], ei];
            face2v[:,Nfaces+3] = elements[[1,2,5,4], ei];
            face2v[:,Nfaces+4] = elements[[1,4,6,3], ei];
            face2v[:,Nfaces+5] = elements[[2,3,6,5], ei];
            face2e[1,Nfaces+1] = ei; #
            face2e[1,Nfaces+2] = ei; #
            face2e[1,Nfaces+3] = ei; #
            face2e[1,Nfaces+4] = ei; #
            face2e[1,Nfaces+5] = ei; #
            e2face[1,ei] = Nfaces+1;
            e2face[2,ei] = Nfaces+2;
            e2face[3,ei] = Nfaces+3;
            e2face[4,ei] = Nfaces+4;
            e2face[5,ei] = Nfaces+5;
            Nfaces += 5;
        elseif etypes[ei] == 7 # 5-pyramid
            face2v[:,Nfaces+1] = elements[[1,4,3,2], ei];
            face2v[1:3,Nfaces+2] = elements[[1,2,5], ei];
            face2v[1:3,Nfaces+3] = elements[[3,4,5], ei];
            face2v[1:3,Nfaces+4] = elements[[2,3,5], ei];
            face2v[1:3,Nfaces+5] = elements[[4,1,5], ei];
            face2e[1,Nfaces+1] = ei; #
            face2e[1,Nfaces+2] = ei; #
            face2e[1,Nfaces+3] = ei; #
            face2e[1,Nfaces+4] = ei; #
            face2e[1,Nfaces+5] = ei; #
            e2face[1,ei] = Nfaces+1;
            e2face[2,ei] = Nfaces+2;
            e2face[3,ei] = Nfaces+3;
            e2face[4,ei] = Nfaces+4;
            e2face[5,ei] = Nfaces+5;
            Nfaces += 5;
        end
    end
    
    # println("before:");
    # println(face2e)
    # println(e2face)
    
    # Remove duplicates
    # This is slow. Checks each face against other faces
    newface2v = zeros(Int, Nfp, Nfaces);
    newface2e = zeros(Int, 2, Nfaces);
    newe2face = zeros(Int, NfacesPerElement, nel);
    newface2v[:,1] = face2v[:,1];
    newface2e[:,1] = face2e[:,1];
    for fk=1:length(e2face[:,face2e[1,1]])
        if e2face[fk,face2e[1,1]] == 1
            newe2face[fk,face2e[1,1]] = 1;
        end
    end
    
    next_ind = 2;
    remove_count = 0;
    found = false;
    both_sides_done = zeros(Bool, Nfaces); # If both sides have been handled, don't need to check anymore.
    for fi=2:Nfaces
        found = false;
        removeinds = face2v[:,fi];
        for fj=1:next_ind-1
            if !both_sides_done[fj]
                keepinds = newface2v[:,fj];
                if shared_face(keepinds, removeinds) # fi is a duplicate. 
                    #newface2e[1,fj] = face2e[1,fj];
                    newface2e[2,fj] = face2e[1,fi];
                    for fk=1:length(e2face[:,face2e[1,fi]])
                        # if e2face[fk,newface2e[1,fj]] == fj
                        #     newe2face[fk,newface2e[1,fj]] = fj;
                        # end
                        if e2face[fk,face2e[1,fi]] == fi
                            newe2face[fk,face2e[1,fi]] = fj;
                        end
                    end
                    
                    remove_count += 1;
                    found = true;
                    both_sides_done[fj] = true;
                    break;
                end
            end
        end
        if !found # fi wasn't a duplicate. give it a new index.
            newface2e[1,next_ind] = face2e[1,fi];
            newface2e[2,next_ind] = 0;
            newface2v[:,next_ind] = face2v[:,fi];
            for fk=1:length(e2face[:,face2e[1,fi]])
                if e2face[fk,face2e[1,fi]] == fi
                    newe2face[fk,face2e[1,fi]] = next_ind;
                end
            end
            
            next_ind += 1;
        end
    end
    
    # cut off the excess
    remaining = Nfaces-remove_count;
    face2e = newface2e[:,1:remaining];
    face2v = newface2v[:,1:remaining];
    e2face = newe2face;
    
    # println("after:");
    # println(face2e)
    # println(e2face)
    
    return (face2v, face2e, e2face);
end

# Compute normal vectors for each face
# Note: this does not account for sign! 
function find_normals(face2v, x)
    dim = size(x,1);
    nfaces = size(face2v,2);
    
    # Special case for dim==1
    if dim == 1
        return ones(1,nfaces);
    end
    
    normals = zeros(dim, nfaces);
    for fi=1:nfaces
        if dim == 2 # faces are lines
            normals[:,fi] = normal2(x[:,face2v[1,fi]], x[:,face2v[2,fi]]);
        elseif dim == 3
            normals[:,fi] = normal3(x[:,face2v[1,fi]], x[:,face2v[2,fi]], x[:,face2v[3,fi]]);
        end
    end
    
    return normals;
end

function normal2(a, b)
    nx = b[2]-a[2];
    ny = a[1]-b[1];
    d = sqrt(nx*nx + ny*ny);
    return [nx/d; ny/d];
end

function normal3(a, b, c)
    v = b.-a;
    w = c.-a;
    nx = v[2]*w[3]-v[3]*w[2];
    ny = v[3]*w[1]-v[1]*w[3];
    nz = v[1]*w[2]-v[3]*w[1];
    d = sqrt(nx*nx + ny*ny + nz*nz);
    return [nx/d; ny/d; nz/d];
end

# Finds faces that have face2e[2,fi] == 0
# Just sets those BIDs to 1
function find_boundaries(face2e)
    nfaces = size(face2e,2);
    bdry = zeros(Int, nfaces);
    for fi=1:nfaces
        if face2e[2,fi] == 0
            bdry[fi] = 1;
        end
    end
    
    return bdry;
end

function shared_face(f1, f2)
    h = 0;
    s = 0;
    for i=1:length(f1)
        if f1[i] > 0
            h += 1;
            for j=1:length(f2)
                if f1[i] == f2[j]
                    s += 1;
                end
            end
        end
    end
    
    return s == h;
end