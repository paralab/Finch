
#=
Makes a graph representation of a mesh for the purpose of partitioning.
Returns the graph::Vector{Vector{Int}} that has lists of connected elements 
for each element.
Also returns the center locations for each element.
=#
function mesh_to_graph(mesh::MeshData)
    dim = size(mesh.nodes,1);
    num_elements = mesh.nel;
    faces_per_el = size(mesh.element2face, 1);
    
    el_centers = zeros(dim, num_elements);
    face_neighbors = fill(zeros(Int,faces_per_el), num_elements);
    
    center = zeros(dim);
    
    for ei=1:num_elements
        # Find center of mass
        center .= 0.0;
        for ni=1:mesh.nv[ei]
            nid = mesh.elements[ni,ei];
            center .+= mesh.nodes[:,nid];
        end
        el_centers[:,ei] .= center ./ mesh.nv[ei];
        
        # Find all face neighbors
        nneighbors = 0;
        for fi=1:faces_per_el
            fid = mesh.element2face[fi,ei];
            if mesh.face2element[1,fid] == ei
                neighbor = mesh.face2element[2,fid]
            else
                neighbor = mesh.face2element[1,fid]
            end
            
            if neighbor > 0
                nneighbors += 1;
                face_neighbors[ei][nneighbors] = neighbor;
            end
        end
        
        # Adjust for boundaries
        if nneighbors < faces_per_el
            face_neighbors[ei] = face_neighbors[ei][1:nneighbors];
        end
    end
    
    return (face_neighbors, el_centers);
end