#=
#
# Writes a .msh file from a MeshData struct
=#
export write_mesh
export write_mesh_graph

function write_mesh(file, format, mesh::MeshData)
    # don't proceed unless the mesh contains things
    if mesh.nx == 0 || mesh.nel == 0
        println("Tried to write a MSH file with an empty mesh.");
        return;
    end
    # select your file format
    if format == MSH_V2
        write_mesh_MSHv2(file, mesh);
    elseif format == MSH_V4
        write_mesh_MSHv4(file, mesh);
    elseif format == "medit"
        write_mesh_medit(file, mesh);
    end
end

# Writes to the file stream using the mesh
function write_mesh_MSHv2(file, mesh::MeshData)
    # Write the format info
    println(file, "\$MeshFormat");
    println(file, "2.2 0 8");
    println(file, "\$EndMeshFormat");
    
    # Write the nodes
    println(file, "\$Nodes");
    println(file, string(mesh.nx));
    tmp = [0.0,0.0,0.0];
    for i=1:mesh.nx
        if finch_state.config.dimension == 1
            tmp[1] = mesh.nodes[1,i];
        elseif finch_state.config.dimension == 2
            tmp[1:2] = mesh.nodes[:,i];
        else
            tmp = mesh.nodes[:,i];
        end
        println(file, string(mesh.indices[i])*" "*string(tmp[1])*" "*string(tmp[2])*" "*string(tmp[3]));
    end
    println(file, "\$EndNodes");
    
    # Print the elements
    println(file, "\$Elements");
    println(file, string(mesh.nel));
    for i=1:mesh.nel
        line = string(i)*" "*string(mesh.etypes[i])*" 0 ";
        for j=1:mesh.nv[i]
            line = line*" "*string(mesh.elements[j,i]);
        end
        println(file, line);
    end
    println(file, "\$EndElements");
end

# Writes to the file stream using the mesh
function write_mesh_MSHv4(file, mesh::MeshData)
    # Write the format info
    println(file, "\$MeshFormat");
    println(file, "4.1 0 8");
    println(file, "\$EndMeshFormat");
    
    # Write the nodes
    minn = mesh.indices[1];
    maxn = minn;
    for i=2:mesh.nx
        minn = min(minn,mesh.indices[i]);
        maxn = max(maxn,mesh.indices[i]);
    end
    dim = 1;
    if mesh.etypes[1] > 1
        dim = 2;
    end
    if mesh.etypes[1] > 3
        dim = 3;
    end
    println(file, "\$Nodes");
    println(file, "1 "*string(mesh.nx)*" "*string(minn)*" "*string(maxn));
    println(file, string(dim)*" 1 0 "*string(mesh.nx));
    for i=1:mesh.nx
        println(file, string(mesh.indices[i]));
    end
    for i=1:mesh.nx
        println(file, string(mesh.nodes[i,1])*" "*string(mesh.nodes[i,2])*" "*string(mesh.nodes[i,3]))
    end
    println(file, "\$EndNodes");
    
    # Print the elements
    println(file, "\$Elements");
    println(file, "1 "*string(mesh.nel)*" 1 "*string(mesh.nel));
    println(file, string(dim)*" 1 "*string(mesh.etypes[1])*" "*string(mesh.nel));
    for i=1:mesh.nel
        line = string(i);
        for j=1:mesh.nv[1]
            line = line*" "*string(mesh.elements[i,j]);
        end
        println(file, line);
    end
    println(file, "\$EndElements");
end

function write_mesh_medit(file, mesh)
    dim = 1;
    if mesh.etypes[1] > 1
        dim = 2;
    end
    if mesh.etypes[1] > 3
        dim = 3;
    end
    
    println(file, "MeshVersionFormatted 1");
    println(file, "# Created by Finch.jl");
    println(file, "");
    println(file, "Dimension");
    println(file, string(dim));
    
    println(file, "Vertices");
    println(file, string(mesh.nx));
    for i=1:mesh.nx
        for j=1:dim
            print(file, string(mesh.nodes[j,i]) * "    ");
        end
        println(file, string(mesh.indices[i]));
    end
    
    # Assume all elements are same type
    if mesh.etypes[1] == 1 # line
    elseif mesh.etypes[1] == 2 # tri
        println(file, "Triangles");
    elseif mesh.etypes[1] == 3 # quad
        println(file, "Quadrilaterals");
    elseif mesh.etypes[1] == 4 # tet
        println(file, "Tetrahedra");
    elseif mesh.etypes[1] == 5 # hex
        println(file, "Hexahedra");
    end
    println(file, string(mesh.nel));
    for i=1:mesh.nel
        for j=1:mesh.nv[1]
            print(file, string(mesh.elements[j,i]) * "    ");
        end
        println(file, "0");
    end
    println(file, "End");
end

# Output the connectivity graph info for a mesh.
# An edge exists between two elements that share a node.
# The weight is how many shared nodes.
# The format is 
# - number of dimensions (2 or 3)
# - number of elements (integer)
# - number of edges (integer)
# - number of boundary elements (integer)
# - element center coordinates (two or three floats each)
# - edge pairs and weight (three integers each)
# - boundary element index and BID (two integers)
function write_mesh_graph(file, mesh::MeshData)
    dim = size(mesh.nodes,1);
    num_elements = mesh.nel;
    total_edges = 0; # to be found
    el_centers = zeros(dim, num_elements);
    faces_per_el = size(mesh.element2face, 1);
    
    # temporary storage
    el_connections = zeros(Int, num_elements * mesh.nv[1]);
    el_edges = zeros(Int, 2, num_elements);
    center = zeros(dim);
    lines = fill("",0);
    bdryLines = fill("",0);
    
    print("mesh to graph progress %0");
    progress_interval = num_elements/20;
    next_progress = progress_interval;
    next_percent = 5;
    dot_interval = num_elements/100;
    next_dot = dot_interval;
    
    for ei=1:num_elements
        # Find center of mass
        center .= 0.0;
        for ni=1:mesh.nv[ei]
            nid = mesh.elements[ni,ei];
            center .+= mesh.nodes[:,nid];
        end
        el_centers[:,ei] .= center ./ mesh.nv[ei];
        
        # find all connections with other elements
        # only check ej<ei to avoid double counting
        num_connections = 0;
        for ej=1:ei-1
            for ni=1:mesh.nv[ei]
                nid = mesh.elements[ni,ei];
                for nj=1:mesh.nv[ej]
                    if nid == mesh.elements[nj,ej]
                        num_connections += 1;
                        el_connections[num_connections] = ej;
                    end
                end
            end
        end
        
        # make pairs and weights
        num_edges = 0;
        for ej=1:num_connections
            found=false;
            for edgei=1:num_edges
                if el_connections[ej] == el_edges[1,edgei]
                    found = true;
                    el_edges[2,edgei] += 1;
                    break;
                end
            end
            if !found
                num_edges += 1;
                el_edges[1, num_edges] = el_connections[ej];
                el_edges[2, num_edges] = 1;
            end
        end
        
        # Write this element's edges to the lines for the file
        tmp_lines = "";
        for edgei=1:num_edges
            tmp_lines *= string(ei) * " " * string(el_edges[1,edgei]) * " " * string(el_edges[2,edgei]) * "\n";
        end
        push!(lines, tmp_lines);
        total_edges += num_edges;
        
        # If this is a boundary element, write its index and BID to bdryLines
        for fi=1:faces_per_el
            fid = mesh.element2face[fi, ei];
            if mesh.bdryID[fid] > 0
                push!(bdryLines, string(ei) * " " * string(mesh.bdryID[fid]) * "\n");
                break;
            end
        end
        
        if ei > next_progress
            print(string(next_percent))
            next_percent = next_percent + 5;
            next_progress += progress_interval;
            
        elseif ei > next_dot
            print(".")
            next_dot += dot_interval;
        end
        
    end
    
    # Now all the numbers have been collected.
    # Print to the file.
    println(file, string(dim));
    println(file, string(num_elements));
    println(file, string(total_edges));
    println(file, string(length(bdryLines)));
    # element center coords
    for ei=1:num_elements
        if dim==2
            println(file, string(el_centers[1,ei]) * " " * string(el_centers[2,ei]));
        else
            println(file, string(el_centers[1,ei]) * " " * string(el_centers[2,ei]) * " " * string(el_centers[3,ei]));
        end
    end
    # graph edges
    for i=1:length(lines)
        print(file, lines[i]);
    end
    # bdry elements and ID
    for i=1:length(bdryLines)
        print(file, bdryLines[i]);
    end
end

