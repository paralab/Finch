#=
#
# Writes a .msh file from a MeshData struct
=#
export write_mesh

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