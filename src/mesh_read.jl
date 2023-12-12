#=
#
# Reads a mesh file and builds a MeshData struct
=#
export read_mesh

# Reads from the file stream
# Returns a MeshData struct
function read_mesh(file)
    # Determine the MSH format version
    mesh_file_type = 0; # 1=msh_v2, 2=msh_v4, 3=medit
    line_number = 0;
    while !eof(file)
        line = readline(file); line_number += 1;
        if occursin("\$MeshFormat", line)
            # Check if the MSH version is old(2) or new(4)
            vals = split(readline(file), " ", keepempty=false);
            if parse_check(Float64, vals[1], line_number) >= 4
                mesh_file_type = 2;
            else
                mesh_file_type = 1;
            end
            break;
        elseif occursin("MeshVersionFormatted", line)
            mesh_file_type = 3;
        end
    end
    
    # use the appropriate reader
    seekstart(file);
    if mesh_file_type == 1
        return read_msh_v2(file);
    elseif mesh_file_type == 2
        return read_msh_v4(file);
    elseif mesh_file_type == 3
        return read_medit(file);
    else
        printerr("Couldn't recognize mesh file type. Use GMSH .msh version 2 or 4, or Medit .mesh")
    end
end

function read_msh_v2(file)
    # numbers of CORNER nodes for first and second order elements as defined by GMSH
    etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5]; # number of vertices for each type
    etypetodim= [1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3]; # dimension of each type
    
    # Find the beginning of the Nodes and Elements sections.
    # Then read in the numbers.
    line_number = 0;
    nodes_done = false;
    elements_done = false;
    nx = 0;
    nel = 0;
    nodes = zeros(0,0);
    indices = zeros(Int,0);
    elements = zeros(Int,0,0);
    etypes = zeros(Int,0);
    nv = zeros(Int,0);
    
    while((!nodes_done || !elements_done) && !eof(file))
        line = readline(file); line_number += 1;
        if occursin("\$Nodes", line)
            # The Nodes section
            line = readline(file); line_number += 1;
            nx = parse_check(Int, split(line, " ", keepempty=false)[1], line_number);
            if nx > 0
                # parse each node's coordinates
                nodes = zeros(3, nx);
                indices = zeros(Int, nx);
                i = 1;
                line = readline(file); line_number += 1;
                while !occursin("\$EndNodes", line) && !eof(file)
                    vals = split(line, " ", keepempty=false);
                    indices[i] = parse_check(Int, vals[1], line_number);
                    nodes[1,i] = parse_check(Float64, vals[2], line_number);
                    nodes[2,i] = parse_check(Float64, vals[3], line_number);
                    nodes[3,i] = parse_check(Float64, vals[4], line_number);
                    i += 1;
                    line = readline(file); line_number += 1;
                end
                
                nodes_done = true;
            end
        elseif occursin("\$Elements", line)
            # The Elements section
            line = readline(file); line_number += 1;
            nel = parse_check(Int, split(line, " ", keepempty=false)[1], line_number);
            if nel > 0
                # parse each element's numbers
                nv = zeros(Int, nel);
                etypes = zeros(Int, nel);
                elements = zeros(Int, 8, nel);
                i = 1;
                line = readline(file); line_number += 1;
                while !occursin("\$EndElements", line) && !eof(file)
                    vals = split(line, " ", keepempty=false);
                    # Skip lower dimensional elements
                    tmptype = parse_check(Int, vals[2], line_number);
                    if etypetodim[tmptype] == finch_state.config.dimension
                        etypes[i] = tmptype;
                        nv[i] = etypetonv[tmptype];
                        offset = parse_check(Int, vals[3], line_number) + 3;
                        for j = 1:nv[i]
                            elements[j, i] = parse_check(Int, vals[offset + j], line_number);
                        end
                        i += 1;
                    end
                    
                    line = readline(file); line_number += 1;
                end
                # adjust to correct number of elements
                nel = i - 1;
                nv = nv[1 : nel];
                etypes = etypes[1 : nel];
                elements = elements[:, 1:nel];
                
                elements_done = true;
            end
        end
    end
    nodes = nodes[1:finch_state.config.dimension, :];
    
    return MeshData(nx, nodes, indices, nel, elements, etypes, nv);
end

function read_msh_v4(file)
    # numbers of CORNER nodes for first and second order elements as defined by GMSH
    etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5]; # number of vertices for each type
    etypetodim= [1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3]; # dimension of each type
    
    # Find the beginning of the Nodes and Elements sections.
    # Then read in the numbers.
    line_number = 0;
    nodes_done = false;
    elements_done = false;
    nx = 0;
    nel = 0;
    nodes = [];
    indices = [];
    elements = [];
    etypes = [];
    nv = [];
    while((!nodes_done || !elements_done) && !eof(file))
        line = readline(file); line_number += 1;
        if occursin("\$Nodes", line)
            # The Nodes section
            #=
            $Nodes
            numEntityBlocks(size_t) numNodes(size_t) minNodeTag(size_t) maxNodeTag(size_t)
            entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodesInBlock(size_t)
                nodeTag(size_t)
                ...
                x(double) y(double) z(double)
                < u(double; if parametric and entityDim >= 1) >
                < v(double; if parametric and entityDim >= 2) >
                < w(double; if parametric and entityDim == 3) >
                ...
            ...
            $EndNodes
            =#
            line = readline(file); line_number += 1;
            nx = parse_check(Int, split(line, " ", keepempty=false)[2], line_number);
            if nx > 0
                # parse node info
                nodes = zeros(3,nx);
                indices = zeros(Int, nx);
                i = 1;
                line = readline(file); line_number += 1;
                while !occursin("\$EndNodes", line) && !eof(file)
                    # process entity blocks 
                    # The first line has entity info. 
                    entinfo = split(line, " ", keepempty=false);
                    entdim = parse_check(Int, entinfo[1], line_number);
                    parametric = parse_check(Int, entinfo[3], line_number);
                    entnx = parse_check(Int, entinfo[4], line_number);
                    # Node indices
                    for ni=1:entnx
                        line = readline(file); line_number += 1;
                        vals = split(line, " ", keepempty=false);
                        indices[i-1 + ni] = parse_check(Int, vals[1], line_number);
                    end
                    # Node coordinates
                    for ni=1:entnx
                        line = readline(file); line_number += 1;
                        vals = split(line, " ", keepempty=false);
                        nodes[1,i-1 + ni] = parse_check(Float64, vals[1], line_number);
                        nodes[2,i-1 + ni] = parse_check(Float64, vals[2], line_number);
                        nodes[3,i-1 + ni] = parse_check(Float64, vals[3], line_number);
                    end
                    # Parametric (ignore)
                    if parametric > 0
                        printerr("Paremetric entity found in mesh file. Ignoring parameters.")
                        for ni=1:entdim
                            line = readline(file); line_number += 1;
                        end
                    end
                    i += entnx;
                    line = readline(file); line_number += 1;
                end
                
                nodes_done = true;
            end
        elseif occursin("\$Elements", line)
            # The Elements section
            #=
            $Elements
            1 2 1 2          1 entity bloc, 2 elements total, min/max element tags: 1 and 2
            2 1 3 2          2D entity (surface) 1, element type 3 (4-node quad), 2 elements
            1 1 2 3 4          quad tag #1, nodes 1 2 3 4
            2 2 5 6 3          quad tag #2, nodes 2 5 6 3
            $EndElements
            =#
            line = readline(file); line_number += 1;
            nel = parse_check(Int, split(line, " ", keepempty=false)[2], line_number);
            if nel > 0
                # parse elements
                nv = zeros(Int, nel);
                etypes = zeros(Int, nel);
                elements = zeros(Int, 8, nel);
                i = 1;
                line = readline(file); line_number += 1;
                while !occursin("\$EndElements", line) && !eof(file)
                    # Process entity blocks
                    # The first line has entity info. 
                    vals = split(line, " ", keepempty=false);
                    entdim = parse_check(Int, vals[1], line_number);
                    enttype = parse_check(Int, vals[3], line_number);
                    entnel = parse_check(Int, vals[4], line_number);
                    
                    # Ignore elements with lower dimension (faces, points etc.)
                    if entdim == finch_state.config.dimension
                        # read element info
                        for ei=1:entnel
                            line = readline(file); line_number += 1;
                            vals = split(line, " ", keepempty=false);
                            etypes[i-1 + ei] = enttype;
                            entnv = etypetonv[enttype];
                            nv[i-1 + ei] = entnv;
                            for j=1:entnv
                                elements[j, i-1 + ei] = parse_check(Int, vals[1 + j], line_number);
                            end
                        end
                        
                        i += entnel;
                        line = readline(file); line_number += 1;
                    else
                        for ei=1:entnel
                            line = readline(file); line_number += 1; # skip lines
                        end
                        line = readline(file); line_number += 1;
                    end
                    
                end
                
                # adjust nel to remove unwanted elements
                nel = i-1;
                elements = elements[:,1:nel];
                etypes = etypes[1:nel];
                nv = nv[1:nel];
                
                elements_done = true;
            end
        end
    end
    nodes = nodes[1:finch_state.config.dimension, :];
    
    return MeshData(nx, nodes, indices, nel, elements, etypes, nv);
end

function read_medit(file)
    # Whitespace is inconsistent, but essentially it will look like
    #
    # MeshVersionFormatted 1
    # Dimension 3
    # Vertices 12345
    # 0.0509745 -0.268407 -0.198911 0      <-- Note: always 4 values (4th is a tag?)
    # ...
    # Tetrahedra 1234
    # 1 2 3 4 0     <-- again an extra tag value
    # ...
    # End
    #
    # There may be lower dimensional elements for surfaces, but it's inconsistent
    line_number = 0;
    nodes_done = false;
    elements_done = false;
    nx = 0;
    nel = 0;
    nodes = [];
    indices = [];
    elements = [];
    etypes = [];
    eflags = [];
    nv = [];
    dim = 1;
    el_type_count = 0;
    el_offset = 0;
    while(!eof(file))
        line = readline(file); line_number += 1;
        if occursin("Dimension", line)
            tokens = split(line, " ", keepempty=false);
            if length(tokens) > 1
                dim = parse_check(Int, tokens[2], line_number);
            else
                line = readline(file); line_number += 1;
                dim = parse_check(Int, line, line_number);
            end
            if finch_state.config.dimension != dim
                # Maybe someone forgot to set dimension.
                # Maybe the mesh file was poorly created.
                # Let's find the highest dmension elements.
                bookmark = position(file);
                highest_dim = 1;
                while !eof(file)
                    line = readline(file); line_number += 1;
                    if (occursin("Triangles", line) || occursin("Quadrilaterals", line))
                        highest_dim = max(highest_dim,2);
                    end
                    if (occursin("Tetrahedra", line) || occursin("Hexahedra", line))
                        highest_dim = max(highest_dim,3);
                    end
                end
                seek(file, bookmark); # rewind
                
                if finch_state.config.dimension != highest_dim
                    printerr("Dimension in mesh file doesn't match config. Updating config to dimension = " * string(highest_dim));
                    finch_state.config.dimension = highest_dim;
                else
                    # It was a poorly made file.
                end
                dim = highest_dim;
            end
            
        elseif occursin("Vertices", line)
            tokens = split(line, " ", keepempty=false);
            if length(tokens) > 1
                nx = parse_check(Int, tokens[2], line_number);
            else
                line = readline(file); line_number += 1;
                nx = parse_check(Int, line, line_number);
            end
            
            # parse vertex info
            nodes = zeros(3,nx);
            flags = zeros(Int, nx);
            for i=1:nx
                line = readline(file); line_number += 1;
                vals = split(line, " ", keepempty=false);
                for j=1:dim
                    nodes[j,i] = parse_check(Float64, vals[j], line_number);
                end
                flags[i] = parse_check(Int, vals[end], line_number);
            end
            
            nodes_done = true;
            
        elseif dim == 1 && (occursin("Edges", line))
            tokens = split(line, " ", keepempty=false);
            if length(tokens) > 1
                nel = parse_check(Int, tokens[2], line_number);
            else
                line = readline(file); line_number += 1;
                nel = parse_check(Int, line, line_number);
            end
            
            nv = fill(2, nel);
            etypes = fill(1, nel);
            elements = zeros(Int, 8, nel);
            eflags = zeros(Int, nel);
            for i=1:nel
                line = readline(file); line_number += 1;
                vals = split(line, " ", keepempty=false);
                elements[1,i] = parse_check(Int, vals[1], line_number);
                elements[2,i] = parse_check(Int, vals[2], line_number);
                eflags[i] = parse_check(Int, vals[4], line_number);
            end
            elements_done = true;
            
        elseif dim == 2 && (occursin("Triangles", line) || occursin("Quadrilaterals", line))
            el_type_count += 1;
            el_offset = nel;
            tokens = split(line, " ", keepempty=false);
            if length(tokens) > 1
                nel = parse_check(Int, tokens[2], line_number);
            else
                line = readline(file); line_number += 1;
                nel = parse_check(Int, line, line_number);
            end
            if tokens[1] == "Triangles"
                if el_type_count > 1
                    append!(nv, fill(3, nel));
                    append!(etypes, fill(2, nel));
                    elements = hcat(elements, zeros(Int, 8, nel));
                    append!(eflags, zeros(Int, nel));
                else
                    nv = fill(3, nel);
                    etypes = fill(2, nel);
                    elements = zeros(Int, 8, nel);
                    eflags = zeros(Int, nel);
                end
                
                for i=(el_offset + 1):(el_offset + nel)
                    line = readline(file); line_number += 1;
                    vals = split(line, " ", keepempty=false);
                    elements[1,i] = parse_check(Int, vals[1], line_number);
                    elements[2,i] = parse_check(Int, vals[2], line_number);
                    elements[3,i] = parse_check(Int, vals[3], line_number);
                    eflags[i] = parse_check(Int, vals[4], line_number);
                end
                
            elseif tokens[1] == "Quadrilaterals"
                if el_type_count > 1
                    append!(nv, fill(4, nel));
                    append!(etypes, fill(3, nel));
                    elements = hcat(elements, zeros(Int, 8, nel));
                    append!(eflags, zeros(Int, nel));
                else
                    nv = fill(4, nel);
                    etypes = fill(3, nel);
                    elements = zeros(Int, 8, nel);
                    eflags = zeros(Int, nel);
                end
                for i=(el_offset + 1):(el_offset + nel)
                    line = readline(file); line_number += 1;
                    vals = split(line, " ", keepempty=false);
                    elements[1,i] = parse_check(Int, vals[1], line_number);
                    elements[2,i] = parse_check(Int, vals[2], line_number);
                    elements[3,i] = parse_check(Int, vals[3], line_number);
                    elements[4,i] = parse_check(Int, vals[4], line_number);
                    eflags[i] = parse_check(Int, vals[5], line_number);
                end
            end
            nel = nel + el_offset;
            elements_done = true;
            
        elseif dim == 3 && (occursin("Tetrahedra", line) || occursin("Hexahedra", line))
            el_type_count += 1;
            el_offset = nel;
            tokens = split(line, " ", keepempty=false);
            if length(tokens) > 1
                nel = parse_check(Int, tokens[2], line_number);
            else
                line = readline(file); line_number += 1;
                nel = parse_check(Int, line, line_number);
            end
            
            if tokens[1] == "Tetrahedra"
                if el_type_count > 1
                    append!(nv, fill(4, nel));
                    append!(etypes, fill(4, nel));
                    elements = hcat(elements, zeros(Int, 8, nel));
                    append!(eflags, zeros(Int, nel));
                else
                    nv = fill(4, nel);
                    etypes = fill(4, nel);
                    elements = zeros(Int, 8, nel);
                    eflags = zeros(Int, nel);
                end
                for i=(el_offset + 1):(el_offset + nel)
                    line = readline(file); line_number += 1;
                    vals = split(line, " ", keepempty=false);
                    elements[1,i] = parse_check(Int, vals[1], line_number);
                    elements[2,i] = parse_check(Int, vals[2], line_number);
                    elements[3,i] = parse_check(Int, vals[3], line_number);
                    elements[4,i] = parse_check(Int, vals[4], line_number);
                    eflags[i] = parse_check(Int, vals[5], line_number);
                end
                
            elseif tokens[1] == "Hexahedra"
                if el_type_count > 1
                    append!(nv, fill(8, nel));
                    append!(etypes, fill(5, nel));
                    elements = hcat(elements, zeros(Int, 8, nel));
                    append!(eflags, zeros(Int, nel));
                else
                    nv = fill(8, nel);
                    etypes = fill(5, nel);
                    elements = zeros(Int, 8, nel);
                    eflags = zeros(Int, nel);
                end
                for i=(el_offset + 1):(el_offset + nel)
                    line = readline(file); line_number += 1;
                    vals = split(line, " ", keepempty=false);
                    elements[1,i] = parse_check(Int, vals[1], line_number);
                    elements[2,i] = parse_check(Int, vals[2], line_number);
                    elements[3,i] = parse_check(Int, vals[3], line_number);
                    elements[4,i] = parse_check(Int, vals[4], line_number);
                    elements[5,i] = parse_check(Int, vals[5], line_number);
                    elements[6,i] = parse_check(Int, vals[6], line_number);
                    elements[7,i] = parse_check(Int, vals[7], line_number);
                    elements[8,i] = parse_check(Int, vals[8], line_number);
                    eflags[i] = parse_check(Int, vals[9], line_number);
                end
            end
            nel = nel + el_offset;
            elements_done = true;
        end
    end
    if !nodes_done || !elements_done
        printerr("Did not successfully read the mesh file.")
        return nothing;
    end
    nodes = nodes[1:finch_state.config.dimension, :];
    indices = Array(1:nx);
    
    return MeshData(nx, nodes, indices, nel, elements, etypes, nv);
end

# Check input numbers for NaN or parse errors.
# If input is bad, print error and exit.
# Input: the type to parse to and the string to parse
# Returns: the parsed value
function parse_check(type, string, line_number)
    try
        result = parse(type, string);
        if isnan(result) || isinf(result)
            printerr(finch_state, "Unexpected value in mesh file on line $(line_number): $string", fatal = true);
        end
        
        return result;
        
    catch
        printerr(finch_state, "Error parsing mesh file on line $(line_number): $string should be of type $type", fatal = true);
    end
end