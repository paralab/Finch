#=
#
# Reads a .msh file and builds a MeshData struct
=#
export read_mesh

# numbers of CORNER nodes for first and second order elements as defined by GMSH
# TODO add higher order types
etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5]; # number of vertices for each type
etypetodim= [1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3]; # dimension of each type

# Reads from the file stream
# Returns a MeshData struct
function read_mesh(file)
    # Determine the MSH format version
    msh_version = 2;
    while !eof(file)
        line = readline(file);
        if occursin("\$MeshFormat", line)
            # Check if the MSH version is old(2) or new(4)
            vals = split(readline(file), " ", keepempty=false);
            if parse(Float64, vals[1]) >= 4
                msh_version = 4;
            end
            break;
        end
    end
    
    # use the appropriate reader
    seekstart(file);
    if msh_version == 2
        return read_msh_v2(file);
    elseif msh_version == 4
        return read_msh_v4(file);
    end
end

function read_msh_v2(file)
    # Find the beginning of the Nodes and Elements sections.
    # Then read in the numbers.
    nodes_done = false;
    elements_done = false;
    msh_version = 2;
    nx = 0;
    nel = 0;
    nodes = [];
    indices = [];
    elements = [];
    etypes = [];
    nv = [];
    while((!nodes_done || !elements_done) && !eof(file))
        line = readline(file);
        if occursin("\$Nodes", line)
            # The Nodes section
            line = readline(file);
            nx = parse(Int, split(line, " ", keepempty=false)[1]);
            if nx > 0
                # parse each node's coordinates
                nodes = zeros(3, nx);
                indices = zeros(Int, nx);
                i = 1;
                line = readline(file);
                while !occursin("\$EndNodes", line) && !eof(file)
                    vals = split(line, " ", keepempty=false);
                    indices[i] = parse(Int, vals[1]);
                    nodes[1,i] = parse(Float64, vals[2]);
                    nodes[2,i] = parse(Float64, vals[3]);
                    nodes[3,i] = parse(Float64, vals[4]);
                    i += 1;
                    line = readline(file);
                end
                
                nodes_done = true;
            end
        elseif occursin("\$Elements", line)
            # The Elements section
            line = readline(file);
            nel = parse(Int, split(line, " ", keepempty=false)[1]);
            if nel > 0
                # parse each element's numbers
                nv = zeros(Int, nel);
                etypes = zeros(Int, nel);
                elements = zeros(Int, 8, nel);
                i = 1;
                line = readline(file);
                while !occursin("\$EndElements", line) && !eof(file)
                    vals = split(line, " ", keepempty=false);
                    # Skip lower dimensional elements
                    tmptype = parse(Int, vals[2]);
                    if etypetodim[tmptype] == config.dimension
                        etypes[i] = tmptype;
                        nv[i] = etypetonv[tmptype];
                        offset = parse(Int, vals[3]) + 3;
                        for j=1:nv[i]
                            elements[j,i] = parse(Int, vals[offset + j]);
                        end
                        i += 1;
                    end
                    
                    line = readline(file);
                end
                # adjust to correct number of elements
                nel = i-1;
                nv = nv[1:nel];
                etypes = etypes[1:nel];
                elements = elements[:, 1:nel];
                
                elements_done = true;
            end
        end
    end
    nodes = nodes[1:config.dimension, :];
    
    return MeshData(nx, nodes, indices, nel, elements, etypes, nv);
end

function read_msh_v4(file)
    # Find the beginning of the Nodes and Elements sections.
    # Then read in the numbers.
    nodes_done = false;
    elements_done = false;
    msh_version = 2;
    nx = 0;
    nel = 0;
    nodes = [];
    indices = [];
    elements = [];
    etypes = [];
    nv = [];
    while((!nodes_done || !elements_done) && !eof(file))
        line = readline(file);
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
            line = readline(file);
            nx = parse(Int, split(line, " ", keepempty=false)[2]);
            if nx > 0
                # parse node info
                nodes = zeros(3,nx);
                indices = zeros(Int, nx);
                i = 1;
                line = readline(file);
                while !occursin("\$EndNodes", line) && !eof(file)
                    # process entity blocks 
                    # The first line has entity info. 
                    entinfo = split(line, " ", keepempty=false);
                    entdim = parse(Int, entinfo[1]);
                    parametric = parse(Int, entinfo[3]);
                    entnx = parse(Int, entinfo[4]);
                    # Node indices
                    for ni=1:entnx
                        line = readline(file);
                        vals = split(line, " ", keepempty=false);
                        indices[i-1 + ni] = parse(Int, vals[1]);
                    end
                    # Node coordinates
                    for ni=1:entnx
                        line = readline(file);
                        vals = split(line, " ", keepempty=false);
                        nodes[1,i-1 + ni] = parse(Float64, vals[1]);
                        nodes[2,i-1 + ni] = parse(Float64, vals[2]);
                        nodes[3,i-1 + ni] = parse(Float64, vals[3]);
                    end
                    # Parametric (ignore)
                    if parametric > 0
                        printerr("Paremetric entity found in mesh file. Ignoring parameters.")
                        for ni=1:entdim
                            line = readline(file);
                        end
                    end
                    i += entnx;
                    line = readline(file);
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
            line = readline(file);
            nel = parse(Int, split(line, " ", keepempty=false)[2]);
            if nel > 0
                # parse elements
                nv = zeros(Int, nel);
                etypes = zeros(Int, nel);
                elements = zeros(Int, 8, nel);
                i = 1;
                line = readline(file);
                while !occursin("\$EndElements", line) && !eof(file)
                    # Process entity blocks
                    # The first line has entity info. 
                    vals = split(line, " ", keepempty=false);
                    entdim = parse(Int, vals[3]);
                    enttype = parse(Int, vals[3]);
                    entnel = parse(Int, vals[4]);
                    
                    # Ignore elements with lower dimension (faces, points etc.)
                    if entdim == config.dimension
                        # read element info
                        for ei=1:entnel
                            line = readline(file);
                            vals = split(line, " ", keepempty=false);
                            etypes[i-1 + ei] = enttype;
                            entnv = etypetonv[enttype];
                            nv[i-1 + ei] = entnv;
                            for j=1:entnv
                                elements[j, i-1 + ei] = parse(Int, vals[1 + j]);
                            end
                        end
                        
                        i += entnel;
                        line = readline(file);
                    else
                        for ei=1:entnel
                            line = readline(file); # skip lines
                        end
                        line = readline(file);
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
    nodes = nodes[1:config.dimension, :];
    
    return MeshData(nx, nodes, indices, nel, elements, etypes, nv);
end
