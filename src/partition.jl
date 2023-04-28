#=
Partition info holds the mesh partition info as well as the indexed variable partition info.
The mesh portion is written at mesh creation.
The index portion is written at assembly loop generation.
=#
mutable struct Partition_Info
    total_partitions::Int # Number of all partitions which ideally is the number of MPI processes
    
    mesh_partitions::Int # Number of mesh partitions which must divide total_partitions
    this_mesh_partition::Int # Index of this proc's mesh partition
    mesh_comm::MPI.Comm # The MPI communicator for this mesh partition
    
    index_partitions::Array{Int,1} # Number of partitions for each indexer
    this_index_partition::Array{Array{Int,1},1} # Range of indices for this proc for each indexer
    
    
    Partition_Info(total, mesh, this_mesh) = new(total, mesh, this_mesh, [], [])
end

# Chooses a partitioning of the mesh and indices such that their product
# is the total number of procs.
# This should really be done by hand, but this is an option.
# Ranges will be an array of integers that are the number of pieces to partition.
# They are listed in order of priority.
# For example, nel=512, ind1=32, ind2=8, num_procs=256
# -> [nel,ind1,ind2] will give [256,1,1]
# -> [ind1,nel,ind2] will give [32,8,1]
# When the numbers do not divide nicely, some procs may be idle, which is not good.
function get_partition_scheme(num_procs, ranges)
    n=length(ranges);
    result = ones(Int, n);
    remaining = num_procs;
    check = 1;
    for i=1:n
        if ranges[i] <= remaining
            result[i] = ranges[i];
            remaining = Int(floor(remaining/ranges[i]));
        elseif remaining > 0
            result[i] = remaining;
            remaining = 0;
        end
        check *= result[i];
    end
    if !(check == num_procs)
        printerr("Number of procs("*string(num_procs)*") could not be nicely divided among "*string(ranges));
        # return ones(Int, n);
    end
    
    return result;
end

#=
Partitions a mesh into np parts using the specified partitioning method.
=#
function get_element_partitions(mesh, np, partitioner)
    if partitioner == METIS
        return get_element_partitions_metis(mesh, np)
    elseif partitioner == FENNEL
        return get_element_partitions_fennel(mesh, np)
    else
        printerr("Unknown partitioning method: $partitioner. Using Metis as default.");
        return get_element_partitions_metis(mesh, np)
    end
end

#=
Uses LocalFennelPartitioning.jl
Returns a list of partition numbers for each element.
=#
function get_element_partitions_fennel(mesh, np)
    (graph, locations) = mesh_to_graph(mesh);
    
    (partitions, map) = local_fennel_sim(graph, locations, np)
    
    return map;
end

#=
Metis is used to partition the mesh.
This uses METIS_jll: libmetis

This is the primary function used to partition a mesh into np partitions.
It returns a list of partition numbers for each element.
=#
function get_element_partitions_metis(mesh, np)
    METIS_NOPTIONS = 40

    # codes returned by metis functions
    METIS_OK = Cint(1)
    METIS_ERROR_INPUT = Cint(-2)
    METIS_ERROR_MEMORY = Cint(-3)
    METIS_ERROR = Cint(-4)

    # location in options vector
    METIS_OPTION_PTYPE     = 1
    METIS_OPTION_OBJTYPE   = 2
    METIS_OPTION_CTYPE     = 3
    METIS_OPTION_IPTYPE    = 4
    METIS_OPTION_RTYPE     = 5
    METIS_OPTION_DBGLVL    = 6
    METIS_OPTION_NITER     = 7
    METIS_OPTION_NCUTS     = 8
    METIS_OPTION_SEED      = 9
    METIS_OPTION_NO2HOP    = 10
    METIS_OPTION_MINCONN   = 11
    METIS_OPTION_CONTIG    = 12
    METIS_OPTION_COMPRESS  = 13
    METIS_OPTION_CCORDER   = 14
    METIS_OPTION_PFACTOR   = 15
    METIS_OPTION_NSEPS     = 16
    METIS_OPTION_UFACTOR   = 17
    METIS_OPTION_NUMBERING = 18
    METIS_OPTION_HELP      = 19
    METIS_OPTION_TPWGTS    = 20
    METIS_OPTION_NCOMMON   = 21
    METIS_OPTION_NOOUTPUT  = 22
    METIS_OPTION_BALANCE   = 23
    METIS_OPTION_GTYPE     = 24
    METIS_OPTION_UBVEC     = 25

    # Partitioning method
    METIS_PTYPE_RB   = Cint(0)
    METIS_PTYPE_KWAY = Cint(1)
    # objectives
    METIS_OBJTYPE_CUT  = Cint(0)
    METIS_OBJTYPE_VOL  = Cint(1)
    METIS_OBJTYPE_NODE = Cint(2)
    # Coarsening method
    METIS_CTYPE_RM   = Cint(0)
    METIS_CTYPE_SHEM = Cint(1)
    # Initial partitioning schemes
    METIS_IPTYPE_GROW    = Cint(0)
    METIS_IPTYPE_RANDOM  = Cint(1)
    METIS_IPTYPE_EDGE    = Cint(2)
    METIS_IPTYPE_NODE    = Cint(3)
    METIS_IPTYPE_METISRB = Cint(4)
    # Refinement schemes
    METIS_RTYPE_FM        = Cint(0)
    METIS_RTYPE_GREEDY    = Cint(1)
    METIS_RTYPE_SEP2SIDED = Cint(2)
    METIS_RTYPE_SEP1SIDED = Cint(3)
    # Debug levels (bit positions)
    METIS_DBG_INFO       = Cint(1)    # Shows various diagnostic messages
    METIS_DBG_TIME       = Cint(2)    # Perform timing analysis
    METIS_DBG_COARSEN    = Cint(4)    # Show the coarsening progress
    METIS_DBG_REFINE     = Cint(8)    # Show the refinement progress
    METIS_DBG_IPART      = Cint(16)   # Show info on initial partitioning
    METIS_DBG_MOVEINFO   = Cint(32)   # Show info on vertex moves during refinement
    METIS_DBG_SEPINFO    = Cint(64)   # Show info on vertex moves during sep refinement
    METIS_DBG_CONNINFO   = Cint(128)  # Show info on minimization of subdomain connectivity
    METIS_DBG_CONTIGINFO = Cint(256)  # Show info on elimination of connected components
    METIS_DBG_MEMORY     = Cint(2048) # Show info related to wspace allocation
    # Graph types
    METIS_GTYPE_DUAL  = Cint(0)
    METIS_GTYPE_NODAL = Cint(1)
    #############################################################################################
    
    # number of common nodes between neighbors
    if finch_state.config.dimension == 1
        common_nodes = 1;
    elseif finch_state.config.dimension == 2
        common_nodes = 2;
    elseif finch_state.config.dimension == 3
        if mesh.nv[1] == 4
            common_nodes = 3;
        elseif mesh.nv[1] == 8
            common_nodes = 4;
        end
    end
    
    # options
    options = fill(Cint(-1), METIS_NOPTIONS);
    # It turns out the default values are all -1 
    # r = ccall((:METIS_SetDefaultOptions, libmetis), Cint, (Ptr{Cint},), options);
    # These can now be set as desired
    
    # Input arguments for metis
    ne = Cint(mesh.nel);# number of elements
    nn = Cint(mesh.nx); # number of nodes
    (eptr, eind) = make_eptr_eind(mesh); # connectivity in a particular format
    vwgt = C_NULL;      # vertex weights (not used)
    vsize = C_NULL;     # vertex size (not used)
    ncommon = Cint(common_nodes);  # number of common nodes to between neighbors
    nparts = Cint(np);  # number of partitions
    tpwgts = C_NULL;    # partition weights (not used)
    # options = C_NULL; # options vector
    # Output arguments
    objval = fill(Cint(0),1);  # objective value
    epart = fill(Cint(0), ne); # partitioned elements
    npart = fill(Cint(0), nn); # partitioned nodes
    
    r = ccall((:METIS_PartMeshDual, libmetis), Cint,
        (Ref{Cint}, Ref{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ref{Cint}, 
         Ref{Cint}, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        ne, nn, eptr, eind, vwgt, vsize, ncommon, nparts, tpwgts, options, objval, epart, npart);
    
    if !(r == METIS_OK)
        #METIS_ERROR_INPUT = Cint(-2)
        #METIS_ERROR_MEMORY = Cint(-3)
        #METIS_ERROR = Cint(-4)
        err_msg = ["Input arugment error","Memory allocation error","Unknown error"];
        printerr("METIS returned error: "*err_msg[-r-1], fatal=true);
    end
    
    return epart;
end

# Writes the element data in a special format for metis.
function make_eptr_eind(mesh)
    nodes_per_element = mesh.nv[1];
    eptr = fill(Cint(0), mesh.nel+1);
    eind = fill(Cint(0), mesh.nel * nodes_per_element);
    
    #eptr: the start index in eind for this element's nodes
    #eind: the indices of nodes
    for ei=1:mesh.nel
        offset = (ei-1)*nodes_per_element; # for 0-based index
        eptr[ei] = Cint(offset);
        for ni=1:nodes_per_element
            eind[offset+ni] = Cint(mesh.elements[ni,ei] - 1);
        end
    end
    eptr[end] = Cint(length(eind));
    
    return (eptr, eind);
end
