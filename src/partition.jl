#=
Metis is used to partition the mesh.
This uses METIS_jll: libmetis
=#
const METIS_NOPTIONS = 40

# codes returned by metis functions
const METIS_OK = Cint(1)
const METIS_ERROR_INPUT = Cint(-2)
const METIS_ERROR_MEMORY = Cint(-3)
const METIS_ERROR = Cint(-4)

# location in options vector
const METIS_OPTION_PTYPE     = 1
const METIS_OPTION_OBJTYPE   = 2
const METIS_OPTION_CTYPE     = 3
const METIS_OPTION_IPTYPE    = 4
const METIS_OPTION_RTYPE     = 5
const METIS_OPTION_DBGLVL    = 6
const METIS_OPTION_NITER     = 7
const METIS_OPTION_NCUTS     = 8
const METIS_OPTION_SEED      = 9
const METIS_OPTION_NO2HOP    = 10
const METIS_OPTION_MINCONN   = 11
const METIS_OPTION_CONTIG    = 12
const METIS_OPTION_COMPRESS  = 13
const METIS_OPTION_CCORDER   = 14
const METIS_OPTION_PFACTOR   = 15
const METIS_OPTION_NSEPS     = 16
const METIS_OPTION_UFACTOR   = 17
const METIS_OPTION_NUMBERING = 18
const METIS_OPTION_HELP      = 19
const METIS_OPTION_TPWGTS    = 20
const METIS_OPTION_NCOMMON   = 21
const METIS_OPTION_NOOUTPUT  = 22
const METIS_OPTION_BALANCE   = 23
const METIS_OPTION_GTYPE     = 24
const METIS_OPTION_UBVEC     = 25

# Partitioning method
const METIS_PTYPE_RB   = Cint(0)
const METIS_PTYPE_KWAY = Cint(1)
# objectives
const METIS_OBJTYPE_CUT  = Cint(0)
const METIS_OBJTYPE_VOL  = Cint(1)
const METIS_OBJTYPE_NODE = Cint(2)
# Coarsening method
const METIS_CTYPE_RM   = Cint(0)
const METIS_CTYPE_SHEM = Cint(1)
# Initial partitioning schemes
const METIS_IPTYPE_GROW    = Cint(0)
const METIS_IPTYPE_RANDOM  = Cint(1)
const METIS_IPTYPE_EDGE    = Cint(2)
const METIS_IPTYPE_NODE    = Cint(3)
const METIS_IPTYPE_METISRB = Cint(4)
# Refinement schemes
const METIS_RTYPE_FM        = Cint(0)
const METIS_RTYPE_GREEDY    = Cint(1)
const METIS_RTYPE_SEP2SIDED = Cint(2)
const METIS_RTYPE_SEP1SIDED = Cint(3)
# Debug levels (bit positions)
const METIS_DBG_INFO       = Cint(1)    # Shows various diagnostic messages
const METIS_DBG_TIME       = Cint(2)    # Perform timing analysis
const METIS_DBG_COARSEN    = Cint(4)    # Show the coarsening progress
const METIS_DBG_REFINE     = Cint(8)    # Show the refinement progress
const METIS_DBG_IPART      = Cint(16)   # Show info on initial partitioning
const METIS_DBG_MOVEINFO   = Cint(32)   # Show info on vertex moves during refinement
const METIS_DBG_SEPINFO    = Cint(64)   # Show info on vertex moves during sep refinement
const METIS_DBG_CONNINFO   = Cint(128)  # Show info on minimization of subdomain connectivity
const METIS_DBG_CONTIGINFO = Cint(256)  # Show info on elimination of connected components
const METIS_DBG_MEMORY     = Cint(2048) # Show info related to wspace allocation
# Graph types
const METIS_GTYPE_DUAL  = Cint(0)
const METIS_GTYPE_NODAL = Cint(1)

#=
This is the primary function used to partition a mesh into np partitions.
It returns a list of partition numbers for each element.
=#
function get_element_partitions(mesh, np)
    # number of common nodes between neighbors
    if config.dimension == 1
        common_nodes = 1;
    elseif config.dimension == 2
        common_nodes = 2;
    elseif config.dimension == 3
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
        #const METIS_ERROR_INPUT = Cint(-2)
        #const METIS_ERROR_MEMORY = Cint(-3)
        #const METIS_ERROR = Cint(-4)
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
