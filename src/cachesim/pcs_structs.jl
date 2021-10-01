#=
Struct defs for using pycachesim
=#
function make_power_of_two(n)
    if n < 2
        return Culong(2);
    end
    if (n & (n - 1)) == 0
        return Culong(n);
    else
        # not a power of two, reduce to next lower power of two
        # n is a Cuint which is 32bit unsigned
        m = 2^31;
        while m>n
            m = m/2;
        end
        return Culong(m);
    end
end

function rep_policy_id(replacement_policy)
    if replacement_policy == "FIFO"
        return Cint(0);
    elseif replacement_policy == "LRU"
        return Cint(1);
    elseif replacement_policy == "MRU"
        return Cint(2);
    elseif replacement_policy == "RR"
        return Cint(3);
    else
        return Cint(1); # default to LRU
    end
end

########################################################################################################

struct cache_entry
    cl_id::Culong;

    dirty::Culong;  # if 0, content is in sync with main memory. if 1, it is not.
                        #used for write-back
    invalid::Culong;# denotes an entry which does not contain a valid cacheline.
                        # it is empty.
    cache_entry() = new(0,1,1);
end

########################################################################################################

struct addr_range 
    # Address range used to communicate consecutive accesses
    # last addr of range is addr+length-1
    addr::Culong;
    length::Culong;
end

########################################################################################################

struct stats
    count::Culonglong;
    byte::Culonglong;
    #cl::Culong; # might be used later
end

########################################################################################################

struct Cache
    name::Cstring;
    sets::Culong;
    ways::Culong;
    cl_size::Culong;
    cl_bits::Culong;
    subblock_size::Culong;
    subblock_bits::Culong;
    replacement_policy_id::Cint;    #// 0 = FIFO, 1 = LRU, 2 = MRU, 3 = RR
                                    #// (state is kept in the ordering)
                                    #// for LFU an additional field would be required to capture state
    write_back::Cint;       # 1 = write-back
                            # 0 = write-through
    write_allocate::Cint;   # 1 = write-allocate,
                            # 0 = non-write-allocate
    write_combining::Cint;  # 1 = this is a write-combining cache
                            # 0 = regular cache_entry
                         
    load_from::Ptr{Cache};
    store_to::Ptr{Cache};
    victims_to::Ptr{Cache};
    
    swap_on_load::Cint;

    placement::Ptr{cache_entry};
    subblock_bitfield::Ptr{UInt8};

    LOAD::stats;
    STORE::stats;
    HIT::stats;
    MISS::stats;
    EVICT::stats;

    verbosity::Cint;
    
    Cache(name, sets, ways, cl_size, replacement_policy, load_from = 0, store_to = 0) = (
        cname = Base.unsafe_convert(Cstring, name);
        replacement_policy_id = rep_policy_id(replacement_policy);
        cl_size2 = make_power_of_two(cl_size);
        load_from = load_from==0 ? Ptr{Cache}() : load_from;
        store_to = store_to==0 ? Ptr{Cache}() : store_to;
        
        new(cname, sets, ways, cl_size2, 0, cl_size2, 0, replacement_policy_id, 1, 1, 0, 
            load_from, store_to, Ptr{Cache}(), 0, Ptr{cache_entry}(), Ptr{UInt8}(), 
            stats(0,0),stats(0,0),stats(0,0),stats(0,0),stats(0,0), 0)
    );
end