#=
A Cache simulator based on pycachesim by https://github.com/RRZE-HPC/pycachesim
pycachesim is licenced under AGPLv3 as detailed on the github page above.

Ported to Julia and modified by Eric Heisler

=#
module CacheSim

# These functions are exported, but all others are accessible
export Cache, cache_load, cache_store, 
       get_cacheSim_from_file, build_cache,
       print_stats, get_stats_string

"""
    CacheEntry

A Cache entry has a line number, dirty bit and invalid bit
"""
mutable struct CacheEntry
    cl_id::Int64  # cache line number
    dirty::Int8   # if 0, content is in sync with main memory. if 1, it is not. used for write-back
    invalid::Int8 # denotes an entry which does not contain a valid cacheline. it is empty.
end

"""
    Stats

Stats keeps track of counts and bytes for misses, hits, etc.
"""
mutable struct Stats
    count::Int64
    byte::Int64
end

"""
    Cache

The cache heirarchy is held in Cache objects. Each Cache represents one 
level and has info about the memory structure, policies, heirarchy, and stats.
"""
mutable struct Cache
    name::String
    sets::Int64
    ways::Int8
    cl_size::Int8
    cl_bits::Int8
    subblock_size::Int64
    subblock_bits::Int64
    
    replacement_policy_id::Int8 # 0 = FIFO, 1 = LRU, 2 = MRU, 3 = RR
    
    write_back::Int8      # 1 = write-back, 0 = write-through
    write_allocate::Int8  # 1 = write-allocate, 0 = non-write-allocate
    write_combining::Int8 # 1 = this is a write-combining cache, 0 = regular cache
    
    load_from::Union{Cache,Nothing}
    store_to::Union{Cache,Nothing}
    victims_to::Union{Cache,Nothing}
    
    swap_on_load::Int
    
    placement::Vector{CacheEntry}
    subblock_bitfield::Vector{UInt8}
    
    LOAD::Stats
    STORE::Stats
    HIT::Stats
    MISS::Stats
    EVICT::Stats
    
    verbosity::Int8
end

"""
    AddressRange

Convenience struct for representing a contiguous set of addresses.
"""
struct AddressRange
    addr::Int64
    length::Int64
end

"""
    cache_load(cache::Cache, addr::Int, bytes::Int)

Load this many bytes from this address.
"""
function cache_load(cache::Cache, addr::Int, bytes::Int)
    Cache__load(cache, AddressRange(addr, bytes));
end

"""
    cache_store(cache::Cache, addr::Int, bytes::Int)

Store this many bytes to this address.
"""
function cache_store(cache::Cache, addr::Int, bytes::Int)
    Cache__store(cache, AddressRange(addr, bytes), 0);
end

"""
    get_cacheSim_from_file(cache_file::String)

Build a set of caches from a definition file.
"""
function get_cacheSim_from_file(cache_file::String)
    #file for log output and errors/warnings (needed because stdout and stderr for some reason cannot be linked when using pin)
    # file  = fopen ("log_cachesim","w");
    # fprintf(file, "Cache* get_cacheSim_from_file(\"%s\"):\n\n", cache_file);
    # fflush(file);
    stream = open(cache_file, "r");
    lines = readlines(stream);
    close(stream);
    
    if length(lines) < 2
        println("Only found "*string(length(lines))*" lines in cache definition file. Need at least two.");
        return nothing;
    end
    
    num_levels = parse(Int, lines[1]);
    if num_levels > (length(lines)-1)
        println("Expected "*string(num_levels)*" cache levels, but only found "*string(length(lines))*" lines in cache definition file.");
        return nothing;
    end

    # buffer for cache objects
    cachelevels = [];
    #buffers to save information about which caches to link at the end
    load_from_buff = fill("",num_levels);
    store_to_buff = fill("",num_levels);
    victims_to_buff = fill("",num_levels);
    counter = 0;
    
    # Each line represents a cache
    # key value pairs seperated by ',' and spaces are ignored (key1=value1, key2=value2, ...)
    for i=2:length(lines)
        if length(lines[i]) == 0 || lines[i][1] == '#' # ignore empty lines or comments with '#'
            continue;
        end
        
        # create an empty cache to be set up
        newcache = make_empty_cache();
        
        kvpairs = split(lines[i], ',', keepempty=false);
        for j=1:length(kvpairs)
            keyvalue = split(kvpairs[j], [' ','='], keepempty=false);
            if length(keyvalue) != 2
                println("Error parsing key/value pair: "*string(kvpairs[j])*"\nExpected pattern is key=value");
                return nothing;
            end
            key = keyvalue[1];
            value = keyvalue[2];
            
            if key == "name"
                newcache.name = value;
            elseif key == "sets"
                newcache.sets = parse(Int, value);
            elseif key == "ways"
                newcache.ways = parse(Int, value);
            elseif key == "cl_size"
                newcache.cl_size = parse(Int, value);
            elseif key == "cl_bits"
                newcache.cl_bits = parse(Int, value);
            elseif key == "subblock_size"
                newcache.subblock_size = parse(Int, value);
            elseif key == "subblock_bits"
                newcache.subblock_bits = parse(Int, value);
            elseif key == "replacement_policy_id"
                newcache.replacement_policy_id = parse(Int, value);
            elseif key == "write_back"
                newcache.write_back = parse(Int, value);
            elseif key == "write_allocate"
                newcache.write_allocate = parse(Int, value);
            elseif key == "write_combining"
                newcache.write_combining = parse(Int, value);
            elseif key == "load_from"
                load_from_buff[counter+1] = value;
            elseif key == "store_to"
                store_to_buff[counter+1] = value;
            elseif key == "victims_to"
                victims_to_buff[counter+1] = value;
            elseif key == "swap_on_load"
                newcache.swap_on_load = parse(Int, value);
            else
                println("ignoring unrecognized parameter: "*key*" in "*kvpairs[j]);
            end
        end
        
        # This is the minimum required info
        if (newcache.name == "empty") println("Cache definition needs \"name\""); return nothing; end
        if (newcache.sets == 0) println("Cache definition needs \"sets\""); return nothing; end
        if (newcache.ways == 0) println("Cache definition needs \"ways\""); return nothing; end
        if (newcache.cl_size == 0) println("Cache definition needs \"cl_size\""); return nothing; end
        
        # anything that wasn't initialized is set to defaults
        if newcache.subblock_size == 0
            newcache.subblock_size = newcache.cl_size;
        end
        # Check if cl_size is a power of two
        if !(isPowerOfTwo(newcache.cl_size))
            println("cl_size is not a power of 2!");
            return nothing;
        end
        # Get number of bits in cacheline adressing
        newcache.cl_bits = log2uint(newcache.cl_size);
        
        # Check if subblock_size is a divisor of cl_size
        if (mod(newcache.cl_size, newcache.subblock_size) != 0)
            println("subblock_size needs to be a devisor of cl_size!");
            return nothing;
        end
        newcache.subblock_bits = Int(newcache.cl_size / newcache.subblock_size);
        
        # subblock_bitfield
        if (newcache.write_combining > 0 && newcache.subblock_size != newcache.cl_size)
            # Subblocking will be used:
            # since UInt8 is used as type, we need upper(subblock_bits/8) per placement
            newcache.subblock_bitfield = zeros(UInt8, BITNSLOTS(newcache.sets*newcache.ways*newcache.subblock_bits));
        end
        
        # init cache
        newcache.placement = Vector{CacheEntry}(undef, newcache.sets * newcache.ways);
        for j = 1:(newcache.sets * newcache.ways)
            newcache.placement[j] = CacheEntry(0,1,0);
        end
        
        push!(cachelevels, newcache);
        counter+=1;
    end
    
    if counter < num_levels
        println("Expected "*string(num_levels)*" cache levels, but only found "*string(counter)*" in cache definition file.");
        return nothing;
    end
    
    # link caches
    for i = 1:num_levels
        for j = 1:num_levels
            if (load_from_buff[i] == cachelevels[j].name)
                cachelevels[i].load_from = cachelevels[j];
            end
            if (store_to_buff[i] == cachelevels[j].name)
                cachelevels[i].store_to = cachelevels[j];
            end
            if (victims_to_buff[i] == cachelevels[j].name)
                cachelevels[i].victims_to = cachelevels[j];
            end
        end
    end

    # find top level cache as interface for the cacheSimulator
    first_level = nothing;
    for i = 1:num_levels
        if (cachelevels[i].load_from == nothing && cachelevels[i].store_to == nothing && cachelevels[i].victims_to == nothing)
            first_level = cachelevels[i];
        end
    end
    if first_level == nothing
        println("Couldn't determine top level cache. Make sure top level does not specify load_from/store_to/victims_to");
        return nothing;
    end
    
    return first_level;
end

"""
    build_cache(name::String, cl_size::Int, sets::Int, ways::Int; [kwargs])

Build a cache with these parameters.
kwargs(default):
 - `subblock_size`(cl_size)
 - `replacement_policy_id`(1=LRU)
 - `write_back`(0=no)
 - `write_allocate`(0=no)
 - `write_combining`(0=no)
 - `load_from`(nothing)
 - `store_to`(nothing)
 - `victims_to`(nothing)
 - `swap_on_load`(0)
 - `verbosity`(0)
"""
function build_cache(name::String, cl_size::Int, sets::Int, ways::Int; 
                     subblock_size=0, replacement_policy_id=1, write_back=0, 
                     write_allocate=0, write_combining=0, load_from=nothing,
                     store_to=nothing, victims_to=nothing, swap_on_load=0, verbosity=0)
    if subblock_size==0
        subblock_size = cl_size;
    end
    
    # To determine: cl_bits, subblock_bits, placement, subblock_bitfield, verbosity
    # Check if cl_size is a power of two
    if !(isPowerOfTwo(cl_size))
        println("cl_size is not a power of 2!");
        return nothing;
    end
    # Get number of bits in cacheline adressing
    cl_bits = log2uint(cl_size);
    
    # Check if subblock_size is a divisor of cl_size
    if (mod(cl_size, subblock_size) != 0)
        println("subblock_size needs to be a devisor of cl_size!");
        return nothing;
    end
    subblock_bits = Int(cl_size / subblock_size);
    
    # subblock_bitfield
    if (write_combining > 0 && subblock_size != cl_size)
        # Subblocking will be used:
        # since UInt8 is used as type, we need upper(subblock_bits/8) per placement
        subblock_bitfield = zeros(UInt8, BITNSLOTS(sets*ways*subblock_bits));
    else
        subblock_bitfield = zeros(UInt8,0);
    end
    
    # init placement
    placement = Vector{CacheEntry}(undef, sets * ways);
    for j = 1:(sets * ways)
        placement[j] = CacheEntry(0,1,0);
    end
    
    newcache = Cache(name, sets, ways, cl_size, cl_bits, subblock_size, subblock_bits,
                     replacement_policy_id, write_back, write_allocate, write_combining,
                     load_from, store_to, victims_to, swap_on_load, placement, 
                     subblock_bitfield, Stats(0,0), Stats(0,0), Stats(0,0), Stats(0,0), 
                     Stats(0,0), verbosity);
    
    return newcache;
end

"""
    print_stats(cache::Cache)

Print all stats for each cache;
"""
function print_stats(cache::Cache)
    println(cache.name*":");
    println("LOAD:  "*string(cache.LOAD.count)*"   size: "*string(cache.LOAD.byte));
    println("STORE: "*string(cache.STORE.count)*"   size: "*string(cache.STORE.byte));
    println("HIT:   "*string(cache.HIT.count)*"   size: "*string(cache.HIT.byte));
    println("MISS:  "*string(cache.MISS.count)*"   size: "*string(cache.MISS.byte));
    println("EVICT: "*string(cache.EVICT.count)*"   size: "*string(cache.EVICT.byte));

    if cache.load_from != nothing
        print_stats(cache.load_from);
    end
    if (cache.store_to != nothing && cache.store_to != cache.load_from)
        print_stats(cache.store_to);
    end
    if (cache.victims_to != nothing && cache.store_to != cache.load_from && cache.store_to != cache.victims_to)
        print_stats(cache.victims_to);
    end
end

"""
    get_stats_string(cache::Cache)

Similar to print_stats, but returns the output string rather than printing it.
"""
function get_stats_string(cache::Cache)
    result = "";
    result *= cache.name*":\n";
    result *= "LOAD:  "*string(cache.LOAD.count)*"   size: "*string(cache.LOAD.byte)*"\n";
    result *= "STORE: "*string(cache.STORE.count)*"   size: "*string(cache.STORE.byte)*"\n";
    result *= "HIT:   "*string(cache.HIT.count)*"   size: "*string(cache.HIT.byte)*"\n";
    result *= "MISS:  "*string(cache.MISS.count)*"   size: "*string(cache.MISS.byte)*"\n";
    result *= "EVICT: "*string(cache.EVICT.count)*"   size: "*string(cache.EVICT.byte)*"\n";

    if cache.load_from != nothing
        result *= get_stats_string(cache.load_from);
    end
    if (cache.store_to != nothing && cache.store_to != cache.load_from)
        result *= get_stats_string(cache.store_to);
    end
    if (cache.victims_to != nothing && cache.store_to != cache.load_from && cache.store_to != cache.victims_to)
        result *= get_stats_string(cache.victims_to);
    end
    return result;
end

function make_empty_cache()
    return Cache(
        "empty",0,0,0,0,0,0,1,0,0,0,
        nothing,nothing,nothing,
        0,[],[],
        Stats(0,0),Stats(0,0),Stats(0,0),Stats(0,0),Stats(0,0),
        0
    );
end

# Array of bits as found in comp.lang.c FAQ Question 20.8: http://c-faq.com/misc/bitsets.html
function BITMASK(b::Int64) return (1 << mod(b,8)); end
function BITSLOT(b::Int64) return Int(floor(b / 8)); end
function BITNSLOTS(nb::Int64) return Int(floor((nb + 7) / 8)); end
function BITSET(a::Vector{UInt8}, b::Int64) a[BITSLOT(b)+1] |= BITMASK(b); end
function BITCLEAR(a::Vector{UInt8}, b::Int64) a[BITSLOT(b)+1] &= ~BITMASK(b); end
function BITTEST(a::Vector{UInt8}, b::Int64) return (a[BITSLOT(b)+1] & BITMASK(b)); end

function log2uint(x::Int)
    ans = 0;
    ux::UInt = x;
    while (ux >> 1) > 0
        ux >>= 1;
        ans+=1;
    end
    return ans;
end

function isPowerOfTwo(x::Int)
    return ((x != 0) && (x & (x - 1))==0);
end

function Cache__get_cacheline_id(self::Cache, addr::Int64)
    return addr >> self.cl_bits;
end

function Cache__get_set_id(self::Cache, cl_id::Int64)
    return cl_id % self.sets;
end

function __range_from_addrs(addr::Int64, last_addr::Int64)
    return AddressRange(addr, last_addr-addr-1);
end

function Cache__get_addr_from_cl_id(self::Cache, cl_id::Int64)
    return cl_id << self.cl_bits;
end

function Cache__get_range_from_cl_id(self::Cache, cl_id::Int64)
    return AddressRange(Cache__get_addr_from_cl_id(self, cl_id), self.cl_size);
end

# Creates a range from a cacheline id and another range.
# The returned range will always be a subset of (or at most the same as) the given range
function Cache__get_range_from_cl_id_and_range(self::Cache, cl_id::Int64, range::AddressRange)
    addr = Cache__get_addr_from_cl_id(self, cl_id);
    addr = addr > range.addr ? addr : range.addr;
    length = (addr + self.cl_size) < (range.addr + range.length) ?
        self.cl_size : range.addr + range.length - addr;
    return AddressRange(addr, length);
end

# Returns the location a cacheline has in a cache
# if cacheline is not present, returns -1
# TODO use sorted data structure for faster searches in case of large number of
# ways or full-associativity?
function Cache__get_location(self::Cache, cl_id::Int64, set_id::Int64)
    for i=1:self.ways
        if (self.placement[set_id*self.ways + i].invalid == 0 &&
            self.placement[set_id*self.ways + i].cl_id == cl_id)
            return i;
        end
    end
    return -1; # Not found
end

# Signals request of addr range by higher level. This handles hits and misses.
function Cache__load(self::Cache, range::AddressRange)
    self.LOAD.count+=1;
    self.LOAD.byte += range.length;
    placement_idx = -1;

    # Handle range:
    last_cl_id = Cache__get_cacheline_id(self, range.addr+range.length-1);
    for cl_id = Cache__get_cacheline_id(self, range.addr):last_cl_id
        set_id = Cache__get_set_id(self, cl_id);
        
        # Check if cl_id is already cached
        location = Cache__get_location(self, cl_id, set_id);
        if location != -1
            # HIT: Found it!
            self.HIT.count+=1;
            # We only add actual bytes that were requested to hit.byte
            self.HIT.byte += self.cl_size < range.length ? self.cl_size : range.length;
            
            entry = self.placement[set_id*self.ways + location + 1];
            
            if (self.replacement_policy_id == 0 || self.replacement_policy_id == 3)
                # FIFO: nothing to do
                # RR: nothing to do
                placement_idx = self.ways-1;
                continue;
            else # if(self.replacement_policy_id == 1 || self.replacement_policy_id == 2) {
                # LRU: Reorder elements to account for access to element
                # MRU: Reorder elements to account for access to element
                if location != 0
                    for j = location:-1:1 # (int j=location; j>0; j--) {
                        self.placement[set_id*self.ways + j + 1] = self.placement[set_id*self.ways + j];
                        
                        # Reorder bitfild in accordance to queue
                        if self.write_combining == 1
                            for i = 0:self.subblock_bits #(long i=0; i<self->subblock_bits; i+=1) {
                                if (BITTEST(self.subblock_bitfield,
                                           set_id*self.ways*self.subblock_bits +
                                           (j-1)*self.subblock_bits + i))
                                    BITSET(self.subblock_bitfield,
                                           set_id*self.ways*self.subblock_bits +
                                           j*self.subblock_bits + i);
                                else
                                    BITCLEAR(self.subblock_bitfield,
                                             set_id*self.ways*self.subblock_bits +
                                             j*self.subblock_bits + i);
                                end
                            end
                        end
                    end
                    self.placement[set_id*self.ways + 1] = entry;
                end
                placement_idx = 0;
                continue;
            end
            # TODO if this is an exclusive cache, swap delivered cacheline with swap_cl_id (here and at end -> DO NOT RETURN)
        end

        # MISS!
        self.MISS.count+=1;
        # We only add actual bytes that were requested to miss.byte
        self.MISS.byte += self.cl_size < range.length ? self.cl_size : range.length;
        
        # Load from lower cachelevel
        # Check victim cache, if available
        victim_hit = 0;
        if self.victims_to != nothing
            # check for hit in victim cache
            victim_set_id = Cache__get_set_id(self.victims_to, cl_id);
            victim_location_victim = Cache__get_location(self.victims_to, cl_id, victim_set_id);
            if victim_location_victim != -1
                # hit in victim cache
                # load data from victim cache
                Cache__load(self.victims_to, Cache__get_range_from_cl_id(self, cl_id));
                # do NOT go onto load_from cache
                victim_hit = 1;
            end
        end
        # If no hit in victim cache, or no victim cache available, go to next cache level
        if (victim_hit==0 && self.load_from != nothing)
            # TODO use replace_entry to inform other cache of swap (in case of exclusive caches)
            Cache__load(self.load_from, Cache__get_range_from_cl_id(self, cl_id));
            # TODO, replace_cl_id);
        end # else last-level-cache
        
        entry = CacheEntry(cl_id, 0, 0);
        
        # Inject new entry into own cache. This also handles replacement.
        placement_idx = Cache__inject(self, entry);

        # TODO if this is an exclusive cache (swap_on_load = True), swap delivered cacheline with swap_cl_id (here and at hit)
    end
    # TODO Does this make sens or multiple cachelines? It is atm only used by write-allocate,
    # which should be fine, because requests are already split into individual cachelines
    return placement_idx;
end

function Cache__store(self::Cache, range::AddressRange, non_temporal::Int)
    self.STORE.count+=1;
    self.STORE.byte += range.length;
    # Handle range:
    last_cl_id = Cache__get_cacheline_id(self, range.addr+range.length-1);
    for cl_id = Cache__get_cacheline_id(self, range.addr):last_cl_id
        set_id = Cache__get_set_id(self, cl_id);
        location = Cache__get_location(self, cl_id, set_id);
        
        if (self.write_allocate == 1 && non_temporal == 0)
            # Write-allocate policy
            # Make sure line is loaded into cache (this will produce HITs and MISSes, iff it is 
            # not present in this cache):
            if location == -1
                # TODO does this also make sens if store with write-allocate and MISS happens on L2?
                # or would this inject byte loads instead of CL loads into the statistic
                # TODO makes no sens if first level is write-through (all byte requests hit L2)
                location = Cache__load(self, Cache__get_range_from_cl_id(self, cl_id));
            end
        elseif (location == -1 && self.write_back == 1)
            # In non-temporal store case, write-combining or write-through:
            # If the cacheline is not yet present, we inject a cachelien without loading it
            entry = CacheEntry(cl_id, 1, 0);
            location = Cache__inject(self, entry);
        end

        # Mark address range as cached in the bitfield
        if self.write_combining == 1
            # If write_combining is active, set the touched bits:
            # Extract local range
            cl_start = Cache__get_addr_from_cl_id(self, cl_id);
            start = range.addr > cl_start ? range.addr : cl_start;
            aend = range.addr+range.length < cl_start+self.cl_size ?
                                range.addr+range.length : cl_start+self.cl_size;
            for i = (start-cl_start):(aend-cl_start-1) #(long i=start-cl_start; i<end-cl_start; i+=1) {
                BITSET(self.subblock_bitfield,
                       set_id*self.ways*self.subblock_bits + location*self.subblock_bits + i);
            end
        end

        if (self.write_back == 1 && location != -1)
            # Write-back policy and cache-line in cache
            # Mark cacheline as dirty for later write-back during eviction
            self.placement[set_id*self.ways+location + 1].dirty = 1;
        else
            # Write-through policy or cache-line not in cache
            # Store to lower cachelevel
            # TODO use Cache__inject
            if self.store_to != nothing
                store_range = Cache__get_range_from_cl_id_and_range(self, cl_id, range);
                self.EVICT.count+=1;
                self.EVICT.byte += store_range.length;
                Cache__store(self.store_to,
                             store_range,
                             non_temporal);
            end # else last-level-cache
        end
    end
end

#=
Injects a cache entry into a cache and handles all side effects that might occur:
    - choose replacement according to policy
    - reorder queues
    - inform victim caches
    - handle write-back on replacement
=#
function Cache__inject(self::Cache, entry::CacheEntry)
    set_id = Cache__get_set_id(self, entry.cl_id);
    
    # Get cacheline id to be replaced according to replacement strategy
    if self.replacement_policy_id == 0 || self.replacement_policy_id == 1
        # FIFO: replace end of queue
        # LRU: replace end of queue
        replace_idx = 0;
        replace_entry = self.placement[set_id*self.ways + self.ways];

        # Reorder queue
        for i = self.ways:-1:2 #(long i=self->ways-1; i>0; i--) {
            self.placement[set_id*self.ways + i] = self.placement[set_id*self.ways + i-1];

            # Reorder bitfild in accordance to queue
            if self.write_combining == 1
                for j=0:self.subblock_bits-1 #(long j=0; j<self->subblock_bits; j+=1) {
                    if BITTEST(self.subblock_bitfield, set_id*self.ways*self.subblock_bits +
                               (i-2)*self.subblock_bits + j)
                        BITSET(self.subblock_bitfield, set_id*self.ways*self.subblock_bits +
                               (i-1)*self.subblock_bits + j);
                    else
                        BITCLEAR(self.subblock_bitfield, set_id*self.ways*self.subblock_bits +
                                 (i-1)*self.subblock_bits + j);
                    end
                end
            end
        end
    elseif self.replacement_policy_id == 2
        # MRU: replace first of queue
        replace_idx = self.ways - 1;
        replace_entry = self.placement[set_id*self.ways + 1];

        # Reorder queue
        for i=1:self.ways #(long i=0; i>self->ways-1; i+=1) {
            self.placement[set_id*self.ways+i] = self.placement[set_id*self.ways+i+1];

            # Reorder bitfild in accordance to queue
            if self.write_combining == 1
                for j=0:self.subblock_bits-1 #(long j=0; j<self->subblock_bits; j+=1) {
                    if (BITTEST(self.subblock_bitfield, set_id*self.ways*self.subblock_bits +
                               (i)*self.subblock_bits + j))
                        BITSET(self.subblock_bitfield, set_id*self.ways*self.subblock_bits +
                               (i-1)*self.subblock_bits + j);
                    else
                        BITCLEAR(self.subblock_bitfield, set_id*self.ways*self.subblock_bits +
                                 (i-1)*self.subblock_bits + j);
                    end
                end
            end
        end
    else # if (self.replacement_policy_id == 3)
        # RR: replace random element
        replace_idx = rand(0:(self.ways-1));
        replace_entry = self.placement[set_id*self.ways+replace_idx + 1];
    end

    # Replace other cacheline according to replacement strategy (using placement order as state)
    self.placement[set_id*self.ways + replace_idx + 1] = entry;

    # ignore invalid cache lines for write-back or victim cache
    if replace_entry.invalid == 0
        # write-back: check for dirty bit of replaced and inform next lower level of store
        if (self.write_back == 1 && replace_entry.dirty == 1)
            self.EVICT.count+=1;
            self.EVICT.byte += self.cl_size;
            if self.store_to != nothing
                non_temporal = 0; # default for non write-combining caches

                if self.write_combining == 1
                    # Check if non-temporal store may be used or write-allocate is necessary
                    non_temporal = 1;
                    for i=0:(self.subblock_bits-1) #(long i=0; i<self->subblock_bits; i+=1) {
                        if (BITTEST(self.subblock_bitfield,
                                set_id*self.ways*self.subblock_bits +
                                replace_idx*self.subblock_bits + i)) 
                            # incomplete cacheline, thus write-allocate is necessary
                            non_temporal = 0;
                        end
                        # Clear bits for future use
                        BITCLEAR(self.subblock_bitfield,
                                set_id*self.ways*self.subblock_bits +
                                replace_idx*self.subblock_bits + i);
                    end
                end
                # TODO addrs vs cl_id is not nicely solved here
                Cache__store(
                    self.store_to,
                    Cache__get_range_from_cl_id(self, replace_entry.cl_id),
                    non_temporal);
                    
            end # else last-level-cache
        elseif (self.victims_to != nothing)
            # Deliver replaced cacheline to victim cache, if neither dirty or already write_back
            # (if it were dirty, it would have been written to store_to if write_back is enabled)
            # Inject into victims_to
            Cache__inject(self.victims_to, replace_entry);
            # Take care to include into evict stats
            self.EVICT.count+=1;
            self.EVICT.byte += self.cl_size;
            victims_to.STORE.count+=1;
            victims_to.STORE.byte += self.cl_size;
        end
    end
    
    return replace_idx;
end

# A simple test just to make sure something works.
# name=L1,sets=64,ways=8,cl_size=64,replacement_policy_id=1,write_back=1,write_allocate=1,subblock_size=64,load_from=L2,store_to=L2
# name=L2,sets=512,ways=8,cl_size=64,replacement_policy_id=1,write_back=1,write_allocate=1,load_from=L3,store_to=L3
# name=L3,sets=9216,ways=16,cl_size=64,replacement_policy_id=1,write_back=1,write_allocate=1
function test_cachesim()
    L3 = build_cache("L3", 64, 9216, 16, write_back=1, write_allocate=1);
    L2 = build_cache("L2", 64, 512, 8, write_back=1, write_allocate=1, load_from=L3, store_to=L3);
    L1 = build_cache("L1", 64, 64, 8, write_back=1, write_allocate=1, load_from=L2, store_to=L2);
    
    cache_load(L1, 2342, 8);
    cache_store(L1, 512, 128);
    cache_load(L1, 512, 16);
    
    print_stats(L1);
end

end # module