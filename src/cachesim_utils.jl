#=
The cache simulator module is in CacheSim.jl

This file has functions for Finch to make use of CacheSim
=#

# A collection of memory chunks to access
mutable struct CachesimArrays
    address::Vector{Int} # start addresses of each array
    length::Vector{Int}  # number of 64-bit numbers in each array (not bytes)
    CachesimArrays() = new(zeros(Int,0),zeros(Int,0));
end

# Add an array to the collection
function add_cachesim_array(arrays::CachesimArrays, size::Int)
    next_index = length(arrays)+1;
    if length(arrays.address) == 0
        next_address = 0;
    else
        next_address = arrays.address[end] + 8 * arrays.length[end];
    end
    
    push!(arrays.address, next_address);
    push!(arrays.length, size);
    
    return next_index;
end

# Load all of the bytes of an array specified in the collection
# The access is sequential and done all at once.
function load_cachesim_array(cache::Cache, arrays::CachesimArrays, ind::Int)
    cache_load(cache, arrays.address[ind], arrays.length[ind]*8);
end

# Similar to above, but for storing
function store_cachesim_array(cache::Cache, arrays::CachesimArrays, ind::Int)
    cache_store(cache, arrays.address[ind], arrays.length[ind]*8);
end

# Load specified elements of an array in the collection
# This is done one at a time.
function load_cachesim_array_indexed(cache::Cache, arrays::CachesimArrays, array_id::Int, indices::Vector{Int})
    for i in indices
        cache_load(cache, arrays.address[ind] + (i-1)*8, 8);
    end
end

# Similar to above, but for storing
function store_cachesim_array_indexed(cache::Cache, arrays::CachesimArrays, array_id::Int, indices::Vector{Int})
    for i in indices
        cache_store(cache, arrays.address[ind] + (i-1)*8, 8);
    end
end
