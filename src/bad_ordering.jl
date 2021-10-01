#=
# A random node reordering.
=#
export reorder_grid_random!, random_order

function random_order(N, seed)
    rng = Random.MersenneTwister(seed); # Use a specific RNG seed to make results consistent
    return randperm(rng, N); # Randomly permute node order
end

function reorder_grid_random!(grid, seed = 17)
    nnodes = size(grid.allnodes,2);
    
    new_order= random_order(nnodes, seed);
    reorder_grid_nodes!(grid, new_order);
    return grid;
end