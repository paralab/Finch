#=
# Cachesim example
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("cachesim");
useLog("cachesimlog")

#################################################
# For cache sim
cachesim(true); # Select cachesim. Normal solve won't happen
l1 = build_cache_level(1, 32, 8, 64, "LRU");    #
l2 = build_cache_level(2, 256, 4, 64, "LRU");   # Set up the cache specs (this is my laptop)
l3 = build_cache_level(3, 12288, 16, 64, "LRU");# (level, sets, ways, line_size, policy)
build_cache([l1, l2, l3]);
#################################################

n = 10;
ord = 2;

# Set up the configuration 
domain(3)
functionSpace(order=ord)

# Specify the problem
mesh(HEXMESH, elsperdim=n)

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
coefficient("f", "-3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)")
weakForm(u, "-dot(grad(u), grad(v)) - f*v")

cachesimSolve(u);

finalize_finch()
