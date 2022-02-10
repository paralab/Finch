#=
# 2D Poisson, Dirichlet bc
# Uses Matlab imported as a custom gen target
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("poisson2dmatlab");

# Try making an optional log
useLog("poisson2dmatlablog")

# Generate for the target in customtarget.jl
generateFor(MATLAB)

# Set up the configuration (order doesn't matter)
domain(2)
functionSpace(order=3)

# Specify the problem
mesh(QUADMESH, elsperdim=30)

u = variable("u")

testSymbol("v")

boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
coefficient("f", "-10*pi*pi*sin(pi*x)*sin(3*pi*y)")
weakForm(u, "-dot(grad(u),grad(v)) - f*v")

solve(u);

finalize_finch()
