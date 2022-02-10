#=
# 2D Poisson, Dirichlet bc
# CG, Linear element
# Simplest test possible
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("poisson2d");

useLog("poisson2dlog", level=3)

# Set up the configuration (order doesn't matter)
domain(2)              # dimension
functionSpace(order=2) # basis function polynomial order

# Specify the problem
mesh(QUADMESH, elsperdim=[20,10], interval=[0,2,0,1])

variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
coefficient("f", "-2*pi*pi*sin(pi*x)*sin(pi*y)")
weakForm(u, "-dot(grad(u), grad(v)) - f*v")

solve(u);

finalize_finch()

# exact solution is sin(pi*x)*sin(pi*y)
# check error
maxerr = 0;
exact(x,y) = sin(pi*x)*sin(pi*y);

for i=1:size(Finch.grid_data.allnodes,2)
    x = Finch.grid_data.allnodes[1,i];
    y = Finch.grid_data.allnodes[2,i];
    err = abs(u.values[i] - exact(x,y));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

# solution is stored in the variable's "values"
# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st=:surface))

# output_values(u, "poisson2d", format="vtk");