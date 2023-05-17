#=
# 2D Poisson, Dirichlet bc
# Using a mixed-element mesh (quads and triangles)
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("mixed2d");

useLog("mixed2dlog", level=3)

domain(2)
# functionSpace(order=2) # basis function polynomial order

mesh("mixedmesh.mesh")

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, 0)

coefficient("f", "-8*pi*pi*sin(2*pi*x)*sin(2*pi*y)")

weakForm(u, "dot(grad(u), grad(v)) + f*v")

exportCode("mixed2dcode");
importCode("mixed2dcodein");

solve(u);

finalizeFinch()

# exact solution is sin(2*pi*x)*sin(2*pi*y)
maxerr = 0;
exact(x,y) = sin(pi*x*2)*sin(pi*y*2);
xy = Finch.finch_state.grid_data.allnodes;

for i=1:size(xy,2)
    x = xy[1,i];
    y = xy[2,i];
    err = abs(u.values[i] - exact(x,y));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

## Uncomment below to plot or output result

using Plots
pyplot();

display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))

# output_values(u, "poisson2d", format="vtk");