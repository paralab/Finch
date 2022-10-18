#=
# 2D Poisson, Dirichlet bc
=#

### If the Finch package has already been added, use this line #########
# using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
##########################################################################

init_finch("poisson2d");

useLog("poisson2dlog", level=3)

domain(2)
functionSpace(order=2) # basis function polynomial order

# mesh
import_mesh=false;
if import_mesh
    # a 0.1 x 0.3 rectangle with triangle elements
    mesh("src/examples/utriangle.msh")
    # same but with unstructured quads
    # mesh("src/examples/uquad.msh")
    # uniform mesh, but considered unstructured
    # mesh("src/examples/squad.msh")
else
    mesh(QUADMESH, elsperdim=[10,30], interval=[0,0.1,0,0.3])
end

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, 0)

coefficient("f", "-(100 + 1/0.3/0.3)*pi*pi*sin(pi*x/0.1)*sin(pi*y/0.3)")

weakForm(u, "dot(grad(u), grad(v)) + f*v")

exportCode("poisson2dcode");

solve(u);

finalize_finch()

# exact solution is sin(pi*x)*sin(pi*y)
maxerr = 0;
exact(x,y) = sin(pi*x/0.1)*sin(pi*y/0.3);

for i=1:size(Finch.grid_data.allnodes,2)
    x = Finch.grid_data.allnodes[1,i];
    y = Finch.grid_data.allnodes[2,i];
    err = abs(u.values[i] - exact(x,y));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

## Uncomment below to plot or output result

# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st=:surface))

# output_values(u, "poisson2d", format="vtk");