#=
# Import a simple tet or hex mesh from a .msh file.
# Then use it with an reaction-diffusion equation.
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("unstruct3d");

useLog("unstruct3dlog", level=3)

domain(3, grid=UNSTRUCTURED)
functionSpace(order=2)

# Unit cube
mesh("src/examples/utet.msh")  # Using tets
#mesh("src/examples/uhex.msh")   # Using hexes

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, "sin(0.5*pi*x)*sin(pi*y)*sin(pi*z)")

coefficient("K", "x*x + 4") # non-constant diffusion
coefficient("C", 10) # reaction
coefficient("f", "x*pi*cos(pi*0.5*x)*sin(pi*y)*sin(pi*z) + (-10-(x*x+4)*pi*pi*2.25)*sin(pi*0.5*x)*sin(pi*y)*sin(pi*z)")

weakForm(u, "K*dot(grad(u), grad(v)) + C*u*v + f*v")

solve(u);

finalizeFinch();

# exact solution is sin(pi*0.5*x)*sin(pi*y)*sin(pi*z)
# check error
maxerr = 0;
for i=1:size(Finch.grid_data.allnodes,2)
    x = Finch.grid_data.allnodes[1,i];
    y = Finch.grid_data.allnodes[2,i];
    z = Finch.grid_data.allnodes[3,i];
    err = abs(u.values[i] - sin(0.5*pi*x)*sin(pi*y)*sin(pi*z));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));
