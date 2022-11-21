#=
# Import a simple triangle or quad mesh from a .msh file.
# Then use it with an reaction-diffusion equation.
=#

### If the Finch package has already been added, use this line #########
# using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
##########################################################################

initFinch("unstruct2dtest");

useLog("unstruct2dlog", level=3)

domain(2, grid=UNSTRUCTURED)
functionSpace(order=2)

# This rectangle covers [0, 0.1]x[0, 0.3]
# Uncomment the desired mesh.
mesh("src/examples/utriangle.msh")  # Using triangles
#mesh("src/examples/uquad.msh")     # Using quads

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, 0)

coefficient("f", "(-10-(x+1)*200*pi*pi)*sin(10*pi*x)*sin(10*pi*y) + 10*pi*cos(10*pi*x)*sin(10*pi*y)")
coefficient("k", "x+1") # non-constant diffusion
coefficient("C", 10) # reaction
weakForm(u, "k*dot(grad(u), grad(v)) + C*u*v+ f*v")

solve(u);

finalizeFinch();

# exact solution is sin(10*pi*x)*sin(10*pi*y)
# check error
maxerr = 0;
for i=1:size(Finch.finch_state.grid_data.allnodes,2)
    x = Finch.finch_state.grid_data.allnodes[1,i];
    y = Finch.finch_state.grid_data.allnodes[2,i];
    err = abs(u.values[i] - sin(10*pi*x)*sin(10*pi*y));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));
