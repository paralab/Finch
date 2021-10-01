#=
# Import a simple tet or hex mesh from a .msh file.
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("unstruct3d");

useLog("unstruct3dlog")

domain(3, grid=UNSTRUCTURED)
functionSpace(order=2)

# Unit cube
mesh("utet.msh")  # Using tets
#mesh("uhex.msh")   # Using hexes

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, "sin(0.5*pi*x)*sin(pi*y)*sin(pi*z)")

coefficient("K", "x*x + 4")
coefficient("C", 10)
coefficient("f", "x*pi*cos(pi*0.5*x)*sin(pi*y)*sin(pi*z) + (-10-(x*x+4)*pi*pi*2.25)*sin(pi*0.5*x)*sin(pi*y)*sin(pi*z)")

# PDE: div(K*grad(u)) -C*u = f
weakForm(u, "K*dot(grad(u), grad(v)) + C*u*v + f*v")

solve(u);

# exact solution is sin(pi*0.5*x)*sin(pi*y)*sin(pi*z)
# check error
maxerr = 0;
allerr = zeros(size(u.values));
exact(x,y,z) = sin(0.5*pi*x)*sin(pi*y)*sin(pi*z);

for i=1:size(Finch.grid_data.allnodes,2)
    x = Finch.grid_data.allnodes[1,i];
    y = Finch.grid_data.allnodes[2,i];
    z = Finch.grid_data.allnodes[3,i];
    err = abs(u.values[i] - exact(x,y,z));
    allerr[i] = err;
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], allerr[:], st=:surface))

finalize_finch();
