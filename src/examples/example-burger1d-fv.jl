#=
1D Burger's equation using 1st order Godunov flux.
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################
init_finch("FVburger1d");

useLog("FVburger1dlog", level=3)

# Configuration setup
domain(1)
solverType(FV)
timeStepper(RK4, cfl=1)

# Mesh
n = 70 # number of elements
mesh(LINEMESH, elsperdim=n)

# Variables and BCs
u = variable("u", SCALAR, CELL)
boundary(u, 1, NO_BC)

v = variable("v", SCALAR, CELL)
boundary(v, 1, NO_BC)

# Time interval and initial condition
T = 1;
timeInterval(T)
initial(v, "sin(pi*2*x)")
initial(u, "x<0.5 ? -1 : 1")

# The flux term of the conservation equation
# See the file fv_ops.jl for the definition of burgerGodunov(u, f)
flux([u,v], ["0.1*burgerGodunov(u,u*u)", "0.1*burgerGodunov(v,v*v)"]) 

# printLatex(u)

# @exportCode("fvburgercode") # uncomment to export generated code to a file
# @importCode("fvburgercode") # uncomment to import code from a file

solve([u,v])

# output_values([u,v], "burger1d", format="vtk");

finalize_finch()

##### Uncomment below to plot #####

x = Finch.fv_info.cellCenters[:]
n = length(x);
u_ic = zeros(n);
v_ic = zeros(n);
for i=1:n
    u_ic[i] = x[i]<0.5 ? -1 : 1
    v_ic[i] = sin(2*pi*x[i])
end

using Plots
pyplot();
pu = plot([x x], [u.values[:] u_ic], markershape=:circle, label=["Godunov" "initial"])
pv = plot([x x], [v.values[:] v_ic], markershape=:circle, label=["Godunov" "initial"])
display(plot(pu, pv, layout=2))
