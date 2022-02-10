#=
# 1D advection matching the example in the book: Nodal Discontinuous Galerkin Method.
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("advection1d");

useLog("dgnewadvection1dlog", level=3)

domain(1)
solverType(DG)
functionSpace(order=3)
timeStepper(LSRK4, cfl=0.01)

mesh(LINEMESH, elsperdim=10, bids=2, interval=[0,2])

u = variable("u")
testSymbol("v")

T = 1;
timeInterval(T)
initial(u, "sin(x)")

boundary(u, 1, DIRICHLET, "-sin(2*pi*t)")
boundary(u, 2, NO_BC)

coefficient("a", 2*pi)
coefficient("beta", 0.2)

###################################################
# Uncomment one of the following two groups of code(or both)

### To generate the Julia code and export it to a set of files, use these two lines.
weakForm(u, "Dt(u*v) + a*grad(u)*v + surface(-a*normal()*u*v + a*normal()*ave(u)*v - 0.5*(1-beta)*a*jump(u)*v)")
exportCode("advec1d")
###

### To import previously generated or modified code, use this line
#importCode("advec1d")
###

###################################################

solve(u);

# Uncomment to plot
# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], u.values[1,:], marker=:circle, reuse=false))

# @finalize()
