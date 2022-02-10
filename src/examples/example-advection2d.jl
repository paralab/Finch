#=
# 2D advection using DG
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("advection2d");

useLog("advection2dlog")

domain(2)
solverType(DG)
functionSpace(order=3)
timeStepper(CRANK_NICHOLSON, cfl=0.1)

mesh(QUADMESH, elsperdim=8)

u = variable("u")
testSymbol("v")

T = 2;
timeInterval(T)
initial(u, "0")

### Add boundary regions and set BCs ###############
addBoundaryID(2, (x,y) -> (x <= 0 || y <= 0))
addBoundaryID(3, (x,y) -> (x < 0.2 && y < 0.2))
boundary(u, 1, NO_BC)
boundary(u, 2, DIRICHLET, 0)
boundary(u, 3, DIRICHLET, "sin(pi*(t - 0.7*(x+y)))^2")
####################################################

coefficient("a", type=VECTOR, [0.5, 0.5])
coefficient("beta", 0)
weakForm(u, "Dt(u*v) + a*grad(u)*v + surface(-dot(a, normal())*u*v + dot(a, normal())*ave(u)*v - 0.5*(1-beta)*dot(a, jump(u))*v)")

solve(u);

# Uncomment to plot
# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st=:surface, reuse=false))

finalize_finch()
