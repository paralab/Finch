#=
# 2D heat equation, Dirichlet bc, CG
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("heat2d");

useLog("heat2dlog", level=3)

# Set up the configuration
domain(2)              # dimension
functionSpace(order=2) # basis polynomial order
timeStepper(RK4)

mesh(QUADMESH, elsperdim=20) # 20x20 elements

u = variable("u")
testSymbol("v")

T = 1;
timeInterval(T) # The time interval is 0 to T
initial(u, "abs(x-0.5)+abs(y-0.5) < 0.2 ? 1 : 0") # initial condition

boundary(u, 1, DIRICHLET, 0)  # boundary condition for BID 1 is Dirichlet with value 0

coefficient("f", "0.5*sin(6*pi*x)*sin(6*pi*y)")

weakForm(u, "Dt(u*v) + 0.01 * dot(grad(u),grad(v)) - f*v")

solve(u);

exportCode("heat2dcode");

finalizeFinch()

## Uncomment below to plot ##

# # solution is stored in the variable's "values"
using Plots
pyplot();
display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st = :surface))


