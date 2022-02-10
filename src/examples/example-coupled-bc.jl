#=
A diffusion equation with different species coupled through boundary conditions only.
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("coupledbc1d");

useLog("coupledbc1dlog")

# Configuration setup
domain(1)
solverType(CG)
functionSpace(order=2)
timeStepper(EULER_IMPLICIT)

# Mesh
n = 40 # number of elements
mesh(LINEMESH, elsperdim=n)

# Variables and BCs
a = variable("a")
b = variable("b")
testSymbol("v")
boundary(a, 1, NEUMANN, 0)
boundary(b, 1, DIRICHLET, "2*a/(1+a*b)")

# Time interval and initial condition
T = 1;
timeInterval(T)
initial(a, "sin(pi*x)^2")
initial(b, 0)

weakForm([a,b], ["Dt(a*v) + 0.01 * dot(grad(a),grad(v))", 
                    "Dt(b*v) + 0.01 * dot(grad(b),grad(v))"])

exportCode("coupledbc1dcode") # uncomment to export generated code to a file
#importCode("coupledbc1dcode") # uncomment to import code from a file

solve([a,b])

finalize_finch()

##### Uncomment below to plot
# x = Finch.grid_data.allnodes[1,:];
# using Plots
# pyplot();
# display(plot([x x], [a.values[:] b.values[:]], marker=:circle))
