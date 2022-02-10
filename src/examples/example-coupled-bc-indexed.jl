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

useLog("coupledbc1dindlog")

# Configuration setup
domain(1, grid=UNSTRUCTURED)
solverType(CG)
functionSpace(order=2)
timeStepper(EULER_IMPLICIT)

# Mesh
n = 40 # number of elements
mesh(LINEMESH, elsperdim=n, bids=2)

# Variables and BCs
flavor = index("flavor", range = [1,3])
c = variable("c", type=VAR_ARRAY, index = flavor)

boundary(c, 1, DIRICHLET,  ["0.5*cos(pi*t)^2", 
                            "0.5 - (c[1]-c[3])/(2+(c[1]-c[3]))",
                            "(c[1] + c[2])/2"])
boundary(c, 2, NEUMANN,  [0,0,0])

testSymbol("v")

# Time interval and initial condition
T = 10.25;
timeInterval(T)
# initial(c, ["sin(pi*x)^2",
#             "sin(pi*2*x)^2",
#             "sin(pi*3*x)^2"])
initial(c, [0.5,0.5,0.5])
# initial(c, ["sin(pi*x)^2",
#             "0"])
# initial(c, "sin(pi*x)^2")

weakForm(c, "Dt(c[flavor]*v) + 0.1 * dot(grad(c[flavor]),grad(v))")

assemblyLoops(c, ["elements", flavor]);

# exportCode("coupledbcind1dcode") # uncomment to export generated code to a file
# importCode("coupledbc1dcode") # uncomment to import code from a file

solve(c)

finalize_finch()

##### Uncomment below to plot
# x = Finch.grid_data.allnodes[1,:];
# using Plots
# pyplot();
# display(plot([x x x], [c.values[1,:] c.values[2,:] c.values[3,:]], marker=:circle))


# display(plot([x x], [c.values[1,:] c.values[2,:]], marker=:circle))
# display(plot([x], [c.values[1,:]], marker=:circle))
