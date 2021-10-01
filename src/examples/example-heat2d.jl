#=
# 2D heat equation, Dirichlet bc, CG
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("heat2d");

# Optionally generate a log
useLog("heat2dlog")

# Set up the configuration (order doesn't matter)
domain(2)                   # dimension
solverType(CG)              # Use CG solver (default)
functionSpace(order=2)      # basis polynomial order
nodeType(LOBATTO)           # GLL elemental node arrangement (default)
timeStepper(RK4)            # time stepper (optional second arg is CFL#)

# Specify the problem
mesh(QUADMESH, elsperdim=10) # 10x10 elements

u = variable("u")            # make a scalar variable with symbol u
testSymbol("v")              # sets the symbol for a test function

T = 1;
timeInterval(T)              # The time interval is 0 to T
initial(u, "abs(x-0.5)+abs(y-0.5) < 0.2 ? 1 : 0") # initial condition

boundary(u, 1, DIRICHLET, 0)  # boundary condition for BID 1 is Dirichlet with value 0

# Write the weak form
coefficient("f", "0.5*sin(6*pi*x)*sin(6*pi*y)")
weakForm(u, "Dt(u*v) + 0.01 * dot(grad(u),grad(v)) - f*v")

solve(u);

# solution is stored in the variable's "values"
# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st = :surface))

# Dump things to the log if desired
log_dump_config();
log_dump_prob();

finalize_finch() # Finish writing and close any files
