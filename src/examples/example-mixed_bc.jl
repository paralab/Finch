#=
# 2D Poisson, mixed bc
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("mixedbc");

# Try making an optional log
@useLog("mixedbclog")

# Set up the configuration (order doesn't matter)
@domain(2)                          # dimension
@solver(CG)                         # Use CG solver
@functionSpace(LEGENDRE, 3)         # basis function, order
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(QUADMESH, 10, 4)               # build uniform mesh. 2nd arg=# of elements, (optional)3rd arg=# of BIDs

@variable(u)                        # same as @variable(u, SCALAR)

@testSymbol(v)                      # sets the symbol for a test function

@boundary(u, 1, NEUMANN, "-2*pi*cos(pi/4) * sin(2*pi*(y+0.125))")       # x=0
@boundary(u, 2, NEUMANN, "2*pi*cos(2*pi*1.125) * sin(2*pi*(y+0.125))")  # x=1
@boundary(u, 3, DIRICHLET, "sin(pi/4) * sin(2*pi*(x+0.125))")           # y=0
@boundary(u, 4, DIRICHLET, "sin(2*pi*1.125) * sin(2*pi*(x+0.125))")     # y=1

# Write the weak form 
@coefficient(f, "-8*pi*pi*sin(2*pi*(x+0.125))*sin(2*pi*(y+0.125))")
@weakForm(u, "-grad(u)*grad(v) - f*v")

solve(u);

# exact solution
# check error
maxerr = 0;
exact(x,y) = sin(2*pi*(x+0.125))*sin(2*pi*(y+0.125));

for i=1:size(Finch.grid_data.allnodes,2)
    x = Finch.grid_data.allnodes[1,i];
    y = Finch.grid_data.allnodes[2,i];
    err = abs(u.values[i] - exact(x,y));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

# solution is stored in the variable's "values"
#using Plots
#pyplot();
#display(plot(Finch.grid_data.allnodes[:,1], Finch.grid_data.allnodes[:,2], u.values[:], st=:surface, reuse=false))

# check
log_dump_config();
log_dump_prob();

@finalize()
