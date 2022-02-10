#=
# 1D Poisson, Dirichlet bc
# CG, Linear element
# Simplest test possible
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("poisson1d");

# Optionally generate a log
useLog("poisson1dlog")

# Set up the configuration (order doesn't matter)
domain(1)                      # dimension
functionSpace(order=3)         # basis function polynomial order

# Specify the problem (mesh comes first)
mesh(LINEMESH, elsperdim=20)   # build uniform LINEMESH with 20 elements

u = variable("u")              # make a scalar variable with symbol u
testSymbol("v")                # sets the symbol for a test function

boundary(u, 1, DIRICHLET, 0)  # boundary condition for BID 1 is Dirichlet with value 0

# Write the weak form 
coefficient("f", "-100*pi*pi*sin(10*pi*x)*sin(pi*x) - pi*pi*sin(10*pi*x)*sin(pi*x) + 20*pi*pi*cos(10*pi*x)*cos(pi*x)")
weakForm(u, "-grad(u)*grad(v) - f*v")

solve(u);

# exact solution is sin(10*pi*x)*sin(pi*x)
# check error
allerr = zeros(size(Finch.grid_data.allnodes,2));
exact(x) = sin(10*pi*x)*sin(pi*x);

for i=1:size(Finch.grid_data.allnodes,2)
    x = Finch.grid_data.allnodes[1,i];
    err = abs(u.values[i] - exact(x));
    allerr[i] = err;
end
maxerr = maximum(abs, allerr);
println("max error = "*string(maxerr));

# solution is stored in the variable's "values"
# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[:], u.values[:], markershape=:circle, legend=false))

# # Dump things to the log if desired
# log_dump_config();
# log_dump_prob();

finalize_finch() # Finish writing and close any files
