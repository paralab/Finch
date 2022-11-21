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

initFinch("poisson1d");

useLog("poisson1dlog", level=3)

# Set up the configuration
domain(1) # dimension

mesh(LINEMESH, elsperdim=200)   # build uniform LINEMESH with 180 elements

u = variable("u")              # make a scalar variable with symbol u
testSymbol("v")                # sets the symbol for a test function

boundary(u, 1, DIRICHLET, 0)  # boundary condition for BID 1 is Dirichlet with value 0

# Write the weak form 
coefficient("f", "-100*pi*pi*sin(10*pi*x)*sin(pi*x) - pi*pi*sin(10*pi*x)*sin(pi*x) + 20*pi*pi*cos(10*pi*x)*cos(pi*x)")
weakForm(u, "-grad(u)*grad(v) - f*v")

exportCode("poisson1dcode");
# importCode("poisson1dcode");

solve(u);

finalizeFinch()

# exact solution is sin(10*pi*x)*sin(pi*x)
# check error
allerr = zeros(size(Finch.finch_state.grid_data.allnodes,2));

for i=1:size(Finch.finch_state.grid_data.allnodes,2)
    x = Finch.finch_state.grid_data.allnodes[1,i];
    exact = sin(10*pi*x)*sin(pi*x);
    allerr[i] = abs(u.values[i] - exact);
end
maxerr = maximum(abs, allerr);
println("max error = "*string(maxerr));

# ### uncomment below to plot ###

# # solution is stored in the variable's "values"
# using Plots
# pyplot();
# display(plot(Finch.finch_state.grid_data.allnodes[:], u.values[:], markershape=:circle, legend=false))
