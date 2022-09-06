#=
# 1D nonlinear problem solved in various ways
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("nltest");

useLog("nltestlog", level=3)

domain(1)

mesh(LINEMESH, elsperdim=100, bids=2)

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, 0)
boundary(u, 2, DIRICHLET, 1)

# Initialize u
# initial(u, 0)
# Finch.eval_initial_conditions();

C = 4 # for convenience
coefficient("C", C)

# Goal - input the nonlinear form
# weakForm(u, "grad(u)*grad(v) + exp(C*u)*v")
# and have it automatically linearized

############################################################################################
# Way 1 - linearized version given as input   f(u) -> f(u_o) + f'(u_o)*(u-u_o)
# uold = variable("uold")
# weakForm(u, "grad(u)*grad(v) + C*exp(C*uold)*u*v - C*uold*exp(C*uold)*v + exp(C*uold)*v")

############################################################################################
# Way 2 - nonlinear term and its derivative are functions and weak form is more general
# uold = variable("uold")
# f(u) = exp(C*u);
# df(u) = C*exp(C*u);
# callbackFunction(f);
# callbackFunction(df);
# weakForm(u, "grad(u)*grad(v) + f(uold)*v + df(uold)*u*v - df(uold)*uold*v")

############################################################################################
# Way 3 - nonlinear term is a function, the derivative is found with AD
# uold = variable("uold")
# using Zygote
# f(u) = exp(C*u);
# df(u) = gradient(f,u)[1];
# callbackFunction(f);
# callbackFunction(df);
# weakForm(u, "grad(u)*grad(v) + f(uold)*v + df(uold)*u*v - df(uold)*uold*v")

############################################################################################
# Ways 1 to 3 require a hand coded iteration
# du = 1;
# tol = 1e-10;
# max_iters = 1000;
# i = 0;
# while i<max_iters && du > tol
#     global i = i+1;
#     solve(u);
#     global du = sqrt(sum((u.values - uold.values).^2));
#     uold.values .= u.values;
#     println("iter: "*string(i)*", du = "*string(du))
# end
# println("total iters = "*string(i)*", du = "*string(du))
############################################################################################

############################################################################################
# Step 4 - Generate the iteration and OLDu automatically by signaling with nonlinear()
# using Zygote
# f(u,c) = exp(c*u);
# df(u,c) = gradient(f,u,c)[1];
# callbackFunction(f);
# callbackFunction(df);
# nonlinear(maxIters=100, relativeTol=1e-8, absoluteTol=1e-8)
# weakForm(u, "grad(u)*grad(v) + f(OLDu,C)*v + df(OLDu,C)*u*v - df(OLDu,C)*OLDu*v")
# solve(u);

############################################################################################
# Step 5 - Detect nonlinear term and linearize automatically
nonlinear(maxIters=100, relativeTol=1e-8, absoluteTol=1e-8)
weakForm(u, "grad(u)*grad(v) + exp(C*u)*v")
solve(u)

# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[:], u.values[:], markershape=:circle, legend=false))

# exportCode("nltestcode");

finalize_finch() # Finish writing and close any files
