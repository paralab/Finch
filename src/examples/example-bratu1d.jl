#=
# 1D Bratu (or LBG) equation
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("bratu1d");
useLog("bratu1dlog", level=3)

domain(1)
mesh(LINEMESH, elsperdim=50)

nlevels = 20; 

ind = index("ind", range = [1,nlevels])

cvals = zeros(nlevels);
# These numbers set an initial guess that will result in
# the upper branch. These were selected empirically.
ivals = [6,5,4,4,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2]
init_str = fill("*sin(pi*x)", nlevels);
dc = 3.49/nlevels;
for i=1:nlevels
    cvals[i] =  dc * i; # coefficient values for each level
    init_str[i] = string(ivals[i]) * init_str[i]; # initial guess for upper branch
end

u = variable("u", type=VAR_ARRAY, index=ind) # upper branch
l = variable("l", type=VAR_ARRAY, index=ind) # lower branch
testSymbol("v")

boundary(u, 1, DIRICHLET, 0)
boundary(l, 1, DIRICHLET, 0)

initial(u, init_str)
initial(l, "0.1*sin(pi*x)")

coefficient("C", cvals, type=VAR_ARRAY)

############################################################################################
#### Select one ############################################################################
# # linearize automatically using AD for the derivative
# nonlinear(maxIters=100, relativeTol=1e-8, absoluteTol=1e-8, derivative="AD")

# linearize automatically using symbolic derivatives
nonlinear(maxIters=100, relativeTol=1e-8, absoluteTol=1e-8, derivative="symbolic")
############################################################################################

weakForm(u, "-grad(u[ind])*grad(v) + C[ind]*exp(u[ind])*v")
weakForm(l, "-grad(l[ind])*grad(v) + C[ind]*exp(l[ind])*v")

exportCode("bratu1dcode");
# importCode("bratu1dcode");

solve(u)
solve(l)

finalizeFinch()

##### Uncomment below to plot ######

# using Plots
# pyplot();
# x = Finch.grid_data.allnodes[1,:];
# # display(plot([x x x x x x], [u.values[1,:] u.values[2,:] u.values[3,:] l.values[1,:] l.values[2,:] l.values[3,:]], markershape=:circle))#, label=["upper branch" "lower branch"]))
# display(plot([x], [u.values' l.values'], markershape=:circle, legend=false))

# outfile = "bratudat.csv";
# x = Finch.grid_data.allnodes[1,:];
# t = 0;
# open(outfile, "w") do f
#     println(f, "x, u, c")
#     for ci=1:nlevels
#         for ni=1:length(x)
#             println(f, string(x[ni]) * ", " * string(u.values[ci,ni]) * ", " * string(cvals[ci]));
#         end
#         for ni=1:length(x)
#             println(f, string(x[ni]) * ", " * string(l.values[ci,ni]) * ", " * string(cvals[ci]));
#         end
#     end
#     close(f)
# end


