#=
# 1D advection-diffusion using FV with indexed equations:
#   flux_ij = a_i.n * u_ij + D_j * grad(u_ij).n
# Indices cover a range of advection speeds and diffusion rates.
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("addiff1dindexed");

useLog("addiff1dindexedlog", level=3)

# Configuration setup
domain(1)
solverType(FV)
timeStepper(EULER_EXPLICIT)

# Mesh
n = 60 # number of elements
mesh(LINEMESH, elsperdim=n, bids=2)

# Set up indices and indexed values
nspeeds = 5
ndiffs = 4
speed = index("speed", range = [1,nspeeds])
diff = index("diff", range = [1,ndiffs])

speeds = zeros(nspeeds);
diffs = zeros(ndiffs);
for j=1:ndiffs
    diffs[j] = 0.03 * (j-1)/(ndiffs-1);
end
for i=1:nspeeds
    speeds[i] = 0.1 + i/nspeeds;
end

# Variables and BCs
u = variable("u", type=VAR_ARRAY, location=CELL, index = [speed, diff])
boundary(u, 1, FLUX, 0)
boundary(u, 2, NO_BC)

# Time interval and initial condition
T = 0.3;
timeInterval(T)
initial(u, "x<=0.3&&x>0.1 ? 1 : 0")

# Arrange the loops as desired. Outermost first.
# assemblyLoops(["elements", speed, diff])
assemblyLoops([speed, diff, "elements"])

# The flux and source terms of the conservation equation
coefficient("a", speeds, type=VAR_ARRAY) # advection velocity
coefficient("d", diffs, type=VAR_ARRAY) # diffusion rate

conservationForm(u, "surface(upwind(a[speed],u[speed, diff]) - d[diff] * dot(grad(u[speed, diff]),normal()))")

exportCode("addiff1dindexedcode") # uncomment to export generated code to a file
# importCode("addiff1dindexedcode") # uncomment to import code from a file

solve(u)

finalizeFinch()

##### Uncomment below to plot solution

# # The exact solution with constant velocity v
# a = 1;
# ini = zeros(n);
# x = Finch.finch_state.fv_info.cellCenters[:]
# for i=1:n
#     if x[i] <= 0.3 && x[i] > 0.1
#         ini[i] = 1
#     end
# end

# using Plots
# pyplot();
# ps = Array{Any,1}(undef, ndiffs);
# for j=1:ndiffs
#     global xs = x;
#     global labels = "initial";
#     global vals = ini;
#     for i=1:nspeeds
#         global xs;
#         global vals;
#         global labels;
#         xs = [xs x];
#         vals = [vals u.values[i+nspeeds*(j-1),:]];
#         labels = [labels "a="*string(speeds[i])*", d="*string(diffs[j])];
#     end
#     ps[j] = plot(xs, vals, markershape=:circle, label=labels, reuse=false)
# end
# display(plot(ps[1], ps[2], ps[3], ps[4], layout=4))