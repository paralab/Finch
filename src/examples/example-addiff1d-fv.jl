#=
1D advection-diffusion using FV
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("FVaddiff1d");

useLog("FVaddiff1dlog", level=3)

# Configuration setup
domain(1)
solverType(FV)
timeStepper(RK4)

# Mesh
n = 40 # number of elements
mesh(LINEMESH, elsperdim=n, bids=2)

# Variables and BCs
u = variable("u", location=CELL)
v = variable("v", location=CELL)
boundary(u, 1, FLUX, 0)
boundary(u, 2, NO_BC)
boundary(v, 1, FLUX, 0)
boundary(v, 2, NO_BC)

# Time interval and initial condition
T = 0.5;
timeInterval(T)
initial(u, "x<0.1 ? 1 : 0")
initial(v, "x<0.1 ? 1 : 0")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
coefficient("a", 1) # advection velocity
# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
coefficient("D", 0.01) # Diffusion rate
flux([u, v], ["upwind(a,u) - D*dot(grad(u),normal())", "upwind(a,v)"]) 
# Note that there is no source() for this problem. 

#@exportCode("fvad1dcode") # uncomment to export generated code to a file
#@importCode("fvad1dcode") # uncomment to import code from a file

solve([u,v])

finalize_finch()

##### Uncomment below to plot

# # The initial condition
# init = zeros(n);
# x = Finch.fv_info.cellCenters[:]
# for i=1:n
#     if x[i] < 0.1
#         init[i] = 1
#     end
# end

# using Plots
# pyplot();
# display(plot([x x x], [init u.values[:] v.values[:]], markershape=:circle, label=["initial" "upwind advection, central diffusion" "advection only"]))
