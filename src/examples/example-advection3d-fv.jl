#=
3D advection using FV
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("FVadvection3d");

useLog("FVadvection3dlog")

# Configuration setup
domain(3)
solverType(FV)
timeStepper(RK4)

# Mesh
n = 10 # number of elements in each direction
mesh(HEXMESH, elsperdim=n, bids=2)

# Variables and BCs
u = variable("u", location=CELL)
boundary(u, 1, FLUX, "sin(pi*y)*sin(pi*z)*sin(2*pi*t)^2") # x=0
boundary(u, 2, NO_BC) # other

# Time interval and initial condition
T = 0.7;
timeInterval(T)
initial(u, "0")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
coefficient("a", [1,0,0], type=VECTOR) # advection velocity
# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
flux(u, "upwind(a,u)") 
# Note that there is no source() for this problem

#exportCode("fvad2dcode") # uncomment to export generated code to a file
#importCode("fvad2dcode") # uncomment to import code from a file

solve(u)

finalize_finch()
