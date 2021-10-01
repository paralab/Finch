if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("FVadvection2d");

useLog("FVadvection2dlog")

# Configuration setup
domain(2)
solverType(FV)
timeStepper(RK4)

# Mesh
n = 15 # number of elements in each direction
mesh(QUADMESH, elsperdim=n, bids=4)

# Variables and BCs
u = variable("u", SCALAR, CELL)
boundary(u, 1, FLUX, "(abs(y-0.2) < 0.11) ? sin(2*pi*t)^2 : 0") # x=0
boundary(u, 2, NO_BC) # x=1
boundary(u, 3, NO_BC) # y=0
boundary(u, 4, NO_BC) # y=1

# Time interval and initial condition
T = 1;
timeInterval(T)
initial(u, "0")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
coefficient("a", ["cos(pi*x/2)","sin(pi*x/2)"], VECTOR) # advection velocity
# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
flux(u, "upwind(a,u)") 
# Note that there is no source() for this problem

exportCode("fvad2dcode") # uncomment to export generated code to a file
#importCode("fvad2dcode") # uncomment to import code from a file

solve(u)

finalize_finch()

##### Uncomment below to plot

xy = Finch.fv_info.cellCenters

using Plots
pyplot();
display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
