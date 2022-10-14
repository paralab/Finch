#=
2D advection using higher order FV or structured or unstructured mesh
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("FVadvection2d");

useLog("FVadvection2dlog", level=3)

# Configuration setup
domain(2)
solverType(FV)

use_tri=false;
if use_tri
    timeStepper(RK4)
    # n = 30
    # mesh(QUADMESH, elsperdim=n, bids=4)
    
    mesh("src/examples/utriangle.msh")  # Using triangles
    
    add_boundary_ID(2, (x,y) -> (x >= 0.1));
    add_boundary_ID(3, (x,y) -> (y <= 0));
    add_boundary_ID(4, (x,y) -> (y >= 0.3));
    finiteVolumeOrder(2);
else
    timeStepper(EULER_EXPLICIT)
    setSteps(0.000278, 1800)
    n = 30
    mesh(QUADMESH, elsperdim=n, bids=4)
    # finiteVolumeOrder(2);
end

# Variables and BCs
u = variable("u", location=CELL)
# boundary(u, 1, FLUX, "(abs(y-0.2) < 0.11) ? sin(2*pi*t)^2 : 0") # x=0
if use_tri
    boundary(u, 1, FLUX, "(abs(y-0.06) < 0.033 && sin(6*pi*t)>0) ? 1 : 0") # x=0
else
    boundary(u, 1, FLUX, "(abs(y-0.2) < 0.11 && sin(6*pi*t)>0) ? 1 : 0") # x=0
end
boundary(u, 2, NO_BC) # x=1
boundary(u, 3, NO_BC) # y=0
boundary(u, 4, NO_BC) # y=1

# Time interval and initial condition
T = 0.5;
timeInterval(T)
initial(u, "0")

# Coefficients
if use_tri
    coefficient("a", ["0.1*cos(pi*x/2/0.1)","0.3*sin(pi*x/2/0.1)"], type=VECTOR) # advection velocity
else
    coefficient("a", ["cos(pi*x/2)","sin(pi*x/2)"], type=VECTOR) # advection velocity
end
coefficient("s", ["sin(pi*x)^4 * sin(pi*y)^4"]) # source

# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
conservationForm(u, "0.1 * s + surface(upwind(a,u))");

exportCode("fvad2dcodeout") # uncomment to export generated code to a file
# importCode("fvad2dcode") # uncomment to import code from a file

solve(u)

# outputValues(u, "fvad2d", format="vtk");

finalizeFinch()

##### Uncomment below to plot

# xy = Finch.fv_info.cellCenters

# using Plots
# pyplot();
# display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
