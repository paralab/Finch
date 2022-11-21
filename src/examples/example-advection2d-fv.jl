#=
2D advection using structured or unstructured mesh
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("FVadvection2d");

useLog("FVadvection2dlog", level=3)

# Configuration setup
domain(2)
solverType(FV)

timeStepper(EULER_IMPLICIT)

use_unstructured=false;
if use_unstructured
    # Using an unstructured mesh of triangles or irregular quads
    # This is a 0.1 x 0.3 rectangle domain
    mesh("src/examples/utriangle.msh")
    # mesh("src/examples/uquad.msh")
    
    addBoundaryID(2, (x,y) -> (x >= 0.1));
    addBoundaryID(3, (x,y) -> (y <= 0));
    addBoundaryID(4, (x,y) -> (y >= 0.3));
    
else
    # a uniform grid of quads on a 0.1 x 0.3 rectangle domain
    mesh(QUADMESH, elsperdim=[15, 45], bids=4, interval=[0, 0.1, 0, 0.3])
end

# Variables and BCs
u = variable("u", location=CELL)
boundary(u, 1, FLUX, "(abs(y-0.06) < 0.033 && sin(3*pi*t)>0) ? 1 : 0") # x=0
boundary(u, 2, NO_BC) # x=0.1
boundary(u, 3, NO_BC) # y=0
boundary(u, 4, NO_BC) # y=0.3

# Time interval and initial condition
T = 1.3;
timeInterval(T)
initial(u, "0")

# Coefficients
coefficient("a", ["0.1*cos(pi*x/2/0.1)","0.3*sin(pi*x/2/0.1)"], type=VECTOR) # advection velocity
coefficient("s", ["0.1 * sin(pi*x)^4 * sin(pi*y)^4"]) # source

# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
conservationForm(u, "s + surface(upwind(a,u))");

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
