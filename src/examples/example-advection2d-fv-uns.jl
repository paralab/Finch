#=
2D advection on unstructured grid using FV
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

useLog("FVadvection2dlog")

# Configuration setup
domain(2)
solverType(FV)
timeStepper(RK4)

# Mesh
# This rectangle covers [0, 0.1]x[0, 0.3]
# Uncomment the desired mesh.
mesh("utriangle.msh")  # Using triangles
#mesh("uquad.msh")     # Using quads

add_boundary_ID(2, (x,y) -> (x >= 0.1));
add_boundary_ID(3, (x,y) -> (y <= 0));
add_boundary_ID(4, (x,y) -> (y >= 0.3));

# Variables and BCs
u = variable("u", location=CELL)
boundary(u, 1, FLUX, 0) # x=0
boundary(u, 2, NO_BC) # x=1
boundary(u, 3, FLUX, "sin(10*pi*x)*sin(2*pi*t)^2") # y=0
boundary(u, 4, NO_BC) # y=1

# Time interval and initial condition
T = 0.7;
timeInterval(T)
initial(u, "0")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
coefficient("a", [0.1, 1], type=VECTOR) # advection velocity
# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
flux(u, "upwind(a,u)") 
# Note that there is no source() for this problem

#exportCode("fvad2dcode") # uncomment to export generated code to a file
#importCode("fvad2dcode") # uncomment to import code from a file

solve(u)

finalize_finch()

##### Uncomment below to plot

# xy = Finch.fv_info.cellCenters

# using Plots
# pyplot();
# display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
