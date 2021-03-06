

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("FVheat2d");

useLog("FVheat2dlog")

# Configuration setup
domain(2)
solverType(FV)
timeStepper(RK4)

# Mesh
n = 20 # number of elements
mesh(QUADMESH, elsperdim=n)

# Variables and BCs
u = variable("u", location=CELL)
boundary(u, 1, FLUX, "0")

# Time interval and initial condition
T = 0.1;
timeInterval(T)
initial(u, "(sin(pi*x)*sin(2*pi*y))^4")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)/A) = int(S dx) - int(F.n ds)
coefficient("D", 0.1) # Diffusion rate
flux(u, "-D*dot(grad(u),normal())") 

printLatex(u)

#exportCode("fvheat2dcode") # uncomment to export generated code to a file
#importCode("fvheat2dcode") # uncomment to import code from a file

solve(u)

# output result to file
output_values(u, "fvheat2d", format="vtk");

@finalize()

##### Uncomment below to plot
# xy = Finch.fv_info.cellCenters

# using Plots
# pyplot();
# display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
