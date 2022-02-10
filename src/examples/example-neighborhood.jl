#=
This example demonstrates the use of callback functions working with a neighborhood struct.
It is a 2D advection-reaction equation.
It is solved in two ways: 
 - u1 using a symbolic upwinding operator for the flux.
 - u2 using a callback function for the flux that uses the neighborhood of a face to do the same upwinding.
 
The execution time is also compared, but doesn't account for compilation/first run time.
Also, the number of flops using the symbolic version is much higher because it is actually doing
    1/2*dot(a,n)*(side1+side2) + (1-alpha)/2*abs(dot(a,n))*(side1-side2)
with a parameter alpha rather than
    dot(a,n) > 0 ? dot(a,n)*side1 : dot(a,n)*side2
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("FVneighborhood");

useLog("FVneighborhoodlog", level=3)

# Configuration setup
domain(2)
solverType(FV)

timeStepper(RK4)
n = 20
mesh(QUADMESH, elsperdim=n, bids=4)
# finiteVolumeOrder(2);

# Variables and BCs
u1 = variable("u1", location=CELL) # will use symbolic upwind operator
u2 = variable("u2", location=CELL) # will use callback with neighborhood

# This boundary provides a pulsing input through a small window on one side
boundary(u1, 1, FLUX, "(abs(y-0.2) < 0.11 && sin(6*pi*t)>0) ? 1 : 0") # x=0
boundary(u1, 2, NO_BC) # x=1
boundary(u1, 3, NO_BC) # y=0
boundary(u1, 4, NO_BC) # y=1
boundary(u2, 1, FLUX, "(abs(y-0.2) < 0.11 && sin(6*pi*t)>0) ? 1 : 0") # x=0
boundary(u2, 2, NO_BC) # x=1
boundary(u2, 3, NO_BC) # y=0
boundary(u2, 4, NO_BC) # y=1

# Time interval and initial condition
T = 1;
timeInterval(T)
initial(u1, "0")
initial(u2, "0")

coefficient("a", ["cos(pi*x/2)","sin(pi*x/2)"], type=VECTOR) # advection velocity
coefficient("s", "sin(pi*x)^4 * sin(pi*y)^4") # source

# Here is the callback function for the flux that simply does first order upwinding.
@callbackFunction(
    function flux_fun(neighborhood, a, normal)
        # neighborhood is a struct containing:
        # - neighborhood.cells -> an array of two arrays of integers like [[L1, L2,...], [R1, R2,...]] 
        #                         where L1 is the nearest neighbor on the left, L2 is next-nearest, etc.
        # - neighborhood.centers -> similar, but holding cell center coordinates like [ [L1x L2x ...  , [R1x R2x ...
        #                                        *note these are 2d arrays->             L1y L2y ...]    R1y R2y ...] ]
        # - neighborhood.values -> holds the values of the requested variable(s) like [[u_L1, u_L2,...], [u_R1, u_R2,...]]
        #                          *vector valued variables or multiple variables will have columns of values
        
        a_dot_n = a[1]*normal[1] + a[2]*normal[2]; # a and normal are 2D vectors
        if a_dot_n > 0
            result = neighborhood.values[1][1] * a_dot_n; # use the first value on the left: values[1][1]
        else
            result = neighborhood.values[2][1] * a_dot_n; # use the first value on the right: values[2][1]
        end
        
        return result;
    end
)

# The PDE expression input
# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
flux(u1, "upwind(a, u1)") 
source(u1, "0.1 * s")

flux(u2, "flux_fun(neighborhood(u2), a, normal())") 
source(u2, "0.1 * s")

exportCode("fvneighborhoodcode") # uncomment to export generated code to a file
#importCode("fvneighborhoodcode") # uncomment to import code from a file

t1 = @elapsed(solve(u1))
t2 = @elapsed(solve(u2))

println("time for u1: "*string(t1))
println("time for u2: "*string(t2))

finalize_finch()

# uncomment to output to file
# output_values([u1,u2], "fvad2d", format="vtk");

##### Uncomment below to plot

# xy = Finch.fv_info.cellCenters

# using Plots
# pyplot();
# p1 = plot(xy[1,:], xy[2,:], u1.values[:], st=:surface)
# p2 = plot(xy[1,:], xy[2,:], u2.values[:], st=:surface)
# display(plot(p1, p2, layout=2))
