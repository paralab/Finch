#=
# 2D Heat, circle, Dirichlet bc for Dendrite
=#

### If the Finch package has already been added, use this line #########
# using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
##########################################################################

initFinch("heat2d");
useLog("heat2d-dendritelog", level=3)

# Add any other parameters to be added to config.txt here.
# config.txt values will override any configuration set
# elsewhere in this script.
generateFor("Dendrite", refineLevel=4, min=[-1.0,-1.0], max=[1.0,1.0], 
                        geometries=[Dict([
                                        (:meshFile, "circle.msh"),
                                        (:meshName, "circle"),
                                        (:boundaryTypes, ["sbm"]),
                                        (:bids, [1])
                                    ])])

domain(2)
functionSpace(order=1)
timeStepper(BDF2)
# Since the mesh is not known, Finch can't automatically choose time steps.
dt = 0.01; nSteps = 100;
setSteps(dt, nSteps) 

u = variable("u")
testSymbol("v")

initial(u, "cos(pi/2*sqrt(x*x+y*y))") # initial condition

# No mesh is created for Dendrite target, so BIDs need to be specified.
addBoundaryID(1, "(x*x + y*y) > (1.0 - 1e-12) && y < 0.0")
addBoundaryID(2, "(x*x + y*y) > (1.0 - 1e-12) && y >= 0.0")
boundary(u, 1, DIRICHLET, 0.0)
boundary(u, 2, DIRICHLET, "sin(pi/2 * y)")

coefficient("f", "0.1 * pi*pi/4 * cos(pi/2*sqrt(x*x+y*y))")

weakForm(u, "Dt(u*v) + 0.1 * dot(grad(u),grad(v)) - f*v")

solve(u);

finalizeFinch()
