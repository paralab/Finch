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
generateFor("Dendrite", baseRefineLevel=4, min=[-1.0,-1.0], max=[1.0,1.0], 
                        refineWhere="level < (sqrt(x*x+y*y) * 8.2) && level < 9",
                        geometries=[Dict([
                                        (:meshFile, "circle.msh"),
                                        (:meshName, "circle"),
                                        (:boundaryTypes, ["sbm"]),
                                        (:bids, [1,2,3,4]),
                                        (:refineLevel, 7)
                                    ])])

domain(2)
functionSpace(order=1)
timeStepper(BDF2)
# Since the mesh is not known, Finch can't automatically choose time steps.
dt = 0.01; nSteps = 100;
setSteps(dt, nSteps) 

u = variable("u")
testSymbol("v")

initial(u, 0.0) # initial condition

# No mesh is created for Dendrite target, so BIDs need to be specified.
# Later ones can overlap. They are tested sequentially and the first one satisfied is used.
addBoundaryID(1, "(x*x + y*y) > (1.0 - 1e-12) && y < 0 && abs(x) < 0.5") # bottom sixth
addBoundaryID(2, "(x*x + y*y) > (1.0 - 1e-12) && y >= 0.0 && x < -0.5") # top left sixth
addBoundaryID(3, "(x*x + y*y) > (1.0 - 1e-12) && y >= 0.0 && x > 0.5") # top right sixth
addBoundaryID(4, "(x*x + y*y) > (1.0 - 1e-12)") # everything else

boundary(u, 1, DIRICHLET, "cos(pi * x)")
boundary(u, 2, DIRICHLET, "sin(pi / 0.866 * y)")
boundary(u, 3, DIRICHLET, "-sin(pi / 0.866 * y)")
boundary(u, 4, DIRICHLET, 0.0)

coefficient("f", "4 * cos(pi/2*sqrt(x*x+y*y))^4")

# Implicit SBM is handled by target
# weakForm(u, "Dt(u*v) + dot(grad(u),grad(v)) - f*v")

# Explicitly written SBM
coefficient("alpha", 200) # penalty
weakForm(u, # volume integral
            "Dt(u*v) + dot(grad(u),grad(v)) - f*v + "*
            # Dirichlet boundary integral
            "dirichletBoundary("*
                "-dot(grad(u), normal()) * v - "*
                "dot(grad(v), normal()) * (u + dot(grad(u), distanceToBoundary()) - dirichletValue()) + "*
                "alpha / elementDiameter() * (u + dot(grad(u), distanceToBoundary()) - dirichletValue()) * (v + dot(grad(v), distanceToBoundary()))"*
            ") + "*
            # Neumann boundary integral
            "neumannBoundary("*
                "dot(normal(), trueNormal()) * (neumannValue() + dot(grad(u), trueNormal())) * v - "*
                "dot(grad(u), normal()) * v"*
            ")")

solve(u);

finalizeFinch()
