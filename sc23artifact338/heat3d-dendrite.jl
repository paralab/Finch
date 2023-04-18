#=
# 3D Heat, bunny mesh, Dirichlet bc for Dendrite
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("heat3d");
# useLog("heat3d-dendritelog", level=3)

# Add any other parameters to be added to config.txt here.
# config.txt values will override any configuration set
# elsewhere in this script.
generateFor("Dendrite", baseRefineLevel=5, min=[0.0,0.0,0.0], max=[1.0,0.775294,0.985516],
                        geometries=[Dict([
                                        (:meshFile, "bunny.stl"),
                                        (:meshName, "bunny"),
                                        (:boundaryTypes, ["sbm"]),
                                        (:bids, [1]),
                                        (:refineLevel, 9)
                                    ])])

domain(3)
functionSpace(order=1)
timeStepper(BDF2)
# Since the mesh is not known, Finch can't automatically choose time steps.
dt = 0.1; nSteps = 10;
setSteps(dt, nSteps) 

u = variable("u")
testSymbol("v")

initial(u, 0.0) # initial condition

# No mesh is created for Dendrite target, so BIDs need to be specified.
# Later ones can overlap. They are tested sequentially and the first one satisfied is used.
addBoundaryID(1, "true") # everywhere

boundary(u, 1, DIRICHLET, "exp(-z*z / 0.04)") # Hot feet

coefficient("f", 0.0)

# PDE with SBM
coefficient("alpha", 200) # penalty
weakForm(u, # volume integral
            "Dt(u*v) + dot(grad(u),grad(v)) - f*v + "*
            # Dirichlet boundary region
            "dirichletBoundary("*
                "-dot(grad(u), normal()) * v - "*
                "dot(grad(v), normal()) * (u + dot(grad(u), distanceToBoundary()) - dirichletValue()) + "*
                "alpha / elementDiameter() * (u + dot(grad(u), distanceToBoundary()) - dirichletValue()) * (v + dot(grad(v), distanceToBoundary()))"*
            ") + "*
            # Neumann boundary region
            "neumannBoundary("*
                "dot(normal(), trueNormal()) * (neumannValue() + dot(grad(u), trueNormal())) * v - "*
                "dot(grad(u), normal()) * v"*
            ")")

solve(u);

finalizeFinch()
