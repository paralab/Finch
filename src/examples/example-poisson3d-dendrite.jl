#=
# 3D Poisson, Dirichlet bc for Dendrite
=#

### If the Finch package has already been added, use this line #########
# using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
##########################################################################

initFinch("poisson3d");
useLog("poisson3d-dendritelog", level=3)

# Add any other parameters to be added to config.txt here.
# config.txt values will override any configuration set
# elsewhere in this script.
generateFor("Dendrite", refineLevel=4)

domain(3)
functionSpace(order=2)

u = variable("u")
testSymbol("v")

# No mesh is created for Dendrite target, so BIDs need to be specified.
addBoundaryID(1, "XMIN || XMAX || YMIN || YMAX || ZMIN || ZMAX")
# Alternatively "(abs(x) < eps()) || (abs(x-1) < eps()) || (abs(y) < eps()) etc..."

boundary(u, 1, DIRICHLET, 0)

coefficient("f", "-14*pi*pi*sin(pi*x)*sin(2*pi*y)*sin(3*pi*z)")

weakForm(u, "-dot(grad(u), grad(v)) - f*v")

solve(u);

finalizeFinch()
