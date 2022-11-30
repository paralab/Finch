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

generateFor("Dendrite")

domain(3)
# functionSpace(order=2)

u = variable("u")
testSymbol("v")

addBoundaryID(1, "XMIN || XMAX || YMIN || YMAX || ZMIN || ZMAX")
# boundary(u, 1, DIRICHLET, 0)
boundary(u, 1, DIRICHLET, "0")

coefficient("f", "-14*pi*pi*sin(pi*x)*sin(2*pi*y)*sin(3*pi*z)")

weakForm(u, "-dot(grad(u), grad(v)) - f*v")

solve(u);

finalizeFinch()
