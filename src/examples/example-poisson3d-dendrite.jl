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

ord = 1;

domain(3)
functionSpace(order=ord)

u = variable("u")

testSymbol("v")

boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
coefficient("f", "-3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)")
weakForm(u, "-dot(grad(u), grad(v)) - f*v")

solve(u);

finalizeFinch()
