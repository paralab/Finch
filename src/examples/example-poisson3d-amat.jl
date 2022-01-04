#=
# 3D Poisson, Dirichlet bc
# Uses aMat target
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("p3damat");

useLog("p3damatlog", level=3)

generateFor("../targets/target_amat_cg.jl")

domain(3)
functionSpace(order=2)

mesh(HEXMESH, elsperdim=20)

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, "0")

# Write the weak form 
coefficient("f", "-29*pi*pi*sin(2*pi*x)*sin(3*pi*y)*sin(4*pi*z)")
weakForm(u, "-dot(grad(u),grad(v)) - f*v")

solve(u);

finalize_finch()
