#=
# 3D heat, Dirichlet bc
# Uses Dendro imported as a custom gen target
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("heat3ddendro");

useLog("heat3ddendrolog", level=3)

# default values (max_depth=6, wavelet_tol = 0.1, partition_tol = 0.3, solve_tol = 1e-6, max_iters = 100)
generateFor("target_dendro_cg.jl", params=(7, 0.01, 0.3, 0.000001, 100))

domain(3)
functionSpace(order=2)
timeStepper(EULER_IMPLICIT)
timeInterval(1)

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, "0")
initial(u, "(sin(pi*x)*sin(pi*y)*sin(pi*z))^4")

# Write the weak form 
coefficient("f", "2*sin(6*pi*x)*sin(pi*x)*sin(6*pi*y)*sin(pi*y)*sin(pi*z)")
weakForm(u, "Dt(u*v) + 0.01 * dot(grad(u),grad(v)) - f*v")

solve(u);

finalize_finch()
