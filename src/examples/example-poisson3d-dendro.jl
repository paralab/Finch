#=
# 3D Poisson, Dirichlet bc
# Uses Dendro imported as a custom gen target
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("poissondendro");

useLog("poissondendrolog", level=3)

# default values (max_depth=6, wavelet_tol = 0.1, partition_tol = 0.3, solve_tol = 1e-6, max_iters = 100)
generateFor(DENDRO, params=(7, 0.01, 0.3, 0.000001, 1000))

domain(3)
functionSpace(order=4)

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, "0")

f = coefficient("f", "-30*pi*pi*sin(5*pi*x)*sin(2*pi*y)*sin(pi*z)")
weakForm(u, "-dot(grad(u),grad(v)) - f*v")

build_octree_with(f);
solve(u);

finalize_finch()
