#=
# 3D Poisson, Dirichlet bc
# Uses Dendro imported as a custom gen target
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("dendrotest");

useLog("dendrotestlog", level=3)

# default values (max_depth=6, wavelet_tol = 0.1, partition_tol = 0.3, solve_tol = 1e-6, max_iters = 100)
generateFor("target_dendro_cg.jl", params=(8, 0.001, 0.3, 0.000001, 500))

domain(3)
functionSpace(order=2)

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, "0")

# Write the weak form 
coefficient("k1", "1.1")
coefficient("k2", "0.1")
coefficient("k3", "0.1")
coefficient("k4", "0.1")
coefficient("k5", "0.1")
f = coefficient("f", "-1.1*14*pi*pi*sin(3*pi*x)*sin(2*pi*y)*sin(pi*z) + 0.1*3*pi*cos(3*pi*x)*sin(2*pi*y)*sin(pi*z) + 0.1*2*pi*sin(3*pi*x)*cos(2*pi*y)*sin(pi*z) + 0.1*pi*sin(3*pi*x)*sin(2*pi*y)*cos(pi*z) - 0.1*sin(3*pi*x)*sin(2*pi*y)*sin(pi*z)")
weakForm(u, "-k1*dot(grad(u),grad(v)) + k2*deriv(u,1)*v + k3*deriv(u,2)*v + k4*deriv(u,3)*v - k5*u*v - f*v")

build_octree_with(f);
solve(u);

finalize_finch()
