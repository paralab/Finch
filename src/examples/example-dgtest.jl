#=
# DG test
# 1D Poisson using SIPG formulation for DG
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("test");

@useLog("dgtestlog")

@domain(1)
@solver(DG)
@functionSpace(LEGENDRE, 2)

@mesh(LINEMESH, 10)

@variable(u)

@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)

@coefficient(f, "-pi*pi*sin(pi*x)")
@coefficient(beta, 1)
@weakForm(u, "dot(grad(u),grad(v)) + f*v - surface( ave_normdotgrad(u) * jump(v)) - surface(jump(u) * ave_normdotgrad(v)) + surface(beta*jump(u)*jump(v))")
@weakForm(u, "dot(grad(u),grad(v)) + f*v - surface( ave_normdotgrad(u) * jump(v)) - surface(jump(u) * ave_normdotgrad(v)) + surface(beta*jump(u)*jump(v))")

solve(u);

# solution is stored in the variable's "values"
using Plots
pyplot();
display(plot(Finch.grid_data.allnodes[1,:], u.values[1,:], reuse=false))

@finalize()
