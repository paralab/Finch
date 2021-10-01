#=
# 2D DG poisson
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("dgpoisson2d");

@useLog("dgpoisson2dlog")

@domain(2)
@solver(DG)
@functionSpace(LEGENDRE, 1)

n = 20;
@mesh(QUADMESH, n)

@variable(u)
@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)

@coefficient(f, "-2*pi*pi*sin(pi*x)*sin(pi*y)")
@coefficient(alpha, 2*n)
@weakForm(u, "-dot(grad(u), grad(v)) - f*v + surface(dot(ave(grad(u)), jump(v)) + dot(ave(grad(v)), jump(u)) - alpha*dot(jump(u), jump(v)))")

solve(u);

# # solution is stored in the variable's "values"
using Plots
pyplot();
display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st=:surface))

@finalize()
