#=
# 2D advection.
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("advection2d");

@useLog("advection2dlog")

@domain(2)
@solver(DG)
@functionSpace(LEGENDRE, 3)
@stepper(CRANK_NICHOLSON, 0.1)

@mesh(QUADMESH, 8)

@variable(u)
@testSymbol(v)

T = 2;
@timeInterval(T)
@initial(u, "0")

### Add boundary regions and set BCs ###############
@addBoundaryID(2, x <= 0 || y <= 0)
@addBoundaryID(3, x < 0.2 && y < 0.2)
@boundary(u, 1, NO_BC)
@boundary(u, 2, DIRICHLET, 0)
@boundary(u, 3, DIRICHLET, "sin(pi*(t - 0.7*(x+y)))^2")
####################################################

@coefficient(a, VECTOR, [0.5, 0.5])
@coefficient(beta, 0)
@weakForm(u, "Dt(u*v) + a*grad(u)*v + surface(-dot(a, normal())*u*v + dot(a, normal())*ave(u)*v - 0.5*(1-beta)*dot(a, jump(u))*v)")

solve(u);

# Uncomment to plot
# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st=:surface, reuse=false))

@finalize()
