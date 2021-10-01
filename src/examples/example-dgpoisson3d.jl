#=
# 2D DG poisson
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("dgpoisson3d");

@useLog("dgpoisson3dlog")

@domain(3)
@solver(DG)
@functionSpace(LEGENDRE, 1)

n = 12;
@mesh(HEXMESH, n)

@variable(u)
@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)

@coefficient(f, "-3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@coefficient(beta, 100000*n)
@weakForm(u, "-dot(grad(u), grad(v)) - f*v + surface(dot(ave(grad(u)), jump(v)) + dot(ave(grad(v)), jump(u)) - beta*dot(jump(u), jump(v)))")

solve(u);

# # solution is stored in the variable's "values"
# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st=:surface))

# exact solution is sin(pi*x)*sin(pi*y)*sin(pi*z)
# check error
maxerr = 0;
exact(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z);

for i=1:size(Finch.grid_data.allnodes,2)
    x = Finch.grid_data.allnodes[1,i];
    y = Finch.grid_data.allnodes[2,i];
    z = Finch.grid_data.allnodes[3,i];
    err = abs(u.values[i] - exact(x,y,z));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

@finalize()
