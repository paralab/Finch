#=
# 2D DG poisson
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("dgpoisson2d");

useLog("dgpoisson2dlog", level=3)

n = 10;
ord = 2;

@domain(2)
@solver(DG)
@functionSpace(LEGENDRE, ord)

@mesh(QUADMESH, n)

@variable(u)
@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)

@coefficient(f, "-2*pi*pi*sin(pi*x)*sin(pi*y)")
@coefficient(beta, 1000*n)
@weakForm(u, "dot(grad(u), grad(v)) + f*v + surface(dot(ave(grad(u)), jump(v)) + dot(ave(grad(v)), jump(u)) - beta*dot(jump(u), jump(v)))")
# @weakForm(u, "dot(grad(u), grad(v)) + f*v + surface(-beta*dot(jump(u), jump(v)))")

solve(u);

# exact solution is sin(pi*x)*sin(pi*y)
# check error
maxerr = 0;
exact(x,y) = sin(pi*x)*sin(pi*y);

for i=1:size(Finch.grid_data.allnodes,2)
    x = Finch.grid_data.allnodes[1,i];
    y = Finch.grid_data.allnodes[2,i];
    err = abs(u.values[i] - exact(x,y));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));


# # solution is stored in the variable's "values"
# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st=:surface, reuse=false))

@finalize()
