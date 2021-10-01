#=
# 3D CG reaction diffusion with tet mesh
# div(K*grad(u)) -C*u = f
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("RDtet");

@useLog("RDtetlog")

@domain(3, SQUARE, UNSTRUCTURED)
@functionSpace(LEGENDRE, 1)

@mesh("cubefine.msh")

@variable(u)
@testSymbol(v)

@boundary(u, 1, DIRICHLET, "sin(pi*0.5*x)*sin(pi*y)*sin(pi*z)")

@coefficient(K, "x*x + 4")
@coefficient(C, 10)
@coefficient(f, "x*pi*cos(pi*0.5*x)*sin(pi*y)*sin(pi*z) + (-10-(x*x+4)*pi*pi*2.25)*sin(pi*0.5*x)*sin(pi*y)*sin(pi*z)")

# PDE: div(K*grad(u)) -C*u = f
@weakForm(u, "K*dot(grad(u), grad(v)) + C*u*v + f*v")

solve(u);

# exact solution is sin(pi*0.5*x)*sin(pi*y)*sin(pi*z)
# check error
maxerr = 0;
exact(x,y,z) = sin(pi*0.5*x)*sin(pi*y)*sin(pi*z);

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
