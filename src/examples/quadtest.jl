#=
# test a simple quad mesh from a file
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("quadtest");

@useLog("quadtestlog")

@domain(2)
@functionSpace(LEGENDRE, 1)

# First make a structured mesh
@mesh(QUADMESH, 10)

# write a msh file for it
@outputMesh("quadmesh.msh")

# read it back in as an unstructured mesh
@mesh("quadmesh.msh")

@variable(u)
@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)

@coefficient(f, "-2*pi*pi*sin(pi*x)*sin(pi*y)")
@weakForm(u, "-dot(grad(u), grad(v)) - f*v")

solve(u);

@finalize()

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


# Plot the mesh to check
# m = Finch.mesh_data;
# using Plots
# pyplot();
# quadx = zeros(4, m.nel);
# quady = zeros(4, m.nel);
# centers = zeros(2, m.nel);
# for i=1:m.nel
#     quadx[:,i] = m.nodes[1,m.elements[1:4,i]];
#     quady[:,i] = m.nodes[2,m.elements[1:4,i]];
#     centers[1,i] = sum(quadx[:,i])/4;
#     centers[2,i] = sum(quady[:,i])/4;
# end
# p1 = plot(quadx, quady, legend=false)
# for i=1:m.nel
#     annotate!(centers[1,i], centers[2,i], string(i));
# end
# display(p1)