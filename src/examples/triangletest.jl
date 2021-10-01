#=
# test a simple triangle mesh from a file
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("tritest");

@useLog("tritestlog")

@domain(2)
@functionSpace(LEGENDRE, 1)

#@mesh("src/examples/square.msh")
@mesh("src/examples/square2.msh")

@finalize()

# Plot the triangles to check
m = Finch.mesh_data;
using Plots
pyplot();
trianglex = zeros(3, m.nel);
triangley = zeros(3, m.nel);
centers = zeros(2, m.nel);
for i=1:m.nel
    trianglex[:,i] = m.nodes[m.elements[1:3,i],1];
    triangley[:,i] = m.nodes[m.elements[1:3,i],2];
    centers[1,i] = sum(trianglex[:,i])/3;
    centers[2,i] = sum(triangley[:,i])/3;
end
p1 = plot(trianglex, triangley, legend=false)
for i=1:m.nel
    annotate!(centers[1,i], centers[2,i], string(i));
end
display(p1)