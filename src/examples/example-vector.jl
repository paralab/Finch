#=
# Tests vector unknown capability with a simple Poisson-like problem.
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("vector");

# Try making an optional log
useLog("vectorlog")

# Set up the configuration
domain(2)
functionSpace(order = 3)

# Specify the problem
mesh(QUADMESH, elsperdim=20)

u = variable("u", VECTOR)
p = variable("p")

testSymbol("v", VECTOR)
testSymbol("w")

boundary(u, 1, DIRICHLET, [0, 0])
boundary(p, 1, DIRICHLET, 0)

# Write the weak form
coefficient("f", ["-25*pi*pi*sin(pi*x)*sin(2*pi*y)", "-125*pi*pi*sin(3*pi*x)*sin(4*pi*y)"], VECTOR)
coefficient("g", "-2*pi*pi*sin(pi*x)*sin(pi*y)")
coefficient("a", 5)

weakForm([u,p], ["-a*inner(grad(u), grad(v)) - dot(f,v)", "-dot(grad(p), grad(w)) - g*w"])

solve([u,p]);

# # exact solution is [sin(pi*x)*sin(2*pi*y), sin(3*pi*x)*sin(4*pi*y)]
# # check error
erroru = zeros(size(u.values));
maxerru = 0
maxerrp = 0
exactu1(x,y) = sin(pi*x)*sin(2*pi*y);
exactu2(x,y) = sin(3*pi*x)*sin(4*pi*y);
exactp(x,y) = sin(pi*x)*sin(pi*y);

for i=1:size(Finch.grid_data.allnodes,2)
    x = Finch.grid_data.allnodes[1,i];
    y = Finch.grid_data.allnodes[2,i];
    exac = [exactu1(x,y), exactu2(x,y)];
    for j=1:size(Finch.grid_data.allnodes,1)
        erroru[j,i] = u.values[j,i] - exac[j];
        global maxerru;
        maxerru = max(abs(erroru[j,i]),maxerru);
    end
    global maxerrp;
    maxerrp = max(abs(p.values[i] - exactp(x,y)),maxerrp);
end
println("u max error = "*string(maxerru));
println("p max error = "*string(maxerrp));

# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[1,:], st=:surface))

# check
# log_dump_config();
# log_dump_prob();

output_values([u,p], "vector2d", format="vtk");

finalize_finch();
