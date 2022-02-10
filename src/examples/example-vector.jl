#=
# Tests vector unknown capability with a simple Poisson-like problem.
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

init_finch("vector");

# Try making an optional log
useLog("vectorlog", level=3)

# Set up the configuration
domain(2)
functionSpace(order = 3)

# Specify the problem
mesh(QUADMESH, elsperdim=20, bids=3)

u = variable("u", type=VECTOR)
p = variable("p")

testSymbol("v", type=VECTOR)
testSymbol("w")

# boundary(u, 2, DIRICHLET, [0.1, 0.1])
# boundary(p, 2, DIRICHLET, 0.1)
boundary(u, 1, NEUMANN, ["-pi*sin(2*pi*y)", "-2*pi*sin(3*pi*y)"])
boundary(p, 1, NEUMANN, "-pi*sin(pi*y)")
boundary(u, 2, DIRICHLET, [0, 0])
boundary(p, 2, DIRICHLET, 0)
boundary(u, 3, DIRICHLET, [0, 0])
boundary(p, 3, DIRICHLET, 0)

# Write the weak form
coefficient("f", ["-25*pi*pi*sin(pi*x)*sin(2*pi*y)", "-65*pi*pi*sin(2*pi*x)*sin(3*pi*y)"], type=VECTOR)
coefficient("g", "-2*pi*pi*sin(pi*x)*sin(pi*y)")
coefficient("a", 5)

weakForm([u,p], ["-a*inner(grad(u), grad(v)) - dot(f,v)", "-dot(grad(p), grad(w)) - g*w"])

solve([u,p]);

# # exact solution is [sin(pi*x)*sin(2*pi*y), sin(2*pi*x)*sin(3*pi*y)]
# # check error
erroru = zeros(size(u.values));
maxerru = 0
maxerrp = 0
exactu1(x,y) = sin(pi*x)*sin(2*pi*y);
exactu2(x,y) = sin(2*pi*x)*sin(3*pi*y);
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

using Plots
pyplot();
display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], erroru[1,:], st=:surface))

# check
# log_dump_config();
# log_dump_prob();

output_values([u,p], "vector2d", format="vtk");

finalize_finch();
