#=
# Tests multiple interdependent scalar variables.
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("multivar");

# Try making an optional log
@useLog("multivarlog")

# Set up the configuration
@domain(1)
@functionSpace(LEGENDRE, 3)

@matrixFree(200,1e-6)

# Build a simple mesh
@mesh(LINEMESH, 30)

# Variables
@variable(u)
@variable(q)

@testSymbol(v)

#Boundary conditions
@boundary(u, 1, DIRICHLET, "0")
@boundary(q, 1, DIRICHLET, "0")

#Equations
@coefficient(f, "-4*pi*pi*sin(2*pi*x)")
@coefficient(g, "-9*pi*pi*sin(3*pi*x) + 0.5*sin(2*pi*x)")

@weakForm([u, q], ["-dot(grad(u), grad(v)) - f*v", "-dot(grad(q), grad(v)) + 0.5*u*v - g*v"])

solve([u, q]);

# exact solution is u = sin(2*pi*x), q = sin(3*pi*x)
# check error
maxerru = 0;
maxerrq = 0;
exactu(x) = sin(2*pi*x);
exactq(x) = sin(3*pi*x);

for i=1:size(Finch.grid_data.allnodes,1)
    x = Finch.grid_data.allnodes[i,1];
    erru = abs(u.values[i] - exactu(x));
    errq = abs(q.values[i] - exactq(x));
    global maxerru;
    global maxerrq;
    maxerru = max(erru,maxerru);
    maxerrq = max(errq,maxerrq);
end
println("max error(u) = "*string(maxerru));
println("max error(q) = "*string(maxerrq));

# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[:], u.values[:], markershape=:circle, reuse=false))
# display(plot(Finch.grid_data.allnodes[:], q.values[:], markershape=:circle, reuse=false))

# check
log_dump_config();
log_dump_prob();

@finalize()
