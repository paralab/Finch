#=
# 1D Poisson, Dirichlet bc
# CG, Linear element
# Simplest test possible
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("poisson2d");

# Try making an optional log
@useLog("poisson2dlog")

# Set up the configuration (order doesn't matter)
@domain(2, SQUARE, UNIFORM_GRID)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 1)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(QUADMESH, [20,10], 1, [0,2,0,1])                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)

@testSymbol(v)                    # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
@coefficient(f, "-2*pi*pi*sin(pi*x)*sin(pi*y)")
@weakForm(u, "-dot(grad(u), grad(v)) - f*v")

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

# solution is stored in the variable's "values"
# using Plots
# pyplot();
# display(plot(Finch.grid_data.allnodes[1,:], Finch.grid_data.allnodes[2,:], u.values[:], st=:surface))

# check
log_dump_config();
log_dump_prob();

@finalize()
