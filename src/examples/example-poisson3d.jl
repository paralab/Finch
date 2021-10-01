#=
# 3D Poisson, Dirichlet bc
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("poisson3d");

# Try making an optional log
@useLog("poisson3dlog")

n = 10;
ord = 3;

# Set up the configuration (order doesn't matter)
@domain(3, SQUARE, UNIFORM_GRID)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, ord)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(HEXMESH, n)                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)

@testSymbol(v)                    # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
@coefficient(f, "-3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@weakForm(u, "-dot(grad(u), grad(v)) - f*v")

solve(u);

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

# solution is stored in the variable's "values"
# using Plots
# pyplot();
# N = n*ord+1;
# half = Int(round(N/2));
# range = (N*N*half+1):(N*N*(half+1));
# display(plot(Finch.grid_data.allnodes[1,range], Finch.grid_data.allnodes[2,range], u.values[range], st=:surface))

# check
log_dump_config();
log_dump_prob();

@finalize()
