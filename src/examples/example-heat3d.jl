#=
# 3D heat eq. Dirichlet bc, CG
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
init_finch("heat3d");

# Try making an optional log
@useLog("heat3dlog")

n = 6;
ord = 2;

# Set up the configuration (order doesn't matter)
@domain(3, SQUARE, UNIFORM_GRID)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, ord)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement
@stepper(EULER_IMPLICIT)            # time stepper (optional second arg is CFL#)

# Specify the problem
@mesh(HEXMESH, n)                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)

@testSymbol(v)                    # sets the symbol for a test function

T = 1;
@timeInterval(T)                    # (start, end) using this sets problem to time dependent
@initial(u, "abs(x-0.5)+abs(y-0.5)+abs(z-0.5) < 0.2 ? 1 : 0")  # initial condition needed if time dependent

@boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
@coefficient(f, "-0.1*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@weakForm(u, "Dt(u*v) + 0.01 * dot(grad(u),grad(v)) - f*v")

solve(u);

# solution is stored in the variable's "values"
#using Plots
#pyplot();
#N = n*ord+1;
#half = Int(round(N/2));
#range = (N*N*half+1):(N*N*(half+1));
#display(plot(Finch.grid_data.allnodes[1,range], Finch.grid_data.allnodes[2,range], u.values[range], st=:surface))

# check
log_dump_config(Finch.config);
log_dump_prob(Finch.prob);

@finalize()
