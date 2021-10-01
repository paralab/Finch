#=
# Test out different node orderings
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
n = 6;
ord = 7;

init_finch("efordering");
@useLog("eforderinglog")

@domain(3)
@functionSpace(LEGENDRE, ord)

@matrixFree(200,1e-6)

@mesh(HEXMESH, n)

@variable(u)

@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)

@coefficient(f, "-3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@weakForm(u, "-dot(grad(u), grad(v)) - f*v")

#
#using Plots
#pyplot();

# exact solution is sin(pi*x)*sin(pi*y)*sin(pi*z)
# check error
function check_error()
    maxerr = 0;
    exact(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z);

    for i=1:size(Finch.grid_data.allnodes,1)
        x = Finch.grid_data.allnodes[i,1];
        y = Finch.grid_data.allnodes[i,2];
        z = Finch.grid_data.allnodes[i,3];
        err = abs(u.values[i] - exact(x,y,z));
        tmp = maxerr;
        maxerr = max(err,maxerr);
    end
    println("max error for u = "*string(maxerr));
end

# First do u (one DOF)
times = 1; # average times times
println("One DOF per node, average of "*string(times)*" runs")
#warm up
solve(u);

# Lex. ordering
regtime = Base.Libc.time();
for iter=1:times
    solve(u);
end
regtime = Base.Libc.time() - regtime;
regtime /= times;
println("Lex. time = "*string(regtime)*"sec.");
check_error();
#display(spy(Finch.CGSolver.Amat, legend=nothing, reuse=false));

# morton ordering
@mesh(HEXMESH, n)
hilbert_elements([n,n,n]);
ef_nodes();

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve(u);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;

println("ef = "*string(tiletime)*"sec.");
check_error()
#display(spy(Finch.CGSolver.Amat, legend=nothing, reuse=false));

@finalize()
