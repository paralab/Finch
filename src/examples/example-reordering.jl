#=
# Test out different node orderings
=#
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
n = 5;
ord = 4;
gd = n*ord + 1; # grid size in each dimension
griddim = [gd,gd,gd];
te = 1; # elemental loop order is also tiled with this width
tds = [ord*te+1]; # tile size in each dimension

init_finch("reordering");
@useLog("reorderinglog")

#cachesim(true);

@domain(3)
@functionSpace(LEGENDRE, ord)

#@matrixFree(200,1e-6)

@mesh(HEXMESH, n)

@variable(u)
@variable(q1)
@variable(q2)
@variable(q3)
@variable(q4)
@variable(q5)
@variable(q6)

@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)
@boundary(q1, 1, DIRICHLET, 0)
@boundary(q2, 1, DIRICHLET, 0)
@boundary(q3, 1, DIRICHLET, 0)
@boundary(q4, 1, DIRICHLET, 0)
@boundary(q5, 1, DIRICHLET, 0)
@boundary(q6, 1, DIRICHLET, 0)

@coefficient(f, "-3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@weakForm(u, "-dot(grad(u), grad(v)) - f*v")

@coefficient(g, "(-3*pi*pi + 5)*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@weakForm([q1,q2,q3,q4,q5,q6], ["-dot(grad(q1), grad(v)) + q2*v + q3*v + q4*v + q5*v + q6*v - g*v",
                                "-dot(grad(q2), grad(v)) + q3*v + q4*v + q5*v + q6*v + q1*v - g*v",
                                "-dot(grad(q3), grad(v)) + q4*v + q5*v + q6*v + q1*v + q2*v - g*v",
                                "-dot(grad(q4), grad(v)) + q5*v + q6*v + q1*v + q2*v + q3*v - g*v",
                                "-dot(grad(q5), grad(v)) + q6*v + q1*v + q2*v + q3*v + q4*v - g*v",
                                "-dot(grad(q6), grad(v)) + q1*v + q2*v + q3*v + q4*v + q5*v - g*v"])

#
#using Plots
#pyplot();

# First do u (one DOF)
times = 4; # average times times
timings = zeros(7);
timings6 = zeros(7);
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
timings[1] = regtime;
println("Lex. time = "*string(regtime)*"sec.");
#display(spy(Finch.CGSolver.Amat, legend=nothing, reuse=false));

# Tiled ordering
for ti=1:length(tds)
    td = tds[ti];
    @mesh(HEXMESH, n)
    tiled_nodes(griddim,(td,td,td));
    tiled_elements((n,n,n),(te,te,te));
    
    #solve(u);
    tiletime = Base.Libc.time();
    for iter=1:times
        solve(u);
    end
    tiletime = Base.Libc.time() - tiletime;
    tiletime /= times;
    timings[2] = tiletime;
    println("tiled("*string(td)*") = "*string(tiletime)*"sec.");
end
#display(spy(Finch.CGSolver.Amat, legend=nothing, reuse=false));

# morton ordering
@mesh(HEXMESH, n)
morton_nodes(griddim);
morton_elements([n,n,n]);

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve(u);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;
timings[3] = tiletime;
println("morton = "*string(tiletime)*"sec.");
#display(spy(Finch.CGSolver.Amat, legend=nothing, reuse=false));

# hilbert ordering
@mesh(HEXMESH, n)
hilbert_nodes(griddim);
hilbert_elements([n,n,n]);

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve(u);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;
timings[4] = tiletime;
println("hilbert = "*string(tiletime)*"sec.");


# Lex. + ef ordering
@mesh(HEXMESH, n)
ef_nodes();
regtime = Base.Libc.time();
for iter=1:times
    solve(u);
end
regtime = Base.Libc.time() - regtime;
regtime /= times;
timings[5] = regtime;
println("Lex + ef time = "*string(regtime)*"sec.");
#display(spy(Finch.CGSolver.Amat, legend=nothing, reuse=false));

# morton + ef ordering
@mesh(HEXMESH, n)
morton_elements([n,n,n]);
ef_nodes();

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve(u);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;
timings[6] = tiletime;
println("morton + ef = "*string(tiletime)*"sec.");
#display(spy(Finch.CGSolver.Amat, legend=nothing, reuse=false));

# hilbert ordering
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
timings[7] = tiletime;
println("hilbert + ef = "*string(tiletime)*"sec.");

# using Plots
# pyplot();
# #display(plot(timings, marker=:circle, reuse=false))
# labels = ["lex", "tiled", "hilb", "mort", "lex+ef", "hilb+ef", "mort+ef"];
# display(bar(labels, timings, legend=false, reuse=false))

###############################################################################
# then do qn (6 DOF)
times = 1;
println("Six DOF per node, averaged "*string(times)*" times")

@mesh(HEXMESH, n)

#warm up
solve([q1,q2,q3,q4,q5,q6]);

# Lex. ordering
regtime = Base.Libc.time();
for iter=1:times
    solve([q1,q2,q3,q4,q5,q6]);
end
regtime = Base.Libc.time() - regtime;
regtime /= times;
timings6[1] = regtime;
println("Lex. time = "*string(regtime)*"sec.");

# Tiled ordering
for ti=1:length(tds)
    td = tds[ti];
    @mesh(HEXMESH, n)
    tiled_nodes((gd,gd,gd),(td,td,td));
    #tiled_elements((n,n,n),(2,2,2));
    
    #solve([q1,q2,q3,q4,q5,q6]);
    tiletime = Base.Libc.time();
    for iter=1:times
        solve([q1,q2,q3,q4,q5,q6]);
    end
    tiletime = Base.Libc.time() - tiletime;
    tiletime /= times;
    timings6[2] = tiletime;
    println("tiled("*string(td)*") = "*string(tiletime)*"sec.");
end

# morton ordering
@mesh(HEXMESH, n)
morton_nodes(griddim);
morton_elements([n,n,n]);

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve([q1,q2,q3,q4,q5,q6]);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;
timings6[3] = tiletime;
println("morton = "*string(tiletime)*"sec.");
#display(spy(Finch.CGSolver.Amat, legend=nothing, reuse=false));

# hilbert ordering
@mesh(HEXMESH, n)
hilbert_nodes(griddim);
hilbert_elements([n,n,n]);

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve([q1,q2,q3,q4,q5,q6]);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;
timings6[4] = tiletime;
println("hilbert = "*string(tiletime)*"sec.");

# Lex. + ef ordering
@mesh(HEXMESH, n)
ef_nodes();
regtime = Base.Libc.time();
for iter=1:times
    solve([q1,q2,q3,q4,q5,q6]);
end
regtime = Base.Libc.time() - regtime;
regtime /= times;
timings6[5] = regtime;
println("Lex + ef time = "*string(regtime)*"sec.");
#display(spy(Finch.CGSolver.Amat, legend=nothing, reuse=false));

# morton + ef ordering
@mesh(HEXMESH, n)
morton_elements([n,n,n]);
ef_nodes();

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve([q1,q2,q3,q4,q5,q6]);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;
timings6[6] = tiletime;
println("morton + ef = "*string(tiletime)*"sec.");
#display(spy(Finch.CGSolver.Amat, legend=nothing, reuse=false));

# hilbert ordering
@mesh(HEXMESH, n)
hilbert_elements([n,n,n]);
ef_nodes();

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve([q1,q2,q3,q4,q5,q6]);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;
timings6[7] = tiletime;
println("hilbert + ef = "*string(tiletime)*"sec.");

using Plots
pyplot();
#display(plot(timings, marker=:circle, reuse=false))
labels = ["lex", "tiled", "hilb", "mort", "lex+ef", "hilb+ef", "mort+ef"];
display(bar(labels, timings, legend=false, reuse=false))
display(bar(labels, timings6, legend=false, reuse=false))

# exact solution is sin(pi*x)*sin(pi*y)*sin(pi*z)
# check error
maxerr = 0;
maxerrq = 0;
exact(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z);
exactq(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z);

for i=1:size(Finch.grid_data.allnodes,2)
    x = Finch.grid_data.allnodes[1,i];
    y = Finch.grid_data.allnodes[2,i];
    z = Finch.grid_data.allnodes[3,i];
    err = abs(u.values[i] - exact(x,y,z));
    errq = abs(q1.values[i] - exactq(x,y,z)) + abs(q2.values[i] - exactq(x,y,z)) + abs(q3.values[i] - exactq(x,y,z)) + 
            abs(q4.values[i] - exactq(x,y,z)) + abs(q5.values[i] - exactq(x,y,z)) + abs(q6.values[i] - exactq(x,y,z));
    global maxerr;
    global maxerrq;
    maxerr = max(err,maxerr);
    maxerrq = max(errq,maxerrq);
end
println("max error for u = "*string(maxerr));
println("max error for q = "*string(maxerrq));

@finalize()
