# hilbert ordering
@mesh(HEXMESH, n)
hilbert_elements([n,n,n]);
ef_nodes();

solve([q1,q2,q3,q4,q5,q6]);

regtime = Base.Libc.time();
for iter=1:times
    solve([q1,q2,q3,q4,q5,q6]);
end
regtime = Base.Libc.time() - regtime;
regtime /= times;
timings[1] = regtime;
println("Finished Hilbet.");