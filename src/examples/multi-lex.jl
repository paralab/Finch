# Lexicographic
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
timings[2] = regtime;
println("Finished Lex.");
