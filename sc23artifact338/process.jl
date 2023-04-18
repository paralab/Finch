#=
Produces a plot similar to Figure 7.
=#
using Plots
pyplot();

# number of processes/GPUs
np = [512,256,128,64,32,16,8,4,2,1];

# solve time read from file
time = zeros(10);
ideal = zeros(10);
file = open("heat3d-test/solvetime.txt", "r")
lines = readlines(file);
for i=1:10
    time[i] = parse(Float64, lines[i]);
end
close(file);

ideal[10] = time[10];
for i=9:-1:1
    ideal[i] = ideal[i+1]/2;
end

# Create plot
p1 = plot(np, [time ideal], shape=[:circle :circle], ms=[6 3], 
            xticks=([1,2,4,8,16,32,64,128,256,512],["1","2","4","8","16","32","64","128","256","512"]), xscale=:log2, yscale=:log2,
            label=["time" "ideal"])
xlabel!("Number of processes")
ylabel!("Execution time(s)")

png(p1, "results.png")
