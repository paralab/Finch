#=
Produces a plot similar to Figure 7.
=#
using Plots
pyplot();

# number of processes/GPUs
cpu_np = [55,40,20,10,5,2,1];
gpu_np = [8,4,2,1];

# solve time read from file
cpu_time = zeros(7);
gpu_time = zeros(4);
file = open("solvetime.txt", "r")
lines = readlines(file);
for i=1:7
    cpu_time[i] = parse(Float64, lines[i])
end
for i=1:4
    gpu_time[i] = parse(Float64, lines[i+7])
end
close(file);

# Create plot
p1 = plot([cpu_np, gpu_np], [cpu_time, gpu_time], shape=[:circle :utriangle], ms=[6 6], 
            xticks=([1,2,4,8,16,32,64],["1","2","4","8","16","32","64"]), xscale=:log2, yscale=:log2,
            label=["CPU-only" "GPU"])
xlabel!("Number of processes/GPUs")
ylabel!("Solve time(s)")

png(p1, "results.png")
