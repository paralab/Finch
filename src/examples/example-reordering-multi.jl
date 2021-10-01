#=
# Test out different node orderings
=#
if !@isdefined(Finch)
    @everywhere include("../Finch.jl");
    @everywhere using .Finch
end

@everywhere include("multi-setup.jl")

println("Six DOF per node, averaged "*string(times)*" times")

@everywhere include("multi-hilb.jl")
@everywhere include("multi-lex.jl")
@everywhere include("multi-rand.jl")

using Plots
pyplot();
labels = ["Hilbert", "Lex.", "Random"];
display(bar(labels, timings, legend=false, reuse=false))

@finalize()
