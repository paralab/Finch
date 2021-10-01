Level1 = true
Level2 = true
Level3 = true
BenchNode = Dict{String,Vector{Float64}}()

function set_benchmark_levels(level1, level2, level3)
    global Level1 = level1 
    global Level2 = level2
    global Level3 = level3
end

#=
macro level_bench(args...)
    return _level_bench(args[1], esc(args[2]))
end

function _level_bench(arg1, arg2)
    if eval(arg1) == true
        t0 = time_ns()
        val = $(arg2)
        println(val)
        t1 = time_ns()
        name = string(arg2.args[1])
        if haskey(BenchNode, name) != true 
            BenchNode[name] = [(t1-t0)/1e9]
        else
            push!(BenchNode[name], (t1-t0)/1e9)
        end
        println("elapsed time in ", name, " ", (t1-t0)/1e9, " seconds")
        return val
    else 
        return arg2
    end 
end
=#

macro level_bench(arg1, arg2)
    name = string(arg2)
    quote
        while false; end
        if $(esc(arg1)) == true
            local t0 = time_ns()
            val = $(esc(arg2))
            local t1 = time_ns()
            name = $(esc(name))
            if haskey(BenchNode, name) != true 
                BenchNode[name] = [(t1-t0)/1e9]
            else
                push!(BenchNode[name], (t1-t0)/1e9)
            end
            #println("elapsed time in ", name, " ", (t1-t0)/1e9, " seconds")
            val
        else 
            $(esc(arg2))
        end 
    end
end