function process_pre_post_args(
    callback_fun::CallbackFunction,
    func_call::GenFunction,
    add_indexers::Vector{Indexer} = Vector{Indexer}(undef,0),
)::Tuple{OrderedDict{String, Any}, OrderedDict{String, String}, Vector{String}}
    args_to_tokens::Dict{String, Any} = OrderedDict()
    variables = OrderedDict{String, String}()
    indexers = String[]
    add_indexers_symbols = [string(indexer.symbol) for indexer in add_indexers]
    # match the function call arguments with the corresponding types
    genfun_expr_args = [string(arg) for arg in func_call.expr.args[2:end]]; # skip the first arg which is the function name
    var_coeff_index_matches = [match(r"^.*(variables|coefficients|indices)\[(\d+)\].*$", arg) for arg in genfun_expr_args];
    for (i, arg) in enumerate(callback_fun.args)
        # If the argument is a variable, coefficient, or index
        if arg in ["x", "y", "z"] # not included
        elseif arg == "normal"
            args_to_tokens[arg] = "mesh_facenormals_gpu"
        elseif var_coeff_index_matches[i] !== nothing
            match = var_coeff_index_matches[i]
            if match.captures[1] == "variables"
                args_to_tokens[arg] = "variables_$(parse(Int, match.captures[2]))_values_gpu"
                variables["variables[$(parse(Int, match.captures[2]))].values"] = args_to_tokens[arg]
            elseif match.captures[1] == "coefficients"
                args_to_tokens[arg] = "coefficients_$(parse(Int, match.captures[2]))_value_gpu"
            elseif match.captures[1] == "indices"
                args_to_tokens["MAX_$(arg)"] = finch_state.indexers[parse(Int, match.captures[2])].range[end]
                pushfirst!(indexers, arg)
            end
        elseif arg in add_indexers_symbols
            args_to_tokens["MAX_$(arg)"] = add_indexers[findfirst(==(arg), add_indexers_symbols)].range[end]
            pushfirst!(indexers, arg)
        else
            # If the argument is a constant or an expression
            call_arg = genfun_expr_args[i]
            int_parsed = tryparse(Int64, call_arg)
            float_parsed = tryparse(Float64, call_arg)
            if int_parsed !== nothing
                args_to_tokens[arg] = int_parsed
            elseif float_parsed !== nothing
                args_to_tokens[arg] = float_parsed
            else 
                args_to_tokens[arg] = call_arg
            end
        end
    end

    return args_to_tokens, variables, indexers
end

function gen_pre_post_kernel_call(
    callback_fun::CallbackFunction,
    args_to_tokens::OrderedDict{String, Any},
    variables::OrderedDict{String, String},
    indexers::Vector{String},
)
    nthreads::Int = 256
    total_threads::Int = 1
    for indexer in indexers
        total_threads *= args_to_tokens["MAX_$(indexer)"]
    end
    nblocks = ceil(Int, total_threads / nthreads)
    call::String = "@CUDA.sync @cuda threads = $(nthreads) blocks = $(nblocks) $(callback_fun.name)_gpu("
    args::Vector{String} = [string(arg) for arg in values(args_to_tokens)]

    call *= join(args, ", ") * ")\n"

    # include copying variables to and from GPU
    copytogpu = ""
    for var in variables.keys
        copytogpu *= "copyto!($(variables[var]), $(string(var)))\n"
    end
    copyfromgpu = ""
    for var in variables.keys
        copyfromgpu *= "copyto!($(string(var)), $(variables[var]))\n"
    end
    call = copytogpu * "\n" * call * "\n" * copyfromgpu

    return call
end

function gen_pre_post_kernel_header(
    calback_fun::CallbackFunction,
    args_to_tokens::OrderedDict{String, Any},
)
    header::String = "function $(calback_fun.name)_gpu("
    args::Vector{String} = [arg for arg in keys(args_to_tokens)]

    header *= join(args, ", ") * ")\n"

    return header
end

function gen_pre_post_kernel_indexers(
    indexers::Vector{String},
)
    indexer_loops::String = "thread_id = threadIdx().x + blockDim().x * (blockIdx().x - 1)\n"
    indexer_maxes = ["MAX_$(indexer)" for indexer in indexers]
    max_thread_id = join(indexer_maxes, " * ")

    indexer_loops *=
"""
if thread_id > $(max_thread_id)
    return
end

"""

    indexer_loops *= "# All index variables\n"
    for (i, indexer) in enumerate(indexers)
        if i == length(indexers)
            thread_ids_per_indexer = "1"
        else
            thread_ids_per_indexer = join(["MAX_$(indexer)" for indexer in indexers[i+1:end]], " * ")
        end

        if i == 1
            indexer_loops *= "$(indexer) = (thread_id - 1) ÷ ($(thread_ids_per_indexer)) + 1\n"
        else
            indexer_loops *= "$(indexer) = ((thread_id - 1) ÷ ($(thread_ids_per_indexer))) % MAX_$(indexer) + 1\n"
        end
    end

    return indexer_loops
end

function gen_pre_post_kernel_body(
    callback_fun::CallbackFunction,
    indexers::Vector{String}
)
    kernel_body::String = gen_pre_post_kernel_indexers(indexers)
    kernel_body *= "\n"
    body_expr::Expr = Meta.parse(callback_fun.body)
    Base.remove_linenums!(body_expr)

    kernel_body *= "# Callback function body\n"
    for part in body_expr.args
        kernel_body *= "$(string(part))\n"
    end

    kernel_body *= "\n"
    kernel_body *= "return nothing\n" # return nothing at the end of the kernel

    return kernel_body
end

function gen_pre_post_kernel(
    calback_fun::CallbackFunction,
    args_to_tokens::OrderedDict{String, Any},
    indexers::Vector{String},
)::String
    kernel_header::String = gen_pre_post_kernel_header(calback_fun, args_to_tokens)
    kernel_body::String = gen_pre_post_kernel_body(calback_fun, indexers)

    kernel_code::String = kernel_header * "\n" * kernel_body * "\nend\n"

    return kernel_code
end