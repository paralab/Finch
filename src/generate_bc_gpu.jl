function process_bc_args(
    callback_fun::CallbackFunction,
    func_call::GenFunction,
    bi::Int,
    var_id::Int,
)::Tuple{OrderedDict{String, Any}, OrderedDict{String, String}, Vector{String}}
    args_to_tokens::Dict{String, Any} = OrderedDict()
    args_to_expr::Dict{String, String} = OrderedDict()
    indexers = String[]
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
            elseif match.captures[1] == "coefficients"
                args_to_tokens[arg] = "coefficients_$(parse(Int, match.captures[2]))_value_gpu"
            elseif match.captures[1] == "indices"
                args_to_tokens["MAX_$(arg)"] = finch_state.indexers[parse(Int, match.captures[2])].range[end]
                pushfirst!(indexers, arg)
            end
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
                args_to_expr[arg] = call_arg
            end
        end
    end

    args_to_tokens["MAX_fi"] = length(finch_state.grid_data.bdryface[bi])
    pushfirst!(indexers, "fi") # fi is the boundary face index
    
    args_to_tokens["mesh_bdryface"] = "mesh_bdryface_$(bi)_gpu"
    args_to_tokens["mesh_face2element"] = "mesh_face2element_gpu"
    args_to_tokens["mesh_bids"] = "mesh_bids_gpu"
    args_to_tokens["geometric_factors_volume"] = "geometric_factors_volume_gpu"
    args_to_tokens["geometric_factors_area"] = "geometric_factors_area_gpu"
    args_to_tokens["dofs_per_node"] = finch_state.variables[var_id].total_components
    args_to_tokens["boundary_flux"] = "boundary_flux_gpu"
    args_to_tokens["boundary_dof_index"] = "boundary_dof_index_gpu"
    args_to_tokens["global_vector"] = "global_vector_gpu"
    if "t" ∉ keys(args_to_tokens)
        args_to_tokens["t"] = "t"
    end
    args_to_tokens["fv_info_faceCenters"] = "fv_info_faceCenters_gpu"

    return args_to_tokens, args_to_expr, indexers
end

function gen_bc_kernel_call(
    calback_fun::CallbackFunction,
    bi::Int,
    var_id::Int,
    args_to_tokens::OrderedDict{String, Any},
    indexers::Vector{String},
)
    nthreads::Int = 256
    total_threads::Int = 1
    for indexer in indexers
        total_threads *= args_to_tokens["MAX_$(indexer)"]
    end
    nblocks = ceil(Int, total_threads / nthreads)
    call::String = "@CUDA.sync @cuda threads = $(nthreads) blocks = $(nblocks) $(calback_fun.name)_bi_$(bi)_var_$(var_id)_gpu("
    args::Vector{String} = [string(arg) for arg in values(args_to_tokens)]

    call *= join(args, ", ") * ")\n"

    return call
end

function gen_bc_kernel(
    calback_fun::CallbackFunction,
    bi::Int,
    var_id::Int,
    args_to_tokens::OrderedDict{String, Any},
    bc_expr_args::OrderedDict{String, String},
    indexers::Vector{String},
)::String
    kernel_header::String = gen_bc_kernel_header(calback_fun, args_to_tokens, bi, var_id)
    kernel_body::String = gen_bc_kernel_body(calback_fun, indexers, bc_expr_args, bi)

    kernel_code::String = kernel_header * "\n" * kernel_body * "\nend\n"

    return kernel_code
end

function gen_bc_kernel_header(
    calback_fun::CallbackFunction,
    args_to_tokens::OrderedDict{String, Any},
    bi::Int,
    var_id::Int,
)::String
    header::String = "function $(calback_fun.name)_bi_$(bi)_var_$(var_id)_gpu("
    args::Vector{String} = [arg for arg in keys(args_to_tokens)]

    header *= join(args, ", ") * ")\n"

    return header
end

"""
This function essentially turns the body of a loop nest into a CUDA kernel. For example:

for a = 1:max_a
    for b = 1:max_b
        for c = 1:max_c
            result = a + b + c
            result_vector[a * (max_b * max_c) + b * (max_c) + c] = result
        end
    end
end

is transformed to a CUDA kernel similar to this:

a = (thread_id - 1) ÷ (max_b * max_c) + 1
b = ((thread_id - 1) ÷ (max_c)) % max_b + 1
c = ((thread_id - 1) ÷ (1)) % max_c + 1
result = a + b + c
result_vector[thread_id] = result
"""
function gen_bc_kernel_indexers(
    indexers::Vector{String},
    bi::Int
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
    indexer_loops *= "bi = $(bi)\n"
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

function gen_bc_kernel_body(
    callback_fun::CallbackFunction,
    indexers::Vector{String},
    bc_expr_args::OrderedDict{String, String},
    bi::Int
)
    kernel_body::String = gen_bc_kernel_indexers(indexers, bi)

    kernel_body::String *=
"""
fid = mesh_bdryface[fi]
eid = mesh_face2element[1, fid]
fbid = mesh_bids[bi]
volume = geometric_factors_volume[eid]
area = geometric_factors_area[fid]
area_over_volume = area / volume
"""
    index_offset_parts::Vector{String} = []
    index_var_maxes::Vector{String} = ["MAX_$(index_var)" for index_var in indexers[2:end]]
    for (i, indexer) in enumerate(indexers[2:end])
        if i == 1
            multiplier = "1"
        else
            multiplier = join(index_var_maxes[begin : i - 1], " * ")
        end
        push!(index_offset_parts, "($(indexer) - 1) * ($(multiplier))")
    end
    index_offset = join(index_offset_parts, " + ")

    kernel_body *=
"""
index_offset = $(index_offset)
row_index = index_offset + 1 + dofs_per_node * (eid - 1)

"""

    kernel_body *= "x = fv_info_faceCenters[1, fid]\n"
    if finch_state.config.dimension >= 2
        kernel_body *= "y = fv_info_faceCenters[2, fid]\n"
    end
    if finch_state.config.dimension == 3
        kernel_body *= "z = fv_info_faceCenters[3, fid]\n"
    end

    for (param, expr) in bc_expr_args
        kernel_body *= "$(param) = $(expr)\n\n"
    end

    body_expr::Expr = Meta.parse(callback_fun.body)
    Base.remove_linenums!(body_expr)
    update_variable_indexing!(body_expr)

    kernel_body *= "# Callback function body\n"
    for part in body_expr.args
        kernel_body *= "$(string(part))\n"
    end
    kernel_body *= "\n"

    kernel_body *=
"""
# Update output
boundary_flux[thread_id] = result * area_over_volume
boundary_dof_index[thread_id] = row_index
global_vector[row_index] += boundary_flux[thread_id]

return nothing
"""

    return kernel_body
end

function update_variable_indexing!(expr::Expr)
    finch_pde_variables = [string(v.symbol) for v in finch_state.variables]

    # Some types of vector accesses need to have an additional index
    if expr.head == :ref
        var = expr.args[1]
        
        if string(var) in finch_pde_variables
            push!(expr.args, :eid)
        elseif var == :normal
            push!(expr.args, :fid)
        end
    end

    # Remove "return result" statement
    filter!(part -> part != :(return result), expr.args)

    for i in eachindex(expr.args)
        part = expr.args[i]
        if typeof(part) == Expr
            update_variable_indexing!(part)
        end
    end
end