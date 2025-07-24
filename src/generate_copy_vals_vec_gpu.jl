copy_bdry_vals_to_vec_kernel::String =
"""
function copy_bdry_val_to_vector_kenel(var_values, vec, MAX_comp, MAX_face, offset, dofs_per_node, mesh_bdryface_bid, mesh_face2element)
    thread_id = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if thread_id > MAX_comp * MAX_face
        return
    end
    comp = (thread_id - 1) ÷ MAX_face + 1
    face = (thread_id - 1) % MAX_face + 1
    fid = mesh_bdryface_bid[face]
    eid = mesh_face2element[1, fid]
    vec[(eid - 1) * dofs_per_node + offset + comp] = var_values[comp, eid]

    return nothing
end
"""

place_vector_in_var_kernel::String =
"""
function place_vector_in_var_kernel(var_values, vec, MAX_comp, MAX_dof, offset, totalcomponents, fv_dofs_partition)
    thread_id = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if thread_id > MAX_comp * MAX_dof
        return
    end
    comp = (thread_id - 1) ÷ MAX_dof + 1
    dof = (thread_id - 1) % MAX_dof + 1
    index = (dof - 1) * totalcomponents + comp + offset
    if index > fv_dofs_partition
        return
    end
    var_values[comp, dof] = vec[index]

    return nothing
end
"""

function gen_copy_bdry_vals_to_vec_calls()
    kernel_calls::String = ""
    offset = 0
    for var_id in 1:size(finch_state.prob.bc_func, 1)
        MAX_comp = length(finch_state.variables[var_id].symvar)
        for bid in 1:size(finch_state.prob.bc_func, 2)
            MAX_face = length(finch_state.grid_data.bdryface[bid])
            if finch_state.prob.bc_type[var_id, bid] == DIRICHLET
                kernel_calls *= "CUDA.@sync @cuda threads = 256 blocks = ceil(Int, $(MAX_comp) * $(MAX_face) / 256) copy_bdry_vals_to_vector_kenel(" *
                    "variables_$(var_id)_values_gpu, solution_gpu, $(MAX_comp), $(MAX_face), $(offset), " *
                    "$(finch_state.variables[var_id].total_components), mesh_bdryface_$(bid)_gpu, mesh_face2element_gpu)\n"
            end
        end
        offset += MAX_comp
    end

    return kernel_calls
end

function gen_place_vector_in_var_calls()
    kernel_calls::String = ""
    totalcomponents = sum(finch_state.variables[var_id].total_components for var_id in 1:length(finch_state.variables))
    offset = 0
    for var_id in 1:length(finch_state.variables)
        MAX_comp = finch_state.variables[var_id].total_components
        MAX_dof = size(finch_state.variables[var_id].values, 2)
        kernel_calls *= "CUDA.@sync @cuda threads = 256 blocks = ceil(Int, $(MAX_comp) * $(MAX_dof) / 256) place_vector_in_var_kernel(" *
            "variables_$(var_id)_values_gpu, solution_gpu, $(MAX_comp), $(MAX_dof), $(offset), $(totalcomponents), fv_dofs_partition)\n"
        offset += MAX_comp
    end

    return kernel_calls
end
