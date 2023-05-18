#=
Generated functions for mixed2d
=#

# begin solve function for u

function generated_solve_function_for_u(var::Vector{Variable{FT}}, mesh::Grid, refel::Vector{Refel{FT}}, geometric_factors::GeometricFactors, config::FinchConfig, coefficients::Vector{Coefficient}, variables::Vector{Variable{FT}}, test_functions::Vector{Coefficient}, ordered_indexers::Vector{Indexer}, prob::FinchProblem, time_stepper::Stepper, buffers::ParallelBuffers, timer_output::TimerOutput, nl_var=nothing) where FT<:AbstractFloat
    
    # User specified data types for int and float
    # int type is Int64
    # float type is Float64
    
    # pre/post step functions if defined
    pre_step_function = prob.pre_step_function;
    post_step_function = prob.post_step_function;
    
    function genfunction_1(x::Union{Float64,Float64},y::Union{Float64,Float64},z::Union{Float64,Float64},t::Union{Float64,Float64},node_index::Int,face_index::Int,indices::Vector{Int}) return (-8*pi*pi*sin(2*pi*x)*sin(2*pi*y)); end
    
    
    # Prepare some useful numbers
    dofs_per_node = 1;
    dofs_per_loop = 1;
    dof_offsets = [0];
    nnodes_partition = size(mesh.allnodes,2);
    nnodes_global = nnodes_partition;
    num_elements = mesh.nel_owned;
    num_elements_global = mesh.nel_global;
    num_elements_ghost = mesh.nel_ghost;
    num_faces = mesh.nface_owned + mesh.nface_ghost;
    
    dofs_global = dofs_per_node * nnodes_global;
    fv_dofs_global = dofs_per_node * num_elements_global;
    dofs_partition = dofs_per_node * nnodes_partition;
    fv_dofs_partition = dofs_per_node * (num_elements + num_elements_ghost);
    num_partitions = config.num_partitions;
    proc_rank = config.proc_rank;
    
    nodes_per_element_multi = [refel[1].Np, refel[2].Np]; ###mixed
    qnodes_per_element_multi = [refel[1].Nqp, refel[2].Nqp]; ###mixed
    faces_per_element_multi = [refel[1].Nfaces, refel[2].Nfaces]; ###mixed
    nodes_per_face_multi = [refel[1].Nfp[1], refel[2].Nfp[1]]; ###mixed
    dofs_per_element_multi = [dofs_per_node * nodes_per_element_multi[1], dofs_per_node * nodes_per_element_multi[2]]; ###mixed
    local_system_size_multi = [dofs_per_loop * nodes_per_element_multi[1], dofs_per_loop * nodes_per_element_multi[2]]; ###mixed
    
    nodes_per_element = maximum(nodes_per_element_multi); ###mixed
    qnodes_per_element = maximum(qnodes_per_element_multi); ###mixed
    dofs_per_element = maximum(dofs_per_element_multi); ###mixed
    local_system_size = maximum(local_system_size_multi); ###mixed
    
    # FEM specific pieces
    Q_multi = [refel[1].Q, refel[2].Q]; ###mixed
    wg_multi = [refel[1].wg, refel[2].wg]; ###mixed
    surf_wg_multi = [refel[1].surf_wg[1], refel[2].surf_wg[2]]; ###mixed
    gness = size(mesh.face2glb,2); # CG=1, DG=2
    
    global_IJV_offset = 0; ###mixed
    
    # For partitioned meshes
    (partitioned_order, partitioned_sizes) = get_partitioned_ordering(dofs_per_node, mesh, config);
    #= Allocate global matrix(IJV form) and vector. =#
    allocated_nonzeros = (num_elements * dofs_per_element * dofs_per_element)
    next_nonzero_index = (allocated_nonzeros + 1)
    global_matrix_I::Vector{Int64} = zeros(Int64, allocated_nonzeros)
    global_matrix_J::Vector{Int64} = zeros(Int64, allocated_nonzeros)
    global_matrix_V::Vector{Float64} = zeros(Float64, allocated_nonzeros)
    global_vector::Vector{Float64} = zeros(Float64, dofs_global)
    global_solution::Vector{Float64} = zeros(Float64, dofs_global)
    if (num_partitions > 1)
        solution::Vector{Float64} = zeros(Float64, dofs_partition)
    else
        solution = global_solution
    end

    #= I and J vectors should init as ones =#
    global_matrix_I .= 1
    global_matrix_J .= 1
    #= Allocate elemental matrix and vector. =#
    element_matrix::Matrix{Float64} = zeros(Float64, local_system_size, local_system_size)
    element_vector::Vector{Float64} = zeros(Float64, local_system_size)
    #= Boundary done flag for each node. =#
    bdry_done::Vector{Int64} = zeros(Int64, nnodes_global)
    #= No indexed variables =#
    index_values::Vector{Int64} = zeros(Int64, 0)
    #= Allocate coefficient vectors. =#
    detj::Vector{Float64} = zeros(Float64, qnodes_per_element)
    #= Allocate for derivative matrices. =#
    RQ1::Matrix{Float64} = zeros(Float64, qnodes_per_element, nodes_per_element)
    RQ2::Matrix{Float64} = zeros(Float64, qnodes_per_element, nodes_per_element)
    NODALvalue__f_1::Vector{Float64} = zeros(Float64, nodes_per_element)
    value__f_1::Vector{Float64} = zeros(Float64, qnodes_per_element)
    t = 0.0
    @timeit timer_output "assembly" begin
        for ei = 1:num_elements
            eid = mesh.elemental_order[ei]
            index_offset = 0
            
            ###mixed
            refel_ind = mesh.refel_ind[eid];
            nodes_per_element = nodes_per_element_multi[refel_ind]; ###mixed
            qnodes_per_element = qnodes_per_element_multi[refel_ind]; ###mixed
            dofs_per_element = dofs_per_element_multi[refel_ind]; ###mixed
            local_system_size = local_system_size_multi[refel_ind]; ###mixed
            
            # FEM specific pieces
            Q = Q_multi[refel_ind]; ###mixed
            wg = wg_multi[refel_ind]; ###mixed
            surf_wg = surf_wg_multi[refel_ind]; ###mixed
            
            build_derivative_matrix(refel[refel_ind], geometric_factors, 1, eid, 0, RQ1) ###mixed
            build_derivative_matrix(refel[refel_ind], geometric_factors, 2, eid, 0, RQ2) ###mixed
            
            #= Prepare derivative matrices. =#
            #= Evaluate coefficients. =#
            for ni = 1:nodes_per_element
                nodeID = mesh.loc2glb[ni, eid]
                x = mesh.allnodes[1, nodeID]
                y = mesh.allnodes[2, nodeID]
                z = 0.0
                NODALvalue__f_1[ni]::Float64 = Float64(genfunction_1(x,y,z,t,nodeID, 0, index_values))
            end
            
            for col = 1:qnodes_per_element
                value__f_1[col] = 0.0
                for row = 1:nodes_per_element
                    value__f_1[col] = (value__f_1[col] + (Q[col, row] * NODALvalue__f_1[row]))
                end

            end
            
            for qnode_i = 1:qnodes_per_element
                detj[qnode_i] = geometric_factors.detJ[qnode_i, eid]
            end
            
            
            @inbounds begin
                for col=1:nodes_per_element
                    for row=1:nodes_per_element
                        element_matrix[row, col] = 0;
                        
                        for i=1:qnodes_per_element
                            element_matrix[row, col] += ((RQ1[i, row] * (wg[i] * detj[i]) * RQ1[i, col]) + (RQ2[i, row] * (wg[i] * detj[i]) * RQ2[i, col]));
                            
                        end

                    end

                end

            end

            
            @inbounds begin
                for row=1:nodes_per_element
                    element_vector[row] = 0;
                    
                    for col=1:qnodes_per_element
                        element_vector[row] += (Q[col, row] * (wg[col] * detj[col] * value__f_1[col] * -1));
                        
                    end

                end

            end
            
            #= No face loop needed. =#
            #= Apply boundary conditions. =#
            apply_boundary_conditions_elemental(var, eid, mesh, refel[refel_ind], geometric_factors, prob, t, element_matrix, element_vector, bdry_done, index_offset, index_values) ###mixed
            #= Place elemental parts in global system. =#
            
            next_ind = global_IJV_offset + 1; ###mixed
            for ni = 1:nodes_per_element
                glb_i = mesh.loc2glb[ni, eid]
                global_vector[glb_i] = (global_vector[glb_i] + element_vector[ni])
                for nj = 1:nodes_per_element
                    glb_j = mesh.loc2glb[nj, eid]
                    global_matrix_I[next_ind] = glb_i
                    global_matrix_J[next_ind] = glb_j
                    global_matrix_V[next_ind] = element_matrix[ni, nj]
                    next_ind = (next_ind + 1)
                end

            end

            global_IJV_offset += nodes_per_element * nodes_per_element; ###mixed
        end

        
    end # timer:assembly

    
    if (num_partitions > 1)
        (global_matrix, global_vector) = gather_system(global_matrix_I, global_matrix_J, global_matrix_V, global_vector, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config, buffers);
        if length(global_vector) > length(global_solution)
            global_solution = zeros(length(global_vector));
        end

    
    else
        global_matrix = sparse(@view(global_matrix_I[1:(next_nonzero_index-1)]), @view(global_matrix_J[1:(next_nonzero_index-1)]), @view(global_matrix_V[1:(next_nonzero_index-1)]));
        
    end

    @timeit timer_output "lin_solve" begin
        linear_system_solve!(global_matrix, global_vector, global_solution, config);
    end # timer:lin_solve

    
    if (num_partitions > 1)
        distribute_solution!(global_solution, solution, nnodes_partition, dofs_per_node, partitioned_order, partitioned_sizes, config, buffers);
    end

    @timeit timer_output "scatter" begin
        place_vector_in_vars(var, solution);
    end # timer:scatter

    
    
    return nothing;
    
end # function



# end solve function for u

