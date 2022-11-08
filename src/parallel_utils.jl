#=
Utilities for distributed parallel.

=#

# FE tools
######################################################################################
# For partitioned meshes
# Returns (b_order, b_sizes) only for proc 0
# b_order is the global indices of each proc's vector
# b_sizes is the size of each proc's vector
function get_partitioned_ordering(dofs_per_node::Int, grid_data::Grid, config::Finch_config)
    if config.num_procs > 1
        nnodes = size(grid_data.allnodes,2);
        # The global ordering of this partition's b vector entries
        b_order_loc = zeros(Int, nnodes*dofs_per_node);
        for ni=1:nnodes
            for di=1:dofs_per_node
                b_order_loc[(ni-1)*dofs_per_node + di] = (grid_data.partition2global[ni]-1)*dofs_per_node + di;
            end
        end
        
        # only proc 0 needs this info for now
        if config.proc_rank == 0
            # The b_sizes
            send_buf = [nnodes*dofs_per_node];
            recv_buf = zeros(Int, config.num_procs);
            MPI.Gather!(send_buf, recv_buf, 0, MPI.COMM_WORLD);
            
            # b_order
            chunk_sizes = recv_buf;
            displacements = zeros(Int, config.num_procs); # for the irregular gatherv
            total_length = chunk_sizes[1];
            for proc_i=2:config.num_procs
                displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
                total_length += chunk_sizes[proc_i];
            end
            full_b_order = zeros(Int, total_length);
            b_order_buf = MPI.VBuffer(full_b_order, chunk_sizes, displacements, MPI.Datatype(Int));
            MPI.Gatherv!(b_order_loc, b_order_buf, 0, MPI.COMM_WORLD);
            
            return (full_b_order, chunk_sizes);
            
        else
            send_buf = [nnodes*dofs_per_node];
            MPI.Gather!(send_buf, nothing, 0, MPI.COMM_WORLD);
            MPI.Gatherv!(b_order_loc, nothing, 0, MPI.COMM_WORLD);
            
            return ([1], [1]);
        end
        
    else
        return (zeros(Int,0), zeros(Int,0));
    end
end

# Share time stepper step size and number so that all procs are the same.
function share_time_step_info(stepper::Stepper, config::Finch_config)
    if config.num_procs > 1
        # gather a list of Nsteps for each proc
        send_buf = [stepper.Nsteps];
        recv_buf = zeros(Int, config.num_procs);
        MPI.Allgather!(send_buf, recv_buf, MPI.COMM_WORLD);
        
        # All procs use the maximal Nsteps and recalculate dt
        max_nsteps = maximum(recv_buf);
        end_time = stepper.dt * stepper.Nsteps;
        stepper.Nsteps = max_nsteps;
        stepper.dt = end_time / max_nsteps;
    end
end

# For multiple processes, gather the system, distribute the solution
# Note: rescatter_b only applies when rhs_only
function gather_system(A, b, nnodes, dofs_per_node, b_order, b_sizes, config; rescatter_b=false)
    rhs_only = (A===nothing);
    if config.num_procs > 1
        if !rhs_only
            # For now just gather all of A in proc 0 to assemble.
            # The row and column indices have to be changed according to partition2global.
            # Also, b must be reordered on proc 0, so send the needed indices as well.
            (AI, AJ, AV) = findnz(A);
            for i=1:length(AI)
                dof = mod(AI[i]-1,dofs_per_node)+1;
                node = Int(floor((AI[i]-dof) / dofs_per_node) + 1);
                AI[i] = (grid_data.partition2global[node] - 1) * dofs_per_node + dof;
                
                dof = mod(AJ[i]-1,dofs_per_node)+1;
                node = Int(floor((AJ[i]-dof) / dofs_per_node) + 1);
                AJ[i] = (grid_data.partition2global[node] - 1) * dofs_per_node + dof;
            end
            
            if config.proc_rank == 0
                # First figure out how long each proc's arrays are
                send_buf = [length(AI)];
                recv_buf = zeros(Int, config.num_procs);
                MPI.Gather!(send_buf, recv_buf, 0, MPI.COMM_WORLD);
                
                # Use gatherv to accumulate A
                chunk_sizes = recv_buf;
                displacements = zeros(Int, config.num_procs); # for the irregular gatherv
                total_length = chunk_sizes[1];
                for proc_i=2:config.num_procs
                    displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
                    total_length += chunk_sizes[proc_i];
                end
                full_AI = zeros(Int, total_length);
                full_AJ = zeros(Int, total_length);
                full_AV = zeros(Float64, total_length);
                AI_buf = MPI.VBuffer(full_AI, chunk_sizes, displacements, MPI.Datatype(Int));
                AJ_buf = MPI.VBuffer(full_AJ, chunk_sizes, displacements, MPI.Datatype(Int));
                AV_buf = MPI.VBuffer(full_AV, chunk_sizes, displacements, MPI.Datatype(Float64));
                
                MPI.Gatherv!(AI, AI_buf, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AJ, AJ_buf, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AV, AV_buf, 0, MPI.COMM_WORLD);
                
                # # Modify AI and AJ to global indices using b_order
                # b_start = 0;
                # for proc_i=1:config.num_procs
                #     for ai=(displacements[proc_i]+1):(displacements[proc_i] + chunk_sizes[proc_i])
                #         full_AI[ai] = b_order[b_start + full_AI[ai]];
                #         full_AJ[ai] = b_order[b_start + full_AJ[ai]];
                #     end
                #     b_start += b_sizes[proc_i];
                # end
                
                # Assemble A
                full_A = sparse(full_AI, full_AJ, full_AV);
                
                # Next gather b
                chunk_sizes = b_sizes;
                displacements = zeros(Int, config.num_procs); # for the irregular gatherv
                total_length = chunk_sizes[1];
                for proc_i=2:config.num_procs
                    displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
                    total_length += chunk_sizes[proc_i];
                end
                full_b = zeros(total_length);
                b_buf = MPI.VBuffer(full_b, chunk_sizes, displacements, MPI.Datatype(Float64));
                MPI.Gatherv!(b, b_buf, 0, MPI.COMM_WORLD);
                
                # Overlapping values will be added.
                new_b = zeros(grid_data.nnodes_global * dofs_per_node);
                for i=1:total_length
                    new_b[b_order[i]] += full_b[i];
                end
                
                return (full_A, new_b);
                
            else # other procs just send their data
                send_buf = [length(AI)];
                MPI.Gather!(send_buf, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AI, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AJ, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AV, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(b, nothing, 0, MPI.COMM_WORLD);
                
                return (ones(1,1),[1]);
            end
            
        else # RHS only\
            if config.proc_rank == 0
                # gather b
                chunk_sizes = b_sizes;
                displacements = zeros(Int, config.num_procs); # for the irregular gatherv
                total_length = chunk_sizes[1];
                for proc_i=2:config.num_procs
                    displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
                    total_length += chunk_sizes[proc_i];
                end
                full_b = zeros(total_length);
                b_buf = MPI.VBuffer(full_b, chunk_sizes, displacements, MPI.Datatype(Float64));
                MPI.Gatherv!(b, b_buf, 0, MPI.COMM_WORLD);
                
                # Overlapping values will be added.
                new_b = zeros(grid_data.nnodes_global * dofs_per_node);
                for i=1:total_length
                    new_b[b_order[i]] += full_b[i];
                end
                
                if rescatter_b
                    return distribute_solution(new_b, nnodes, dofs_per_node, b_order, b_sizes, config)
                else
                    return new_b;
                end
                
            else # other procs just send their data
                MPI.Gatherv!(b, nothing, 0, MPI.COMM_WORLD);
                
                if rescatter_b
                    return distribute_solution(nothing, nnodes, dofs_per_node, b_order, b_sizes, config)
                else
                    return [1];
                end
            end
        end
        
    else # one process
        if rhs_only
            return b;
        else
            return (A, b);
        end
    end
end

function distribute_solution(sol::Vector, nnodes::Int, dofs_per_node::Int, b_order::Vector, b_sizes::Vector, config::Finch_config)
    if config.num_procs > 1
        my_sol = zeros(nnodes*dofs_per_node);
        if config.proc_rank == 0 # 0 has the full sol
            # Need to reorder b and put in a larger array according to b_order
            total_length = sum(b_sizes);
            full_sol = zeros(total_length);
            for i=1:total_length
                full_sol[i] = sol[b_order[i]];
            end
            
            # scatter b
            chunk_sizes = b_sizes;
            displacements = zeros(Int, config.num_procs); # for the irregular scatterv
            for proc_i=2:config.num_procs
                displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
            end
            sol_buf = MPI.VBuffer(full_sol, chunk_sizes, displacements, MPI.Datatype(Float64));
            
            # Scatter it amongst the little ones
            MPI.Scatterv!(sol_buf, my_sol, 0, MPI.COMM_WORLD);
            
        else # Others don't have sol
            MPI.Scatterv!(nothing, my_sol, 0, MPI.COMM_WORLD);
        end
        return my_sol;
        
    else
        return sol;
    end
end

# Simply does a reduction.
# Works for scalar values or vectors, but vectors are not themselves reduced.
# 1, 2, 3 -> 6
# [1,10,100], [2,20,200], [3,30,300] -> [6,60,600]
function combine_values(val; combine_op = +)
    if config.num_procs > 1
        rval = MPI.Allreduce(val, combine_op, MPI.COMM_WORLD);
        return rval;
        
    else
        return val;
    end
end

# This reduces a global vector to one value
# [1,10,100], [2,20,200], [3,30,300] -> 666
function reduce_vector(vec::Array)
    if config.num_procs > 1
        rval = MPI.Allreduce(sum(vec), +, MPI.COMM_WORLD);
        return rval;
        
    else
        return sum(vec);
    end
end

# Exchange variable values so that each partition has values for neighboring cells.
function exchange_ghosts_fv(var::Vector, grid::Grid, dofs_per_node::Int, tag::Int)
    n_neighbors = grid.num_neighbor_partitions;
    if n_neighbors < 1 # do nothing if no neighbors
        return;
    end
    
    # How big does each buffer need to be?
    buffer_size = zeros(Int, n_neighbors);
    send_arrays = fill(zeros(0), n_neighbors);
    recv_arrays = fill(zeros(0), n_neighbors);
    requests = Vector{MPI.Request}(undef, n_neighbors * 2);
    
    for ni=1:n_neighbors
        buffer_size[ni] = dofs_per_node * grid.ghost_counts[ni]; # number of values to send to neighbor ni
        send_arrays[ni] = zeros(buffer_size[ni]);
        recv_arrays[ni] = zeros(buffer_size[ni]);
        
        # Fill the send buffer
        next_ind = 1;
        for vi=1:length(var)
            for ei=1:grid.ghost_counts[ni]
                for comp_i=1:var[vi].total_components
                    send_arrays[ni][next_ind] = var[vi].values[comp_i,grid.ghost_index[ni][2,ei]];
                    next_ind += 1;
                end
            end
        end
    end
    
    # send ghosts
    for ni=1:n_neighbors
        requests[ni*2-1] = MPI.Isend(send_arrays[ni], grid.neighboring_partitions[ni], tag, MPI.COMM_WORLD);
    end
    
    # receive ghosts
    for ni=1:n_neighbors
        requests[ni*2] = MPI.Irecv!(recv_arrays[ni], grid.neighboring_partitions[ni], tag, MPI.COMM_WORLD);
    end
    
    # wait for it
    stats = MPI.Waitall!(requests);
    
    
    # Put the received ghost values in place
    for ni=1:n_neighbors
        next_ind = 1;
        for vi=1:length(var)
            for ei=1:grid.ghost_counts[ni]
                for comp_i=1:var[vi].total_components
                    var[vi].values[comp_i,grid.ghost_index[ni][1,ei]] = recv_arrays[ni][next_ind];
                    next_ind += 1;
                end
            end
        end
    end
end