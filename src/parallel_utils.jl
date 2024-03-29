#=
Utilities for distributed parallel.
=#

# # Buffers for the global system
# mutable struct ParallelBuffers
#     full_AI::Vector
#     full_AJ::Vector
#     full_AV::Vector
#     full_b::Vector
    
#     vec_b::Vector
# end

# FE tools
######################################################################################
# For partitioned meshes
# Returns (b_order, b_sizes) only for proc 0
# b_order is the global indices of each proc's vector
# b_sizes is the size of each proc's vector
function get_partitioned_ordering(dofs_per_node::Int, grid_data::Grid, config::FinchConfig)
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
function share_time_step_info(stepper::Stepper, config::FinchConfig)
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

# For multiple processes, gather the system directly using MPI
# Note: rescatter_b only applies when rhs_only
function gather_system_assembled(A::Union{Nothing, SparseMatrixCSC}, b::Vector, nnodes::Int, dofs_per_node::Int, 
                        b_order::Vector, b_sizes::Vector, config::FinchConfig, buffers::ParallelBuffers; rescatter_b::Bool=false)
    rhs_only = (A===nothing);
    if config.num_procs > 1
        if !isbitstype(config.float_type) || !isbitstype(config.index_type)
            printerr("MPI can only work with types for which isbitstype == true.\nCurrent types are: int="*
                        string(config.index_type)*", float="*string(config.float_type), fatal=true);
        end
        grid_data = finch_state.grid_data;
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
                # If buffers aren't allocated yet, do so
                if length(buffers.full_AI) == 0
                    buffers.full_AI = zeros(config.index_type, total_length);
                    buffers.full_AJ = zeros(config.index_type, total_length);
                    buffers.full_AV = zeros(config.float_type, total_length);
                end
                
                AI_buf = MPI.VBuffer(buffers.full_AI, chunk_sizes, displacements, MPI.Datatype(config.index_type));
                AJ_buf = MPI.VBuffer(buffers.full_AJ, chunk_sizes, displacements, MPI.Datatype(config.index_type));
                AV_buf = MPI.VBuffer(buffers.full_AV, chunk_sizes, displacements, MPI.Datatype(config.float_type));
                
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
                full_A = sparse(buffers.full_AI, buffers.full_AJ, buffers.full_AV);
                
                # Next gather b
                chunk_sizes = b_sizes;
                displacements = zeros(Int, config.num_procs); # for the irregular gatherv
                total_length = chunk_sizes[1];
                for proc_i=2:config.num_procs
                    displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
                    total_length += chunk_sizes[proc_i];
                end
                if length(buffers.full_b) == 0
                    buffers.full_b = zeros(config.float_type, total_length);
                    buffers.vec_b = zeros(config.float_type, grid_data.nnodes_global * dofs_per_node);
                end
                b_buf = MPI.VBuffer(buffers.full_b, chunk_sizes, displacements, MPI.Datatype(config.float_type));
                MPI.Gatherv!(b, b_buf, 0, MPI.COMM_WORLD);
                
                # Overlapping values will be added.
                for i=1:total_length
                    buffers.vec_b[b_order[i]] += buffers.full_b[i];
                end
                
                return (full_A, buffers.vec_b);
                
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
                if length(buffers.full_b) == 0
                    buffers.full_b = zeros(config.float_type, total_length);
                    buffers.vec_b = zeros(config.float_type, grid_data.nnodes_global * dofs_per_node);
                end
                
                b_buf = MPI.VBuffer(buffers.full_b, chunk_sizes, displacements, MPI.Datatype(config.float_type));
                MPI.Gatherv!(b, b_buf, 0, MPI.COMM_WORLD);
                
                # Overlapping values will be added.
                for i=1:total_length
                    buffers.vec_b[b_order[i]] += buffers.full_b[i];
                end
                
                if rescatter_b
                    distribute_solution!(buffers.vec_b, b, nnodes, dofs_per_node, b_order, b_sizes, config, buffers);
                end
                return buffers.vec_b;
                
            else # other procs just send their data
                MPI.Gatherv!(b, nothing, 0, MPI.COMM_WORLD);
                
                if rescatter_b
                    distribute_solution!([], b, nnodes, dofs_per_node, b_order, b_sizes, config, buffers);
                end
                return [1];
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

# For multiple processes, gather the system directly using MPI
# Note: rescatter_b only applies when rhs_only
function gather_system(AI::Union{Nothing, Vector{Int}}, AJ::Union{Nothing, Vector{Int}}, AV::Union{Nothing, Vector}, 
                        b::Vector, nnodes::Int, dofs_per_node::Int, b_order::Vector{Int}, b_sizes::Vector{Int}, 
                        config::FinchConfig, buffers::ParallelBuffers; rescatter_b::Bool=false, return_assembled::Bool=true)
    rhs_only = (AI===nothing);
    if config.num_procs > 1
        if !isbitstype(config.float_type)
            printerr("MPI can only work with types for which isbitstype == true.\nCurrent type is float="*string(config.float_type), fatal=true);
        end
        grid_data = finch_state.grid_data;
        if !rhs_only
            # For now just gather all of A in proc 0 to assemble.
            # The row and column indices have to be changed according to partition2global.
            # Also, b must be reordered on proc 0, so send the needed indices as well.
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
                # If buffers aren't allocated yet, do so
                if length(buffers.full_AI) == 0
                    buffers.full_AI = zeros(config.index_type, total_length);
                    buffers.full_AJ = zeros(config.index_type, total_length);
                    buffers.full_AV = zeros(config.float_type, total_length);
                end
                
                AI_buf = MPI.VBuffer(buffers.full_AI, chunk_sizes, displacements, MPI.Datatype(config.index_type));
                AJ_buf = MPI.VBuffer(buffers.full_AJ, chunk_sizes, displacements, MPI.Datatype(config.index_type));
                AV_buf = MPI.VBuffer(buffers.full_AV, chunk_sizes, displacements, MPI.Datatype(config.float_type));
                
                MPI.Gatherv!(AI, AI_buf, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AJ, AJ_buf, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AV, AV_buf, 0, MPI.COMM_WORLD);
                
                # # Assemble A
                # full_A = sparse(buffers.full_AI, buffers.full_AJ, buffers.full_AV);
                
                # Next gather b
                chunk_sizes = b_sizes;
                displacements = zeros(Int, config.num_procs); # for the irregular gatherv
                total_length = chunk_sizes[1];
                for proc_i=2:config.num_procs
                    displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
                    total_length += chunk_sizes[proc_i];
                end
                if length(buffers.full_b) == 0
                    buffers.full_b = zeros(config.float_type, total_length);
                    buffers.vec_b = zeros(config.float_type, grid_data.nnodes_global * dofs_per_node);
                end
                b_buf = MPI.VBuffer(buffers.full_b, chunk_sizes, displacements, MPI.Datatype(config.float_type));
                MPI.Gatherv!(b, b_buf, 0, MPI.COMM_WORLD);
                
                # Overlapping values will be added.
                buffers.vec_b .= 0.0;
                for i=1:total_length
                    buffers.vec_b[b_order[i]] += buffers.full_b[i];
                end
                
                if return_assembled
                    return (sparse(buffers.full_AI, buffers.full_AJ, buffers.full_AV), buffers.vec_b);
                else
                    return (buffers.full_AI, buffers.full_AJ, buffers.full_AV, buffers.vec_b);
                end
                
            else # other procs just send their data
                send_buf = [length(AI)];
                MPI.Gather!(send_buf, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AI, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AJ, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AV, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(b, nothing, 0, MPI.COMM_WORLD);
                
                if return_assembled
                    return (sparse(AI, AJ, AV), b);
                else
                    return (AI, AJ, AV, b);
                end
                
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
                if length(buffers.full_b) == 0
                    buffers.full_b = zeros(config.float_type, total_length);
                    buffers.vec_b = zeros(config.float_type, grid_data.nnodes_global * dofs_per_node);
                end
                
                b_buf = MPI.VBuffer(buffers.full_b, chunk_sizes, displacements, MPI.Datatype(config.float_type));
                MPI.Gatherv!(b, b_buf, 0, MPI.COMM_WORLD);
                
                # Overlapping values will be added.
                buffers.vec_b .= 0.0;
                for i=1:total_length
                    buffers.vec_b[b_order[i]] += buffers.full_b[i];
                end
                
                if rescatter_b
                    distribute_solution!(buffers.vec_b, b, nnodes, dofs_per_node, b_order, b_sizes, config, buffers);
                end
                return buffers.vec_b;
                
            else # other procs just send their data
                MPI.Gatherv!(b, nothing, 0, MPI.COMM_WORLD);
                
                if rescatter_b
                    distribute_solution!([], b, nnodes, dofs_per_node, b_order, b_sizes, config, buffers);
                end
                return b;
            end
        end
        
    else # one process
        if rhs_only
            return b;
        else
            if return_assembled
                return (sparse(AI, AJ, AV), b);
            else
                return (AI, AJ, AV, b);
            end
        end
    end
end

# For multiple processes, gather the system directly using MPI
# Note: rescatter_b only applies when rhs_only
function gather_system_FV(AI::Union{Nothing, Vector{Int}}, AJ::Union{Nothing, Vector{Int}}, AV::Union{Nothing, Vector}, 
                        b::Vector, nel::Int, dofs_per_node::Int, 
                        config::FinchConfig, buffers::ParallelBuffers; rescatter_b::Bool=false)
    rhs_only = (AI===nothing);
    if config.num_procs > 1
        if !isbitstype(config.float_type)
            printerr("MPI can only work with types for which isbitstype == true.\nCurrent type is float="*string(config.float_type), fatal=true);
        end
        grid_data::Grid = finch_state.grid_data;
        if !rhs_only
            # For now just gather all of A in proc 0 to assemble.
            # The row and column indices have to be changed according to partition2global.
            # Also, b must be reordered on proc 0, so send the needed indices as well.
            nzlength = length(AI);
            for i=1:length(AI)
                if AI[i] == 0
                    nzlength = i-1;
                    break;
                end
                dof = mod(AI[i]-1,dofs_per_node)+1;
                el = Int((AI[i]-dof) / dofs_per_node + 1);
                AI[i] = (grid_data.partition2global_element[el] - 1) * dofs_per_node + dof;
                
                dof = mod(AJ[i]-1,dofs_per_node)+1;
                el = Int((AJ[i]-dof) / dofs_per_node + 1);
                AJ[i] = (grid_data.partition2global_element[el] - 1) * dofs_per_node + dof;
            end
            if nzlength < length(AI)
                # trim them to nonzero parts
                AI = AI[1:nzlength];
                AJ = AJ[1:nzlength];
                AV = AV[1:nzlength];
            end
            
            # Send things to proc 0 to solve the system.
            # I know this is not how to do it, but let's just get it working for now.
            if config.proc_rank == 0
                # First figure out how long each proc's arrays are
                send_buf = [nzlength];
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
                # If buffers aren't allocated yet, do so
                if length(buffers.full_AI) == 0
                    buffers.full_AI = zeros(config.index_type, total_length);
                    buffers.full_AJ = zeros(config.index_type, total_length);
                    buffers.full_AV = zeros(config.float_type, total_length);
                end
                
                AI_buf = MPI.VBuffer(buffers.full_AI, chunk_sizes, displacements, MPI.Datatype(config.index_type));
                AJ_buf = MPI.VBuffer(buffers.full_AJ, chunk_sizes, displacements, MPI.Datatype(config.index_type));
                AV_buf = MPI.VBuffer(buffers.full_AV, chunk_sizes, displacements, MPI.Datatype(config.float_type));
                
                MPI.Gatherv!(AI, AI_buf, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AJ, AJ_buf, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AV, AV_buf, 0, MPI.COMM_WORLD);
                
                # Next gather b
                owned_dofs = dofs_per_node * nel;
                send_buf = [owned_dofs];
                recv_buf = zeros(Int, config.num_procs);
                MPI.Gather!(send_buf, recv_buf, 0, MPI.COMM_WORLD);
                
                chunk_sizes = recv_buf;
                displacements = zeros(Int, config.num_procs); # for the irregular gatherv
                total_length = chunk_sizes[1];
                for proc_i=2:config.num_procs
                    displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
                    total_length += chunk_sizes[proc_i];
                end
                if length(buffers.full_b) == 0
                    buffers.full_b = zeros(config.float_type, total_length);
                    buffers.vec_b = zeros(config.float_type, grid_data.nel_global * dofs_per_node);
                    buffers.b_order = zeros(config.float_type, total_length);
                end
                b_buf = MPI.VBuffer(buffers.full_b, chunk_sizes, displacements, MPI.Datatype(config.float_type));
                MPI.Gatherv!(b[1:owned_dofs], b_buf, 0, MPI.COMM_WORLD);
                
                # get the ordering for b
                bord_buf = MPI.VBuffer(buffers.b_order, chunk_sizes, displacements, MPI.Datatype(Int));
                MPI.Gatherv!(grid_data.partition2global_element[1:owned_dofs], bord_buf, 0, MPI.COMM_WORLD);
                
                # Overlapping values will be added.
                for i=1:total_length
                    buffers.vec_b[buffers.b_order[i]] = buffers.full_b[i]; 
                end
                
                return (sparse(buffers.full_AI, buffers.full_AJ, buffers.full_AV), buffers.vec_b)
                # return (buffers.full_AI, buffers.full_AJ, buffers.full_AV, buffers.vec_b);
                
            else # other procs just send their data
                send_buf = [nzlength];
                owned_dofs = dofs_per_node * nel;
                MPI.Gather!(send_buf, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AI, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AJ, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(AV, nothing, 0, MPI.COMM_WORLD);
                MPI.Gather!([owned_dofs], nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(b[1:owned_dofs], nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(grid_data.partition2global_element[1:owned_dofs], nothing, 0, MPI.COMM_WORLD);
                
                return (sparse(AI,AJ,AV),b);
                # return (AI,AJ,AV,b);
            end
            
        else # RHS only\
            if config.proc_rank == 0
                owned_dofs = dofs_per_node * nel;
                send_buf = [owned_dofs];
                recv_buf = zeros(Int, config.num_procs);
                MPI.Gather!(send_buf, recv_buf, 0, MPI.COMM_WORLD);
                
                # gather b
                chunk_sizes = recv_buf;
                displacements = zeros(Int, config.num_procs); # for the irregular gatherv
                total_length = chunk_sizes[1];
                for proc_i=2:config.num_procs
                    displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
                    total_length += chunk_sizes[proc_i];
                end
                if length(buffers.full_b) == 0
                    buffers.full_b = zeros(config.float_type, total_length);
                    buffers.vec_b = zeros(config.float_type, grid_data.nel_global * dofs_per_node);
                    buffers.b_order = zeros(config.float_type, total_length);
                end
                
                b_buf = MPI.VBuffer(buffers.full_b, chunk_sizes, displacements, MPI.Datatype(config.float_type));
                MPI.Gatherv!(b[1:owned_dofs], b_buf, 0, MPI.COMM_WORLD);
                
                # get the ordering for b
                bord_buf = MPI.VBuffer(buffers.b_order, chunk_sizes, displacements, MPI.Datatype(Int));
                MPI.Gatherv!(grid_data.partition2global_element[1:owned_dofs], bord_buf, 0, MPI.COMM_WORLD);
                
                # Overlapping values will be added.
                for i=1:total_length
                    buffers.vec_b[buffers.b_order[i]] = buffers.full_b[i]; 
                end
                
                if rescatter_b
                    distribute_solution_FV!(buffers.vec_b, b, nel, dofs_per_node, config, buffers);
                end
                return buffers.vec_b;
                
            else # other procs just send their data
                MPI.Gather!([length(b)], nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(b, nothing, 0, MPI.COMM_WORLD);
                MPI.Gatherv!(grid_data.partition2global_element, nothing, 0, MPI.COMM_WORLD);
                
                if rescatter_b
                    distribute_solution_FV!([], b, nel, dofs_per_node, config, buffers);
                end
                return b;
            end
        end
        
    else # one process
        if rhs_only
            return b;
        else
            return (AI, AJ, AV, b);
        end
    end
end

function distribute_solution!(sol::Vector, local_sol::Vector, nnodes::Int, dofs_per_node::Int, 
                            b_order::Vector, b_sizes::Vector, config::FinchConfig, buffers::ParallelBuffers)
    if config.num_procs > 1
        if !isbitstype(config.float_type) || !isbitstype(config.index_type)
            printerr("MPI can only work with types for which isbitstype == true.\nCurrent types are: int="*
                        string(config.index_type)*", float="*string(config.float_type), fatal=true);
        end
        
        if config.proc_rank == 0 # 0 has the full sol
            # Need to reorder b and put in a larger array according to b_order
            total_length = sum(b_sizes);
            for i=1:total_length
                buffers.full_b[i] = sol[b_order[i]];
            end
            
            # scatter b
            chunk_sizes = b_sizes;
            displacements = zeros(Int, config.num_procs); # for the irregular scatterv
            for proc_i=2:config.num_procs
                displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
            end
            sol_buf = MPI.VBuffer(buffers.full_b, chunk_sizes, displacements, MPI.Datatype(config.float_type));
            
            # Scatter it amongst the little ones
            MPI.Scatterv!(sol_buf, local_sol, 0, MPI.COMM_WORLD);
            
        else # Others don't have sol
            MPI.Scatterv!(nothing, local_sol, 0, MPI.COMM_WORLD);
        end
    end
end

function distribute_solution_FV!(sol::Vector, local_sol::Vector, nel::Int, dofs_per_node::Int, 
                                config::FinchConfig, buffers::ParallelBuffers)
    if config.num_procs > 1
        if !isbitstype(config.float_type) || !isbitstype(config.index_type)
            printerr("MPI can only work with types for which isbitstype == true.\nCurrent types are: int="*
                        string(config.index_type)*", float="*string(config.float_type), fatal=true);
        end
        
        if config.proc_rank == 0 # 0 has the full sol
            send_buf = [nel*dofs_per_node];
            recv_buf = zeros(Int, config.num_procs);
            MPI.Gather!(send_buf, recv_buf, 0, MPI.COMM_WORLD);
            
            # Need to reorder b and put in a larger array according to b_order
            total_length = sum(recv_buf);
            for i=1:total_length
                buffers.full_b[i] = sol[buffers.b_order[i]];
            end
            
            # scatter b
            chunk_sizes = recv_buf;
            displacements = zeros(Int, config.num_procs); # for the irregular scatterv
            for proc_i=2:config.num_procs
                displacements[proc_i] = displacements[proc_i-1] + chunk_sizes[proc_i-1];
            end
            sol_buf = MPI.VBuffer(buffers.full_b, chunk_sizes, displacements, MPI.Datatype(config.float_type));
            
            # Scatter it amongst the little ones
            MPI.Scatterv!(sol_buf, local_sol, 0, MPI.COMM_WORLD);
            
        else # Others don't have sol
            send_buf = [nel*dofs_per_node];
            recv_buf = zeros(Int, config.num_procs);
            MPI.Gather!(send_buf, recv_buf, 0, MPI.COMM_WORLD);
            
            MPI.Scatterv!(nothing, local_sol, 0, MPI.COMM_WORLD);
        end
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
function exchange_ghosts_fv(var::Vector, grid::Grid, dofs_per_node::Int, tag::Int, config::FinchConfig)
    if !isbitstype(config.float_type) || !isbitstype(config.index_type)
        printerr("MPI can only work with types for which isbitstype == true.\nCurrent types are: int="*
                    string(config.index_type)*", float="*string(config.float_type), fatal=true);
    end
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