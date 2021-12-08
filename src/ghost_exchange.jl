#=
Functions related to exchanging ghost data.
=#

#=
Exchange ghosts.
This does both send and receive.
Values are sent/received directly from the variable.values arrays

NOTE: This is currently only set up for the case where each process has a distinct partition.
For other strategies, this needs to be modified.
=#
function exchange_ghosts(var, grid, tag=0)
    tag = config.num_procs + Int(round(tag));
    if config.solver_type == FV
        return exchange_ghosts_fv(var, grid, tag);
    else
        printerr("Ghost exchange not ready for "*config.solver_type, true)
    end
end

function exchange_ghosts_fv(var, grid, tag)
    n_neighbors = grid.num_neighbor_partitions;
    # How big does each buffer need to be?
    buffer_size = zeros(Int, n_neighbors);
    send_arrays = fill(zeros(0), n_neighbors);
    recv_arrays = fill(zeros(0), n_neighbors);
    requests = Vector{MPI.Request}(undef, n_neighbors * 2);
    if typeof(var) <: Array
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        for vi=1:length(var)
            dofs_per_node += var[vi].total_components;
        end
    else
        # one variable
        dofs_per_node = var.total_components;
    end
    for ni=1:n_neighbors
        buffer_size[ni] = dofs_per_node * grid.ghost_counts[ni]; # number of values to send to neighbor ni
        send_arrays[ni] = zeros(buffer_size[ni]);
        recv_arrays[ni] = zeros(buffer_size[ni]);
        
        # Fill the send buffer
        next_ind = 1;
        if typeof(var) <: Array
            for vi=1:length(var)
                for ei=1:grid.ghost_counts[ni]
                    for comp_i=1:var[vi].total_components
                        send_arrays[ni][next_ind] = var[vi].values[comp_i,grid.ghost_index[ni][2,ei]];
                        next_ind += 1;
                    end
                end
            end
        else
            for ei=1:grid.ghost_counts[ni]
                for comp_i=1:var.total_components
                    send_arrays[ni][next_ind] = var.values[comp_i,grid.ghost_index[ni][2,ei]];
                    next_ind += 1;
                end
            end
        end
    end
    
    # send ghosts
    if config.use_mpi
        for ni=1:grid.num_neighbor_partitions
            requests[ni*2-1] = MPI.Isend(send_arrays[ni], grid.neighboring_partitions[ni], tag, MPI.COMM_WORLD);
        end
    end
    # receive ghosts
    if config.use_mpi
        for ni=1:grid.num_neighbor_partitions
            requests[ni*2] = MPI.Irecv!(recv_arrays[ni], grid.neighboring_partitions[ni], tag, MPI.COMM_WORLD);
        end
    end
    # wait for it
    stats = MPI.Waitall!(requests);
    
    # Put the received ghost values in place
    for ni=1:n_neighbors
        next_ind = 1;
        if typeof(var) <: Array
            for vi=1:length(var)
                for ei=1:grid.ghost_counts[ni]
                    for comp_i=1:var[vi].total_components
                        var[vi].values[comp_i,grid.ghost_index[ni][1,ei]] = recv_arrays[ni][next_ind];
                        next_ind += 1;
                    end
                end
            end
        else
            for ei=1:grid.ghost_counts[ni]
                for comp_i=1:var.total_components
                    var.values[comp_i,grid.ghost_index[ni][1,ei]] = recv_arrays[ni][next_ind];
                    next_ind += 1;
                end
            end
        end
    end
end
