#=
Functions related to exchanging ghost data.
=#

"""
    ExchangeGroup

A description of information that needs to be communicated between processes during
a ghost exchange. It includes the ID of the involved processes, the included elements,
and the variable and component indices. The ghost exchange function will use this
to set up the data and path to communicate.
"""
mutable struct ExchangeGroup
    members::Array{Int,1} # ranks of the members of the group
    
    ghost_elements_only::Bool   # True if the exchange is for ghost element pairs
    exchange_all_data::Bool     # True if all variables/components are exchanged
    
    elements::Array{Int,2}  # 1xN list of global indices of elements to share. 2xN Pairs if ghost_elements_only.
    variables::Array{Int,1} # indices of variables to share (empty if exchange_all_data=true)
    components::Array{Array{Int,1},1} # components of each variable to share (empty if exchange_all_data=true)
end

"""

"""
function build_exchange_groups(grid; members=[], variables=[], components=[], ghost_elements_only=true, exchange_all_data=true)
    # This list of exchange groups will be returned
    ex_groups = [];
    
    if ghost_elements_only
        if exchange_all_data
            # If ghost_elements_only and exchange_all_data, this matches the ghost info from the grid
            ngroups = grid.num_neighbor_partitions;
            for gi=1:ngroups
                neighbor = grid.neighboring_partitions[gi];
                
                push!(ex_groups, ExchangeGroup([neighbor], ghost_elements_only, exchange_all_data, grid.ghost_index[gi], [], []));
            end
        else
            # If exchange_all_data = false, this matches the ghost info from the grid, but specifies some components
            ngroups = grid.num_neighbor_partitions;
            for gi=1:ngroups
                neighbor = grid.neighboring_partitions[gi];
                
                push!(ex_groups, ExchangeGroup([neighbor], ghost_elements_only, exchange_all_data, grid.ghost_index[gi], variables, components));
            end
        end
    else # all local elements are shared with ranks that have same partition
        if exchange_all_data
            
        else
            
        end
    end
    
    return ex_groups;
end

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
        printerr("Ghost exchange not ready for "*config.solver_type, fatal=true)
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
        
        # wait for it
        stats = MPI.Waitall!(requests);
    end
    
    
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
