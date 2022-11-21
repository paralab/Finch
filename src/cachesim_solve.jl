# cachsim output
export linear_solve_cachesim
# import ..Finch: init_cachesimout, add_cachesim_array, cachesim_load, cachesim_store, cachesim_load_range, cachesim_store_range

function linear_solve_cachesim(var, bilinear, linear, stepper=nothing)
    # # if config.linalg_matrixfree
    # #     return solve_matrix_free_sym(var, bilinear, linear, stepper);
    # #     #return solve_matrix_free_asym(var, bilinear, linear, stepper);
    # # end
    # # If more than one variable
    # if typeof(var) <: Array
    #     # multiple variables being solved for simultaneously
    #     dofs_per_node = 0;
    #     dofs_per_loop = 0;
    #     for vi=1:length(var)
    #         dofs_per_loop += length(var[vi].symvar);
    #         dofs_per_node += var[vi].total_components;
    #     end
    # else
    #     # one variable
    #     dofs_per_loop = length(var.symvar);
    #     dofs_per_node = var.total_components;
    # end
    # N1 = size(finch_state.grid_data.allnodes,2);
    # Nn = dofs_per_node * N1;
    # Np = finch_state.refel.Np;
    # nel = finch_state.mesh_data.nel;
    
    # init_cachesimout(N1, finch_state.refel, nel, dofs_per_node, finch_state.variables);
    
    # if finch_state.prob.time_dependent && !(finch_state.stepper === nothing)
    #     #TODO time dependent coefficients
    #     assemble_t = @elapsed((A, b) = assemble_cachesim(var, bilinear, linear, 0, stepper.dt));
    #     log_entry("Assembly took "*string(assemble_t)*" seconds");

    #     log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
    #     t = 0;
    #     sol = [];
    #     start_t = Base.Libc.time();
    #     for i=1:stepper.Nsteps
    #         b = assemble_rhs_only_cachesim(var, linear, t, stepper.dt);
    #         #sol = A\b;
    #         cachesim_load_range(1);
    #         cachesim_load_range(2);
            
    #         t += stepper.dt;
    #     end
    #     end_t = Base.Libc.time();

    #     log_entry("Cachesim solve took "*string(end_t-start_t)*" seconds");
    #     return 0;

    # else
    #     assemble_t = @elapsed((A, b) = assemble_cachesim(var, bilinear, linear));
    #     #sol_t = @elapsed(sol = A\b);
    #     cachesim_load_range(1);
    #     cachesim_load_range(2);

    #     log_entry("Assembly took "*string(assemble_t)*" seconds");
        
    #     return 0;
    # end
end

function assemble_cachesim(var, bilinear, linear, t=0.0, dt=0.0)
    # Np = refel.Np;
    # nel = mesh_data.nel;
    # N1 = size(grid_data.allnodes,2);
    # multivar = typeof(var) <: Array;
    # if multivar
    #     # multiple variables being solved for simultaneously
    #     dofs_per_node = 0;
    #     dofs_per_loop = 0;
    #     for vi=1:length(var)
    #         dofs_per_loop += length(var[vi].symvar);
    #         dofs_per_node += var[vi].total_components;
    #     end
    # else
    #     # one variable
    #     dofs_per_loop = length(var.symvar);
    #     dofs_per_node = var.total_components;
    # end
    # Nn = dofs_per_node * N1;
    
    # b = 0;
    # A = 0;
    # AI = 0;
    # AJ = 0;
    # AV = 0;
    
    # allnodes_id = add_cachesim_array(size(grid_data.allnodes),8);
    # ak_id = add_cachesim_array(Np*Np*dofs_per_node*dofs_per_node,8);
    # aktmp_id = add_cachesim_array(Np*Np*dofs_per_node*dofs_per_node,8);
    # bk_id = add_cachesim_array(Np*dofs_per_node,8);
    # bktmp_id = add_cachesim_array(Np*dofs_per_node,8);
    
    # #  Elemental loop follows elemental ordering
    # for e=elemental_order;
    #     gis = grid_data.glbvertex[:,e];
    #     vx = grid_data.allnodes[:,gis];         # coordinates of element's vertices
    #     glb = grid_data.loc2glb[:,e];                 # global indices of this element's nodes for extracting values from var arrays
    #     xe = grid_data.allnodes[:,glb[:]];  # coordinates of this element's nodes for evaluating coefficient functions
    #     cachesim_load_range(allnodes_id, gis[:], 1:refel.dim);
    #     cachesim_load_range(allnodes_id, glb[:], 1:refel.dim);
        
    #     cachesim_load_range(aktmp_id);
    #     cachesim_load_range(bktmp_id);
        
    #     Astart = (e-1)*Np*dofs_per_node*Np*dofs_per_node + 1; # The segment of AI, AJ, AV for this element
        
    #     # The linear part. Compute the elemental linear part for each dof
    #     rhsargs = (var, e, 0, grid_data, geo_factors, refel, t, dt, 0, 0);
    #     lhsargs = (var, e, 0, grid_data, geo_factors, refel, t, dt, 0, 0);
    #     if dofs_per_node == 1
    #         linchunk = linear.func(rhsargs);  # get the elemental linear part
    #         #b[glb] .+= linchunk;
    #         cachesim_load_range(2, glb);
    #         cachesim_load_range(bk_id);
    #         cachesim_store_range(2, glb);

    #         bilinchunk = bilinear.func(lhsargs); # the elemental bilinear part
    #         #A[glb, glb] .+= bilinchunk;         # This will be very inefficient for sparse A
    #         Arange = Astart:(Astart + Np*Np*dofs_per_node*dofs_per_node);
    #         cachesim_load_range(1, Arange, 1:3);
    #         cachesim_load_range(ak_id);
    #         cachesim_store_range(1, Arange, 1:3);
            
    #     elseif typeof(var) == Variable
    #         # only one variable, but more than one dof
    #         linchunk = linear.func(rhsargs);
    #         insert_linear_cachesim!(bk_id, linchunk, glb, 1:dofs_per_node, dofs_per_node);

    #         bilinchunk = bilinear.func(lhsargs);
    #         insert_bilinear_cachesim!(ak_id, Astart, bilinchunk, glb, 1:dofs_per_node, dofs_per_node);
    #     else
    #         linchunk = linear.func(rhsargs);
    #         insert_linear_cachesim!(bk_id, linchunk, glb, 1:dofs_per_node, dofs_per_node);

    #         bilinchunk = bilinear.func(lhsargs);
    #         insert_bilinear_cachesim!(ak_id, Astart, bilinchunk, glb, 1:dofs_per_node, dofs_per_node);
    #     end
    # end
    
    # cachesim_load_range(1);
    
    # # Boundary conditions are not considered

    # return (A, b);
end

# assembles the A and b in Au=b
function assemble_rhs_only_cachesim(var, linear, t=0.0, dt=0.0)
    # Np = refel.Np;
    # nel = mesh_data.nel;
    # N1 = size(grid_data.allnodes,2);
    # multivar = typeof(var) <: Array;
    # if multivar
    #     # multiple variables being solved for simultaneously
    #     dofs_per_node = 0;
    #     var_to_dofs = [];
    #     for vi=1:length(var)
    #         tmp = dofs_per_node;
    #         dofs_per_node += length(var[vi].symvar);
    #         push!(var_to_dofs, (tmp+1):dofs_per_node);
    #     end
    # else
    #     # one variable
    #     dofs_per_node = length(var.symvar);
    # end
    # Nn = dofs_per_node * N1;

    # #b = zeros(Nn);

    # for e=1:nel;
    #     glb = grid_data.loc2glb[:,e];                 # global indices of this element's nodes for extracting values from var arrays
    #     xe = grid_data.allnodes[:,glb[:]];  # coordinates of this element's nodes for evaluating coefficient functions

    #     rhsargs = (var, e, 0, grid_data, geo_factors, refel, t, dt, 0, 0);

    #     #linchunk = linear.func(args);  # get the elemental linear part
    #     if dofs_per_node == 1
    #         linchunk = linear.func(rhsargs);  # get the elemental linear part
    #         #b[glb] .+= linchunk;

    #     elseif typeof(var) == Variable
    #         # only one variable, but more than one dof
    #         linchunk = linear.func(rhsargs);
    #         insert_linear_cachesim!(bk_id, linchunk, glb, 1:dofs_per_node, dofs_per_node);

    #     else
    #         linchunk = linear.func(rhsargs);
    #         insert_linear_cachesim!(bk_id, linchunk, glb, 1:dofs_per_node, dofs_per_node);

    #     end
    # end
    
    # return b;
end

# Inset the single dof into the greater construct
function insert_linear_cachesim!(bk_id, bel, glb, dof, Ndofs)
    # # group nodal dofs
    # for d=1:length(dof)
    #     ind = glb.*Ndofs .- (Ndofs-dof[d]);
    #     ind2 = ((d-1)*length(glb)+1):(d*length(glb));

    #     #b[ind] = b[ind] + bel[ind2];
        
    #     cachesim_load_range(2,ind);
    #     cachesim_load_range(bk_id,ind2);
    #     cachesim_store_range(2,ind);
    # end
end

function insert_bilinear_cachesim!(ak_id, Astart, ael, glb, dof, Ndofs)
    # Np = length(glb);
    # # group nodal dofs
    # for dj=1:length(dof)
    #     indj = glb.*Ndofs .- (Ndofs-dof[dj]);
    #     indj2 = ((dj-1)*Np+1):(dj*Np);
    #     for di=1:length(dof)
    #         indi = glb.*Ndofs .- (Ndofs-dof[di]);
    #         indi2 = ((di-1)*Np+1):(di*Np);
            
    #         for jj=1:Np
    #             offset = Astart + (jj-1 + Np*(dj-1))*Np*Ndofs + Np*(di-1) - 1;
    #             Arange = (offset+1):(offset + Np);
    #             cachesim_load_range(1, Arange, 1:3);
    #             cachesim_load_range(ak_id, indi2, indj2);
    #             cachesim_store_range(1, Arange, 1:3);
    #         end
    #     end
    # end
    
end
