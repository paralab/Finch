#=
# FV solver
=#
module FVSolver

export solve, nonlinear_solve

# See finch_import_symbols.jl for a list of all imported symbols.
import ..Finch: @import_finch_symbols
@import_finch_symbols()

using LinearAlgebra, SparseArrays

include("fv_boundary.jl");

# Things to do before and after each step(or stage) ########
function default_pre_step() end
function default_post_step() end

pre_step_function = default_pre_step;
post_step_function = default_post_step;

function set_pre_step(fun) global pre_step_function = fun; end
function set_post_step(fun) global post_step_function = fun; end
############################################################

function init_solver()
    global pre_step_function = default_pre_step;
    global post_step_function = default_post_step;
end

function linear_solve(var, source_lhs, source_rhs, flux_lhs, flux_rhs, stepper=nothing, assemble_func=nothing)
    # Treat explicit and implicit solvers differently
    if stepper === nothing
        printerr("FV assumes a time dependent problem. Set time stepper and initial conditions.", fatal=true)
    elseif stepper.implicit
        return linear_solve_implicit(var, source_lhs, source_rhs, flux_lhs, flux_rhs, stepper, assemble_func)
    else # explicit
        return linear_solve_explicit(var, source_lhs, source_rhs, flux_lhs, flux_rhs, stepper, assemble_func)
    end
end

function linear_solve_explicit(var, source_lhs, source_rhs, flux_lhs, flux_rhs, stepper=nothing, assemble_func=nothing)
    # If more than one variable
    if typeof(var) <: Array
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        dofs_per_loop = 0;
        for vi=1:length(var)
            dofs_per_loop += length(var[vi].symvar);
            dofs_per_node += var[vi].total_components;
        end
    else
        # one variable
        dofs_per_loop = length(var.symvar);
        dofs_per_node = var.total_components;
    end
    if fv_grid === nothing
        grid = grid_data;
        nel = grid_data.nel_owned;
        nfaces = size(grid_data.face2element, 2);
    else
        grid = fv_grid;
        nel = fv_grid.nel_owned;
        nfaces = size(fv_grid.face2element, 2);
    end
    
    Nn = dofs_per_node * nel;
    Nf = dofs_per_node * nfaces
    
    # Allocate arrays that will be used by assemble
    # These vectors will hold the integrated values(one per cell).
    # They will later be combined.
    sourcevec = zeros(Nn);
    fluxvec = zeros(Nn);
    facefluxvec = zeros(Nf);
    face_done = zeros(nfaces); # Increment when the corresponding flux value is computed.
    allocated_vecs = [sourcevec, fluxvec, facefluxvec, face_done];
    
    if prob.time_dependent && !(stepper === nothing)
        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        if config.num_partitions > 1 exchange_ghosts(var, grid, 0); end
        sol = get_var_vals(var);
        
        # allocate storage used by steppers
        if stepper.type == LSRK4
            tmppi = zeros(size(sol));
            tmpki = zeros(size(sol));
        elseif stepper.type == RK4
            tmpvals = zeros(size(sol));
            tmpki = zeros(length(sol), stepper.stages);
        end
        
        start_t = Base.Libc.time();
        last2update = 0;
        last10update = 0;
        if config.proc_rank == 0 print("Time stepping progress(%): 0"); end
        for i=1:stepper.Nsteps
            if stepper.stages > 1
                # LSRK4 is a special case, low storage
                if stepper.type == LSRK4
                    # Low storage RK4: 
                    # p0 = u
                    #   ki = ai*k(i-1) + dt*f(p(i-1), t+ci*dt)
                    #   pi = p(i-1) + bi*ki
                    # u = p5
                    
                    tmppi = get_var_vals(var, tmppi);
                    tmpki = zeros(size(sol));
                    for rki=1:stepper.stages
                        rktime = t + stepper.c[rki]*stepper.dt;
                        # p(i-1) is currently in u
                        
                        if config.num_partitions > 1 exchange_ghosts(var, grid, i); end
                        pre_step_function();
                        
                        sol = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, rktime, stepper.dt, assemble_loops=assemble_func);
                        
                        
                        if rki == 1 # because a1 == 0
                            tmpki = stepper.dt .* sol;
                        else
                            tmpki = stepper.a[rki].*tmpki + stepper.dt.*sol;
                        end
                        tmppi = tmppi + stepper.b[rki].*tmpki
                        
                        FV_copy_bdry_vals_to_vector(var, tmppi, grid, dofs_per_node);
                        place_sol_in_vars(var, tmppi, stepper);
                        
                        post_step_function();
                    end
                    
                else
                    # Explicit multi-stage methods: 
                    # x = x + dt*sum(bi*ki)
                    # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
                    
                    # will hold the final result
                    # sol = get_var_vals(var, sol);
                    # will be placed in var.values for each stage
                    for stage=1:stepper.stages
                        stime = t + stepper.c[stage]*stepper.dt;
                        
                        # Update the values in vars to be used in this stage
                        if stage > 1
                            initialized_tmpvals = false;
                            for j=1:stage
                                if stepper.a[stage, j] > 0
                                    if !initialized_tmpvals
                                        initialized_tmpvals = true;
                                        for k=1:length(sol)
                                            tmpvals[k] = sol[k] + stepper.dt * stepper.a[stage, j] * tmpki[k,j];
                                        end
                                    else
                                        for k=1:length(sol)
                                            tmpvals[k] += stepper.dt * stepper.a[stage, j] * tmpki[k,j];
                                        end
                                    end
                                end
                            end
                            
                            if initialized_tmpvals
                                FV_copy_bdry_vals_to_vector(var, tmpvals, grid, dofs_per_node);
                                place_sol_in_vars(var, tmpvals, stepper);
                            end
                            post_step_function(); # seems weird, but imagine this is happening after stage-1
                        end
                        
                        if config.num_partitions > 1 exchange_ghosts(var, grid, i); end
                        pre_step_function();
                        
                        tmpki[:,stage] = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, stime, stepper.dt, assemble_loops=assemble_func);
                        
                    end
                    # Stages are done. Assemble the final result for this step
                    for stage=1:stepper.stages
                        sol += stepper.dt * stepper.b[stage] .* tmpki[:, stage];
                    end
                    FV_copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node);
                    place_sol_in_vars(var, sol, stepper);
                    
                    post_step_function();
                end
                
            elseif stepper.type == EULER_EXPLICIT
                if config.num_partitions > 1 exchange_ghosts(var, grid, i); end
                pre_step_function();
                
                sol = sol .+ stepper.dt .* assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t, stepper.dt, assemble_loops=assemble_func);
                
                FV_copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node);
                place_sol_in_vars(var, sol, stepper);
                
                post_step_function();
                
            elseif stepper.type == PECE
                # Predictor (explicit Euler)
                if config.num_partitions > 1 exchange_ghosts(var, grid, i); end
                pre_step_function();
                tmpsol = sol .+ stepper.dt .* assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t, stepper.dt, assemble_loops=assemble_func);
                
                FV_copy_bdry_vals_to_vector(var, tmpsol, grid, dofs_per_node);
                place_sol_in_vars(var, tmpsol, stepper);
                post_step_function();
                
                # Corrector (implicit Euler)
                if config.num_partitions > 1 exchange_ghosts(var, grid, i); end
                pre_step_function();
                sol = sol .+ stepper.dt .* assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t+stepper.dt, stepper.dt, assemble_loops=assemble_func);
                
                FV_copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node);
                place_sol_in_vars(var, sol, stepper);
                post_step_function();
                
            else
                printerr("Only explicit time steppers for FV are ready. TODO")
                return sol;
            end
            
            # ########## uncomment to return after one time step
            # return sol
            # ############
            
            t += stepper.dt;
            
            if config.proc_rank == 0
                progressPercent = Int(floor(i*100.0/stepper.Nsteps));
                if progressPercent - last2update >= 2
                    last2update = progressPercent;
                    if progressPercent - last10update >= 10
                        print(progressPercent);
                        last10update = progressPercent;
                    else
                        print(".");
                    end
                end
            end
        end
        if config.proc_rank == 0 println(""); end
        if config.num_partitions > 1 exchange_ghosts(var, grid, 0); end
        end_t = Base.Libc.time();
        
        log_entry("Stepping took "*string(end_t-start_t)*" seconds.");
        return sol;

    else
        # Does it make any sense to do this for time-independent problems?
        return sol;
    end
end

function linear_solve_implicit(var, source_lhs, source_rhs, flux_lhs, flux_rhs, stepper=nothing, assemble_func=nothing)
    # If more than one variable
    if typeof(var) <: Array
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        dofs_per_loop = 0;
        for vi=1:length(var)
            dofs_per_loop += length(var[vi].symvar);
            dofs_per_node += var[vi].total_components;
        end
    else
        # one variable
        dofs_per_loop = length(var.symvar);
        dofs_per_node = var.total_components;
    end
    if fv_grid === nothing
        grid = grid_data;
        nel = grid_data.nel_owned;
        nfaces = size(grid_data.face2element, 2);
    else
        grid = fv_grid;
        nel = fv_grid.nel_owned;
        nfaces = size(fv_grid.face2element, 2);
    end
    
    Nn = dofs_per_node * nel;
    Nf = dofs_per_node * nfaces
    
    # Allocate arrays that will be used by assemble
    face_done = zeros(nfaces); # true when the corresponding flux value is computed.
    facefluxvec = zeros(Nf);
    rhsvec = zeros(Nn);
    lhsmatI = zeros(Int, Nn+Nf);
    lhsmatJ = zeros(Int, Nn+Nf);
    lhsmatV = zeros(Nn+Nf);
    allocated_vecs = [rhsvec, lhsmatI, lhsmatJ, lhsmatV, facefluxvec, face_done];
    
    if prob.time_dependent && !(stepper === nothing)
        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        if config.num_partitions > 1 exchange_ghosts(var, grid, 0); end
        sol = get_var_vals(var);
        
        # allocate storage used by steppers
        if stepper.type == LSRK4
            tmppi = zeros(size(sol));
            tmpki = zeros(size(sol));
        elseif stepper.type == RK4
            tmpki = zeros(length(sol), stepper.stages);
        end
        
        start_t = Base.Libc.time();
        last2update = 0;
        last10update = 0;
        if config.proc_rank == 0 print("Time stepping progress(%): 0"); end
        for i=1:stepper.Nsteps
            if stepper.type == EULER_IMPLICIT || stepper.type == CRANK_NICHOLSON
                if config.num_partitions > 1 exchange_ghosts(var, grid, i); end
                pre_step_function();
                
                (A, b) = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t, stepper.dt, assemble_loops=assemble_func);
                sol = sol .+ stepper.dt .* (A\b);
                
                FV_copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node);
                place_sol_in_vars(var, sol, stepper);
                
                post_step_function();
                
            else
                printerr("Problem with stepper: "*string(stepper), fatal=true)
            end
            
            # ########## uncomment to return after one time step
            # return sol
            # ############
            
            t += stepper.dt;
            
            if config.proc_rank == 0
                progressPercent = Int(floor(i*100.0/stepper.Nsteps));
                if progressPercent - last2update >= 2
                    last2update = progressPercent;
                    if progressPercent - last10update >= 10
                        print(progressPercent);
                        last10update = progressPercent;
                    else
                        print(".");
                    end
                end
            end
        end
        if config.proc_rank == 0 println(""); end
        if config.num_partitions > 1 exchange_ghosts(var, grid, 0); end
        end_t = Base.Libc.time();
        
        log_entry("Stepping took "*string(end_t-start_t)*" seconds.");
        return sol;

    else
        # Does it make any sense to do this for time-independent problems?
        return sol;
    end
end

#
function assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0; assemble_loops=nothing)
    # If an assembly loop function was provided, use it
    if !(assemble_loops === nothing)
        return assemble_loops.func(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t, dt);
    end
    # If parent maps were created for high order flux, use a slightly different function
    if !(parent_maps === nothing)
        return assemble_using_parent_child(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t, dt);
    end
    
    nel = fv_grid.nel_owned;
    
    # Label things that were allocated externally
    sourcevec = allocated_vecs[1];
    fluxvec = allocated_vecs[2];
    facefluxvec = allocated_vecs[3];
    face_done = allocated_vecs[4];
    
    face_done .= 0;
    
    # Elemental loop
    for ei=1:nel
        eid = elemental_order[ei]; # The index of this element
        # Zero the result vectors for this element
        sourcevec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] .= 0;
        fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] .= 0;
        
        ##### Source integrated over the cell #####
        # Compute RHS volume integral
        if !(source_rhs === nothing)
            #sourceargs = prepare_args(var, eid, 0, RHS, "volume", t, dt); # (var, e, nodex, loc2glb, refel, detj, J, t, dt)
            sourceargs = (var, eid, 0, fv_grid, fv_geo_factors, fv_info, refel, t, dt);
            # source = source_rhs.func(sourceargs) ./ fv_geo_factors.volume[eid];
            source = source_rhs.func(sourceargs);
            # Add to global source vector
            sourcevec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] = source;
        end
        
        ##### Flux integrated over the faces #####
        # Loop over this element's faces.
        for i=1:refel.Nfaces
            fid = fv_grid.element2face[i, eid];
            # Only one element on either side is available here. For more use parent/child version.
            (leftel, rightel) = fv_grid.face2element[:,fid];
            if rightel == 0
                neighborhood = [[leftel],[]];
            else
                neighborhood = [[leftel],[rightel]];
            end
            
            if !(flux_rhs === nothing)
                if face_done[fid] == 0
                    face_done[fid] = 1; # Possible race condition, but in the worst case it will be computed twice.
                    
                    fluxargs = (var, eid, fid, neighborhood, fv_grid, fv_geo_factors, fv_info, refel, t, dt);
                    flux = flux_rhs.func(fluxargs) .* fv_geo_factors.area[fid];
                    # Add to global flux vector for faces
                    facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] = flux;
                    # Combine all flux for this element
                    fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] += flux ./ fv_geo_factors.volume[eid];
                    
                else
                    # This flux has either been computed or is being computed by another thread.
                    # The state will need to be known before paralellizing, but for now assume it's complete.
                    fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] -= facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] ./ fv_geo_factors.volume[eid];
                end
            end
            
            # Boundary conditions are applied to flux
            fbid = fv_grid.facebid[fid]; # BID of this face
            if fbid > 0
                facex = fv_grid.allnodes[:, fv_grid.face2glb[:,1,fid]];  # face node coordinates
                
                if typeof(var) <: Array
                    dofind = 0;
                    for vi=1:length(var)
                        for compo=1:length(var[vi].symvar)
                            dofind = dofind + 1;
                            if prob.bc_type[var[vi].index, fbid] == NO_BC
                                # do nothing
                            elseif prob.bc_type[var[vi].index, fbid] == FLUX
                                # compute the value and add it to the flux directly
                                # Qvec = (refel.surf_wg[fv_grid.faceRefelInd[1,fid]] .* fv_geo_factors.face_detJ[fid])' * (refel.surf_Q[fv_grid.faceRefelInd[1,fid]])[:, refel.face2local[fv_grid.faceRefelInd[1,fid]]]
                                # Qvec = Qvec ./ fv_geo_factors.area[fid];
                                # bflux = FV_flux_bc_rhs_only(prob.bc_func[var[vi].index, fbid][compo], facex, Qvec, t, dofind, dofs_per_node) .* fv_geo_factors.area[fid];
                                bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var[vi].index, fbid][compo], fid, t) .* fv_geo_factors.area[fid];
                                
                                fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ fv_geo_factors.volume[eid];
                                facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                            elseif prob.bc_type[var[vi].index, fbid] == DIRICHLET
                                # Set variable array and handle after the face loop
                                var[vi].values[compo,eid] = FV_evaluate_bc(prob.bc_func[var[vi].index, fbid][compo], eid, fid, t);
                            else
                                printerr("Unsupported boundary condition type: "*prob.bc_type[var[vi].index, fbid]);
                            end
                        end
                    end
                else
                    for d=1:dofs_per_node
                        dofind = d;
                        if prob.bc_type[var.index, fbid] == NO_BC
                            # do nothing
                        elseif prob.bc_type[var.index, fbid] == FLUX
                            # compute the value and add it to the flux directly
                            # Qvec = (refel.surf_wg[fv_grid.faceRefelInd[1,fid]] .* fv_geo_factors.face_detJ[fid])' * (refel.surf_Q[fv_grid.faceRefelInd[1,fid]])[:, refel.face2local[fv_grid.faceRefelInd[1,fid]]]
                            # Qvec = Qvec ./ fv_geo_factors.area[fid];
                            # bflux = FV_flux_bc_rhs_only(prob.bc_func[var.index, fbid][d], facex, Qvec, t, dofind, dofs_per_node) .* fv_geo_factors.area[fid];
                            bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var.index, fbid][d], fid, t) .* fv_geo_factors.area[fid];
                            
                            fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ fv_geo_factors.volume[eid];
                            facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                        elseif prob.bc_type[var.index, fbid] == DIRICHLET
                            # Set variable array and handle after the face loop
                            var.values[compo,eid] = FV_evaluate_bc(prob.bc_func[var.index, fbid][d], eid, fid, t);
                        else
                            printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, fbid]);
                        end
                    end
                end
            end# BCs
            
        end# face loop
        
    end# element loop
    
    return sourcevec + fluxvec;
end

function assemble_implicit(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0; assemble_loops=nothing, rhs_only=false)
    # If an assembly loop function was provided, use it
    if !(assemble_loops === nothing)
        return assemble_loops.func(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t, dt);
    end
    # If parent maps were created for high order flux, use a slightly different function
    if !(parent_maps === nothing)
        return assemble_using_parent_child(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t, dt);
    end
    
    nel = fv_grid.nel_owned;
    
    # Label things that were allocated externally
    b = allocated_vecs[1];
    if !rhs_only
        AI = allocated_vecs[2];
        AJ = allocated_vecs[3];
        AV = allocated_vecs[4];
    end
    facefluxvec = allocated_vecs[5]
    face_done = allocated_vecs[6];
    # zero b
    b .= 0;
    face_done .= 0;
    
    # Elemental loop
    for ei=1:nel
        eid = elemental_order[ei]; # The index of this element
        # Each element owns (1+nfaces)*dofs_per_node elements of the A matrix
        nfaces = size(grid.element2face,1); # eventually this should be modified different element types
        first_index = (eid-1)*(nfaces+1)*dofs_per_node + 1;
        last_index = eid*(nfaces+1)*dofs_per_node;
        # Zero A for this element
        AI[first_index:last_index] .= 0;
        AJ[first_index:last_index] .= 0;
        AV[first_index:last_index] .= 0;
        
        ##### Source
        sourceargs = (var, eid, 0, fv_grid, fv_geo_factors, fv_info, refel, t, dt);
        if !(source_lhs === nothing)
            source = source_lhs.func(sourceargs);
            # These are the diagonal elements of A
            AI[first_index:first_index+dofs_per_node] = ((eid-1)*dofs_per_node + 1):eid*dofs_per_node;
            AJ[first_index:first_index+dofs_per_node] = ((eid-1)*dofs_per_node + 1):eid*dofs_per_node;
            AV[first_index:first_index+dofs_per_node] = source;
        end
        if !(source_rhs === nothing)
            source = source_rhs.func(sourceargs);
            # RHS vector b
            b[((eid-1)*dofs_per_node+1):eid*dofs_per_node] = source;
        end
        
        ##### Flux
        # Loop over this element's faces.
        for i=1:refel.Nfaces
            fid = fv_grid.element2face[i, eid];
            fbid = fv_grid.facebid[fid]; # BID of this face
            
            # Only one element on either side is available here. For more use parent/child version.
            (leftel, rightel) = fv_grid.face2element[:,fid];
            if rightel == 0
                neighborhood = [[leftel],[]];
                neighbor_id = eid;
            else
                neighborhood = [[leftel],[rightel]];
                if leftel == eid
                    neighbor_id = rightel;
                else
                    neighbor_id = leftel;
                end
            end
            
            if face_done[fid] == 0
                face_done[fid] = 1;
                
                fluxargs = (var, eid, fid, neighborhood, fv_grid, fv_geo_factors, fv_info, refel, t, dt);
                if !(flux_rhs === nothing)
                    flux = flux_rhs.func(fluxargs) .* fv_geo_factors.area[fid];
                    # Add to global flux vector for faces
                    facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] = flux;
                    # Combine all flux for this element
                    b[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] += flux ./ fv_geo_factors.volume[eid];
                end
                if !(flux_lhs === nothing)
                    flux = flux_lhs.func(fluxargs) .* fv_geo_factors.area[fid];
                    # These are the off-diagonal elements of A
                    AI[first_index:first_index+dofs_per_node] = ((eid-1)*dofs_per_node + 1):eid*dofs_per_node; # row for this element
                    AJ[first_index:first_index+dofs_per_node] = ((neighbor_id-1)*dofs_per_node + 1):neighbor_id*dofs_per_node; # col for neighbor
                    AV[first_index:first_index+dofs_per_node] = flux ./ fv_geo_factors.volume[eid];
                end
                
            else # This flux has either been computed or is being computed by another thread.
                 # The state will need to be known before paralellizing, but for now assume it's complete.
                if !(flux_rhs === nothing)
                    b[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] -= facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] ./ fv_geo_factors.volume[eid];
                end
                if !(flux_lhs === nothing)
                    fluxargs = (var, eid, fid, neighborhood, fv_grid, fv_geo_factors, fv_info, refel, t, dt);
                    flux = flux_lhs.func(fluxargs) .* fv_geo_factors.area[fid];
                    # These are the off-diagonal elements of A
                    AI[first_index:first_index+dofs_per_node] = ((eid-1)*dofs_per_node + 1):eid*dofs_per_node; # row for this element
                    AJ[first_index:first_index+dofs_per_node] = ((neighbor_id-1)*dofs_per_node + 1):neighbor_id*dofs_per_node; # col for neighbor
                    AV[first_index:first_index+dofs_per_node] = flux ./ fv_geo_factors.volume[eid];
                end
            end
            
            # Boundary conditions are applied to flux
            if fbid > 0
                if typeof(var) <: Array
                    dofind = 0;
                    for vi=1:length(var)
                        for compo=1:length(var[vi].symvar)
                            dofind = dofind + 1;
                            if prob.bc_type[var[vi].index, fbid] == NO_BC
                                # do nothing
                            elseif prob.bc_type[var[vi].index, fbid] == FLUX
                                # compute the value and add it to the flux directly
                                bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var[vi].index, fbid][compo], fid, t) .* fv_geo_factors.area[fid];
                                
                                b[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ fv_geo_factors.volume[eid];
                                facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                                # TODO How to modify LHS
                            elseif prob.bc_type[var[vi].index, fbid] == DIRICHLET
                                # Set variable array and handle after the face loop
                                var[vi].values[compo,eid] = FV_evaluate_bc(prob.bc_func[var[vi].index, fbid][compo], eid, fid, t);
                            else
                                printerr("Unsupported boundary condition type: "*prob.bc_type[var[vi].index, fbid]);
                            end
                        end
                    end
                else
                    for d=1:dofs_per_node
                        dofind = d;
                        if prob.bc_type[var.index, fbid] == NO_BC
                            # do nothing
                        elseif prob.bc_type[var.index, fbid] == FLUX
                            # compute the value and add it to the flux directly
                            # Qvec = (refel.surf_wg[fv_grid.faceRefelInd[1,fid]] .* fv_geo_factors.face_detJ[fid])' * (refel.surf_Q[fv_grid.faceRefelInd[1,fid]])[:, refel.face2local[fv_grid.faceRefelInd[1,fid]]]
                            # Qvec = Qvec ./ fv_geo_factors.area[fid];
                            # bflux = FV_flux_bc_rhs_only(prob.bc_func[var.index, fbid][d], facex, Qvec, t, dofind, dofs_per_node) .* fv_geo_factors.area[fid];
                            bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var.index, fbid][d], fid, t) .* fv_geo_factors.area[fid];
                            
                            b[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ fv_geo_factors.volume[eid];
                            facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                        elseif prob.bc_type[var.index, fbid] == DIRICHLET
                            # Set variable array and handle after the face loop
                            var.values[compo,eid] = FV_evaluate_bc(prob.bc_func[var.index, fbid][d], eid, fid, t);
                        else
                            printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, fbid]);
                        end
                    end
                end
            end# BCs
            
        end# face loop
        
    end# element loop
    
    # Construct the sparse matrix A
    A = Sparse(AI, AJ, AV);
    
    return (A, b);
end

function assemble_using_parent_child(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0)
    debug = false;
    
    # Label things that were allocated externally
    sourcevec = allocated_vecs[1];
    fluxvec = allocated_vecs[2];
    facefluxvec = allocated_vecs[3];
    face_done = allocated_vecs[4];
    
    face_done .= 0;
    fluxvec .= 0; # need to do this here
    
    # Loop over parent elements
    # Since flux depends only on this and nearest neighbor parent cells, it is safe 
    # to partition between parents using only the nearest parents as ghosts.
    Nparents = size(parent_maps.parent2child,2);
    nchildren = size(parent_maps.parent2child,1); # (assumes one element type)
    npfaces = size(fv_grid.element2face,1); # Outer faces of parent. Parents are elements in fv_grid (assumes one element type)
    ncfaces = size(parent_maps.parent2face,1); # All faces of children in a parent (assumes one element type)
    dim = size(fv_grid.allnodes,1);
    
    for parentid=1:Nparents
        # Only work on owned parents
        if !fv_grid.is_subgrid || fv_grid.element_owner[parent_maps.parent2child[1,parentid]] < 0
            # Get the cells in the local patch
            patch_cells = parent_maps.patches[:, parentid];
            
            # Split the work into two parts: 
            # - source terms that don't depend on neighbors and can be done independently (loop over elements)
            # - flux terms that depend on neighbors (loop over faces)
            
            # The loop over elements for source terms.
            for childid = 1:nchildren
                eid = parent_maps.parent2child[childid, parentid];
                # Zero the result vectors for this element
                sourcevec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] .= 0;
                
                # Use source_rhs for RHS volume integral
                if !(source_rhs === nothing)
                    sourceargs = (var, eid, 0, fv_grid, fv_geo_factors, fv_info, fv_refel, t, dt);
                    source = source_rhs.func(sourceargs);
                    # source = source_rhs.func(sourceargs) ./ fv_geo_factors.volume[eid];
                    # Add to global source vector
                    sourcevec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] = source;
                end
            end
            
            # Loop over this parent's faces for flux terms.
            for pfaceid = 1:ncfaces
                fid = parent_maps.parent2face[pfaceid, parentid];
                els = fv_grid.face2element[:,fid];
                eid = els[1];
                neighbor = els[2];
                
                if debug && parentid < 3 && t<=dt
                    println("parent: "*string(parentid)*", fid: "*string(fid)*"("*string(pfaceid)*"), eid/neigh: "*string(eid)*"/"*string(neighbor));
                    println("patch: "*string(patch_cells));
                end
                
                if !(flux_rhs === nothing)
                    if face_done[fid] == 0
                        face_done[fid] = 1; # Possible race condition, but in the worst case it will be computed twice.\
                        
                        (left_cells, right_cells) = get_left_right_cells(patch_cells, pfaceid, parent_maps, dim, fv_info.fluxOrder);
                        
                        # leftness and rightness depends on the normal vector direction
                        # normal points left to right.
                        if !(eid == left_cells[1])
                            #oops, swap them
                            tmp = left_cells;
                            left_cells = right_cells;
                            right_cells = tmp;
                        end
                        
                        # remove zero values
                        left_cells = remove_zero_cells(left_cells);
                        right_cells = remove_zero_cells(right_cells);
                        neighborhood = [left_cells, right_cells];
                        
                        fluxargs = (var, eid, fid, neighborhood, fv_grid, fv_geo_factors, fv_info, fv_refel, t, dt);
                        flux = flux_rhs.func(fluxargs) .* fv_geo_factors.area[fid];
                        
                        # Add to global flux vector for faces
                        facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] = flux;
                        # Contribute to flux for elements on both sides
                        if fv_grid.is_subgrid
                            if fv_grid.element_owner[eid] < 0
                                fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] += flux ./ fv_geo_factors.volume[eid];
                            end
                            if neighbor > 0 && fv_grid.element_owner[neighbor] < 0
                                fluxvec[((neighbor-1)*dofs_per_node + 1):(neighbor*dofs_per_node)] -= flux ./ fv_geo_factors.volume[neighbor];
                            end
                            
                        else
                            fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] += flux ./ fv_geo_factors.volume[eid];
                            if neighbor > 0 # not a boundary, so the other side exists
                                fluxvec[((neighbor-1)*dofs_per_node + 1):(neighbor*dofs_per_node)] -= flux ./ fv_geo_factors.volume[neighbor];
                            end
                        end
                        
                        if debug && parentid < 3 && t<=dt
                            println("cells left: "*string(left_cells)*", right: "*string(right_cells));
                            println("flux: "*string(flux));
                            println("contribute "*string(flux ./ fv_geo_factors.volume[eid])*" to "*string(((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)));
                            if neighbor > 0
                                println("contribute "*string(-flux ./ fv_geo_factors.volume[neighbor])*" to "*string(((neighbor-1)*dofs_per_node + 1):(neighbor*dofs_per_node)));
                            end 
                            println("fluxvec "*string(fluxvec[1:8]));
                            println("");
                        end
                        
                    else
                        # This flux has either been computed or is being computed by another thread.
                        # The state will need to be known before paralellizing, but for now assume it's complete.
                        # fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] .-= facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] ./ fv_geo_factors.volume[eid];
                    end
                end
                
                # Boundary conditions are applied to flux
                fbid = fv_grid.facebid[fid]; # BID of this face
                if fbid > 0
                    if debug && parentid < 3 && t==0
                        println("applying bdry cond "*string(fbid)*" to face "*string(fid));
                    end
                    # facex = fv_grid.allnodes[:, fv_grid.face2glb[:,1,fid]];  # face node coordinates
                    
                    if typeof(var) <: Array
                        dofind = 0;
                        for vi=1:length(var)
                            for compo=1:length(var[vi].symvar)
                                dofind = dofind + 1;
                                if prob.bc_type[var[vi].index, fbid] == NO_BC
                                    # do nothing
                                elseif prob.bc_type[var[vi].index, fbid] == FLUX
                                    # compute the value and add it to the flux directly
                                    # Qvec = (fv_refel.surf_wg[fv_grid.faceRefelInd[1,fid]] .* fv_geo_factors.face_detJ[fid])' * (fv_refel.surf_Q[fv_grid.faceRefelInd[1,fid]])[:, fv_refel.face2local[fv_grid.faceRefelInd[1,fid]]]
                                    # Qvec = Qvec ./ fv_geo_factors.area[fid];
                                    # bflux = FV_flux_bc_rhs_only(prob.bc_func[var[vi].index, fbid][compo], facex, Qvec, t, dofind, dofs_per_node) .* fv_geo_factors.area[fid];
                                    bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var[vi].index, fbid][compo], fid, t) .* fv_geo_factors.area[fid];
                                    
                                    fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ fv_geo_factors.volume[eid];
                                    facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                                elseif prob.bc_type[var[vi].index, fbid] == DIRICHLET
                                    # Set variable array and handle after the face loop
                                    var[vi].values[compo,eid] = FV_evaluate_bc(prob.bc_func[var[vi].index, fbid][compo], eid, fid, t);
                                else
                                    printerr("Unsupported boundary condition type: "*prob.bc_type[var[vi].index, fbid]);
                                end
                            end
                        end
                    else
                        for d=1:dofs_per_node
                            dofind = d;
                            if prob.bc_type[var.index, fbid] == NO_BC
                                # do nothing
                            elseif prob.bc_type[var.index, fbid] == FLUX
                                # compute the value and add it to the flux directly
                                # Qvec = (fv_refel.surf_wg[fv_grid.faceRefelInd[1,fid]] .* fv_geo_factors.face_detJ[fid])' * (fv_refel.surf_Q[fv_grid.faceRefelInd[1,fid]])[:, fv_refel.face2local[fv_grid.faceRefelInd[1,fid]]]
                                # Qvec = Qvec ./ fv_geo_factors.area[fid];
                                # bflux = FV_flux_bc_rhs_only(prob.bc_func[var.index, fbid][d], facex, Qvec, t, dofind, dofs_per_node) .* fv_geo_factors.area[fid];
                                bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var.index, fbid][d], fid, t) .* fv_geo_factors.area[fid];
                                
                                fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ fv_geo_factors.volume[eid];
                                facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                            elseif prob.bc_type[var.index, fbid] == DIRICHLET
                                # Set variable array and handle after the face loop
                                var.values[compo,eid] = FV_evaluate_bc(prob.bc_func[var.index, fbid][d], eid, fid, t);
                            else
                                printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, fbid]);
                            end
                        end
                    end
                end# BCs
            end# face loop
        end # is owned?
        
    end# element loop
    if debug && t<=dt
        println("fluxvec: "*string(fluxvec[1:20]));
    end
    return sourcevec + fluxvec;
end

# Returns a copy of a with zeros removed.
function remove_zero_cells(a)
    b = similar(a);
    nnz = 0;
    for ai in a
        if !(ai == 0)
            nnz += 1;
            b[nnz] = ai;
        end
    end
    b = b[1:nnz];
    
    return b;
end

# Returns arrays of left and right cells
function get_left_right_cells(patch, face, maps, dim, order)
    # Provide all of the cells that can be used for this face in a particular order.
    # Make two arrays, one for left one for right, starting from the nearest neighbor.
    left_cells = [];
    right_cells = [];
    if dim == 1
        left_cell_table_1d = [
            [ # children=1
                [2],
                [1, 2]
            ],
            [ # children=2
                [4, 3],
                [1, 4, 3],
                [2, 1, 4, 3]
            ],
            [ # children=3
                [6, 5, 4],
                [1, 6, 5, 4],
                [2, 1, 6, 5, 4],
                [3, 2, 1, 6, 5, 4]
            ],
            [ # children=4
                [8, 7, 6, 5],
                [1, 8, 7, 6, 5],
                [2, 1, 8, 7, 6, 5],
                [3, 2, 1, 8, 7, 6, 5],
                [4, 3, 2, 1, 8, 7, 6, 5],
            ]
        ]
        right_cell_table_1d = [
            [ # children=1
                [1, 3],
                [3]
            ],
            [ # children=2
                [1, 2, 5, 6],
                [2, 5, 6],
                [5, 6]
            ],
            [ # children=3
                [1, 2, 3, 7, 8, 9],
                [2, 3, 7, 8, 9],
                [3, 7, 8, 9],
                [7, 8, 9]
            ],
            [ # children=4
                [1, 2, 3, 4, 9, 10, 11, 12],
                [2, 3, 4, 9, 10, 11, 12],
                [3, 4, 9, 10, 11, 12],
                [4, 9, 10, 11, 12],
                [9, 10, 11, 12]
            ]
        ]
        nchildren = size(maps.parent2child,1);
        left_cells = left_cell_table_1d[nchildren][face];
        right_cells = right_cell_table_1d[nchildren][face];
        
    elseif dim == 2
        if size(maps.parent2neighbor,1) == 3 # triangles
            # Triangle parents have 9 faces, patches have 16 cells
            # Here Left means toward the center of the central parent
            left_cell_table_triangle = [
                [1, 4, 3, 15, 16, 13, 2],
                [2, 4, 3, 9, 12, 11, 1],
                [2, 4, 1, 7, 8, 5, 3],
                [3, 4, 1, 13, 16, 15, 2],
                [3, 4, 2, 11, 12, 9, 1],
                [1, 4, 2, 5, 8, 7, 3],
                [4, 2, 3, 9, 11, 12],
                [4, 1, 3, 15, 13, 16],
                [4, 1, 2, 5, 7, 8]
            ]
            right_cell_table_triangle = [
                [5, 8, 6, 7],
                [7, 8, 6, 5],
                [9, 12, 10, 11],
                [11, 12, 10, 9],
                [13, 16, 14, 15],
                [15, 16, 14, 13],
                [1, 5, 15, 8, 16],
                [2, 7, 9, 8, 12],
                [3, 11, 13, 12, 16]
            ]
            left_cells = left_cell_table_triangle[face];
            right_cells = right_cell_table_triangle[face];
        else # quads
            # Quad parents have 12 faces, patches have 20 cells
            # Here Left means toward the center of the central parent
            left_cell_table_quad = [
                [1, 4, 2, 16, 15, 3, 13, 14],
                [2, 3, 1, 13, 14, 4, 16, 15],
                [2, 1, 3, 20, 19, 4, 17, 18],
                [3, 4, 2, 17, 18, 1, 20, 19],
                [3, 2, 4, 8, 7, 1, 5, 6],
                [4, 1, 3, 5, 6, 2, 8, 7],
                [4, 3, 1, 12, 11, 2, 9, 10],
                [1, 2, 4, 9, 10, 3, 12, 11],
                [1, 20, 4, 19, 17, 18, 5],
                [2, 8, 1, 7, 5, 6, 9],
                [3, 12, 2, 11, 9, 10, 13],
                [4, 16, 3, 15, 13, 14, 17]
            ]
            right_cell_table_quad = [
                [5, 6, 8, 7],
                [8, 7, 5, 6],
                [9, 10, 12, 11],
                [12, 11, 9, 10],
                [13, 14, 16, 15],
                [16, 15, 13, 14],
                [17, 18, 20, 19],
                [20, 19, 17, 18],
                [2, 9, 3, 10, 12, 11],
                [3, 13, 4, 14, 16, 15],
                [4, 17, 1, 18, 20, 19],
                [1, 5, 2, 6, 8, 7]
            ]
            left_cells = left_cell_table_quad[face];
            right_cells = right_cell_table_quad[face];
        end
        
    elseif dim == 3
        #TODO
    end
    
    if length(left_cells) > order
        left_cells = left_cells[1:order];
    end
    if length(right_cells) > order
        right_cells = right_cells[1:order];
    end
    
    return (patch[left_cells], patch[right_cells]);
end

# Insert the single dof into the greater vector
function insert_linear!(b, bel, glb, dof, Ndofs)
    # group nodal dofs
    for d=1:length(dof)
        ind = glb.*Ndofs .- (Ndofs-dof[d]);
        ind2 = ((d-1)*length(glb)+1):(d*length(glb));
        
        b[ind] = b[ind] + bel[ind2];
    end
end

# place the values in the variable value arrays
function place_sol_in_vars(var, sol, stepper)
    nel = fv_grid.nel_owned;
    if typeof(var) <: Array
        tmp = 0;
        totalcomponents = 0;
        for vi=1:length(var)
            totalcomponents = totalcomponents + var[vi].total_components;
        end
        for vi=1:length(var)
            components = var[vi].total_components;
            for compi=1:components
                var[vi].values[compi,1:nel] = sol[(compi+tmp):totalcomponents:end];
                tmp = tmp + 1;
            end
        end
    else
        components = var.total_components;
        for compi=1:components
            var.values[compi,1:nel] = sol[compi:components:end];
        end
    end
end

# place the variable values in a vector
function get_var_vals(var, vect=nothing)
    nel = fv_grid.nel_owned;
    if typeof(var) <: Array
        tmp = 0;
        totalcomponents = 0;
        for vi=1:length(var)
            totalcomponents = totalcomponents + var[vi].total_components;
        end
        if vect === nothing
            vect = zeros(totalcomponents * nel);
        end
        
        for vi=1:length(var)
            components = var[vi].total_components;
            for compi=1:components
                vect[(compi+tmp):totalcomponents:end] = var[vi].values[compi,1:nel];
                tmp = tmp + 1;
            end
        end
    else
        components = var.total_components;
        if vect === nothing
            vect = zeros(components * nel);
        end
        for compi=1:components
            vect[compi:components:end] = var.values[compi,1:nel];
        end
    end
    
    return vect;
end

end #module
