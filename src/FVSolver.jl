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

function linear_solve(var, source_lhs, source_rhs, flux_lhs, flux_rhs, stepper=nothing, assemble_func=nothing)
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
        nel = size(grid_data.loc2glb, 2);
        nfaces = size(grid_data.face2element, 2);
    else
        grid = fv_grid;
        nel = size(fv_grid.loc2glb, 2);
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
    face_done = zeros(Int, nfaces); # Increment when the corresponding flux value is computed.
    allocated_vecs = [sourcevec, fluxvec, facefluxvec, face_done];
    
    if prob.time_dependent && !(stepper === nothing)
        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
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
        print("Time stepping progress(%): 0");
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
                        
                        pre_step_function();
                        
                        sol = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, stime, stepper.dt, assemble_loops=assemble_func);
                        
                        
                        if rki == 1 # because a1 == 0
                            tmpki = stepper.dt .* sol;
                        else
                            tmpki = stepper.a[rki].*tmpki + stepper.dt.*sol;
                        end
                        tmppi = tmppi + stepper.b[rki].*tmpki
                        
                        copy_bdry_vals_to_vector(var, tmppi, grid, dofs_per_node);
                        place_sol_in_vars(var, tmppi, stepper);
                        
                        post_step_function();
                    end
                    
                else
                    # Explicit multi-stage methods: 
                    # x = x + dt*sum(bi*ki)
                    # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
                    
                    # will hold the final result
                    sol = get_var_vals(var, sol);
                    # will be placed in var.values for each stage
                    tmpvals = sol;
                    for stage=1:stepper.stages
                        stime = t + stepper.c[stage]*stepper.dt;
                        
                        pre_step_function();
                        
                        tmpki[:,stage] = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, stime, stepper.dt, assemble_loops=assemble_func);
                        
                        tmpvals = sol;
                        for j=1:(stage-1)
                            if stepper.a[stage, j] > 0
                                tmpvals += stepper.dt * stepper.a[stage, j] .* tmpki[:,j];
                            end
                        end
                        
                        copy_bdry_vals_to_vector(var, tmpvals, grid, dofs_per_node);
                        place_sol_in_vars(var, tmpvals, stepper);
                        
                        post_step_function();
                    end
                    for stage=1:stepper.stages
                        sol += stepper.dt * stepper.b[stage] .* tmpki[:, stage];
                    end
                    copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node);
                    place_sol_in_vars(var, sol, stepper);
                    
                    post_step_function();
                end
                
            elseif stepper.type == EULER_EXPLICIT
                pre_step_function();
                
                sol = sol .+ stepper.dt .* assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t, stepper.dt, assemble_loops=assemble_func);
                
                copy_bdry_vals_to_vector(var, sol, grid, dofs_per_node);
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
        println("");
        end_t = Base.Libc.time();
        
        log_entry("Stepping took "*string(end_t-start_t)*" seconds.");
        #log_entry("Stepping took "*string(end_t-start_t)*" seconds. ("*string(assemble_t)*" for assembly, "*string(linsolve_t)*" for linear solve)");
        #display(sol);
		# outfile = "linear_sol.txt"
		# open(outfile, "w") do f
  		# 	for ii in sol
    	# 		println(f, ii)
  		# 	end
		# end # the file f is automatically closed after this block finishes
        return sol;

    else
        # Does it make any sense to do this for time-independent problems?
        return sol;
    end
end

# macro loops(indices, ranges, content)
#     n = length(indices.args);
#     if n == 1
#         this_ind = indices.args[1];
#         this_range = ranges.args[1];
#         return esc(quote
#             for $this_ind in $this_range
#                 $content
#             end
#         end)
#     else
#         this_ind = indices.args[1];
#         this_range = ranges.args[1];
#         sub_ind = copy(indices);
#         sub_ind.args = sub_ind.args[2:end];
#         sub_range = copy(ranges);
#         sub_range.args = sub_range.args[2:end];
#         return esc(:(@loops([$this_ind], [$this_range], @loops($sub_ind, $sub_range, $content))));
#     end
# end

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
    
    nel = size(grid_data.loc2glb, 2);
    
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
            sourceargs = (var, eid, 0, grid_data, geo_factors, fv_info, refel, t, dt);
            source = source_rhs.func(sourceargs) ./ geo_factors.volume[eid];
            # Add to global source vector
            sourcevec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] = source;
        end
        
        ##### Flux integrated over the faces #####
        # Loop over this element's faces.
        for i=1:refel.Nfaces
            fid = grid_data.element2face[i, eid];
            # Only one element on either side is available here. For more use parent/child version.
            (leftel, rightel) = grid_data.face2element[:,fid];
            if rightel == 0
                neighborhood = [[leftel],[]];
            else
                neighborhood = [[leftel],[rightel]];
            end
            
            if !(flux_rhs === nothing)
                if face_done[fid] == 0
                    face_done[fid] = 1; # Possible race condition, but in the worst case it will be computed twice.
                    
                    #fluxargs = prepare_args(var, eid, fid, RHS, "surface", t, dt); #(var, (e, neighbor), refel, vol_loc2glb, nodex, cellx, frefelind, facex, face2glb, normal, fdetj, face_area, (J, vol_J_neighbor), t, dt);
                    fluxargs = (var, eid, fid, neighborhood, grid_data, geo_factors, fv_info, refel, t, dt);
                    flux = flux_rhs.func(fluxargs) .* geo_factors.area[fid];
                    # Add to global flux vector for faces
                    facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] = flux;
                    # Combine all flux for this element
                    fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] += flux ./ geo_factors.volume[eid];
                    
                else
                    # This flux has either been computed or is being computed by another thread.
                    # The state will need to be known before paralellizing, but for now assume it's complete.
                    fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] -= facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] ./ geo_factors.volume[eid];
                end
            end
            
            # Boundary conditions are applied to flux
            fbid = grid_data.facebid[fid]; # BID of this face
            if fbid > 0
                facex = grid_data.allnodes[:, grid_data.face2glb[:,1,fid]];  # face node coordinates
                
                if typeof(var) <: Array
                    dofind = 0;
                    for vi=1:length(var)
                        for compo=1:length(var[vi].symvar)
                            dofind = dofind + 1;
                            if prob.bc_type[var[vi].index, fbid] == NO_BC
                                # do nothing
                            elseif prob.bc_type[var[vi].index, fbid] == FLUX
                                # compute the value and add it to the flux directly
                                # Qvec = (refel.surf_wg[grid_data.faceRefelInd[1,fid]] .* geo_factors.face_detJ[fid])' * (refel.surf_Q[grid_data.faceRefelInd[1,fid]])[:, refel.face2local[grid_data.faceRefelInd[1,fid]]]
                                # Qvec = Qvec ./ geo_factors.area[fid];
                                # bflux = FV_flux_bc_rhs_only(prob.bc_func[var[vi].index, fbid][compo], facex, Qvec, t, dofind, dofs_per_node) .* geo_factors.area[fid];
                                bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var[vi].index, fbid][compo], fid, t) .* geo_factors.area[fid];
                                
                                fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ geo_factors.volume[eid];
                                facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                            elseif prob.bc_type[var[vi].index, fbid] == DIRICHLET
                                # Set variable array and handle after the face loop
                                var[vi].values[compo,eid] = evaluate_bc(prob.bc_func[var[vi].index, fbid][compo], eid, fid, t);
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
                            # Qvec = (refel.surf_wg[grid_data.faceRefelInd[1,fid]] .* geo_factors.face_detJ[fid])' * (refel.surf_Q[grid_data.faceRefelInd[1,fid]])[:, refel.face2local[grid_data.faceRefelInd[1,fid]]]
                            # Qvec = Qvec ./ geo_factors.area[fid];
                            # bflux = FV_flux_bc_rhs_only(prob.bc_func[var.index, fbid][d], facex, Qvec, t, dofind, dofs_per_node) .* geo_factors.area[fid];
                            bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var.index, fbid][d], fid, t) .* geo_factors.area[fid];
                            
                            fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ geo_factors.volume[eid];
                            facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                        elseif prob.bc_type[var.index, fbid] == DIRICHLET
                            # Set variable array and handle after the face loop
                            var.values[compo,eid] = evaluate_bc(prob.bc_func[var.index, fbid][d], eid, fid, t);
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
    npfaces = size(grid_data.element2face,1); # Outer faces of parent. Parents are elements in grid_data (assumes one element type)
    ncfaces = size(parent_maps.parent2face,1); # All faces of children in a parent (assumes one element type)
    dim = size(grid_data.allnodes,1);
    
    for parentid=1:Nparents
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
                source = source_rhs.func(sourceargs; patch_cells=patch_cells) ./ fv_geo_factors.volume[eid];
                # Add to global source vector
                sourcevec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] = source;
            end
        end
        
        # Loop over this parent's faces for flux terms.
        for pfaceid = 1:ncfaces
            fid = parent_maps.parent2face[pfaceid, parentid];
            els = fv_grid.face2element[:,fid];
            # # The eid passed to the flux func. will be one inside the parent.
            # if els[1] in parent_maps.parent2child[:,parentid]
            #     eid = els[1];
            #     neighbor = els[2];
            # else
            #     eid = els[2];
            #     neighbor = els[1];
            # end
            eid = els[1];
            neighbor = els[2];
            
            if debug && parentid < 3
                println("parent: "*string(parentid)*", fid: "*string(fid)*"("*string(pfaceid)*"), eid/neigh: "*string(eid)*"/"*string(neighbor));
                println("patch: "*string(patch_cells));
            end
            
            if !(flux_rhs === nothing)
                if face_done[fid] == 0
                    face_done[fid] = 1; # Possible race condition, but in the worst case it will be computed twice.\
                    
                    (left_cells, right_cells) = get_left_right_cells(patch_cells, pfaceid, parent_maps, dim);
                    
                    # leftness and rightness depends on the normal vector direction
                    # normal points left to right.
                    leftcell = fv_grid.face2element[1,fid];
                    if leftcell == right_cells[1]
                        #oops, swap them
                        tmp = left_cells;
                        left_cells = right_cells;
                        right_cells = tmp;
                    end
                    
                    # remove zero values
                    left_cells = remove_zero_cells(left_cells);
                    right_cells = remove_zero_cells(right_cells);
                    neighborhood = [left_cells, right_cells];
                    
                    #fluxargs = prepare_args(var, eid, fid, RHS, "surface", t, dt); #(var, (e, neighbor), refel, vol_loc2glb, nodex, cellx, frefelind, facex, face2glb, normal, fdetj, face_area, (J, vol_J_neighbor), t, dt);
                    fluxargs = (var, eid, fid, neighborhood, fv_grid, fv_geo_factors, fv_info, fv_refel, t, dt);
                    flux = flux_rhs.func(fluxargs) .* fv_geo_factors.area[fid];
                    
                    # Add to global flux vector for faces
                    facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] = flux;
                    # Contribute to flux for elements on both sides
                    fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] += flux ./ fv_geo_factors.volume[eid];
                    if neighbor > 0 # not a boundary, so the other side exists
                        fluxvec[((neighbor-1)*dofs_per_node + 1):(neighbor*dofs_per_node)] -= flux ./ fv_geo_factors.volume[neighbor];
                    end
                    
                    if debug && parentid < 3
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
                if debug 
                    println("applying bdry cond"*string(fbid)*" to "*string(fid));
                end
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
                                # Qvec = (fv_refel.surf_wg[fv_grid.faceRefelInd[1,fid]] .* fv_geo_factors.face_detJ[fid])' * (fv_refel.surf_Q[fv_grid.faceRefelInd[1,fid]])[:, fv_refel.face2local[fv_grid.faceRefelInd[1,fid]]]
                                # Qvec = Qvec ./ fv_geo_factors.area[fid];
                                # bflux = FV_flux_bc_rhs_only(prob.bc_func[var[vi].index, fbid][compo], facex, Qvec, t, dofind, dofs_per_node) .* fv_geo_factors.area[fid];
                                bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var[vi].index, fbid][compo], fid, t) .* geo_factors.area[fid];
                                
                                fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ fv_geo_factors.volume[eid];
                                facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                            elseif prob.bc_type[var[vi].index, fbid] == DIRICHLET
                                # Set variable array and handle after the face loop
                                var[vi].values[compo,eid] = evaluate_bc(prob.bc_func[var[vi].index, fbid][compo], eid, fid, t);
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
                            bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var.index, fbid][d], fid, t) .* geo_factors.area[fid];
                            
                            fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ fv_geo_factors.volume[eid];
                            facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                        elseif prob.bc_type[var.index, fbid] == DIRICHLET
                            # Set variable array and handle after the face loop
                            var.values[compo,eid] = evaluate_bc(prob.bc_func[var.index, fbid][d], eid, fid, t);
                        else
                            printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, fbid]);
                        end
                    end
                end
            end# BCs
        end# face loop
    end# element loop
    if debug
        println("fluxvec: "*string(fluxvec[1:10]));
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
function get_left_right_cells(patch, face, maps, dim)
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
                [1, 4, 3, 15, 16, 13, 2, ],
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
            #TODO
        end
        
    elseif dim == 3
        #TODO
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

function place_sol_in_vars(var, sol, stepper)
    # place the values in the variable value arrays
    if typeof(var) <: Array
        tmp = 0;
        totalcomponents = 0;
        for vi=1:length(var)
            totalcomponents = totalcomponents + var[vi].total_components;
        end
        for vi=1:length(var)
            components = var[vi].total_components;
            for compi=1:components
                if stepper.type == EULER_EXPLICIT
                    var[vi].values[compi,:] = sol[(compi+tmp):totalcomponents:end];
                else 
                    var[vi].values[compi,:] = sol[(compi+tmp):totalcomponents:end];
                end
                tmp = tmp + 1;
            end
        end
    else
        components = var.total_components;
        for compi=1:components
            if stepper.type == EULER_EXPLICIT
                var.values[compi,:] = sol[compi:components:end];
            else
                var.values[compi,:] = sol[compi:components:end];
            end
        end
    end
end

function get_var_vals(var, vect=nothing)
    # place the variable values in a vector
    if typeof(var) <: Array
        tmp = 0;
        totalcomponents = 0;
        for vi=1:length(var)
            totalcomponents = totalcomponents + var[vi].total_components;
        end
        if vect === nothing
            vect = zeros(totalcomponents * length(var[1].values[1,:]));
        end
        
        for vi=1:length(var)
            components = var[vi].total_components;
            for compi=1:components
                vect[(compi+tmp):totalcomponents:end] = var[vi].values[compi,:];
                tmp = tmp + 1;
            end
        end
    else
        components = var.total_components;
        if vect === nothing
            vect = zeros(components * length(var.values[1,:]));
        end
        for compi=1:components
            vect[compi:components:end] = var.values[compi,:];
        end
    end
    
    return vect;
end

end #module
