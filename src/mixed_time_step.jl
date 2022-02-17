#=
Time steps for the MixedSolver module.
Since the FE and FV portions may use different steppers, this must work for both,
but they can be done separately by only passing one set of variables.
For all of them:
input
- fe_var
- fv_var
- vol_lhs
- vol_rhs
- surf_lhs
- surf_rhs
- step_info = (t, dt, dofs_per_node_fe, dofs_per_loop_fe, assemble_func_fe, dofs_per_node_fv, dofs_per_loop_fv, assemble_func_fv)
- fe_storage = (A, b, fe_sol, allocated_vecs, b_order, b_sizes)
- fv_storage = (fv_sol, allocated_vecs_fv)
=#

function mixed_euler_explicit_step(fe_var, fv_var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, t, step_info, fe_storage, fv_storage)
    fe_stepper = step_info[1];
    fv_stepper = step_info[2];
    dofs_per_node_fe = step_info[3];
    dofs_per_loop_fe = step_info[4];
    fe_assemble_func = step_info[5];
    dofs_per_node_fv = step_info[6];
    dofs_per_loop_fv = step_info[7];
    fv_assemble_func = step_info[8];
    dt = fe_stepper.dt;
    assemble_t = 0;
    linsolve_t = 0;
    
    if !(fv_var === nothing) && config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, i); end
    pre_step_function();
    
    ### First step the FE part #######################################################################################
    if !(fe_var === nothing)
        # Unpack things
        A = fe_storage[1];
        b = fe_storage[2];
        fe_sol = fe_storage[3];
        allocated_vecs = fe_storage[4];
        b_order = fe_storage[5];
        b_sizes = fe_storage[6];
        
        N1 = size(grid_data.allnodes,2);
        
        assemble_t += @elapsed begin
                                    b = assemble(fe_var, nothing, vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, t, dt; rhs_only = true, assemble_loops=fe_assemble_func);
                                    b = gather_system(nothing, b, N1, dofs_per_node_fe, b_order, b_sizes);
                                end
        
        linsolve_t += @elapsed(tmpvec = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node_fe, b_order, b_sizes));
        
        # At this point tmpvec holds the boundary values
        # directly write them to the variable values and zero sol.
        copy_bdry_vals_to_variables(fe_var, tmpvec, grid_data, dofs_per_node_fe, zero_vals=true);
        
        for i=1:length(fe_sol)
            fe_sol[i] = fe_sol[i] + dt * tmpvec[i];
        end
        
        copy_bdry_vals_to_vector(fe_var, fe_sol, grid_data, dofs_per_node_fe);
        place_sol_in_vars(fe_var, fe_sol);
        
    end
    
    ### Then step the FV part ######################################################################################
    if !(fv_var === nothing)
        fv_sol = fv_storage[1];
        fv_allocated_vecs = fv_storage[2];
        
        assemble_t += @elapsed(tmpvec = fv_assemble(fv_var, vol_lhs[2], vol_rhs[2], surf_lhs[2], surf_rhs[2], fv_allocated_vecs, dofs_per_node_fv, dofs_per_loop_fv, t, dt, assemble_loops=fv_assemble_func));
        
        for i=1:length(fv_sol)
            fv_sol[i] = fv_sol[i] + dt * tmpvec[i];
        end
        
        FV_copy_bdry_vals_to_vector(fv_var, fv_sol, fv_grid, dofs_per_node_fv);
        place_sol_in_vars(fv_var, fv_sol);
    end
    
    post_step_function();
    
    return (assemble_t, linsolve_t);
end

function mixed_euler_implicit_step(fe_var, fv_var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, t, step_info, fe_storage, fv_storage)
    fe_stepper = step_info[1];
    fv_stepper = step_info[2];
    dofs_per_node_fe = step_info[3];
    dofs_per_loop_fe = step_info[4];
    fe_assemble_func = step_info[5];
    dofs_per_node_fv = step_info[6];
    dofs_per_loop_fv = step_info[7];
    fv_assemble_func = step_info[8];
    dt = fe_stepper.dt;
    assemble_t = 0;
    linsolve_t = 0;
    
    if !(fv_var === nothing) && config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, i); end
    pre_step_function();
    
    ### First step the FE part #######################################################################################
    if !(fe_var === nothing)
        # Unpack things
        A = fe_storage[1];
        b = fe_storage[2];
        fe_sol = fe_storage[3];
        allocated_vecs = fe_storage[4];
        b_order = fe_storage[5];
        b_sizes = fe_storage[6];
        
        N1 = size(grid_data.allnodes,2);
        
        assemble_t += @elapsed begin
                                    b = assemble(fe_var, nothing, vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, t, dt; rhs_only = true, assemble_loops=fe_assemble_func);
                                    b = gather_system(nothing, b, N1, dofs_per_node_fe, b_order, b_sizes);
                                end
        
        linsolve_t += @elapsed(tmpvec = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node_fe, b_order, b_sizes));
        
        # At this point tmpvec holds the boundary values
        # directly write them to the variable values and zero sol.
        copy_bdry_vals_to_variables(fe_var, tmpvec, grid_data, dofs_per_node_fe, zero_vals=true);
        
        for i=1:length(fe_sol)
            fe_sol[i] = tmpvec[i];
        end
        
        copy_bdry_vals_to_vector(fe_var, fe_sol, grid_data, dofs_per_node_fe);
        place_sol_in_vars(fe_var, fe_sol);
    end
    
    ### Then step the FV part ######################################################################################
    if !(fv_var === nothing)
        fv_sol = fv_storage[1];
        fv_allocated_vecs = fv_storage[2];
        
        ## TODO
        
        FV_copy_bdry_vals_to_vector(fv_var, fv_sol, fv_grid, dofs_per_node_fv);
        place_sol_in_vars(fv_var, fv_sol);
    end
    
    post_step_function();
    
    return (assemble_t, linsolve_t);
end

function mixed_crank_nicholson_step(fe_var, fv_var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, t, step_info, fe_storage, fv_storage)
    fe_stepper = step_info[1];
    fv_stepper = step_info[2];
    dofs_per_node_fe = step_info[3];
    dofs_per_loop_fe = step_info[4];
    fe_assemble_func = step_info[5];
    dofs_per_node_fv = step_info[6];
    dofs_per_loop_fv = step_info[7];
    fv_assemble_func = step_info[8];
    dt = fe_stepper.dt;
    assemble_t = 0;
    linsolve_t = 0;
    
    if !(fv_var === nothing) && config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, i); end
    pre_step_function();
    
    ### First step the FE part #######################################################################################
    if !(fe_var === nothing)
        # Unpack things
        A = fe_storage[1];
        b = fe_storage[2];
        fe_sol = fe_storage[3];
        allocated_vecs = fe_storage[4];
        b_order = fe_storage[5];
        b_sizes = fe_storage[6];
        
        N1 = size(grid_data.allnodes,2);
        
        assemble_t += @elapsed begin
            b = assemble(fe_var, nothing, vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, t, dt; rhs_only = true, assemble_loops=fe_assemble_func);
            b = gather_system(nothing, b, N1, dofs_per_node_fe, b_order, b_sizes);
        end

        linsolve_t += @elapsed(tmpvec = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node_fe, b_order, b_sizes));

        # At this point tmpvec holds the boundary values
        # directly write them to the variable values and zero sol.
        copy_bdry_vals_to_variables(fe_var, tmpvec, grid_data, dofs_per_node_fe, zero_vals=true);

        for i=1:length(fe_sol)
        fe_sol[i] = tmpvec[i];
        end
        
        copy_bdry_vals_to_vector(fe_var, fe_sol, grid_data, dofs_per_node_fe);
        place_sol_in_vars(fe_var, fe_sol);
    end
    
    ### Then step the FV part ######################################################################################
    if !(fv_var === nothing)
        fv_sol = fv_storage[1];
        fv_allocated_vecs = fv_storage[2];
        
        ## TODO
        
        FV_copy_bdry_vals_to_vector(fv_var, fv_sol, fv_grid, dofs_per_node_fv);
        place_sol_in_vars(fv_var, fv_sol);
    end
    
    post_step_function();
    
    return (assemble_t, linsolve_t);
end

function mixed_lsrk4_step(fe_var, fv_var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, t, step_info, fe_storage, fv_storage)
    fe_stepper = step_info[1];
    fv_stepper = step_info[2];
    dofs_per_node_fe = step_info[3];
    dofs_per_loop_fe = step_info[4];
    fe_assemble_func = step_info[5];
    dofs_per_node_fv = step_info[6];
    dofs_per_loop_fv = step_info[7];
    fv_assemble_func = step_info[8];
    dt = fe_stepper.dt;
    assemble_t = 0;
    linsolve_t = 0;
    
    ### First step the FE part #######################################################################################
    if !(fe_var === nothing)
        # Unpack things
        A = fe_storage[1];
        b = fe_storage[2];
        fe_sol = fe_storage[3];
        allocated_vecs = fe_storage[4];
        b_order = fe_storage[5];
        b_sizes = fe_storage[6];
        tmppi = fe_storage[7];
        tmpki = fe_storage[8];
        
        N1 = size(grid_data.allnodes,2);
        
        # Low storage RK4: 
        # p0 = u
        #   ki = ai*k(i-1) + dt*f(p(i-1), t+ci*dt)
        #   pi = p(i-1) + bi*ki
        # u = p5
        
        tmppi = get_var_vals(fe_var, tmppi);
        for rki=1:fe_stepper.stages
            rktime = t + fe_stepper.c[rki]*fe_stepper.dt;
            # p(i-1) is currently in u
            pre_step_function();
            
            ### First step the FE part #######################################################################################
            assemble_t += @elapsed begin
                                b = assemble(fe_var, nothing, vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, rktime, fe_stepper.dt; rhs_only = true, assemble_loops=fe_assemble_func);
                                b = gather_system(nothing, b, N1, dofs_per_node_fe, b_order, b_sizes);
                            end
            
            linsolve_t += @elapsed(fe_sol = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node_fe, b_order, b_sizes));
            
            # At this point sol holds the boundary values
            # directly write them to the variable values and zero sol.
            copy_bdry_vals_to_variables(fe_var, fe_sol, grid_data, dofs_per_node_fe, zero_vals=true);
            
            if rki == 1 # because a1 == 0
                tmpki = fe_stepper.dt .* fe_sol;
            else
                tmpki = fe_stepper.a[rki].*tmpki + fe_stepper.dt.*fe_sol;
            end
            tmppi = tmppi + fe_stepper.b[rki].*tmpki
            
            copy_bdry_vals_to_vector(fe_var, tmppi, grid_data, dofs_per_node_fe);
            place_sol_in_vars(fe_var, tmppi);
            
            post_step_function();
        end
        
    end
    
    ### Then step the FV part ######################################################################################
    if !(fv_var === nothing)
        fv_sol = fv_storage[1];
        fv_allocated_vecs = fv_storage[2];
        fv_tmppi = fv_storage[3];
        fv_tmpki = fv_storage[4];
        
        fv_tmppi = get_var_vals(fv_var, fv_tmppi);
        for rki=1:fv_stepper.stages
            rktime = t + fv_stepper.c[rki]*fv_stepper.dt;
            # p(i-1) is currently in u
            if config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, i); end
            pre_step_function();
            
            assemble_t += @elapsed(fv_sol = fv_assemble(fv_var, vol_lhs[2], vol_rhs[2], surf_lhs[2], surf_rhs[2], fv_allocated_vecs, dofs_per_node_fv, dofs_per_loop_fv, rktime, dt, assemble_loops=fv_assemble_func));
                
            if rki == 1 # because a1 == 0
                fv_tmpki = fv_stepper.dt .* fv_sol;
            else
                fv_tmpki = fv_stepper.a[rki].*fv_tmpki + fv_stepper.dt.*fv_sol;
            end
            fv_tmppi = fv_tmppi + fv_stepper.b[rki].*fv_tmpki
            
            FV_copy_bdry_vals_to_vector(fv_var, fv_tmppi, fv_grid, dofs_per_node_fv);
            place_sol_in_vars(fv_var, fv_tmppi);
            
            post_step_function();
        end
    end
    
    return (assemble_t, linsolve_t);
end

# This may be removed.
function mixed_rk4_step(fe_var, fv_var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, t, step_info, fe_storage, fv_storage)
    return mixed_multistage_step(fe_var, fv_var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, t, step_info, fe_storage, fv_storage);
end

function mixed_multistage_step(fe_var, fv_var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, t, step_info, fe_storage, fv_storage)
    fe_stepper = step_info[1];
    fv_stepper = step_info[2];
    dofs_per_node_fe = step_info[3];
    dofs_per_loop_fe = step_info[4];
    fe_assemble_func = step_info[5];
    dofs_per_node_fv = step_info[6];
    dofs_per_loop_fv = step_info[7];
    fv_assemble_func = step_info[8];
    dt = fe_stepper.dt;
    assemble_t = 0;
    linsolve_t = 0;
    
    ### First step the FE part #######################################################################################
    if !(fe_var === nothing)
        # Unpack things
        A = fe_storage[1];
        b = fe_storage[2];
        fe_sol = fe_storage[3];
        allocated_vecs = fe_storage[4];
        b_order = fe_storage[5];
        b_sizes = fe_storage[6];
        tmpresult = fe_storage[7];
        tmpki = fe_storage[8];
        
        N1 = size(grid_data.allnodes,2);
        
        # solution will be placed in var.values for each stage
        for stage=1:fe_stepper.stages
            stime = t + fe_stepper.c[stage]*fe_stepper.dt;
            
            # Update the values in vars to be used in this stage
            if stage > 1
                initialized_tmpresult = false;
                for j=1:stage
                    if fe_stepper.a[stage, j] > 0
                        if !initialized_tmpresult
                            initialized_tmpresult = true;
                            for k=1:length(fe_sol)
                                tmpresult[k] = fe_sol[k] + fe_stepper.dt * fe_stepper.a[stage, j] * tmpki[k,j];
                            end
                        else
                            for k=1:length(fe_sol)
                                tmpresult[k] += fe_stepper.dt * fe_stepper.a[stage, j] * tmpki[k,j];
                            end
                        end
                    end
                end
                
                if initialized_tmpresult
                    copy_bdry_vals_to_vector(fe_var, tmpresult, grid_data, dofs_per_node_fe);
                    place_sol_in_vars(fe_var, tmpresult);
                end
                post_step_function(); # seems weird, but imagine this is happening after stage-1
            end
            
            pre_step_function();
            
            ### First step the FE part #######################################################################################
            assemble_t += @elapsed begin
                                b = assemble(fe_var, nothing, vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, stime, fe_stepper.dt; rhs_only = true, assemble_loops=fe_assemble_func);
                                b = gather_system(nothing, b, N1, dofs_per_node_fe, b_order, b_sizes);
                            end
            
            linsolve_t += @elapsed(tmpki[:,stage] = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node_fe, b_order, b_sizes));
            
            # At this point tmpki[:,stage] holds the boundary values
            # directly write them to the variable values and zero sol.
            copy_bdry_vals_to_variables(fe_var, tmpki[:,stage], grid_data, dofs_per_node_fe, zero_vals=true);
            
        end
        for i=1:length(fe_sol)
            for stage=1:fe_stepper.stages
                fe_sol[i] += fe_stepper.dt * fe_stepper.b[stage] * tmpki[i, stage];
            end
        end
        
        copy_bdry_vals_to_vector(fe_var, fe_sol, grid_data, dofs_per_node_fe);
        place_sol_in_vars(fe_var, fe_sol);
        
        post_step_function();
    end
    
    ### Then step the FV part ######################################################################################
    if !(fv_var === nothing)
        fv_sol = fv_storage[1];
        fv_allocated_vecs = fv_storage[2];
        fv_tmpresult = fv_storage[3];
        fv_tmpki = fv_storage[4];
        
        # solution will be placed in var.values for each stage
        for stage=1:fv_stepper.stages
            stime = t + fv_stepper.c[stage]*fv_stepper.dt;
            
            # Update the values in vars to be used in this stage
            if stage > 1
                initialized_tmpresult = false;
                for j=1:stage
                    if fv_stepper.a[stage, j] > 0
                        if !initialized_tmpresult
                            initialized_tmpresult = true;
                            for k=1:length(fv_sol)
                                fv_tmpresult[k] = fv_sol[k] + fv_stepper.dt * fv_stepper.a[stage, j] * fv_tmpki[k,j];
                            end
                        else
                            for k=1:length(fv_sol)
                                fv_tmpresult[k] += fv_stepper.dt * fv_stepper.a[stage, j] * fv_tmpki[k,j];
                            end
                        end
                    end
                end
                
                if initialized_tmpresult
                    FV_copy_bdry_vals_to_vector(fv_var, fv_tmpresult, fv_grid, dofs_per_node_fv);
                    place_sol_in_vars(fv_var, fv_tmpresult);
                end
                post_step_function(); # seems weird, but imagine this is happening after stage-1
            end
            
            if config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, i); end
            pre_step_function();
            
            tmpki[:,stage] = fv_assemble(fv_var, vol_lhs[2], vol_rhs[2], surf_lhs[2], surf_rhs[2], fv_allocated_vecs, dofs_per_node_fv, dofs_per_loop_fv, stime, fv_stepper.dt, assemble_loops=fv_assemble_func);
        end
        for i=1:length(fv_sol)
            for stage=1:fv_stepper.stages
                fv_sol[i] += fv_stepper.dt * fv_stepper.b[stage] * fv_tmpki[i, stage];
            end
        end
        FV_copy_bdry_vals_to_vector(fv_var, fv_sol, fv_grid, dofs_per_node_fv);
        place_sol_in_vars(fv_var, fv_sol);
        
        post_step_function();
    end
    
    return (assemble_t, linsolve_t);
end