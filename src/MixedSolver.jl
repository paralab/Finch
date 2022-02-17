#=
# Mixed FEM/FVM solver
=#
module MixedSolver

export solve, nonlinear_solve

using LinearAlgebra, SparseArrays

# See finch_import_symbols.jl for a list of all imported symbols.
import ..Finch: @import_finch_symbols
@import_finch_symbols()

include("mixed_time_step.jl");
include("fe_boundary.jl");
include("fv_boundary.jl");
include("cg_matrixfree.jl");
include("level_benchmark.jl");

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

# The input should be in arrays like [fe_version, fv_version]
function linear_solve(var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, stepper=nothing, assemble_func=nothing)
    # if config.linalg_matrixfree
    #     return solve_matrix_free_sym(var, bilinear, linear, stepper, assemble_func=assemble_func);
    #     #return solve_matrix_free_asym(var, bilinear, linear, stepper);
    # end
    
    if stepper === nothing
        printerr("FV assumes a time dependent problem. Set time stepper and initial conditions.", fatal=true)
    end
    
    # Separate FE and FV variables
    # These should both have at least one var
    fe_var = [];
    fv_var = [];
    
    # If more than one variable(should always be the case)
    if typeof(var) <: Array
        # multiple variables being solved for simultaneously
        # separate the FE and FV ones
        dofs_per_node_fe = 0;
        dofs_per_loop_fe = 0;
        dofs_per_node_fv = 0;
        dofs_per_loop_fv = 0;
        for vi=1:length(var)
            if var[vi].discretization == FV
                dofs_per_loop_fv += length(var[vi].symvar);
                dofs_per_node_fv += var[vi].total_components;
                push!(fv_var, var[vi]);
            elseif var[vi].discretization == DG
                printerr("Mixed solver is not ready for DG. Sorry.", fatal=true)
            else
                dofs_per_loop_fe += length(var[vi].symvar);
                dofs_per_node_fe += var[vi].total_components;
                push!(fe_var, var[vi]);
            end
        end
    else
        # one variable
        # In this case it doesn't make any sense to use a mixed solver.
        # Just redirect it to the appropriate one.
        # Note that this does not initialize the solver, so be careful.
        if var.discretization == CG
            return CGSolver.linear_solve(var, vol_lhs, vol_rhs, stepper, assemble_func);
        elseif var.discretization == DG
            return DGSolver.linear_solve(var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, stepper, assemble_func);
        elseif var.discretization == FV
            return FVSolver.linear_solve(var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, stepper, assemble_func);
        end
    end
    
    if typeof(stepper) <: Array
        fe_stepper = stepper[1];
        fv_stepper = stepper[2];
    else
        fe_stepper = stepper;
        fv_stepper = stepper;
    end
    
    if typeof(assemble_func) <: Array
        fe_assemble_func = assemble_func[1];
        fv_assemble_func = assemble_func[2];
    else
        fe_assemble_func = assemble_func;
        fv_assemble_func = assemble_func;
    end
    
    #########################################
    # FE parts
    N1 = size(grid_data.allnodes,2);
    Nn_fe = dofs_per_node_fe * N1;
    Np = refel.Np;
    nel = size(grid_data.loc2glb,2);
    
    # For partitioned meshes, keep these numbers handy
    (b_order, b_sizes) = get_partitioned_ordering(N1, dofs_per_node_fe);
    
    # Allocate arrays that will be used by assemble
    Nsparse = nel*dofs_per_node_fe*Np*dofs_per_node_fe*Np; # reserve space for sparse matrix IJV
    rhsvec = zeros(Nn_fe);
    lhsmatI = zeros(Int, Nsparse);
    lhsmatJ = zeros(Int, Nsparse);
    lhsmatV = zeros(Nsparse);
    allocated_vecs = [rhsvec, lhsmatI, lhsmatJ, lhsmatV];
    
    #########################################
    # FV parts\
    nel_fv = fv_grid.nel_owned;
    nfaces = size(fv_grid.face2element, 2);
    Nn_fv = dofs_per_node_fv * nel_fv;
    Nf_fv = dofs_per_node_fv * nfaces
    
    # Allocate arrays that will be used by assemble
    # These vectors will hold the integrated values(one per cell).
    # They will later be combined.
    sourcevec = zeros(Nn_fv);
    fluxvec = zeros(Nn_fv);
    facefluxvec = zeros(Nf_fv);
    face_done = zeros(nfaces); # Increment when the corresponding flux value is computed.
    fv_allocated_vecs = [sourcevec, fluxvec, facefluxvec, face_done];
    
    #########################################
    
    # First assemble LHS
    assemble_t = @elapsed begin
                    (A, b) = assemble(fe_var, vol_lhs[1], vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, 0, fe_stepper.dt, assemble_loops=fe_assemble_func);
                    (A, b) = gather_system(A, b, N1, dofs_per_node_fe, b_order, b_sizes);
                end
    
    log_entry("Initial assembly took "*string(assemble_t)*" seconds");

    log_entry("Beginning "*string(fe_stepper.Nsteps)*" time steps.");
    t = 0;
    fe_sol = get_var_vals(fe_var);
    if config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, 0); end
    fv_sol = get_var_vals(fv_var);
    
    start_t = Base.Libc.time();
    assemble_t = 0;
    linsolve_t = 0;
    last2update = 0;
    last10update = 0;
    
    # allocate any temporary storage needed by steppers
    if fe_stepper.type == LSRK4
        tmppi = zeros(size(b));
        tmpki = zeros(size(b));
        fe_storage = (A, b, fe_sol, allocated_vecs, b_order, b_sizes, tmppi, tmpki);
        perform_fe_time_step = mixed_lsrk4_step;
        
    elseif fe_stepper.type == RK4
        tmpresult = zeros(size(b));
        tmpki = zeros(length(b), fe_stepper.stages);
        fe_storage = (A, b, fe_sol, allocated_vecs, b_order, b_sizes, tmpresult, tmpki);
        perform_fe_time_step = mixed_rk4_step;
        
    elseif fe_stepper.type == EULER_EXPLICIT
        fe_storage = (A, b, fe_sol, allocated_vecs, b_order, b_sizes);
        perform_fe_time_step = mixed_euler_explicit_step;
        
    elseif fe_stepper.type == EULER_IMPLICIT
        fe_storage = (A, b, fe_sol, allocated_vecs, b_order, b_sizes);
        perform_fe_time_step = mixed_euler_implicit_step;
        
    elseif fe_stepper.type == CRANK_NICHOLSON
        fe_storage = (A, b, fe_sol, allocated_vecs, b_order, b_sizes);
        perform_fe_time_step = mixed_crank_nicholson_step;
        
    elseif fe_stepper.stages > 1
        tmpresult = zeros(size(b));
        tmpki = zeros(length(b), fe_stepper.stages);
        fe_storage = (A, b, fe_sol, allocated_vecs, b_order, b_sizes, tmpresult, tmpki);
        perform_fe_time_step = mixed_multistage_step;
    end
    
    # same for FV
    if fv_stepper.type == LSRK4
        fv_tmppi = zeros(size(fv_sol));
        fv_tmpki = zeros(size(fv_sol));
        fv_storage = (fv_sol, fv_allocated_vecs, fv_tmppi, fv_tmpki);
        perform_fv_time_step = mixed_lsrk4_step;
        
    elseif fv_stepper.type == RK4
        fv_tmpvals = zeros(size(fv_sol));
        fv_tmpki = zeros(length(fv_sol), fv_stepper.stages);
        fv_storage = (fv_sol, fv_allocated_vecs, fv_tmpvals, fv_tmpki);
        perform_fv_time_step = mixed_rk4_step;
        
    elseif fv_stepper.type == EULER_EXPLICIT
        fv_storage = (fv_sol, fv_allocated_vecs);
        perform_fv_time_step = mixed_euler_explicit_step;
        
    elseif fv_stepper.type == EULER_IMPLICIT
        # perform_fv_time_step = mixed_euler_implicit_step;
        printerr("implicit stepper not ready for FV", fatal=true)
        
    elseif fv_stepper.type == CRANK_NICHOLSON
        # perform_fv_time_step = mixed_crank_nicholson_step;
        printerr("implicit stepper not ready for FV", fatal=true)
        
    elseif fv_stepper.stages > 1
        fv_tmpvals = zeros(size(fv_sol));
        fv_tmpki = zeros(length(fv_sol), fv_stepper.stages);
        fv_storage = (fv_sol, fv_allocated_vecs, fv_tmpvals, fv_tmpki);
        perform_fv_time_step = mixed_multistage_step;
    end
    
    # Info to be passed to step functions
    step_info = (fe_stepper, fv_stepper, dofs_per_node_fe, dofs_per_loop_fe, fe_assemble_func, dofs_per_node_fv, dofs_per_loop_fv, fv_assemble_func);
    
    # The time step loop
    print("Time stepping progress(%): 0");
    for i=1:fe_stepper.Nsteps
        
        (at, lt) = perform_fe_time_step(fe_var, nothing, vol_lhs, vol_rhs, surf_lhs, surf_rhs, t, step_info, fe_storage, nothing);
        assemble_t += at;
        linsolve_t += lt;
        
        (at, lt) = perform_fv_time_step(nothing, fv_var, vol_lhs, vol_rhs, surf_lhs, surf_rhs, t, step_info, nothing, fv_storage);
        assemble_t += at;
        linsolve_t += lt;
        
        
        
        # This will be removed eventually, but keep it here for now.
        
        # if stepper.stages > 1
        #     # LSRK4 is a special case, low storage
        #     if stepper.type == LSRK4
        #         # Low storage RK4: 
        #         # p0 = u
        #         #   ki = ai*k(i-1) + dt*f(p(i-1), t+ci*dt)
        #         #   pi = p(i-1) + bi*ki
        #         # u = p5
                
        #         tmppi = get_var_vals(fe_var, tmppi);
        #         fv_tmppi = get_var_vals(fv_var, fv_tmppi);
        #         for rki=1:stepper.stages
        #             rktime = t + stepper.c[rki]*stepper.dt;
        #             # p(i-1) is currently in u
        #             if config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, i); end
        #             pre_step_function();
                    
        #             ### First step the FE part #######################################################################################
        #             assemble_t += @elapsed begin
        #                                 b = assemble(fe_var, nothing, vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, rktime, stepper.dt; rhs_only = true, assemble_loops=assemble_func);
        #                                 b = gather_system(nothing, b, N1, dofs_per_node_fe, b_order, b_sizes);
        #                             end
                    
        #             linsolve_t += @elapsed(fe_sol = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node_fe, b_order, b_sizes));
                    
        #             # At this point sol holds the boundary values
        #             # directly write them to the variable values and zero sol.
        #             copy_bdry_vals_to_variables(fe_var, fe_sol, grid_data, dofs_per_node_fe, zero_vals=true);
                    
        #             if rki == 1 # because a1 == 0
        #                 tmpki = stepper.dt .* fe_sol;
        #             else
        #                 tmpki = stepper.a[rki].*tmpki + stepper.dt.*fe_sol;
        #             end
        #             tmppi = tmppi + stepper.b[rki].*tmpki
                    
        #             copy_bdry_vals_to_vector(fe_var, tmppi, grid_data, dofs_per_node_fe);
        #             place_sol_in_vars(fe_var, tmppi, stepper);
                    
        #             ### Then step the FV part ######################################################################################
        #             fv_sol = fv_assemble(fv_var, vol_lhs[2], vol_rhs[2], surf_lhs[2], surf_rhs[2], fv_allocated_vecs, dofs_per_node_fv, dofs_per_loop_fv, rktime, stepper.dt, assemble_loops=assemble_func);
                        
        #             if rki == 1 # because a1 == 0
        #                 fv_tmpki = stepper.dt .* fv_sol;
        #             else
        #                 fv_tmpki = stepper.a[rki].*fv_tmpki + stepper.dt.*fv_sol;
        #             end
        #             fv_tmppi = fv_tmppi + stepper.b[rki].*fv_tmpki
                    
        #             FV_copy_bdry_vals_to_vector(fv_var, fv_tmppi, fv_grid, dofs_per_node_fv);
        #             place_sol_in_vars(fv_var, fv_tmppi, stepper);
                    
        #             post_step_function();
        #         end
                
        #     else
        #         # Explicit multi-stage methods: 
        #         # x = x + dt*sum(bi*ki)
        #         # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
                
        #         # solution will be placed in var.values for each stage
        #         for stage=1:stepper.stages
        #             stime = t + stepper.c[stage]*stepper.dt;
                    
        #             # Update the values in vars to be used in this stage
        #             if stage > 1
        #                 initialized_tmpresult = false;
        #                 for j=1:stage
        #                     if stepper.a[stage, j] > 0
        #                         if !initialized_tmpresult
        #                             initialized_tmpresult = true;
        #                             for k=1:length(fe_sol)
        #                                 tmpresult[k] = fe_sol[k] + stepper.dt * stepper.a[stage, j] * tmpki[k,j];
        #                             end
        #                             for k=1:length(fv_sol)
        #                                 fv_tmpresult[k] = fv_sol[k] + stepper.dt * stepper.a[stage, j] * fv_tmpki[k,j];
        #                             end
        #                         else
        #                             for k=1:length(fe_sol)
        #                                 tmpresult[k] += stepper.dt * stepper.a[stage, j] * tmpki[k,j];
        #                             end
        #                             for k=1:length(fv_sol)
        #                                 fv_tmpresult[k] += stepper.dt * stepper.a[stage, j] * fv_tmpki[k,j];
        #                             end
        #                         end
        #                     end
        #                 end
                        
        #                 if initialized_tmpresult
        #                     copy_bdry_vals_to_vector(fe_var, tmpresult, grid_data, dofs_per_node_fe);
        #                     place_sol_in_vars(fe_var, tmpresult, stepper);
        #                     FV_copy_bdry_vals_to_vector(fv_var, fv_tmpresult, fv_grid, dofs_per_node_fv);
        #                     place_sol_in_vars(fv_var, fv_tmpresult, stepper);
        #                 end
        #                 post_step_function(); # seems weird, but imagine this is happening after stage-1
        #             end
                    
        #             if config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, i); end
        #             pre_step_function();
                    
        #             ### First step the FE part #######################################################################################
        #             assemble_t += @elapsed begin
        #                                 b = assemble(fe_var, nothing, vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, stime, stepper.dt; rhs_only = true, assemble_loops=assemble_func);
        #                                 b = gather_system(nothing, b, N1, dofs_per_node_fe, b_order, b_sizes);
        #                             end
                    
        #             linsolve_t += @elapsed(tmpki[:,stage] = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node_fe, b_order, b_sizes));
                    
        #             # At this point tmpki[:,stage] holds the boundary values
        #             # directly write them to the variable values and zero sol.
        #             copy_bdry_vals_to_variables(fe_var, tmpki[:,stage], grid_data, dofs_per_node_fe, zero_vals=true);
                    
        #             ### Then step the FV part ######################################################################################
        #             tmpki[:,stage] = fv_assemble(fv_var, vol_lhs[2], vol_rhs[2], surf_lhs[2], surf_rhs[2], fv_allocated_vecs, dofs_per_node_fv, dofs_per_loop_fv, stime, stepper.dt, assemble_loops=assemble_func);
                    
        #         end
        #         for stage=1:stepper.stages
        #             fe_sol += stepper.dt * stepper.b[stage] .* tmpki[:, stage];
        #             fv_sol += stepper.dt * stepper.b[stage] .* fv_tmpki[:, stage];
        #         end
        #         copy_bdry_vals_to_vector(fe_var, fe_sol, grid_data, dofs_per_node_fe);
        #         FV_copy_bdry_vals_to_vector(fv_var, fv_sol, fv_grid, dofs_per_node_fv);
        #         place_sol_in_vars(fe_var, fe_sol, stepper);
        #         place_sol_in_vars(fv_var, fv_sol, stepper);
                
        #         post_step_function();
        #     end
            
        # elseif stepper.type == EULER_EXPLICIT
        #     if config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, i); end
        #     pre_step_function();
            
        #     ### First step the FE part #######################################################################################
        #     assemble_t += @elapsed begin
        #                                 b = assemble(fe_var, nothing, vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, t, stepper.dt; rhs_only = true, assemble_loops=assemble_func);
        #                                 b = gather_system(nothing, b, N1, dofs_per_node_fe, b_order, b_sizes);
        #                             end
            
        #     linsolve_t += @elapsed(tmpvec = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node_fe, b_order, b_sizes));
            
        #     # At this point tmpvec holds the boundary values
        #     # directly write them to the variable values and zero sol.
        #     copy_bdry_vals_to_variables(fe_var, tmpvec, grid_data, dofs_per_node_fe, zero_vals=true);
            
        #     fe_sol = fe_sol .+ stepper.dt .* tmpvec;
            
        #     copy_bdry_vals_to_vector(fe_var, fe_sol, grid_data, dofs_per_node_fe);
        #     place_sol_in_vars(fe_var, fe_sol, stepper);
            
        #     ### Then step the FV part ######################################################################################
        #     fv_sol = fv_sol .+ stepper.dt .*fv_assemble(fv_var, vol_lhs[2], vol_rhs[2], surf_lhs[2], surf_rhs[2], fv_allocated_vecs, dofs_per_node_fv, dofs_per_loop_fv, t, stepper.dt, assemble_loops=assemble_func);
            
        #     FV_copy_bdry_vals_to_vector(fv_var, fv_sol, fv_grid, dofs_per_node_fv);
        #     place_sol_in_vars(fv_var, fv_sol, stepper);
            
        #     post_step_function();
            
        # elseif stepper.type == PECE
        #     # Predictor (explicit Euler)
        #     if config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, i); end
        #     pre_step_function();
            
        #     ### First step the FE part #######################################################################################
        #     assemble_t += @elapsed begin
        #                                 b = assemble(fe_var, nothing, vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, t, stepper.dt; rhs_only = true, assemble_loops=assemble_func);
        #                                 b = gather_system(nothing, b, N1, dofs_per_node_fe, b_order, b_sizes);
        #                             end
            
        #     linsolve_t += @elapsed(fe_tmpsol = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node_fe, b_order, b_sizes));
            
        #     # At this point fe_tmpsol holds the boundary values
        #     # directly write them to the variable values and zero sol.
        #     copy_bdry_vals_to_variables(fe_var, fe_tmpsol, grid_data, dofs_per_node_fe, zero_vals=true);
            
        #     fe_tmpsol = fe_sol .+ stepper.dt .* fe_tmpsol;
            
        #     copy_bdry_vals_to_vector(fe_var, fe_tmpsol, grid_data, dofs_per_node_fe);
        #     place_sol_in_vars(fe_var, fe_tmpsol, stepper);
            
        #     ### Then step the FV part ######################################################################################
        #     fv_tmpsol = fv_sol .+ stepper.dt .*fv_assemble(fv_var, vol_lhs[2], vol_rhs[2], surf_lhs[2], surf_rhs[2], fv_allocated_vecs, dofs_per_node_fv, dofs_per_loop_fv, t, stepper.dt, assemble_loops=assemble_func);
            
        #     FV_copy_bdry_vals_to_vector(fv_var, fv_tmpsol, fv_grid, dofs_per_node_fv);
        #     place_sol_in_vars(fv_var, fv_tmpsol, stepper);
            
        #     post_step_function();
            
        #     # Corrector (implicit Euler)
        #     if config.num_partitions > 1 exchange_ghosts(fv_var, fv_grid, i); end
        #     pre_step_function();
            
        #     ### First step the FE part #######################################################################################
        #     assemble_t += @elapsed begin
        #         b = assemble(fe_var, nothing, vol_rhs[1], allocated_vecs, dofs_per_node_fe, dofs_per_loop_fe, t+stepper.dt, stepper.dt; rhs_only = true, assemble_loops=assemble_func);
        #         b = gather_system(nothing, b, N1, dofs_per_node_fe, b_order, b_sizes);
        #     end

        #     linsolve_t += @elapsed(fe_tmpsol = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node_fe, b_order, b_sizes));

        #     # At this point fe_tmpsol holds the boundary values
        #     # directly write them to the variable values and zero sol.
        #     copy_bdry_vals_to_variables(fe_var, fe_tmpsol, grid_data, dofs_per_node_fe, zero_vals=true);

        #     fe_sol = fe_sol .+ stepper.dt .* fe_tmpsol;

        #     copy_bdry_vals_to_vector(fe_var, fe_sol, grid_data, dofs_per_node_fe);
        #     place_sol_in_vars(fe_var, fe_sol, stepper);

        #     ### Then step the FV part ######################################################################################
        #     fv_sol = fv_sol .+ stepper.dt .*fv_assemble(fv_var, vol_lhs[2], vol_rhs[2], surf_lhs[2], surf_rhs[2], fv_allocated_vecs, dofs_per_node_fv, dofs_per_loop_fv, t+stepper.dt, stepper.dt, assemble_loops=assemble_func);

        #     FV_copy_bdry_vals_to_vector(fv_var, fv_sol, fv_grid, dofs_per_node_fv);
        #     place_sol_in_vars(fv_var, fv_sol, stepper);

        #     post_step_function();
            
        # else
        #     printerr("Unknown explicit stepper: "*string(stepper.type))
        #     return sol;
        # end
        
        t += fe_stepper.dt;
        
        progressPercent = Int(floor(i*100.0/fe_stepper.Nsteps));
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

    log_entry("Stepping took "*string(end_t-start_t)*" seconds. ("*string(assemble_t)*" for assembly, "*string(linsolve_t)*" for linear solve)");
    
    return fe_sol;
end

function nonlinear_solve(var, nlvar, bilinear, linear, stepper=nothing; assemble_loops=nothing)
    printerr("nonlinear solve not ready for mixed solver.", fatal=true)
end

# assembles the A and b in Au=b
function assemble(var, bilinear, linear, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0; assemble_loops=nothing, rhs_only = false)
    # If an assembly loop function was provided, use it
    if !(assemble_loops === nothing)
        return assemble_loops.func(var, bilinear, linear, allocated_vecs, dofs_per_node, dofs_per_loop, t, dt; rhs_only=rhs_only);
    end
    
    # Label things that were allocated externally
    b = allocated_vecs[1];
    if !rhs_only
        AI = allocated_vecs[2];
        AJ = allocated_vecs[3];
        AV = allocated_vecs[4];
    end
    # zero b
    b .= 0;
    
    Np = refel.Np;
    nel = size(grid_data.loc2glb,2);

    # Stiffness and mass are precomputed for uniform grid meshes
    precomputed_mass_stiffness = config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
    if precomputed_mass_stiffness
        wdetj = refel.wg .* geo_factors.detJ[1];
        J = geo_factors.J[1];
        
        if config.dimension == 1
            (RQ1, RD1) = @level_bench Level1 build_deriv_matrix(refel, J);
            TRQ1 = RQ1';
            stiffness = [(TRQ1 * diagm(wdetj) * RQ1)];
        elseif config.dimension == 2
            (RQ1, RQ2, RD1, RD2) = build_deriv_matrix(refel, J);
            (TRQ1, TRQ2) = (RQ1', RQ2');
            stiffness = [(TRQ1 * diagm(wdetj) * RQ1) , (TRQ2 * diagm(wdetj) * RQ2)];
        else
            (RQ1, RQ2, RQ3, RD1, RD2, RD3) = build_deriv_matrix(refel, J);
            (TRQ1, TRQ2, TRQ3) = (RQ1', RQ2', RQ3');
            stiffness = [(TRQ1 * diagm(wdetj) * RQ1) , (TRQ2 * diagm(wdetj) * RQ2) , (TRQ3 * diagm(wdetj) * RQ3)];
        end
        mass = (refel.Q)' * diagm(wdetj) * refel.Q;
    else
        stiffness = 0;
        mass = 0;
    end
    
    loop_time = Base.Libc.time();
    # Elemental loop follows elemental ordering
    for ei=1:nel
        eid = elemental_order[ei];
        loc2glb = grid_data.loc2glb[:,eid]; # global indices of this element's nodes
        
        if !rhs_only
            Astart = (eid-1)*Np*dofs_per_node*Np*dofs_per_node + 1; # The segment of AI, AJ, AV for this element
        end
        
        volargs = (var, eid, 0, grid_data, geo_factors, refel, t, dt, stiffness, mass);
        
        if dofs_per_node == 1
            linchunk = @level_bench Level1 linear.func(volargs);  # get the elemental linear part
            b[loc2glb] += linchunk;
            
            if !rhs_only
                bilinchunk = bilinear.func(volargs); # the elemental bilinear part
                #A[glb, glb] .+= bilinchunk;         # This will be very inefficient for sparse A
                for jj=1:Np
                    offset = Astart - 1 + (jj-1)*Np;
                    for ii=1:Np
                        AI[offset + ii] = loc2glb[ii];
                        AJ[offset + ii] = loc2glb[jj];
                        AV[offset + ii] = bilinchunk[ii, jj];
                    end
                end
            end
            
        else # more than one DOF per node
            linchunk = linear.func(volargs);
            insert_linear!(b, linchunk, loc2glb, 1:dofs_per_node, dofs_per_node);
            
            if !rhs_only
                bilinchunk = bilinear.func(volargs);
                @level_bench Level1  insert_bilinear!(AI, AJ, AV, Astart, bilinchunk, loc2glb, 1:dofs_per_node, dofs_per_node);
            end
        end
    end
    loop_time = Base.Libc.time() - loop_time;
    
    if !rhs_only
        # Build the sparse A. Uses default + to combine overlaps
        A = @level_bench Level1 sparse(AI, AJ, AV);
    end
    
    # Boundary conditions
    bc_time = Base.Libc.time();
    if rhs_only
        b =@level_bench Level1  apply_boundary_conditions_rhs_only(var, b, t);
    else
        (A, b) = @level_bench Level1 apply_boundary_conditions_lhs_rhs(var, A, b, t);
    end
    bc_time = Base.Libc.time() - bc_time;
    
    if rhs_only
        return b;
        
    else
        log_entry("Elemental loop time:     "*string(loop_time));
        log_entry("Boundary condition time: "*string(bc_time));
        return (A, b);
    end
end

function assemble_rhs_only(var, linear, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0; assemble_loops=nothing)
    return assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, t, dt; assemble_loops=assemble_loops, rhs_only=true)
end

function fv_assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, dofs_per_loop=1, t=0, dt=0; assemble_loops=nothing)
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

# ##########################################################

# Inset the single dof into the greater construct
function insert_linear!(b, bel, glb, dof, Ndofs)
    # group nodal dofs
    for d=1:length(dof)
        ind = glb.*Ndofs .- (Ndofs-dof[d]);
        ind2 = ((d-1)*length(glb)+1):(d*length(glb));

        b[ind] = b[ind] + bel[ind2];
    end
end

function insert_bilinear!(AI, AJ, AV, Astart, ael, glb, dof, Ndofs)
    Np = length(glb);
    # group nodal dofs
    for dj=1:length(dof)
        indj = glb.*Ndofs .- (Ndofs-dof[dj]);
        indj2 = ((dj-1)*Np+1):(dj*Np);
        for di=1:length(dof)
            indi = glb.*Ndofs .- (Ndofs-dof[di]);
            indi2 = ((di-1)*Np+1):(di*Np);
            
            #a[indi, indj] = a[indi, indj] + ael[indi2, indj2];
            
            for jj=1:Np
                offset = Astart + (jj-1 + Np*(dj-1))*Np*Ndofs + Np*(di-1) - 1;
                for ii=1:Np
                    AI[offset + ii] = indi[ii];
                    AJ[offset + ii] = indj[jj];
                    AV[offset + ii] = ael[indi2[ii], indj2[jj]];
                end
            end
        end
    end
    
end

function place_sol_in_vars(var, sol, stepper=nothing)
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
                var[vi].values[compi,:] = sol[(compi+tmp):totalcomponents:end];
                tmp = tmp + 1;
            end
        end
    else
        components = var.total_components;
        for compi=1:components
            var.values[compi,:] = sol[compi:components:end];
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

######################################################################################
# For partitioned meshes
function get_partitioned_ordering(nnodes, dofs_per_node)
    if config.num_procs > 1
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
        return (nothing, nothing);
    end
end

# For multiple processes, gather the system, distribute the solution
# Note: rescatter_b only applies when rhs_only
function gather_system(A, b, nnodes, dofs_per_node, b_order, b_sizes; rescatter_b=false)
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
                    return distribute_solution(new_b, nnodes, dofs_per_node, b_order, b_sizes)
                else
                    return new_b;
                end
                
            else # other procs just send their data
                MPI.Gatherv!(b, nothing, 0, MPI.COMM_WORLD);
                
                if rescatter_b
                    return distribute_solution(nothing, nnodes, dofs_per_node, b_order, b_sizes)
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

function distribute_solution(sol, nnodes, dofs_per_node, b_order, b_sizes)
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

################################################################
## Options for solving the assembled linear system
################################################################
function linear_system_solve(A,b)
    if config.num_procs > 1 && config.proc_rank > 0
        # Other procs don't have the full A, just return b
        return b;
    end
    
    if config.linalg_backend == DEFAULT_SOLVER
        return A\b;
    elseif config.linalg_backend == PETSC_SOLVER
        # For now try this simple setup.
        # There will be many options that need to be available to the user. TODO
        # The matrix should really be constructed from the begninning as a PETSc one. TODO
        
        petsclib = PETSc.petsclibs[1]           #
        inttype = PETSc.inttype(petsclib);      # These should maybe be kept in config
        scalartype = PETSc.scalartype(petsclib);#
        
        (I, J, V) = findnz(A);
        (n,n) = size(A);
        nnz = zeros(inttype, n); # number of non-zeros per row
        for i=1:length(I)
            nnz[I[i]] += 1;
        end
        
        petscA = PETSc.MatSeqAIJ{scalartype}(n,n,nnz);
        for i=1:length(I)
            petscA[I[i],J[i]] = V[i];
        end
        PETSc.assemble(petscA);
        
        ksp = PETSc.KSP(petscA; ksp_rtol=1e-8, pc_type="jacobi", ksp_monitor=false); # Options should be available to user
        
        return ksp\b;
    elseif config.linalg_backend == CUDA_SOLVER
        elty = typeof(b[1])
        x = zeros(elty ,length(b))
        tol = convert(real(elty),1e-6)
        x = CUDA.CUSOLVER.csrlsvlu!(A, b, x, tol ,one(Cint),'O')
        return x
    end
end

end #module
