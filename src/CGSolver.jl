#=
# CG solver
=#
module CGSolver

export solve, nonlinear_solve

using LinearAlgebra, SparseArrays, CUDA

# See finch_import_symbols.jl for a list of all imported symbols.
import ..Finch: @import_finch_symbols
@import_finch_symbols()


include("fe_boundary.jl");
include("nonlinear.jl")
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

function linear_solve(var, func, stepper=nothing)
    args = (var, grid_data, refel, geo_factors, config, coefficients, variables, test_functions, indexers, prob);
    return TimerOutputs.@timeit timer_output "cg_colve" func.func(args);
end



function old_linear_solve(var, bilinear, linear, stepper=nothing; assemble_func=nothing)
    if config.linalg_matrixfree
        return solve_matrix_free_sym(var, bilinear, linear, stepper, assemble_func=assemble_func);
        #return solve_matrix_free_asym(var, bilinear, linear, stepper);
    end
    
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
    N1 = size(grid_data.allnodes,2);
    Nn = dofs_per_node * N1;
    Np = refel.Np;
    nel = size(grid_data.loc2glb,2);
    
    # For partitioned meshes, keep these numbers handy
    (b_order, b_sizes) = get_partitioned_ordering(N1, dofs_per_node);
    
    # Allocate arrays that will be used by assemble
    rhsvec = zeros(Nn);
    lhsmatI = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
    lhsmatJ = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
    lhsmatV = zeros(nel*dofs_per_node*Np*dofs_per_node*Np);
    allocated_vecs = [rhsvec, lhsmatI, lhsmatJ, lhsmatV];
    
    if prob.time_dependent && !(stepper === nothing)
        # First assemble both lhs and rhs
        assemble_t = @elapsed begin
                        (A, b) = assemble(var, bilinear, linear, allocated_vecs, dofs_per_node, dofs_per_loop, 0, stepper.dt, assemble_loops=assemble_func);
                        (A, b) = gather_system(A, b, N1, dofs_per_node, b_order, b_sizes);
                    end
        
        log_entry("Initial assembly took "*string(assemble_t)*" seconds");

        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        sol = get_var_vals(var);
        
        start_t = Base.Libc.time();
        assemble_t = 0;
        linsolve_t = 0;
        last2update = 0;
        last10update = 0;
        
        # allocate any temporary storage needed by steppers
        if stepper.type == LSRK4
            tmppi = zeros(size(b));
            tmpki = zeros(size(b));
        elseif stepper.stages > 1
            last_result = zeros(size(b));
            tmpresult = zeros(size(b));
            tmpki = zeros(length(b), stepper.stages);
        end
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
                    for rki=1:stepper.stages
                        rktime = t + stepper.c[rki]*stepper.dt;
                        # p(i-1) is currently in u
                        pre_step_function();
                        
                        assemble_t += @elapsed begin
                                            b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, rktime, stepper.dt; rhs_only = true, assemble_loops=assemble_func);
                                            b = gather_system(nothing, b, N1, dofs_per_node, b_order, b_sizes);
                                        end
                        
                        linsolve_t += @elapsed(sol = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node, b_order, b_sizes));
                        
                        # At this point sol holds the boundary values
                        # directly write them to the variable values and zero sol.
                        copy_bdry_vals_to_variables(var, sol, grid_data, dofs_per_node, zero_vals=true);
                        
                        if rki == 1 # because a1 == 0
                            tmpki = stepper.dt .* sol;
                        else
                            tmpki = stepper.a[rki].*tmpki + stepper.dt.*sol;
                        end
                        tmppi = tmppi + stepper.b[rki].*tmpki
                        
                        copy_bdry_vals_to_vector(var, tmppi, grid_data, dofs_per_node);
                        place_sol_in_vars(var, tmppi, stepper);
                        
                        post_step_function();
                    end
                    
                else
                    # Explicit multi-stage methods: 
                    # x = x + dt*sum(bi*ki)
                    # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
                    
                    # will hold the final result
                    last_result = get_var_vals(var, last_result);
                    # will be placed in var.values for each stage
                    # Storage for each stage
                    # tmpki = zeros(length(b), stepper.stages);
                    for stage=1:stepper.stages
                        stime = t + stepper.c[stage]*stepper.dt;
                        
                        # Update the values in vars to be used in this stage
                        if stage > 1
                            initialized_tmpresult = false;
                            for j=1:stage
                                if stepper.a[stage, j] > 0
                                    if !initialized_tmpresult
                                        initialized_tmpresult = true;
                                        for k=1:length(last_result)
                                            tmpresult[k] = last_result[k] + stepper.dt * stepper.a[stage, j] * tmpki[k,j];
                                        end
                                    else
                                        for k=1:length(last_result)
                                            tmpresult[k] += stepper.dt * stepper.a[stage, j] * tmpki[k,j];
                                        end
                                    end
                                end
                            end
                            
                            if initialized_tmpresult
                                copy_bdry_vals_to_vector(var, tmpresult, grid_data, dofs_per_node);
                                place_sol_in_vars(var, tmpresult, stepper);
                            end
                            post_step_function(); # seems weird, but imagine this is happening after stage-1
                        end
                        
                        pre_step_function();
                        
                        assemble_t += @elapsed begin
                                            b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, stime, stepper.dt; rhs_only = true, assemble_loops=assemble_func)
                                            b = gather_system(nothing, b, N1, dofs_per_node, b_order, b_sizes);
                                        end
                        
                        linsolve_t += @elapsed(tmpki[:,stage] = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node, b_order, b_sizes));
                        
                        # At this point tmpki[:,stage] holds the boundary values
                        # directly write them to the variable values and zero sol.
                        copy_bdry_vals_to_variables(var, tmpki[:,stage], grid_data, dofs_per_node, zero_vals=true);
                        
                    end
                    for stage=1:stepper.stages
                        last_result += stepper.dt * stepper.b[stage] .* tmpki[:, stage];
                    end
                    copy_bdry_vals_to_vector(var, last_result, grid_data, dofs_per_node);
                    place_sol_in_vars(var, last_result, stepper);
                    post_step_function();
                end
                
            elseif stepper.type == EULER_EXPLICIT
                pre_step_function();
                assemble_t += @elapsed begin
                                            b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, t, stepper.dt; rhs_only = true, assemble_loops=assemble_func)
                                            b = gather_system(nothing, b, N1, dofs_per_node, b_order, b_sizes);
                                        end
                
                linsolve_t += @elapsed(tmpvec = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node, b_order, b_sizes));
                
                # At this point tmpvec holds the boundary values
                # directly write them to the variable values and zero sol.
                copy_bdry_vals_to_variables(var, tmpvec, grid_data, dofs_per_node, zero_vals=true);
                
                sol = sol .+ stepper.dt .* tmpvec;
                
                copy_bdry_vals_to_vector(var, sol, grid_data, dofs_per_node);
                place_sol_in_vars(var, sol, stepper);
                post_step_function();
                
            elseif stepper.type == PECE
                # Predictor (explicit Euler)
                pre_step_function();
                assemble_t += @elapsed begin
                                            b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, t, stepper.dt; rhs_only = true, assemble_loops=assemble_func)
                                            b = gather_system(nothing, b, N1, dofs_per_node, b_order, b_sizes);
                                        end
                
                linsolve_t += @elapsed(tmpvec = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node, b_order, b_sizes));
                
                # At this point tmpvec holds the boundary values
                # directly write them to the variable values and zero sol.
                copy_bdry_vals_to_variables(var, tmpvec, grid_data, dofs_per_node, zero_vals=true);
                
                tmpvec = sol .+ stepper.dt .* tmpvec;
                
                copy_bdry_vals_to_vector(var, tmpvec, grid_data, dofs_per_node);
                place_sol_in_vars(var, tmpvec, stepper);
                post_step_function();
                
                # Corrector (implicit Euler)
                pre_step_function();
                assemble_t += @elapsed begin
                                            b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, t+stepper.dt, stepper.dt; rhs_only = true, assemble_loops=assemble_func)
                                            b = gather_system(nothing, b, N1, dofs_per_node, b_order, b_sizes);
                                        end
                
                linsolve_t += @elapsed(tmpvec = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node, b_order, b_sizes));
                
                # At this point tmpvec holds the boundary values
                # directly write them to the variable values and zero sol.
                copy_bdry_vals_to_variables(var, tmpvec, grid_data, dofs_per_node, zero_vals=true);
                
                sol = sol .+ stepper.dt .* tmpvec;
                
                copy_bdry_vals_to_vector(var, sol, grid_data, dofs_per_node);
                place_sol_in_vars(var, sol, stepper);
                post_step_function();
                
            else # implicit methods include the dt part 
                pre_step_function();
                assemble_t += @elapsed begin
                                            b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, t, stepper.dt; rhs_only = true, assemble_loops=assemble_func)
                                            b = gather_system(nothing, b, N1, dofs_per_node, b_order, b_sizes);
                                        end
                
                linsolve_t += @elapsed(sol = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node, b_order, b_sizes));
                
                place_sol_in_vars(var, sol, stepper);
                post_step_function();
            end
            
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

        log_entry("Stepping took "*string(end_t-start_t)*" seconds. ("*string(assemble_t)*" for assembly, "*string(linsolve_t)*" for linear solve)");
        #display(sol);
		# outfile = "linear_sol.txt"
		# open(outfile, "w") do f
  		# 	for ii in sol
    	# 		println(f, ii)
  		# 	end
		# end # the file f is automatically closed after this block finishes
        return sol;

    else
        assemble_t = @elapsed begin
                        (A, b) = assemble(var, bilinear, linear, allocated_vecs, dofs_per_node, dofs_per_loop, 0, 0, assemble_loops=assemble_func);
                        (A, b) = gather_system(A, b, N1, dofs_per_node, b_order, b_sizes);
                    end
        # uncomment to look at A
        # global Amat = A;
        
        sol_t = @elapsed(sol = distribute_solution( linear_system_solve(A,b) , N1, dofs_per_node, b_order, b_sizes));
        
        log_entry("Assembly took "*string(assemble_t)*" seconds");
        log_entry("Linear solve took "*string(sol_t)*" seconds");
        # display(A);
        # display(b);
        # display(sol);
        return sol;
    end
end

function nonlinear_solve(var, nlvar, bilinear, linear, stepper=nothing; assemble_loops=nothing)
    if config.num_pertitions > 1
        printerr("nonlinear solver is not ready for partitioned meshes. sorry.", fatal=true);
    end
    if prob.time_dependent && !(stepper === nothing)
        #TODO time dependent coefficients
        
        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        start_t = Base.Libc.time();
        nl = nonlinear(100, 1e-9, 1e-9);
        init_nonlinear(nl, var, nlvar, bilinear, linear);
        for i=1:stepper.Nsteps
			println("solve for time step ", i);
			newton(nl, assemble, assemble_rhs_only, nlvar, t, stepper.dt);
            t += stepper.dt;
			nlvar[4].values = copy(nlvar[1].values);
			nlvar[5].values = copy(nlvar[2].values);
			if (i % 10 ==0)
				outfile = string("LD_u_",i);
				open(outfile, "w") do f
					for ii=1:length(nlvar[1].values)
						println(f, nlvar[1].values[ii])
					end
					close(f)
				end
				outfile = string("LD_v_",i);
				open(outfile, "w") do f
					for ii=1:length(nlvar[2].values)
						println(f, nlvar[2].values[ii])
					end
					close(f)
				end
			end
        end
        end_t = Base.Libc.time();

        log_entry("Solve took "*string(end_t-start_t)*" seconds");
        #display(sol);
        #return nlvar.values;
        return [];

    else
        start_t = Base.Libc.time();
        nl = nonlinear(100, 1e-12, 1e-12);
        init_nonlinear(nl, var, nlvar, bilinear, linear);
        newton(nl, assemble, assemble_rhs_only, nlvar);
        end_t = Base.Libc.time();

        log_entry("Solve took "*string(end_t-start_t)*" seconds");
        
        #return nlvar.values;
        return [];
    end
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

# ######################################################
# # To be romoved. Fix nonlinear solve first.
# # assembles the A and b in Au=b
# function assemble_rhs_only(var, linear, t=0.0, dt=0.0; keep_geometric_factors = false, saved_geometric_factors = nothing)
#     Np = refel.Np;
#     nel = size(grid_data.loc2glb,2);
#     N1 = size(grid_data.allnodes,2);
#     multivar = typeof(var) <: Array;
#     maxvarindex = 0;
#     if multivar
#         # multiple variables being solved for simultaneously
#         dofs_per_node = 0;
#         var_to_dofs = [];
#         for vi=1:length(var)
#             tmp = dofs_per_node;
#             dofs_per_node += length(var[vi].symvar);
#             push!(var_to_dofs, (tmp+1):dofs_per_node);
#             maxvarindex = max(maxvarindex,var[vi].index);
#         end
#     else
#         # one variable
#         dofs_per_node = length(var.symvar);
#         maxvarindex = var.index;
#     end
#     Nn = dofs_per_node * N1;

#     b = zeros(Nn);
    
#     # Save geometric factors if desired
#     if keep_geometric_factors
#         saved_wdetj = [];
#         saved_J = [];
#     end
    
#     # Stiffness and mass are precomputed for uniform grid meshes
#     precomputed_mass_stiffness = config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
#     if precomputed_mass_stiffness
#         glb = grid_data.loc2glb[:,1];
#         xe = grid_data.allnodes[:,glb[:]];
#         (detJ, J) = geometric_factors(refel, xe);
#         wdetj = refel.wg .* detJ;
#         if keep_geometric_factors
#             saved_wdetj = wdetj;
#             saved_J = J;
#         end
#         if config.dimension == 1
#             (RQ1, RD1) = build_deriv_matrix(refel, J);
#             TRQ1 = RQ1';
#             stiffness = [(TRQ1 * diagm(wdetj) * RQ1)];
#         elseif config.dimension == 2
#             (RQ1, RQ2, RD1, RD2) = build_deriv_matrix(refel, J);
#             (TRQ1, TRQ2) = (RQ1', RQ2');
#             stiffness = [(TRQ1 * diagm(wdetj) * RQ1) , (TRQ2 * diagm(wdetj) * RQ2)];
#         else
#             (RQ1, RQ2, RQ3, RD1, RD2, RD3) = build_deriv_matrix(refel, J);
#             (TRQ1, TRQ2, TRQ3) = (RQ1', RQ2', RQ3');
#             stiffness = [(TRQ1 * diagm(wdetj) * RQ1) , (TRQ2 * diagm(wdetj) * RQ2) , (TRQ3 * diagm(wdetj) * RQ3)];
#         end
#         mass = (refel.Q)' * diagm(wdetj) * refel.Q;
#     else
#         stiffness = 0;
#         mass = 0;
#     end
    
#     # Elemental loop follows elemental ordering
#     for e=elemental_order;
#         glb = grid_data.loc2glb[:,e];       # global indices of this element's nodes for extracting values from var arrays
#         xe = grid_data.allnodes[:,glb[:]];  # coordinates of this element's nodes for evaluating coefficient functions
        
#         if !precomputed_mass_stiffness
#             if saved_geometric_factors === nothing
#                 (detJ, J) = geometric_factors(refel, xe);
#                 wdetj = refel.wg .* detJ;
                
#                 if keep_geometric_factors
#                     push!(saved_wdetj, wdetj);
#                     push!(saved_J, J);
#                 end
#             else
#                 wdetj = saved_geometric_factors[1][e];
#                 J = saved_geometric_factors[2][e];
#             end
#         end
        
#         rhsargs = (var, xe, glb, refel, wdetj, J, RHS, t, dt, stiffness, mass);

#         #linchunk = linear.func(args);  # get the elemental linear part
#         if dofs_per_node == 1
#             linchunk = linear.func(rhsargs);  # get the elemental linear part
#             b[glb] .+= linchunk;

#         elseif typeof(var) == Variable
#             # only one variable, but more than one dof
#             linchunk = linear.func(rhsargs);
#             insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

#         else
#             linchunk = linear.func(rhsargs);
#             insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

#         end
#     end
    
#     # Boundary conditions
#     bidcount = length(grid_data.bids); # the number of BIDs
#     if dofs_per_node > 1
#         if multivar
#             dofind = 0;
#             for vi=1:length(var)
#                 for compo=1:length(var[vi].symvar)
#                     dofind = dofind + 1;
#                     for bid=1:bidcount
#                         if prob.bc_type[var[vi].index, bid] == NO_BC
#                             # do nothing
#                         else
#                             b = dirichlet_bc_rhs_only(b, prob.bc_func[var[vi].index, bid][compo], grid_data.bdry[bid], t, dofind, dofs_per_node);
#                         end
#                     end
#                 end
#             end
#         else
#             for d=1:dofs_per_node
#                 #rows = ((d-1)*length(glb)+1):(d*length(glb));
#                 dofind = d;
#                 for bid=1:bidcount
#                     if prob.bc_type[var.index, bid] == NO_BC
#                         # do nothing
#                     else
#                         b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, bid][d], grid_data.bdry[bid], t, d, dofs_per_node);
#                     end
#                 end
#             end
#         end
#     else
#         for bid=1:bidcount
#             if prob.bc_type[var.index, bid] == NO_BC
#                 # do nothing
#             else
#                 b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, bid][1], grid_data.bdry[bid], t);
#             end
#         end
#     end
    
#     # Reference points
#     if size(prob.ref_point,1) >= maxvarindex
#         if multivar
#             posind = zeros(Int,0);
#             vals = zeros(0);
#             for vi=1:length(var)
#                 if prob.ref_point[var[vi].index,1]
#                     eii = prob.ref_point[var[vi].index, 2];
#                     tmp = (grid_data.glbvertex[eii[1], eii[2]] - 1)*dofs_per_node + var_to_dofs[vi][1];
#                     if length(prob.ref_point[var[vi].index, 3]) > 1
#                         tmp = tmp:(tmp+length(prob.ref_point[var[vi].index, 3])-1);
#                     end
#                     posind = [posind; tmp];
#                     vals = [vals; prob.ref_point[var[vi].index, 3]];
#                 end
#             end
#             if length(vals) > 0
#                 b[posind] = vals;
#             end
            
#         else
#             if prob.ref_point[var.index,1]
#                 eii = prob.ref_point[var.index, 2];
#                 posind = (grid_data.glbvertex[eii[1], eii[2]] - 1)*dofs_per_node + 1;
#                 if length(prob.ref_point[var.index, 3]) > 1
#                     posind = posind:(posind+length(prob.ref_point[var[vi].index, 3])-1);
#                 end
#                 b[posind] = prob.ref_point[var.index, 3];
#             end
#         end
#     end

#     return b;
# end
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
                var[vi].values[compi,:] .= sol[(compi+tmp):totalcomponents:end];
            end
            tmp = tmp + components;
        end
    else
        components = var.total_components;
        for compi=1:components
            var.values[compi,:] .= sol[compi:components:end];
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
