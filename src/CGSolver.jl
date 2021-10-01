#=
# CG solver
=#
module CGSolver

export solve, nonlinear_solve

import ..Finch: JULIA, CPP, MATLAB, DENDRO, HOMG, CUSTOM_GEN_TARGET,
        SQUARE, IRREGULAR, UNIFORM_GRID, TREE, UNSTRUCTURED, 
        CG, DG, HDG, FV,
        NODAL, MODAL, CELL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, 
        NONLINEAR_NEWTON, NONLINEAR_SOMETHING, 
        EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICHOLSON, RK4, LSRK4, ABM4, 
        DEFAULT_SOLVER, PETSC, 
        VTK, RAW_OUTPUT, CUSTOM_OUTPUT, 
        DIRICHLET, NEUMANN, ROBIN, NO_BC, FLUX,
        MSH_V2, MSH_V4,
        SCALAR, VECTOR, TENSOR, SYM_TENSOR, VAR_ARRAY,
        LHS, RHS,
        LINEMESH, QUADMESH, HEXMESH
        
import ..Finch: log_entry, printerr
import ..Finch: config, prob, variables, mesh_data, grid_data, refel, time_stepper, elemental_order, genfunctions
import ..Finch: Variable, Coefficient, GenFunction
import ..Finch: GeometricFactors, geo_factors, geometric_factors, build_deriv_matrix

using LinearAlgebra, SparseArrays

include("fe_boundary.jl");
include("nonlinear.jl")
include("cg_matrixfree.jl");
include("level_benchmark.jl");

function linear_solve(var, bilinear, linear, stepper=nothing; assemble_func=nothing)
    if config.linalg_matrixfree
        return solve_matrix_free_sym(var, bilinear, linear, stepper);
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
    nel = mesh_data.nel;
    
    # Allocate arrays that will be used by assemble
    # These vectors will hold the integrated values(one per cell).
    # They will later be combined.
    rhsvec = zeros(Nn);
    lhsmatI = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
    lhsmatJ = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
    lhsmatV = zeros(nel*dofs_per_node*Np*dofs_per_node*Np);
    allocated_vecs = [rhsvec, lhsmatI, lhsmatJ, lhsmatV];
    
    if prob.time_dependent && !(stepper === nothing)
        # First assemble both lhs and rhs
        assemble_t = @elapsed((A, b) = assemble(var, bilinear, linear, allocated_vecs, dofs_per_node, dofs_per_loop, 0, stepper.dt, assemble_loops=assemble_func));
        
        log_entry("Initial assembly took "*string(assemble_t)*" seconds");

        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        sol = [];
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
                        assemble_t += @elapsed(b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, rktime, stepper.dt; rhs_only = true, assemble_loops=assemble_func));
                        
                        linsolve_t += @elapsed(sol = A\b);
                        if rki == 1 # because a1 == 0
                            tmpki = stepper.dt .* sol;
                        else
                            tmpki = stepper.a[rki].*tmpki + stepper.dt.*sol;
                        end
                        tmppi = tmppi + stepper.b[rki].*tmpki
                        
                        place_sol_in_vars(var, tmppi, stepper);
                    end
                    
                else
                    # Explicit multi-stage methods: 
                    # x = x + dt*sum(bi*ki)
                    # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
                    
                    # will hold the final result
                    last_result = get_var_vals(var, last_result);
                    # will be placed in var.values for each stage
                    # tmpvals = tmpresult;
                    # Storage for each stage
                    # tmpki = zeros(length(b), stepper.stages);
                    for stage=1:stepper.stages
                        stime = t + stepper.c[stage]*stepper.dt;
                        
                        assemble_t += @elapsed(b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, stime, stepper.dt; rhs_only = true, assemble_loops=assemble_func));
                        
                        linsolve_t += @elapsed(tmpki[:,stage] = A\b);
                        
                        tmpresult = get_var_vals(var, tmpresult);
                        for j=1:(stage-1)
                            if stepper.a[stage, j] > 0
                                tmpresult += stepper.dt * stepper.a[stage, j] .* tmpki[:,j];
                            end
                        end
                        if stage < stepper.stages
                            place_sol_in_vars(var, tmpresult, stepper);
                        end
                    end
                    for stage=1:stepper.stages
                        last_result += stepper.dt * stepper.b[stage] .* tmpki[:, stage];
                    end
                    place_sol_in_vars(var, last_result, stepper);
                end
                
            else # single stage methods such as Euler
                assemble_t += @elapsed(b = assemble(var, nothing, linear, allocated_vecs, dofs_per_node, dofs_per_loop, t, stepper.dt; rhs_only = true, assemble_loops=assemble_func));
                
                linsolve_t += @elapsed(sol = A\b);
                
                place_sol_in_vars(var, sol, stepper);
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
        assemble_t = @elapsed((A, b) = assemble(var, bilinear, linear, allocated_vecs, dofs_per_node, dofs_per_loop, assemble_loops=assemble_func));
        # uncomment to look at A
        # global Amat = A;
        
        sol_t = @elapsed(sol = A\b);
        
        log_entry("Assembly took "*string(assemble_t)*" seconds");
        log_entry("Linear solve took "*string(sol_t)*" seconds");
        # display(A);
        # display(b);
        # display(sol);
        return sol;
    end
end

function nonlinear_solve(var, nlvar, bilinear, linear, stepper=nothing; assemble_loops=nothing)
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
    nel = mesh_data.nel;

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
#     nel = mesh_data.nel;
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
            vect = zeros(totalcomponents * length(var.values[1,:]));
        end
        for compi=1:components
            vect[compi:components:end] = var.values[compi,:];
        end
    end
    
    return vect;
end

end #module
