#=
# DG solver
=#
module DGSolver

export solve, nonlinear_solve

# See finch_import_symbols.jl for a list of all imported symbols.
import ..Finch: @import_finch_symbols
@import_finch_symbols()

using LinearAlgebra, SparseArrays

include("fe_boundary.jl"); # Can we use the CG versions here? 
include("nonlinear.jl");
include("cg_matrixfree.jl");

# Things to do before and after each step(or stage) ########
function default_pre_step() end
function default_post_step() end

pre_step_function = default_pre_step;
post_step_function = default_post_step;

function set_pre_step(fun) global pre_step_function = fun; end
function set_post_step(fun) global post_step_function = fun; end
############################################################

function linear_solve(var, bilinear, linear, face_bilinear, face_linear, stepper=nothing)
    if config.dimension > 2
        printerr("DG solver only available for 1D and 2D. Surface quadrature under construction.")
        #return;
    end
    
    if config.linalg_matrixfree
        return solve_matrix_free_sym(var, bilinear, linear, stepper);
        #return solve_matrix_free_asym(var, bilinear, linear, stepper);
    end
    if prob.time_dependent && !(stepper === nothing)
        assemble_t = @elapsed((A, b) = assemble(var, bilinear, linear, face_bilinear, face_linear, 0, stepper.dt));
        log_entry("Initial assembly took "*string(assemble_t)*" seconds");
        # global Amat = A;
        # global bvec = b;
        # display(A);
        # display(b);

        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        sol = sol = get_var_vals(var);
        
        start_t = Base.Libc.time();
        assemble_t = 0;
        linsolve_t = 0;
        last2update = 0;
        last10update = 0;
        print("Time stepping progress(%): 0");
        for i=1:stepper.Nsteps
            stupidjulia=sol;
            if stepper.stages > 1
                # LSRK4 is a special case, low storage
                if stepper.type == LSRK4
                    # Low storage RK4: 
                    # p0 = u
                    #   ki = ai*k(i-1) + dt*f(p(i-1), t+ci*dt)
                    #   pi = p(i-1) + bi*ki
                    # u = p5
                    
                    # store the original values, x, in a temp vector, so they aren't lost.
                    tmppi = get_var_vals(var);
                    tmpki = zeros(size(b));
                    for rki=1:stepper.stages
                        rktime = t + stepper.c[rki]*stepper.dt;
                        # p(i-1) is currently in u
                        pre_step_function();
                        
                        assemble_t += @elapsed(b = assemble(var, nothing, linear, nothing, face_linear, rktime, stepper.dt; rhs_only = true));
                        
                        linsolve_t += @elapsed(sol = A\b);
                        if rki == 1 # because a1 == 0
                            tmpki = stepper.dt .* sol;
                        else
                            tmpki = stepper.a[rki].*tmpki + stepper.dt .* sol;
                        end
                        tmppi = tmppi + stepper.b[rki].*tmpki
                        
                        # Need to apply boundary conditions now to tmppi
                        tmppi = apply_boundary_conditions_rhs_only(var, tmppi, rktime);
                        
                        place_sol_in_vars(var, tmppi, stepper);
                        post_step_function();
                    end
                    
                else
                    # Explicit multi-stage methods: 
                    # x = x + dt*sum(bi*ki)
                    # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
                    
                    # will hold the final result
                    result = get_var_vals(var);
                    # will be placed in var.values for each stage
                    tmpvals = result;
                    # Storage for each stage
                    tmpki = zeros(length(b), stepper.stages);
                    for stage=1:stepper.stages
                        stime = t + stepper.c[stage]*stepper.dt;
                        pre_step_function();
                        
                        assemble_t += @elapsed(b = assemble(var, nothing, linear, nothing, face_linear, stime, stepper.dt; rhs_only = true));
                        
                        linsolve_t += @elapsed(tmpki[:,stage] = A\b);
                        
                        tmpvals = result;
                        for j=1:(stage-1)
                            if stepper.a[stage, j] > 0
                                tmpvals += stepper.dt * stepper.a[stage, j] .* tmpki[:,j];
                            end
                        end
                        
                        # Need to apply boundary conditions now
                        tmpvals = apply_boundary_conditions_rhs_only(var, tmpvals, stime);
                        
                        place_sol_in_vars(var, tmpvals, stepper);
                        post_step_function();
                    end
                    for stage=1:stepper.stages
                        result += stepper.dt * stepper.b[stage] .* tmpki[:, stage];
                    end
                    
                    result = apply_boundary_conditions_rhs_only(var, result, t + stepper.dt);
                    place_sol_in_vars(var, result, stepper);
                    post_step_function();
                end
            else
                ##### TODO ###### Update needed
                pre_step_function();
                assemble_t += @elapsed(b = assemble(var, nothing, linear, nothing, face_linear, t, stepper.dt; rhs_only = true));
                linsolve_t += @elapsed(sol = A\b);
                
                place_sol_in_vars(var, sol, stepper);
                post_step_function();
            end
            
            #println("b(bdry): "*string(b[1])*", "*string(b[end]));
            
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
        assemble_t = @elapsed((A, b) = assemble(var, bilinear, linear, face_bilinear, face_linear));
        # uncomment to look at A
        # global Amat = A;
        
        sol_t = @elapsed(sol = A\b);
        
        log_entry("Assembly took "*string(assemble_t)*" seconds");
        log_entry("Linear solve took "*string(sol_t)*" seconds");
        #display(A);
        #display(b);
        #display(sol);
        return sol;
    end
end

function nonlinear_solve(var, nlvar, bilinear, linear, stepper=nothing)
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

# ########################
# # tmp for testing
# function the_desired(args)
    
# end
# #########################

# assembles the A and b in Au=b
function assemble(var, bilinear, linear, face_bilinear, face_linear, t=0.0, dt=0.0; rhs_only = false)
    Np = refel.Np;
    nel = mesh_data.nel;
    N1 = size(grid_data.allnodes,2);
    multivar = typeof(var) <: Array;
    maxvarindex = 0;
    if multivar
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        var_to_dofs = [];
        for vi=1:length(var)
            tmp = dofs_per_node;
            dofs_per_node += length(var[vi].symvar);
            push!(var_to_dofs, (tmp+1):dofs_per_node);
            maxvarindex = max(maxvarindex,var[vi].index);
        end
    else
        # one variable
        dofs_per_node = length(var.symvar);
        maxvarindex = var.index;
    end
    Nn = dofs_per_node * N1;

    b = zeros(Nn);
    
    if !rhs_only
        AI = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
        AJ = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
        AV = zeros(nel*dofs_per_node*Np*dofs_per_node*Np);
    end
    
    # Stiffness and mass are precomputed for uniform grid meshes
    precomputed_mass_stiffness = config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
    if precomputed_mass_stiffness
        wdetj = refel.wg .* geo_factors.detJ[1];
        J = geo_factors.J[1];
        if config.dimension == 1
            (RQ1, RD1) = build_deriv_matrix(refel, J);
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
        loc2glb = grid_data.loc2glb[:,eid];           # global indices of this element's nodes for extracting values from var arrays
        
        if !rhs_only
            Astart = (eid-1)*Np*dofs_per_node*Np*dofs_per_node + 1; # The segment of AI, AJ, AV for this element
            
        end
        
        volargs = (var, eid, 0, grid_data, geo_factors, refel, t, dt, stiffness, mass);
        
        if dofs_per_node == 1
            linchunk = linear.func(volargs);  # get the elemental linear part
            b[loc2glb] .+= linchunk;
            
            if !rhs_only
                bilinchunk = bilinear.func(volargs); # the elemental bilinear part
                #A[loc2glb, loc2glb] .+= bilinchunk;         # This will be very inefficient for sparse A
                for jj=1:Np
                    offset = Astart - 1 + (jj-1)*Np;
                    for ii=1:Np
                        AI[offset + ii] = loc2glb[ii];
                        AJ[offset + ii] = loc2glb[jj];
                        AV[offset + ii] = bilinchunk[ii, jj];
                    end
                end
            end
            
        elseif typeof(var) == Variable
            # only one variable, but more than one dof
            linchunk = linear.func(volargs);
            insert_linear!(b, linchunk, loc2glb, 1:dofs_per_node, dofs_per_node);
            
            if !rhs_only
                bilinchunk = bilinear.func(volargs);
                insert_bilinear!(AI, AJ, AV, Astart, bilinchunk, loc2glb, 1:dofs_per_node, dofs_per_node);
            end
            
        else
            linchunk = linear.func(volargs);
            insert_linear!(b, linchunk, loc2glb, 1:dofs_per_node, dofs_per_node);
            
            if !rhs_only
                bilinchunk = bilinear.func(volargs);
                insert_bilinear!(AI, AJ, AV, Astart, bilinchunk, loc2glb, 1:dofs_per_node, dofs_per_node);
            end
        end
    end
    
    # surface assembly for dg
    if !(face_linear===nothing) && (rhs_only || !(face_bilinear===nothing))
        for fid = 1:size(grid_data.face2glb,3)
            frefelind = [grid_data.faceRefelInd[1,fid], grid_data.faceRefelInd[2,fid]]; # refel based index of face in both elements
            
            face2glb = grid_data.face2glb[:,:,fid]; # local to global for face
            if frefelind[2] > 0 # not a boundary face
                flocal = [refel.face2local[frefelind[1]], refel.face2local[frefelind[2]]]; # local indices of face in both elements
            else # a boundary face
                face2glb[:,2] = face2glb[:,1]; # same global index for inside and out
                frefelind[2] = frefelind[1];
                flocal = [refel.face2local[frefelind[1]], refel.face2local[frefelind[1]]]; # local indices of face in both elements
            end
            eid = grid_data.face2element[1,fid]; # element on side 1
            
            faceBID = mesh_data.bdryID[fid]; # BID of face (0 for interior)
            
            surfargs = (var, eid, fid, grid_data, geo_factors, refel, t, dt);
            
            if dofs_per_node == 1
                linchunk = face_linear.func(surfargs);  # get the elemental linear part
                #linchunk = the_desired(rhsargs);
                # linchunk[1] = linchunk[1][flocal[1]];
                # linchunk[2] = linchunk[2][flocal[2]];
                if faceBID == 0
                    b[face2glb[:,1]] .+= linchunk[1];
                    b[face2glb[:,2]] .+= linchunk[2];
                elseif grid_data.face2element[1,fid] == 0
                    b[face2glb[:,2]] .+= linchunk[2];
                else
                    b[face2glb[:,1]] .+= linchunk[1];
                end
                
                if !rhs_only
                    bilinchunk = face_bilinear.func(surfargs); # the elemental bilinear part
                    
                    if typeof(bilinchunk[1]) <: Number
                        bilinchunk = [[bilinchunk[1]], [bilinchunk[2]], [bilinchunk[3]], [bilinchunk[4]]];
                    else
                        # bilinchunk[1] = bilinchunk[1][flocal[1],flocal[1]];
                        # bilinchunk[2] = bilinchunk[2][flocal[1],flocal[2]];
                        # bilinchunk[3] = bilinchunk[3][flocal[2],flocal[1]];
                        # bilinchunk[4] = bilinchunk[4][flocal[2],flocal[2]];
                    end
                    #bilinchunk = temporary_lhs_func(lhsargs);
                    #println(bilinchunk)
                    nfp = size(face2glb,1);
                    for jj=1:nfp
                        append!(AI, face2glb[:,1]);
                        append!(AJ, face2glb[jj,1]*ones(Int, nfp));
                        append!(AV, bilinchunk[1][:,jj]);
                        
                        append!(AI, face2glb[:,1]);
                        append!(AJ, face2glb[jj,2]*ones(Int, nfp));
                        append!(AV, bilinchunk[2][:,jj]);
                        
                        append!(AI, face2glb[:,2]);
                        append!(AJ, face2glb[jj,1]*ones(Int, nfp));
                        append!(AV, bilinchunk[3][:,jj]);
                        
                        append!(AI, face2glb[:,2]);
                        append!(AJ, face2glb[jj,2]*ones(Int, nfp));
                        append!(AV, bilinchunk[4][:,jj]);
                    end
                end
            end	 	
        end
    end
    
    loop_time = Base.Libc.time() - loop_time;
    
    if !rhs_only
        # Build the sparse A. Uses default + to combine overlaps
        A = sparse(AI, AJ, AV);
    end
    
    # Boundary conditions
    bc_time = Base.Libc.time();
    if rhs_only
        b = apply_boundary_conditions_rhs_only(var, b, t);
    else
        (A, b) = apply_boundary_conditions_lhs_rhs(var, A, b, t);
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
            totalcomponents = totalcomponents + length(var[vi].symvar);
        end
        for vi=1:length(var)
            components = length(var[vi].symvar);
            for compi=1:components
                if stepper.type == EULER_EXPLICIT
                    var[vi].values[compi,:] += sol[(compi+tmp):totalcomponents:end];
                else
                    var[vi].values[compi,:] = sol[(compi+tmp):totalcomponents:end];
                end
                tmp = tmp + 1;
            end
        end
    else
        components = length(var.symvar);
        for compi=1:components
            if stepper.type == EULER_EXPLICIT
                var.values[compi,:] += sol[compi:components:end];
            else
                var.values[compi,:] = sol[compi:components:end];
            end
        end
    end
end

function get_var_vals(var)
    # place the variable values in a vector
    if typeof(var) <: Array
        tmp = 0;
        totalcomponents = 0;
        for vi=1:length(var)
            totalcomponents = totalcomponents + length(var[vi].symvar);
        end
        vect = zeros(totalcomponents * length(var[1].values[1,:]));
        for vi=1:length(var)
            components = length(var[vi].symvar);
            for compi=1:components
                vect[(compi+tmp):totalcomponents:end] = var[vi].values[compi,:];
                tmp = tmp + 1;
            end
        end
    else
        components = length(var.symvar);
        vect = zeros(components * length(var.values[1,:]));
        for compi=1:components
            vect[compi:components:end] = var.values[compi,:];
        end
    end
    
    return vect;
end

end #module