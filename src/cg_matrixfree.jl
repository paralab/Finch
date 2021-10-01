#=
Functions used by the CGSolver for matrix free solutions.
=#

# Note: This uses the conjugate gradient iterative method,
# which assumes an SPD matrix.
function solve_matrix_free_sym(var, bilinear, linear, stepper=nothing, t=0, dt=0)
    start_time = time_ns();
    tol = config.linalg_matfree_tol;
    maxiters = config.linalg_matfree_max;
    Np = refel.Np;
    nel = mesh_data.nel;
    N1 = size(grid_data.allnodes,2);
    multivar = typeof(var) <: Array;
    if multivar
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        var_to_dofs = [];
        for vi=1:length(var)
            tmp = dofs_per_node;
            dofs_per_node += length(var[vi].symvar);
            push!(var_to_dofs, (tmp+1):dofs_per_node);
        end
    else
        # one variable
        dofs_per_node = length(var.symvar);
    end
    Nn = dofs_per_node * N1;
    
    # if dofs_per_node > 1
    #     println("Oops. Still working on multi DOF matrix-free");
    #     return zeros(Nn);
    # end
    
    if prob.time_dependent && !(stepper === nothing)
        #TODO
    else
        # Use regular rhs assembly
        b = assemble_rhs_only(var, linear, t, dt);
        normb = norm(b, Inf);
        
        x = zeros(Nn);
        
        # Do initial matvec
        #Ax = elem_matvec(x, bilinear, dofs_per_node, var);
        # nevermind, this will just be zeros
        
        r0 = copy(b); # = b - Ax;
        p = copy(r0);
        
        iter = 0;
        err = 1;
        while iter < maxiters && err > tol
            iter = iter+1;
            
            Ap = elem_matvec(p,bilinear, dofs_per_node, var, t, dt);
            
            alpha = dot(r0,r0) / dot(p,Ap);
            
            x = x .+ alpha.*p;
            r1 = r0 .- alpha.*Ap;
            
            beta = dot(r1,r1) / dot(r0, r0);
            p = r1 .+ beta.*p;
            
            r0 = copy(r1);
            err = norm(r0, Inf)/normb;
            
            if iter%50 == 0
                println("iteration "*string(iter)*": res = "*string(err));
            end
        end
        
        #println("Converged to "*string(err)*" in "*string(iter)*" iterations");
        log_entry("Converged to "*string(err)*" in "*string(iter)*" iterations");
        
        total_time = time_ns() - start_time;
        log_entry("Matrix-free solve took "*string(total_time/1e9)*" seconds");
        
        return x;
    end
end

# Note: This uses the stabilized biconjugate gradient method which works for nonsymmetric matrices
function solve_matrix_free_asym(var, bilinear, linear, stepper=nothing, t=0, dt=0)
    start_time = time_ns();
    tol = config.linalg_matfree_tol;
    maxiters = config.linalg_matfree_max;
    Np = refel.Np;
    nel = mesh_data.nel;
    N1 = size(grid_data.allnodes,2);
    multivar = typeof(var) <: Array;
    if multivar
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        var_to_dofs = [];
        for vi=1:length(var)
            tmp = dofs_per_node;
            dofs_per_node += length(var[vi].symvar);
            push!(var_to_dofs, (tmp+1):dofs_per_node);
        end
    else
        # one variable
        dofs_per_node = length(var.symvar);
    end
    Nn = dofs_per_node * N1;
    
    if prob.time_dependent && !(stepper === nothing)
        #TODO
    else
        # Use regular rhs assembly
        b = assemble_rhs_only(var, linear, t, dt);
        normb = norm(b, Inf);
        
        x = zeros(Nn);
        
        r0 = copy(b); # = b - Ax; but x is initially zeros
        rs = copy(b);
        rho0 = 1;
        alpha = 1;
        w0 = 1;
        p = zeros(Nn);
        Ap = zeros(Nn);
        As = zeros(Nn);
        
        iter = 0;
        err = 1;
        while iter < maxiters && err > tol
            iter = iter+1;
            
            rho1 = dot(rs,r0);
            beta = (rho1*alpha) / (rho0*w0);
            p = r0 + beta*(p - w0*Ap);
            
            Ap = elem_matvec(p,bilinear, dofs_per_node, var, t, dt);
            
            alpha = rho1/dot(rs,Ap);
            s = r0 - alpha*Ap;
            
            As = elem_matvec(s,bilinear, dofs_per_node, var, t, dt);
            
            w1 = dot(s,As) / dot(As,As);
            x = x + alpha*p + w1*s;
            r0 = s - w1*As;
            
            err = norm(r0, Inf)/normb;
            
            if iter%50 == 0
                println("iteration "*string(iter)*": res = "*string(err));
            end
        end
        
        println("Converged to "*string(err)*" in "*string(iter)*" iterations");
        log_entry("Converged to "*string(err)*" in "*string(iter)*" iterations");
        
        total_time = time_ns() - start_time;
        log_entry("Matrix free solve took "*string(total_time/1e9)*" seconds");
        
        return x;
    end
end

# Does the elemental matvec Ax=b
# x is input, b is output, A is from bilinear
function elem_matvec(x, bilinear, dofs_per_node, var, t = 0.0, dt = 0.0)
    Np = refel.Np;
    nel = mesh_data.nel;
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
            
            # # check for neumann bcs
            # for bi=1:length(prob.bc_type[vi,:])
            #     if prob.bc_type[vi,bi] == NEUMANN
                    
            #     end
            # end
        end
        
    else
        # one variable
        dofs_per_node = length(var.symvar);
        maxvarindex = var.index;
    end
    
    Ax = zeros(size(x));
    
    # Stiffness and mass are precomputed for uniform grid meshes
    if config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
        glb = grid_data.loc2glb[:,1];
        xe = grid_data.allnodes[:,glb[:]];
        (detJ, J) = geometric_factors(refel, xe);
        wgdetj = refel.wg .* detJ;
        if config.dimension == 1
            (RQ1, RD1) = build_deriv_matrix(refel, J);
            TRQ1 = RQ1';
            stiffness = [(TRQ1 * diagm(wgdetj) * RQ1)];
        elseif config.dimension == 2
            (RQ1, RQ2, RD1, RD2) = build_deriv_matrix(refel, J);
            (TRQ1, TRQ2) = (RQ1', RQ2');
            stiffness = [(TRQ1 * diagm(wgdetj) * RQ1) , (TRQ2 * diagm(wgdetj) * RQ2)];
        else
            (RQ1, RQ2, RQ3, RD1, RD2, RD3) = build_deriv_matrix(refel, J);
            (TRQ1, TRQ2, TRQ3) = (RQ1', RQ2', RQ3');
            stiffness = [(TRQ1 * diagm(wgdetj) * RQ1) , (TRQ2 * diagm(wgdetj) * RQ2) , (TRQ3 * diagm(wgdetj) * RQ3)];
        end
        mass = (refel.Q)' * diagm(wgdetj) * refel.Q;
    else
        stiffness = 0;
        mass = 0;
    end
    
    #Elemental loop follows elemental ordering
    for e=elemental_order;
        glb = grid_data.loc2glb[:,e];                 # global indices of this element's nodes for extracting values from var arrays
        xe = grid_data.allnodes[:,glb[:]];  # coordinates of this element's nodes for evaluating coefficient functions
        
        subx = extract_linear(x, glb, dofs_per_node);
        
        lhsargs = (var, xe, glb, refel, LHS, t, dt, stiffness, mass);
        bilinchunk = bilinear.func(lhsargs); # the elemental bilinear part
        # Neumann bcs need to be applied to this
        
        
        if dofs_per_node == 1
            #bilinchunk = bilinear.func(lhsargs); # the elemental bilinear part
            Ax[glb] = Ax[glb] + bilinchunk * subx;
        elseif typeof(var) == Variable
            # only one variable, but more than one dof
            #bilinchunk = bilinear.func(lhsargs);
            insert_linear_matfree!(Ax, bilinchunk*subx, glb, 1:dofs_per_node, dofs_per_node);
        else
            #bilinchunk = bilinear.func(lhsargs);
            insert_linear_matfree!(Ax, bilinchunk*subx, glb, 1:dofs_per_node, dofs_per_node);
        end
    end
    
    # Apply boudary conditions
    bidcount = length(grid_data.bids); # the number of BIDs
    if dofs_per_node > 1
        if multivar
            d = 0;
            for vi=1:length(var)
                for compo=1:length(var[vi].symvar)
                    d = d + 1;
                    for bid=1:bidcount
                        if prob.bc_type[var[vi].index, bid] == NO_BC
                            # do nothing
                        else
                            Ax = dirichlet_bc_matfree(Ax, x, grid_data.bdry[bid], d, dofs_per_node);
                        end
                    end
                end
            end
        else
            for bid=1:bidcount
                if prob.bc_type[var.index, bid] == NO_BC
                    # do nothing
                else
                    Ax = dirichlet_bc_matfree(Ax, x, grid_data.bdry[bid], 1:dofs_per_node, dofs_per_node);
                end
            end
        end
    else
        for bid=1:bidcount
            if prob.bc_type[var.index, bid] == NO_BC
                # do nothing
            else
                Ax = dirichlet_bc_matfree(Ax, x, grid_data.bdry[bid]);
            end
        end
    end
    
    # Reference points
    if size(prob.ref_point,1) >= maxvarindex
        if multivar
            posind = zeros(Int,0);
            vals = zeros(0);
            for vi=1:length(var)
                if prob.ref_point[var[vi].index,1]
                    eii = prob.ref_point[var[vi].index, 2];
                    tmp = (grid_data.glbvertex[eii[1], eii[2]] - 1)*dofs_per_node + var_to_dofs[vi][1];
                    if length(prob.ref_point[var[vi].index, 3]) > 1
                        tmp = tmp:(tmp+length(prob.ref_point[var[vi].index, 3])-1);
                    end
                    posind = [posind; tmp];
                    vals = [vals; prob.ref_point[var[vi].index, 3]];
                end
            end
            if length(vals) > 0
                Ax[posind] = vals;
            end
            
        else
            if prob.ref_point[var.index,1]
                eii = prob.ref_point[var.index, 2];
                posind = (grid_data.glbvertex[eii[1], eii[2]] - 1)*dofs_per_node + 1;
                if length(prob.ref_point[var.index, 3]) > 1
                    posind = posind:(posind+length(prob.ref_point[var[vi].index, 3])-1);
                end
                Ax[posind] = prob.ref_point[var.index, 3];
            end
        end
    end
    
    return Ax;
end

# sets b[bdry] = x[bdry]
function dirichlet_bc_matfree(b, x, bdryind, dofind=1, totaldofs=1)
    if totaldofs > 1
        for d=dofind
            ind2 = totaldofs.*(bdryind .- 1) .+ d;
            b[ind2] = x[ind2];
        end
    else
        b[bdryind] = x[bdryind];
    end
    
    return b;
end

function extract_linear(b, glb, dofs)
    if dofs == 1
        return b[glb];
    else
        np = length(glb);
        part = zeros(np*dofs);
        
        for d=1:dofs
            ind1 = ((d-1)*np+1):(d*np);
            ind2 = dofs.*(glb .- 1) .+ d;
            part[ind1] = b[ind2];
        end
        
        return part;
    end
end

# Inset the single dof into the greater construct
function insert_linear_matfree!(b, bel, glb, dof, Ndofs)
    # group nodal dofs
    np = length(glb);
    for d=1:length(dof)
        ind1 = ((d-1)*np+1):(d*np);
        ind2 = Ndofs.*(glb .- 1) .+ d;
        
        b[ind2] = b[ind2] + bel[ind1];
    end
end

function block_to_interlace(x, Np, dofs)
    newx = similar(x);
    for d=1:dofs
        newx[d:dofs:end] = x[((d-1)*Np+1):(d*Np)];
    end
    return newx;
end
function interlace_to_block(x, Np, dofs)
    newx = similar(x);
    for d=1:dofs
        newx[((d-1)*Np+1):(d*Np)] = x[d:dofs:end];
    end
    return newx;
end