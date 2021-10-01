# Apply boundary conditions for FV

function apply_boundary_conditions_to_face_rhs(var, fid, facefluxvec, t)
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
    
    # Boundary conditions are applied to flux
    fbid = grid_data.facebid[fid]; # BID of this face
    if fbid > 0
        facex = grid_data.allnodes[:, grid_data.face2glb[:,1,fid]];  # face node coordinates
        
        if multivar
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
                        
                        bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var[vi].index, fbid][compo], facex, t) .* geo_factors.area[fid];
                        
                        fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ geo_factors.volume[eid];
                        facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
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
                    bflux = FV_flux_bc_rhs_only_simple(prob.bc_func[var.index, fbid][d], facex, t) .* geo_factors.area[fid];
                    
                    fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ geo_factors.volume[eid];
                    facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                else
                    printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, fbid]);
                end
            end
        end
    end# BCs
    
    
end



function FV_flux_bc_rhs_only(val, facex, Qvec, t=0, dofind = 1, totaldofs = 1)
    if typeof(val) <: Number
        return val;
        
    elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
        bvals = zeros(size(facex,2));
        if config.dimension == 1
            for i=1:length(bvals)
                bvals[i]=val.value[1].func(facex[1,i],0,0,t);
            end
        elseif config.dimension == 2
            for i=1:length(bvals)
                bvals[i]=val.value[1].func(facex[1,i],facex[2,i],0,t);
            end
        else
            for i=1:length(bvals)
                bvals[i]=val.value[1].func(facex[1,i],facex[2,i],facex[3,i],t);
            end
        end
        
    elseif typeof(val) == GenFunction
        bvals = zeros(size(facex,2));
        if config.dimension == 1
            for i=1:length(bvals)
                bvals[i]=val.func(facex[1,i],0,0,t);
            end
        elseif config.dimension == 2
            for i=1:length(bvals)
                bvals[i]=val.func(facex[1,i],facex[2,i],0,t);
            end
        else
            for i=1:length(bvals)
                bvals[i]=val.func(facex[1,i],facex[2,i],facex[3,i],t);
            end
        end
    end
    
    # Do quadrature over the face
    b = Qvec * bvals;
    
    return b[1];
end

function FV_flux_bc_rhs_only_simple(val, bctype, facex, t=0)
    if typeof(val) <: Number
        return val;
        
    elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
        if config.dimension == 1
            bval = val.value[1].func(facex[1],0,0,t);
        elseif config.dimension == 2
            bval = val.value[1].func(facex[1],facex[2],0,t);
        else
            bval = val.value[1].func(facex[1],facex[2],facex[3],t);
        end
        
    elseif typeof(val) == GenFunction
        if config.dimension == 1
            bval = val.func(facex[1],0,0,t);
        elseif config.dimension == 2
            bval = val.func(facex[1],facex[2],0,t);
        else
            bval = val.func(facex[1],facex[2],facex[3],t);
        end
    end
    
    return bval;
end