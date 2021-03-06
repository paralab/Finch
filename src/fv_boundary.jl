# Apply boundary conditions for FV

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

function FV_flux_bc_rhs_only_simple(val, fid, t=0)
    # if typeof(val) <: Number
    #     return val;
        
    # elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
    #     if config.dimension == 1
    #         bval = val.value[1].func(facex[1],0,0,t);
    #     elseif config.dimension == 2
    #         bval = val.value[1].func(facex[1],facex[2],0,t);
    #     else
    #         bval = val.value[1].func(facex[1],facex[2],facex[3],t);
    #     end
        
    # elseif typeof(val) == GenFunction
    #     if config.dimension == 1
    #         bval = val.func(facex[1],0,0,t);
    #     elseif config.dimension == 2
    #         bval = val.func(facex[1],facex[2],0,t);
    #     else
    #         bval = val.func(facex[1],facex[2],facex[3],t);
    #     end
    # end
    eid = fv_grid.face2element[1,fid];
    
    return FV_evaluate_bc(val, eid, fid, t);
end

function FV_copy_bdry_vals_to_vector(var, vec, grid, dofs_per_node)
    if typeof(var) <: Array
        dofind = 0;
        for vi=1:length(var)
            for compo=1:length(var[vi].symvar)
                dofind = dofind + 1;
                for bid=1:size(prob.bc_type,2)
                    if prob.bc_type[var[vi].index, bid] == DIRICHLET
                        for i = 1:length(grid.bdryface[bid]) # loop over faces with this BID
                            fid = grid.bdryface[bid][i];
                            eid = grid.face2element[1,fid];
                            vec[(eid-1)*dofs_per_node + dofind] = var[vi].values[compo, eid];
                        end
                    end
                end
            end
        end
    else
        for d=1:dofs_per_node
            dofind = d;
            for bid=1:size(prob.bc_type,2)
                if prob.bc_type[var.index, bid] == DIRICHLET
                    for i = 1:length(grid.bdryface[bid]) # loop over faces with this BID
                        fid = grid.bdryface[bid][i];
                        eid = grid.face2element[1,fid];
                        vec[(eid-1)*dofs_per_node + dofind] = var.values[d, eid];
                    end
                end
            end
        end
    end
end

# This evaluates the BC at a specific node.
# That could mean:
# - the value of constant BCs
# - evaluate a genfunction.func for BCs defined by strings->genfunctions
# - evaluate a function for BCs defined by callback functions
function FV_evaluate_bc(val, eid, fid, t)
    dim = config.dimension;
    if typeof(val) <: Number
        result = val;
        
    elseif typeof(val) == Coefficient 
        facex = fv_info.faceCenters[:,fid];
        result = evaluate_coefficient(val, 1, facex, t, eid, fid);
        
    elseif typeof(val) == GenFunction
        facex = fv_info.faceCenters[:,fid];
        if dim == 1
            result=val.func(facex[1],0,0,t,eid,fid);
        elseif dim == 2
            result=val.func(facex[1],facex[2],0,t,eid,fid);
        else
            result=val.func(facex[1],facex[2],facex[3],t,eid,fid);
        end
        
    elseif typeof(val) == CallbackFunction
        facex = fv_info.faceCenters[:,fid];
        #form a dict for the arguments. x,y,z,t are always included
        arg_list = [];
        append!(arg_list, [("x", facex[1]),
                    ("y", dim>1 ? facex[2] : 0),
                    ("z", dim>2 ? facex[3] : 0),
                    ("t", t)]);
        # Add in the requires list
        for r in val.args
            foundit = false;
            # is it an entity?
            if typeof(r) == Variable
                foundit = true;
                if size(r.values,1) == 1
                    push!(arg_list, (string(r.symbol), r.values[1,eid]));
                else
                    push!(arg_list, (string(r.symbol), r.values[:,eid]));
                end
            elseif typeof(r) == Coefficient
                foundit = true;
                cvals = zeros(size(r.value));
                facex = fv_info.faceCenters[:,fid];
                for i=1:length(cvals)
                    cvals[i] = evaluate_coefficient(r, i, facex, t, eid, fid);
                end
                push!(arg_list, (string(r.symbol), cvals));
            elseif typeof(r) == Indexer
                foundit = true;
                push!(arg_list, (string(r.symbol), r.value));
            elseif typeof(r) == String
                # This could also be a variable, coefficient, or special thing like normal, 
                if r in ["x","y","z","t"]
                    foundit = true;
                    # These are already included
                elseif r == "normal"
                    foundit = true;
                    push!(arg_list, ("normal", fv_grid.facenormals[:,fid]));
                else
                    for v in variables
                        if string(v.symbol) == r
                            foundit = true;
                            if size(v.values,1) == 1
                                push!(arg_list, (r, v.values[1,eid]));
                            else
                                push!(arg_list, (r, v.values[:,eid]));
                            end
                            break;
                        end
                    end
                    for c in coefficients
                        if string(c.symbol) == r
                            foundit = true;
                            push!(arg_list, (r, FV_evaluate_bc(c, eid, fid, t)));
                            break;
                        end
                    end
                end
            else
                # What else could it be?
            end

            if !foundit
                # Didn't figure this thing out.
                printerr("Unknown requirement for callback function "*string(fun)*" : "*string(requires[i]));
            end
        end
        
        # Build the dict
        args = Dict(arg_list);
        
        # call the function
        result = val.func(args);
        
    end
    return result;
end