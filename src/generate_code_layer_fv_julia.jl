#=
Code generation functions for FV for Julia.
=#

# Extract the input args. This must match the arguments passed by the solver.
function handle_input_args_fv_julia(lorr, vors)
    # args = (var, eid, fid, grid_data, geo_factors, fv_info, t, dt)
    code =
"var =       args[1];
eid =       args[2];
fid =       args[3];
grid =      args[4];
geo_facs =  args[5];
fv_data =   args[6];
refel =     args[7];
time =      args[8];
dt =        args[9];
"
    return code;
end

# Allocate, compute, or fetch all needed values
function prepare_needed_values_fv_julia(entities, var, lorr, vors)
    # Only gather the pieces that are used. See the end of this function.
    if vors == "volume"
        # el, nodex, loc2glb, detJ, J
        piece_needed = zeros(Bool, 5);
    else # surface
        # els, nodex, loc2glb, cellx, frefelind, facex, face2glb, normal, face_detJ, area, vol_J
        piece_needed = zeros(Bool, 11);
    end
    need_deriv_matrix = false;
    need_deriv_dist = false;
    
    code = "";
    
    # First label the index values
    for i=1:length(indexers)
        if indexers[i].symbol === :elements || indexers[i].symbol === :cells
            # not needed
        else
            code *= "INDEX_VAL_"*string(indexers[i].symbol)*" = kwargs["*string(i)*"]\n";
        end
    end
    
    for i=1:length(entities)
        cname = make_coef_name(entities[i]);
        if is_unknown_var(entities[i], var) && lorr == LHS
            # TODO
            
        else
            # Is coefficient(number or function) or variable(array)?
            (ctype, cval) = get_coef_val(entities[i]);
            if ctype == -1
                # It was a special symbol like dt or FACENORMAL
                if entities[i].name == "FACENORMAL1"
                    code *= cname * " = normal["*string(entities[i].index)*"]; # normal vector component\n"
                    piece_needed[8] = true;
                else entities[i].name == "FACENORMAL2"
                    code *= cname * " = -normal["*string(entities[i].index)*"]; # reverse normal vector component\n"
                    piece_needed[8] = true;
                end
            elseif ctype == 0
                # It was a number, do nothing?
            elseif ctype == 1 # a constant wrapped in a coefficient
                # This generates something like: coef_k_i = 4;
                if length(entities[i].derivs) > 0
                    code *= cname * " = 0; # NOTE: derivative applied to constant coefficient = 0\n";
                else
                    code *= cname * " = " * string(cval) * ";\n";
                end
                
            elseif ctype == 2 # a coefficient function
                # This generates something like:
                ######################################
                # coef_n_i = zeros(refel.Np); # allocate
                # for coefi = 1:refel.Np
                #     coef_k_i[coefi] = (Finch.genfunctions[cval]).func(x[1,coefi], x[2,coefi],x[3,coefi],time); # evaluate at nodes
                # end
                ######################################
                if vors == "surface" && length(entities[i].derivs) == 0
                    nodesymbol = "facex"
                    piece_needed[6] = true;
                else
                    nodesymbol = "nodex"
                    piece_needed[2] = true;
                end
                # cargs = "("*nodesymbol*"[coefi], 0, 0, time)";
                # if config.dimension == 2
                #     cargs = "("*nodesymbol*"[1, coefi], "*nodesymbol*"[2, coefi], 0, time)";
                # elseif config.dimension == 3
                #     cargs = "("*nodesymbol*"[1, coefi], "*nodesymbol*"[2, coefi], "*nodesymbol*"[3, coefi], time)";
                # end
                cargs = "("*nodesymbol*"[:,coefi], time, eid, fid)";
                coef_index = get_coef_index(entities[i]);
                
                if vors == "volume"
                    code *= cname * " = zeros(refel.Np);\n";
                    # code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargs * " end\n";
                    code *= "for coefi = 1:refel.Np " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", "*nodesymbol*"[:,coefi], time, eid, fid) end\n";
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RD"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                        need_deriv_matrix = true;
                        piece_needed[5] = true;
                    end
                    # integrate over cell
                    code *= cname * " = (refel.wg .* detj)' * refel.Q * " * cname * "; # integrate over cell\n";
                    piece_needed[4] = true;
                    
                else # surface
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        code *= cname * " = zeros(refel.Np);\n";
                        # code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargs * " end\n";
                        code *= "for coefi = 1:refel.Np " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", "*nodesymbol*"[:,coefi], time, eid, fid) end\n";
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RD"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                        code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]]; # extract face values only.";
                        need_deriv_matrix = true;
                        piece_needed[[4,5,8]] .= true;
                        
                    else # no derivatives, only need surface nodes
                        code *= cname * " = zeros(refel.Nfp[frefelind[1]]);\n";
                        # code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargs * " end\n";
                        code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", "*nodesymbol*"[:,coefi], time, eid, fid) end\n";
                        piece_needed[5] = true;
                    end
                    # integrate over face
                    if config.dimension == 1
                        # in 1d there is only one face node
                        code *= cname * " = " * cname * "[1]\n";
                    else
                        code *= cname * " = (refel.surf_wg[frefelind[1]] .* face_detJ)' * refel.surf_Q[frefelind[1]][:, refel.face2local[frefelind[1]]] * " * cname * " / area; # integrate over face\n";
                        piece_needed[9] = true;
                        piece_needed[10] = true;
                    end
                end
                
            elseif ctype == 3 # a known variable value
                piece_needed[1] = true;
                # This generates something like: coef_u_1 = copy((Finch.variables[1]).values[1, gbl])
                cellside = 0; # 0 means no side flag
                for flagi=1:length(entities[i].flags)
                    if occursin("DGSIDE1", entities[i].flags[flagi]) || occursin("CELL1", entities[i].flags[flagi])
                        cellside = 1;
                    elseif occursin("DGSIDE2", entities[i].flags[flagi]) || occursin("CELL2", entities[i].flags[flagi])
                        cellside = 2;
                    end
                end
                if vors == "surface"
                    l2gsymbol = "els[1]"
                else
                    l2gsymbol = "el"
                end
                if cellside == 1
                    l2gsymbol = "els[1]"
                elseif cellside == 2
                    l2gsymbol = "els[2]"
                end
                
                if typeof(entities[i].index) <: Array
                    # It is an indexed variable
                    if length(entities[i].index) == 1
                        indstr = "INDEX_VAL_"*entities[i].index[1];
                    else
                        # There is more than one index. Need to form an expression for it.
                        indstr = "(INDEX_VAL_"*entities[i].index[1];
                        indices = variables[cval].indexer;
                        for indi=2:length(entities[i].index)
                            indstr *= " + ("*string(length(indices[indi-1].range))*"*(INDEX_VAL_"*entities[i].index[indi]*"-1)";
                        end
                        for indi=1:length(entities[i].index)
                            indstr *= ")";
                        end
                    end
                    
                else
                    indstr = string(entities[i].index);
                end
                    
                if vors == "volume"
                    # Apply any needed derivative operators.
                    if length(entities[i].derivs) > 0
                        # Need a derivative here...
                        # TODO
                    else
                        code *= cname * " = Finch.variables["*string(cval)*"].values["*indstr*", "*l2gsymbol*"];\n";
                    end
                else
                    if length(entities[i].derivs) > 0
                        code *= cname * " = Finch.variables["*string(cval)*"].values["*indstr*", els[2]] - Finch.variables["*string(cval)*"].values["*indstr*", els[1]];\n";
                        code *= cname * " = (els[1] != els[2] && abs(normal["*string(entities[i].derivs[1])*"]) > 1e-10) ? "*cname*" ./ dxyz["*string(entities[i].derivs[1])*"]  : 0\n"
                        need_deriv_dist = true;
                        piece_needed[4] = true; # cellx
                        piece_needed[8] = true; # normal
                    else
                        if cellside == 0
                            # No side was specified, so use the average
                            code *= cname * " = 0.5 * (Finch.variables["*string(cval)*"].values["*indstr*", els[1]] + Finch.variables["*string(cval)*"].values["*indstr*", els[2]]);\n";
                        else
                            code *= cname * " = Finch.variables["*string(cval)*"].values["*indstr*", "*l2gsymbol*"];\n";
                        end
                    end
                end
                
            elseif ctype == 4 # an indexed coefficient
                # This generates something like:
                ######################################
                # coef_n_1 = zeros(refel.Np); # allocate
                # for coefi = 1:refel.Np
                #     coef_k_1[coefi] = (Finch.coefficients[cval]).value[INDEX_VAL_i].func(x[1,coefi], x[2,coefi],x[3,coefi],time); # evaluate at nodes
                # end
                ######################################
                if vors == "surface" && length(entities[i].derivs) == 0
                    nodesymbol = "facex"
                    piece_needed[6] = true;
                else
                    nodesymbol = "nodex"
                    piece_needed[2] = true;
                end
                # cargs = "("*nodesymbol*"[coefi], 0, 0, time)";
                # if config.dimension == 2
                #     cargs = "("*nodesymbol*"[1, coefi], "*nodesymbol*"[2, coefi], 0, time)";
                # elseif config.dimension == 3
                #     cargs = "("*nodesymbol*"[1, coefi], "*nodesymbol*"[2, coefi], "*nodesymbol*"[3, coefi], time)";
                # end
                
                indstr = "";
                for indi=1:length(entities[i].index)
                    if indi>1
                        indstr *= ",";
                    end
                    indstr *= "INDEX_VAL_"*entities[i].index[indi];
                end
                
                if vors == "volume"
                    code *= cname * " = zeros(refel.Np);\n";
                    # code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.coefficients["*string(cval)*"]).value["*indstr*"].func" * cargs * " end\n";
                    code *= "for coefi = 1:refel.Np " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(cval)*"], ["*indstr*"], "*nodesymbol*"[:,coefi], time, eid, fid) end\n";
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RD"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                        need_deriv_matrix = true;
                        piece_needed[5] = true;
                    end
                    # integrate over cell
                    code *= cname * " = (refel.wg .* detj)' * refel.Q * " * cname * "; # integrate over cell\n";
                    piece_needed[4] = true;
                    
                else # surface
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        code *= cname * " = zeros(refel.Np);\n";
                        # code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.coefficients["*string(cval)*"]).value["*indstr*"].func" * cargs * " end\n";
                        code *= "for coefi = 1:refel.Np " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(cval)*"], ["*indstr*"], "*nodesymbol*"[:,coefi], time, eid, fid) end\n";
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RD"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                        code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]]; # extract face values only.";
                        need_deriv_matrix = true;
                        piece_needed[[4,5,8]] .= true;
                        
                    else # no derivatives, only need surface nodes
                        code *= cname * " = zeros(refel.Nfp[frefelind[1]]);\n";
                        # code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = (Finch.coefficients["*string(cval)*"]).value["*indstr*"].func" * cargs * " end\n";
                        code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(cval)*"], ["*indstr*"], "*nodesymbol*"[:,coefi], time, eid, fid) end\n";
                        piece_needed[5] = true;
                    end
                    # integrate over face
                    if config.dimension == 1
                        # in 1d there is only one face node
                        code *= cname * " = " * cname * "[1]\n";
                    else
                        code *= cname * " = (refel.surf_wg[frefelind[1]] .* face_detJ)' * refel.surf_Q[frefelind[1]][:, refel.face2local[frefelind[1]]] * " * cname * " / area; # integrate over face\n";
                        piece_needed[9] = true;
                        piece_needed[10] = true;
                    end
                end
                
            end
        end # if coefficient
    end # entity loop
    
    
    if vors == "volume"
        # el, nodex, loc2glb, detJ, J
        needed_pieces = "";
        if piece_needed[1]
            needed_pieces *= "el = eid;\n"
        end
        if piece_needed[2] || piece_needed[3]
            needed_pieces *= "loc2glb = grid.loc2glb[:,eid]; # local to global map\n"
        end
        if piece_needed[2]
            needed_pieces *= "nodex = grid.allnodes[:,loc2glb]; # node coordinates\n"
        end
        if piece_needed[4]
            needed_pieces *= "detj = geo_factors.detJ[eid]; # geometric factors\n"
        end
        if piece_needed[5]
            needed_pieces *= "J = geo_factors.J[eid]; # geometric factors\n"
        end
        
        if need_deriv_matrix
            needed_pieces *= 
"
# Note on derivative matrices:
# RQn are vandermond matrices for the derivatives of the basis functions
# with Jacobian factors. They are made like this.
# |RQ1|   | rx sx tx || Qx |
# |RQ2| = | ry sy ty || Qy |
# |RQ3|   | rz sz tz || Qz |

"
            if config.dimension == 1
                needed_pieces *= "(RQ1, RD1) = build_deriv_matrix(refel, J);\n";
                needed_pieces *= "TRQ1 = RQ1';\n"
                
            elseif config.dimension == 2
                needed_pieces *= "(RQ1, RQ2, RD1, RD2) = build_deriv_matrix(refel, J);\n";
                needed_pieces *= "(TRQ1, TRQ2) = (RQ1', RQ2');\n"
                
            elseif config.dimension == 3
                needed_pieces *= "(RQ1, RQ2, RQ3, RD1, RD2, RD3) = build_deriv_matrix(refel, J);\n";
                needed_pieces *= "(TRQ1, TRQ2, TRQ3) = (RQ1', RQ2', RQ3');\n"
            end
        end
        if need_deriv_dist
            needed_pieces *= "TODO: neighboring cell-based derivative for volume. generate_code_layer_fv_julia.jl - prepare_needed_values_fv_julia()"
        end
        
    else # surface
        needed_pieces = 
"if grid.face2element[1, fid] == eid # The normal points out of e"*(piece_needed[8] ? "\n\tnormal = grid.facenormals[:, fid];" : "")*"
    neighbor = grid.face2element[2, fid];"*(piece_needed[5] ? "\n\tfrefelind = [grid.faceRefelInd[1,fid], grid.faceRefelInd[2,fid]]; # refel based index of face in both elements" : "")*"
else # The normal points into e"*(piece_needed[8] ? "\n\tnormal = -grid.facenormals[:, fid];" : "")*"
    neighbor = grid.face2element[1, fid];"*(piece_needed[5] ? "\n\tfrefelind = [grid.faceRefelInd[2,fid], grid.faceRefelInd[1,fid]]; # refel based index of face in both elements" : "")*"
end
if neighbor == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
    neighbor = eid;
end\n
";
        # els, nodex, loc2glb, cellx, frefelind, facex, face2glb, normal, face_detJ, area, vol_J
        if piece_needed[1]
            needed_pieces *= "els = (eid, neighbor); # indices of elements on both sides\n"
        end
        if piece_needed[2] || piece_needed[3]
            needed_pieces *= "loc2glb = (grid.loc2glb[:,eid], grid.loc2glb[:, neighbor]); # volume local to global\n"
        end
        if piece_needed[2]
            needed_pieces *= "nodex = (grid.allnodes[:,loc2glb[1][:]], grid.allnodes[:,loc2glb[2][:]]); # volume node coordinates\n"
        end
        if piece_needed[4]
            needed_pieces *= "cellx = (fv_data.cellCenters[:, eid], fv_data.cellCenters[:, neighbor]); # cell center coordinates\n"
        end
        if piece_needed[5] || piece_needed[6]
            needed_pieces *= "face2glb = grid.face2glb[:,:,fid];         # global index for face nodes for each side of each face\n"
        end
        if piece_needed[6]
            needed_pieces *= "facex = grid.allnodes[:, face2glb[:, 1]];  # face node coordinates\n"
        end
        if piece_needed[9]
            needed_pieces *= "face_detJ = geo_facs.face_detJ[fid];       # detJ on face\n"
        end
        if piece_needed[10]
            needed_pieces *= "area = geo_facs.area[fid];                 # area of face\n"
        end
        if piece_needed[11]
            needed_pieces *= "vol_J = (geo_facs.J[eid], geo_facs.J[neighbor]);\n"
        end
        
        if need_deriv_matrix
            needed_pieces *= "TODO: derivative matrices for surface. generate_code_layer_fv_julia.jl - prepare_needed_values_fv_julia()"
        end
        if need_deriv_dist
            needed_pieces *= "dxyz = norm(cellx[2] - cellx[1]) .* normal; # normal scaled by distance between cell centers\n"
        end
    end
    
    code = needed_pieces * "\n" * code;
    
    return code;
end

function make_elemental_computation_fv_julia(terms, var, dofsper, offset_ind, lorr, vors)
    # Here is where I make some assumption about the form of the expression.
    # Since it was expanded by the parser it should look like a series of terms: t1 + t2 + t3...
    # Where each term, if LHS, involves one unknown component and possibly some coefficients.
    code = "";
    
    # If there were no terms, 
    if length(terms) < 1
        code = "return nothing # There were no terms to compute here";
    end
    
    # Allocate the vector or matrix to be returned if needed
    if dofsper > 1
        if lorr == RHS
            code *= "cell_average = zeros("*string(dofsper)*"); # Allocate for returned values.\n"
        else
            if vors == "volume"
                code *= "cell_matrix = zeros(refel.Np * "*string(dofsper)*", "*string(dofsper)*"); # Allocate for returned matrix.\n"
            else
                code *= "cell_matrix = zeros(refel.Nfp[frefelind[1]] * "*string(dofsper)*", "*string(dofsper)*"); # Allocate for returned matrix.\n"
            end
        end
    end
    
    # Separate the factors of each term and form the calculation
    if dofsper > 1
        # Subvectors for each component
        if lorr == LHS
            subvector = Array{String, 2}(undef, dofsper, dofsper);
        else # RHS
            subvector = Array{String, 1}(undef, dofsper);
        end
        for smi=1:length(subvector)
            subvector[smi] = "";
        end
        
        if typeof(var) <: Array
            for vi=1:length(var) # variables
                # Process the terms for this variable
                for ci=1:length(terms[vi]) # components
                    for i=1:length(terms[vi][ci])
                        (term_result, var_ind) = generate_term_calculation_fv_julia(terms[vi][ci][i], var, lorr);
                        
                        # println(terms)
                        # println(terms[vi])
                        # println(terms[vi][ci])
                        # println(terms[vi][ci][i])
                        # println(term_result * " : "*string(test_ind)*", "*string(trial_ind))
                        
                        # Find the appropriate subvector for this term
                        subveci = offset_ind[vi] + ci;
                        subvecj = var_ind;
                        if lorr == LHS
                            subvec_ind = subveci + dofsper * (subvecj-1);
                        else
                            subvec_ind = subveci;
                        end
                        
                        if length(subvector[subvec_ind]) > 1
                            subvector[subvec_ind] *= " .+ " * term_result;
                        else
                            subvector[subvec_ind] = term_result;
                        end
                    end
                end
                
            end # vi
            
        else # only one variable
            # Process the terms for this variable
            for ci=1:length(terms) # components
                for i=1:length(terms[ci])
                    (term_result, var_ind) = generate_term_calculation_fv_julia(terms[ci][i], var, lorr);
                    
                    # Find the appropriate submatrix for this term
                    subvec_ind = ci;
                    
                    if length(subvector[subvec_ind]) > 1
                        subvector[subvec_ind] *= " .+ " * term_result;
                    else
                        subvector[subvec_ind] = term_result;
                    end
                end
            end
            
        end
        
        # Put the subvector together into cell_average or cell_matrix
        if lorr == LHS
            for emi=1:dofsper
                for emj=1:dofsper
                    if length(subvector[emi]) > 1
                        if vors == "volume"
                            rangei = "("*string(emi-1)*"*refel.Np + 1):("*string(emi)*"*refel.Np)";
                        else
                            rangei = "("*string(emi-1)*"*refel.Nfp[frefelind[1]] + 1):("*string(emi)*"*refel.Nfp[frefelind[1]])";
                        end
                        
                        code *= "cell_matrix["*rangei*"] = " * subvector[emi, emj] * "\n";
                    end
                end
            end
            code *= "return cell_matrix;\n"
            
        else # RHS
            for emi=1:dofsper
                if length(subvector[emi]) > 1
                    code *= "cell_average["*string(emi)*"] = " * subvector[emi] * "\n";
                end
            end
            code *= "return cell_average;\n"
        end
        
        
    else # one dof
        terms = terms[1];
        if lorr == LHS
            if vors == "volume"
                result = "zeros(refel.Np)";
            else
                result = "zeros(refel.Nfp[frefelind[1]])";
            end
            
        else
            result = "[0]";
        end
        
        #process each term
        first_term = true;
        for i=1:length(terms)
            (term_result, var_ind) = generate_term_calculation_fv_julia(terms[i], var, lorr);
            
            if !(term_result == "0")
                if first_term
                    result = term_result;
                    first_term = false;
                else
                    result *= " .+ " * term_result;
                end
            end
        end
        code *= "cell_average = [" * result * "];\n";
        code *= "return cell_average;\n"
    end
    
    return code;
end

function generate_term_calculation_fv_julia(term, var, lorr)
    result = "";
    # Note: separate_factors return test and trial info, but FV will not have any test functions.
    # trial_part refers to the unknown variable part.
    # example:
    # 0.1*D1__u_1*_FACENORMAL1_1    ->  ???               on LHS
    #                             ->  0.1 * (coef_D1xu_1 .* normal[1])     on RHS
    if lorr == LHS
        (test_part, var_part, coef_part, test_ind, var_ind) = separate_factors(term, var);
        # # LHS: ??
        # TODO
        #
        # Seriously, need to do
        
    else
        (test_part, var_part, coef_part, test_ind, var_ind) = separate_factors(term);
        # RHS: coef_part
        if !(coef_part === nothing)
            result = string(replace_entities_with_symbols(coef_part));
        else
            result = "0";
        end
    end
    
    return (result, var_ind);
end

# Generate the assembly loop structures and insert the content
function generate_assembly_loop_fv_julia(indices)
    # Each of the indices must be passed to the functions in a named tuple.
    # Pass all defined indexers.
    index_args = "";
    for i=1:length(indexers)
        if i > 1
            index_args *= ", ";
        end
        index_args *= "index_val_"*string(indexers[i].symbol)*"=indexing_variable_"*string(indexers[i].symbol);
    end
    if length(indexers) == 0
        index_args = "empty_index = 1";
    end
    
    # Make sure the elemental loop is included
    elements_included = false;
    for i=1:length(indices)
        if indices[i] == "elements" || indices[i] == "cells"
            elements_included = true;
        end
    end
    # If elements were not included, make them the outermost loop
    if !elements_included
        indices = ["elements"; indices];
    end
    
    code = 
"# Label things that were allocated externally
sourcevec = allocated_vecs[1];
fluxvec = allocated_vecs[2];
facefluxvec = allocated_vecs[3];
face_done = allocated_vecs[4];

"
    # generate the loop structures
    code *= "# Loops\n"
    loop_start = ""
    loop_end = ""
    ind_offset = ""
    prev_ind = ""
    for i=1:length(indices)
        if indices[i] == "elements" || indices[i] == "cells"
            loop_start *= "for eid in elemental_order\n";
            loop_start *= "face_done .= 0; # Reset face_done in elemental loop\n";
            loop_end *= "end # loop for elements\n";
        else
            loop_start *= "for indexing_variable_"*string(indices[i].symbol)*" in "*string(indices[i].range)*"\n";
            loop_end *= "end # loop for "*string(indices[i].symbol)*"\n";
            if length(prev_ind) == 0
                ind_offset *= "(indexing_variable_"*string(indices[i].symbol);
                prev_ind = string(length(indices[i].range));
            else
                ind_offset *= " + "*prev_ind*" * (indexing_variable_"*string(indices[i].symbol)*" - 1";
                prev_ind = string(length(indices[i].range));
            end
        end
    end
    for i=2:length(indices)
        ind_offset *= ")";
    end
    index_offset = "dofs_per_loop * ("*ind_offset*" - 1 + dofs_per_node * (eid - 1)) + 1"
    face_index_offset = "dofs_per_loop * ("*ind_offset*" - 1 + dofs_per_node * (fid - 1)) + 1"
    
    code *= loop_start;
    
    # insert the content
    code *= "
    # determine index in global vector
    index_offset = "*index_offset*";
    
    # Zero the result vectors for this element
    sourcevec[index_offset:(index_offset + dofs_per_loop-1)] .= 0;
    fluxvec[index_offset:(index_offset + dofs_per_loop-1)] .= 0;
    
    ##### Source integrated over the cell #####
    # Compute RHS volume integral
    if !(source_rhs === nothing)
        
        sourceargs = (var, eid, 0, grid_data, geo_factors, fv_info, refel, t, dt);
        source = source_rhs.func(sourceargs; "*index_args*") ./ geo_factors.volume[eid];
        
        # Add to global source vector
        sourcevec[index_offset:(index_offset + dofs_per_loop-1)] .= source;
    end

    ##### Flux integrated over the faces #####
    # Loop over this element's faces.
    for i=1:refel.Nfaces
        fid = grid_data.element2face[i, eid];
        
        # determine index in global face vector
        face_index_offset = "*face_index_offset*";
        
        if !(flux_rhs === nothing)
            if face_done[fid] < dofs_per_node
                face_done[fid] += dofs_per_loop;
                
                fluxargs = (var, eid, fid, grid_data, geo_factors, fv_info, refel, t, dt);
                flux = flux_rhs.func(fluxargs; "*index_args*") .* geo_factors.area[fid];
                
                # Add to global flux vector for faces
                facefluxvec[face_index_offset:(face_index_offset + dofs_per_loop-1)] .= flux;
                # Combine all flux for this element
                fluxvec[index_offset:(index_offset + dofs_per_loop-1)] .+= flux ./ geo_factors.volume[eid];
                
            else
                # This flux has either been computed or is being computed by another thread.
                # The state will need to be known before paralellizing, but for now assume it's complete.
                flux = -facefluxvec[face_index_offset:(face_index_offset + dofs_per_loop-1)];
                fluxvec[index_offset:(index_offset + dofs_per_loop-1)] .-= facefluxvec[face_index_offset:(face_index_offset + dofs_per_loop-1)] ./ geo_factors.volume[eid];
            end
        end
        
        # Boundary conditions are applied to flux
        fbid = grid_data.facebid[fid]; # BID of this face
        if fbid > 0
            facex = grid_data.allnodes[:, grid_data.face2glb[:,1,fid]];  # face node coordinates
            
            if typeof(var) <: Array
                dofind = 0;
                for vi=1:length(var)
                    compo = "*ind_offset*";
                    dofind = dofind + 1;
                    if prob.bc_type[var[vi].index, fbid] == NO_BC
                        # do nothing
                    elseif prob.bc_type[var[vi].index, fbid] == FLUX
                        # compute the value and add it to the flux directly
                        # Qvec = (refel.surf_wg[grid_data.faceRefelInd[1,fid]] .* geo_factors.face_detJ[fid])' * (refel.surf_Q[grid_data.faceRefelInd[1,fid]])[:, refel.face2local[grid_data.faceRefelInd[1,fid]]]
                        # Qvec = Qvec ./ geo_factors.area[fid];
                        # bflux = FV_flux_bc_rhs_only(prob.bc_func[var[vi].index, fbid][compo], facex, Qvec, t, dofind, dofs_per_node) .* geo_factors.area[fid];
                        
                        bflux = FVSolver.FV_flux_bc_rhs_only_simple(prob.bc_func[var[vi].index, fbid][compo], fid, t);
                        
                        fluxvec[index_offset + dofind-1] += (bflux - facefluxvec[face_index_offset + dofind-1]) ./ geo_factors.volume[eid];
                        facefluxvec[face_index_offset + dofind-1] = bflux;
                    else
                        printerr(\"Unsupported boundary condition type: \"*prob.bc_type[var[vi].index, fbid]);
                    end
                end
                
            else # one variable
                d="*ind_offset*";
                if prob.bc_type[var.index, fbid] == NO_BC
                    # do nothing
                elseif prob.bc_type[var.index, fbid] == FLUX
                    # compute the value and add it to the flux directly
                    # Qvec = (refel.surf_wg[grid_data.faceRefelInd[1,fid]] .* geo_factors.face_detJ[fid])' * (refel.surf_Q[grid_data.faceRefelInd[1,fid]])[:, refel.face2local[grid_data.faceRefelInd[1,fid]]]
                    # Qvec = Qvec ./ geo_factors.area[fid];
                    # bflux = FV_flux_bc_rhs_only(prob.bc_func[var.index, fbid][d], facex, Qvec, t, dofind, dofs_per_node) .* geo_factors.area[fid];
                    
                    bflux = FVSolver.FV_flux_bc_rhs_only_simple(prob.bc_func[var.index, fbid][d], fid, t);
                    
                    
                    fluxvec[index_offset] += (bflux - facefluxvec[face_index_offset]) ./ geo_factors.volume[eid];
                    facefluxvec[face_index_offset] = bflux;
                else
                    printerr(\"Unsupported boundary condition type: \"*prob.bc_type[var.index, fbid]);
                end
            end
        end# BCs
        
    end# face loop
"
    # close loops
    code *= loop_end;
    
    # finish
    code *=
"
if !(source_rhs === nothing)
    return sourcevec + fluxvec;
else
    return fluxvec;
end
"
    # println(code)
    return code;
end