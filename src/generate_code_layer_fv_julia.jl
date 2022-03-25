#=
Code generation functions for FV for Julia.
=#

# Extract the input args. This must match the arguments passed by the solver.
function handle_input_args_fv_julia(lorr, vors)
    if vors == "volume"
        # args = (var, eid, fid, grid_data, geo_factors, fv_info, t, dt)
        code =
"var =       args[1]; # variable list
eid =       args[2]; # element ID
fid =       args[3]; # NA
grid =      args[4]; # grid data struct
geo_facs =  args[5]; # geometric factors
fv_data =   args[6]; # FV specific info
refel =     args[7]; # reference element
time =      args[8]; # time
dt =        args[9]; # dt
"
    else # surface
        # args = (var, eid, fid, els, grid_data, geo_factors, fv_info, t, dt)
        code =
"var =       args[1]; # variable list
eid =       args[2]; # element this flux is applied to
fid =       args[3]; # face ID
els =       args[4]; # neighborhood element IDs [[left], [right]]
grid =      args[5]; # grid data struct
geo_facs =  args[6]; # geometric factors
fv_data =   args[7]; # FV specific info
refel =     args[8]; # reference element
time =      args[9]; # time
dt =        args[10]; # dt
"
    end
    return code;
end

# Allocate, compute, or fetch all needed values
function prepare_needed_values_fv_julia(entities, var, lorr, vors)
    # If there are no entities, do nothing
    if length(entities) < 1
        return "";
    end
    
    # Options
    average_coefficients = false; # true to integrate and average coefficients over cells and faces.
    
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
    fv_order = fv_info.fluxOrder;
    
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
        cname = make_entity_name(entities[i]);
        if is_unknown_var(entities[i], var) && lorr == LHS
            # The unknown variable will not be included as it is being solved for,
            # but it should be checked to make sure it is just a linear dependence. TODO
            
        else
            # Is coefficient(number or function) or variable(array)?
            (ctype, cval) = get_coef_val(entities[i]);
            if ctype == -1
                # It was a special symbol like dt or FACENORMAL
                if entities[i].name == "FACENORMAL1"
                    code *= cname * " = normal["*string(entities[i].index)*"]; # normal vector component\n"
                    piece_needed[8] = true;
                elseif entities[i].name == "FACENORMAL2"
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
                else
                    nodesymbol = "nodex"
                end
                coef_index = get_coef_index(entities[i]);
                
                if vors == "volume"
                    if average_coefficients
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
                        piece_needed[2] = true;
                    else # just evaluate coefficient at center
                        code *= cname * " = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", fv_data.cellCenters[:,eid], time, eid, fid);\n";
                    end
                    
                    
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
                        # integrate over face
                        if config.dimension == 1
                            # in 1d there is only one face node
                            code *= cname * " = " * cname * "[1]\n";
                        else
                            code *= cname * " = (refel.surf_wg[frefelind[1]] .* face_detJ)' * refel.surf_Q[frefelind[1]][:, refel.face2local[frefelind[1]]] * " * cname * " / area; # integrate over face\n";
                            piece_needed[9] = true;
                            piece_needed[10] = true;
                        end
                        need_deriv_matrix = true;
                        piece_needed[[4,5,8]] .= true;
                        
                    else # no derivatives, only need surface nodes
                        if average_coefficients
                            code *= cname * " = zeros(refel.Nfp[frefelind[1]]);\n";
                            # code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargs * " end\n";
                            code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", "*nodesymbol*"[:,coefi], time, eid, fid) end\n";
                            piece_needed[5] = true;
                            piece_needed[6] = true;
                            
                            # integrate over face
                            if config.dimension == 1
                                # in 1d there is only one face node
                                code *= cname * " = " * cname * "[1]\n";
                            else
                                code *= cname * " = (refel.surf_wg[frefelind[1]] .* face_detJ)' * refel.surf_Q[frefelind[1]][:, refel.face2local[frefelind[1]]] * " * cname * " / area; # integrate over face\n";
                                piece_needed[9] = true;
                                piece_needed[10] = true;
                            end
                        else # just evaluate coefficient at center
                            code *= cname * " = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", fv_data.faceCenters[:,fid], time, eid, fid);\n";
                        end
                    end
                end
                
            elseif ctype == 3 # a known variable value
                piece_needed[1] = true;
                # This generates something like: coef_u_1 = copy((Finch.variables[1]).values[1, gbl])
                # cellside = 0; # 0 means no side flag
                # for flagi=1:length(entities[i].flags)
                #     if occursin("DGSIDE1", entities[i].flags[flagi]) || occursin("CELL1", entities[i].flags[flagi])
                #         cellside = 1;
                #     elseif occursin("DGSIDE2", entities[i].flags[flagi]) || occursin("CELL2", entities[i].flags[flagi])
                #         cellside = 2;
                #     elseif occursin("CENTRAL", entities[i].flags[flagi])
                #         cellside = 3;
                #     elseif occursin("NEIGHBORHOOD", entities[i].flags[flagi])
                #         cellside = 4;
                #     end
                # end
                cellside = get_face_side_info(entities[i]); # in code_generator_utils.jl
                if vors == "surface"
                    l2gsymbol = "els1"
                else
                    l2gsymbol = "el"
                end
                if cellside == 1
                    l2gsymbol = "els1"
                elseif cellside == 2
                    l2gsymbol = "els2"
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
                        if variables[cval].discretization == FV
                            code *= cname * " = Finch.variables["*string(cval)*"].values["*indstr*", "*l2gsymbol*"];\n";
                        else
                            code *= cname * " = sum(Finch.variables["*string(cval)*"].values["*indstr*", Finch.grid_data.loc2glb[:,eid]]) / size(Finch.grid_data.loc2glb,1);\n";
                        end
                    end
                    
                else # surface
                    if variables[cval].discretization == FV
                        # Need to reconstruct it using neighboring cells
                        if fv_order == 1
                            if length(entities[i].derivs) > 0
                                code *= cname * " = Finch.variables["*string(cval)*"].values["*indstr*", els2] - Finch.variables["*string(cval)*"].values["*indstr*", els1];\n";
                                code *= cname * " = (els1 != els2 && abs(normal["*string(entities[i].derivs[1])*"]) > 1e-10) ? "*cname*" ./ dxyz["*string(entities[i].derivs[1])*"]  : 0\n"
                                need_deriv_dist = true;
                                piece_needed[4] = true; # cellx
                                piece_needed[8] = true; # normal
                            else
                                if cellside == 0
                                    # No side was specified, so use the average
                                    code *= cname * " = 0.5 * (Finch.variables["*string(cval)*"].values["*indstr*", els1] + Finch.variables["*string(cval)*"].values["*indstr*", els2]);\n";
                                elseif cellside < 3
                                    code *= cname * " = Finch.variables["*string(cval)*"].values["*indstr*", "*l2gsymbol*"];\n";
                                elseif cellside == 3 # central
                                    # same as 0 for this case
                                    code *= cname * " = 0.5 * (Finch.variables["*string(cval)*"].values["*indstr*", els1] + Finch.variables["*string(cval)*"].values["*indstr*", els2]);\n";
                                elseif cellside == 4 # neighborhood
                                    # This is a special case only for callback functions.
                                    # Rather than representing a value, it constructs a Neighborhood object to be passed.
                                    code *= cname * " = Finch.Neighborhood(els, cellx, [Finch.variables["*string(cval)*"].values["*indstr*", els[1]], Finch.variables["*string(cval)*"].values["*indstr*", els[2]]]);\n";
                                    piece_needed[4] = true; # cellx
                                end
                            end
                            
                        else # order > 1
                            # collect the cell info
                            piece_needed[6] = true; # face center coords
                            piece_needed[4] = true; # cell centers
                            if length(entities[i].derivs) > 0
                                code *= cname * " = Finch.variables["*string(cval)*"].values["*indstr*", els[2][1]] - Finch.variables["*string(cval)*"].values["*indstr*", els[1][1]];\n";
                                code *= cname * " = (els[1][1] != els[2][1] && abs(normal["*string(entities[i].derivs[1])*"]) > 1e-10) ? "*cname*" ./ dxyz["*string(entities[i].derivs[1])*"]  : 0\n"
                                need_deriv_dist = true;
                                piece_needed[4] = true; # cellx
                                piece_needed[8] = true; # normal
                            else
                                if cellside == 0 || cellside == 3
                                    # Use the full neighborhood
                                    code *= cname * " = Finch.FV_reconstruct_value([cellx[1] cellx[2]], [Finch.variables["*string(cval)*"].values["*indstr*", els[1]] ; Finch.variables["*string(cval)*"].values["*indstr*", els[2]]], facex);\n";
                                elseif cellside == 1 # left
                                    code *= cname * " = Finch.FV_reconstruct_value_left_right(cellx[1], cellx[2], Finch.variables["*string(cval)*"].values["*indstr*", els[1]], Finch.variables["*string(cval)*"].values["*indstr*", els[2]], facex, limiter=\"vanleer\")[1];\n";
                                elseif cellside == 2 # right
                                    code *= cname * " = Finch.FV_reconstruct_value_left_right(cellx[1], cellx[2], Finch.variables["*string(cval)*"].values["*indstr*", els[1]], Finch.variables["*string(cval)*"].values["*indstr*", els[2]], facex, limiter=\"vanleer\")[2];\n";
                                elseif cellside == 2 # right
                                    code *= cname * " = Finch.Neighborhood(els, cellx, [Finch.variables["*string(cval)*"].values["*indstr*", els[1]], Finch.variables["*string(cval)*"].values["*indstr*", els[2]]]);\n";
                                    piece_needed[4] = true; # cellx
                                end
                            end
                        end
                    else # A nodal variable that has a value specified at the surface
                        if length(entities[i].derivs) > 0
                            code *= cname * "_TMP1 = sum(Finch.variables["*string(cval)*"].values["*indstr*", Finch.grid_data.loc2glb[:,els1]]) / size(Finch.grid_data.loc2glb,1);\n";
                            code *= cname * "_TMP2 = sum(Finch.variables["*string(cval)*"].values["*indstr*", Finch.grid_data.loc2glb[:,els2]]) / size(Finch.grid_data.loc2glb,1);\n";
                            code *= cname * " = "*cname*"_TMP2 - "*cname*"_TMP1;\n";
                            code *= cname * " = (els1 != els2 && abs(normal["*string(entities[i].derivs[1])*"]) > 1e-10) ? "*cname*" ./ dxyz["*string(entities[i].derivs[1])*"]  : 0\n"
                            need_deriv_dist = true;
                            piece_needed[4] = true; # cellx
                            piece_needed[8] = true; # normal
                        else
                            if cellside < 4
                                code *= cname * " = sum(Finch.variables["*string(cval)*"].values["*indstr*", Finch.grid_data.face2glb[:,1,fid]]) / size(Finch.grid_data.face2glb, 1);\n";
                            elseif cellside == 4 # neighborhood
                                # This is a special case only for callback functions.
                                # Rather than representing a value, it constructs a Neighborhood object to be passed.
                                code *= cname * " = sum(Finch.variables["*string(cval)*"].values["*indstr*", Finch.grid_data.face2glb[:,1,fid]]) / size(Finch.grid_data.face2glb, 1);\n";
                                code *= cname * " = Finch.Neighborhood(els, cellx, ["*cname*", "*cname*"]);\n";
                                piece_needed[4] = true; # cellx
                            end
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
                else
                    nodesymbol = "nodex"
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
                    if average_coefficients
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
                        piece_needed[2] = true;
                    else # just evaluate coefficient at center
                        code *= cname * " = evaluate_coefficient(Finch.coefficients["*string(cval)*"], ["*indstr*"], fv_data.cellCenters[:,eid], time, eid, fid);\n";
                    end
                    
                    
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
                        # integrate over face
                        if config.dimension == 1
                            # in 1d there is only one face node
                            code *= cname * " = " * cname * "[1]\n";
                        else
                            code *= cname * " = (refel.surf_wg[frefelind[1]] .* face_detJ)' * refel.surf_Q[frefelind[1]][:, refel.face2local[frefelind[1]]] * " * cname * " / area; # integrate over face\n";
                            piece_needed[9] = true;
                            piece_needed[10] = true;
                        end
                        need_deriv_matrix = true;
                        piece_needed[[4,5,8]] .= true;
                        
                    else # no derivatives, only need surface nodes
                        if average_coefficients
                            code *= cname * " = zeros(refel.Nfp[frefelind[1]]);\n";
                            # code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = (Finch.coefficients["*string(cval)*"]).value["*indstr*"].func" * cargs * " end\n";
                            code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(cval)*"], ["*indstr*"], "*nodesymbol*"[:,coefi], time, eid, fid) end\n";
                            # integrate over face
                            if config.dimension == 1
                                # in 1d there is only one face node
                                code *= cname * " = " * cname * "[1]\n";
                            else
                                code *= cname * " = (refel.surf_wg[frefelind[1]] .* face_detJ)' * refel.surf_Q[frefelind[1]][:, refel.face2local[frefelind[1]]] * " * cname * " / area; # integrate over face\n";
                                piece_needed[9] = true;
                                piece_needed[10] = true;
                            end
                            piece_needed[5] = true;
                            piece_needed[6] = true;
                        else # just evaluate coefficient at center
                            code *= cname * " = evaluate_coefficient(Finch.coefficients["*string(cval)*"], ["*indstr*"], fv_data.faceCenters[:,fid], time, eid, fid);\n";
                        end
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
            needed_pieces *= "detj = geo_facs.detJ[eid]; # geometric factors\n"
        end
        if piece_needed[5]
            needed_pieces *= "J = geo_facs.J[eid]; # geometric factors\n"
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
"if length(els[2]) == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
    els[2] = [eid];
end
";
        if fv_order == 1
            needed_pieces *= "els1 = els[1][1]; els2 = els[2][1]; # For first order only.\n";
        end
        # els, nodex, loc2glb, cellx, frefelind, facex, face2glb, normal, face_detJ, area, vol_J
        if piece_needed[5]
            needed_pieces *= "frefelind = [grid.faceRefelInd[1,fid], grid.faceRefelInd[2,fid]]; # refel based index of face in both elements\n"
        end
        if piece_needed[8] || true # probably always need
            needed_pieces *= "normal = grid.facenormals[:, fid];\n"
        end
        if piece_needed[2] || piece_needed[3]
            needed_pieces *= "loc2glb = [grid.loc2glb[:,els[1][1]], grid.loc2glb[:, els[2][1]]]; # volume local to global\n"
        end
        if piece_needed[2]
            needed_pieces *= "nodex = [grid.allnodes[:,loc2glb[1][:]], grid.allnodes[:,loc2glb[2][:]]]; # volume node coordinates\n"
        end
        if piece_needed[4]
            needed_pieces *= "cellx = [fv_data.cellCenters[:, els[1]], fv_data.cellCenters[:, els[2]]]; # cell center coordinates\n"
        end
        if piece_needed[5] || piece_needed[6]
            needed_pieces *= "face2glb = grid.face2glb[:,:,fid];         # global index for face nodes for each side of each face\n"
        end
        if piece_needed[6]
            # needed_pieces *= "facex = fv_data.faceCenters[:, fid];  # face node coordinates\n"
            needed_pieces *= "facex = grid.allnodes[:, grid.face2glb[:,1,fid]];  # face node coordinates\n"
        end
        if piece_needed[9]
            needed_pieces *= "face_detJ = geo_facs.face_detJ[fid];       # detJ on face\n"
        end
        if piece_needed[10]
            needed_pieces *= "area = geo_facs.area[fid];                 # area of face\n"
        end
        if piece_needed[11]
            needed_pieces *= "vol_J = (geo_facs.J[els[1][1]], geo_facs.J[els[2][1]]);\n"
        end
        
        if need_deriv_matrix
            needed_pieces *= "TODO: derivative matrices for surface. generate_code_layer_fv_julia.jl - prepare_needed_values_fv_julia()"
        end
        if need_deriv_dist
            needed_pieces *= "dxyz = norm(cellx[2][:,1] - cellx[1][:,1]) .* normal; # normal scaled by distance between cell centers\n"
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
    how_many_terms = 0;
    for i=1:length(terms)
        how_many_terms = how_many_terms + length(terms[i]);
    end
    if how_many_terms < 1
        if lorr == LHS && vors == "surface"
            code = "return [0   0] # There were no terms to compute here";
        else
            code = "return [0] # There were no terms to compute here";
        end
        return code;
    end
    
    # Allocate the vector or matrix to be returned if needed
    if dofsper > 1
        if lorr == RHS
            code *= "result = zeros("*string(dofsper)*"); # Allocate for returned values.\n"
        else
            if vors == "volume"
                code *= "result = zeros("*string(dofsper)*", "*string(dofsper)*"); # Allocate for returned values.\n"
            else
                code *= "result = zeros(2 * "*string(dofsper)*", "*string(dofsper)*"); # Allocate for returned values.\n"
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
                        (term_result, var_ind) = generate_term_calculation_fv_julia(terms[vi][ci][i], var, lorr, vors);
                        
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
                    (term_result, var_ind) = generate_term_calculation_fv_julia(terms[ci][i], var, lorr, vors);
                    
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
        
        # Put the subvector together into result or cell_matrix
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
                    code *= "result["*string(emi)*"] = " * subvector[emi] * "\n";
                end
            end
            if vors == "surface"
                code *= "
if els[2][1] == eid && els[1][1] != els[2][1]
    result = -result; # Since this flux is applied to element eid, make sure it's going in the right direction.
end
return result;
"
            end
            
        end
        
        
    else # one dof
        terms = terms[1];
        if lorr == LHS
            if vors == "volume"
                result = "0";
            else
                result1 = "0";
                result2 = "0";
            end
            
        else
            result = "0";
        end
        
        #process each term
        first_term1 = true;
        first_term2 = true;
        for i=1:length(terms)
            (term_result1, term_result2, var_ind) = generate_term_calculation_fv_julia(terms[i], var, lorr, vors);
            if lorr == LHS && vors == "surface"
                # There are two results, one for each side of the face
                if !(term_result1 == "0" || term_result1 == "")
                    if first_term1
                        result1 = term_result1;
                        first_term1 = false;
                    else
                        result1 *= " .+ " * term_result1;
                    end
                end
                if !(term_result2 == "0" || term_result2 == "")
                    if first_term2
                        result2 = term_result2;
                        first_term2 = false;
                    else
                        result2 *= " .+ " * term_result2;
                    end
                end
            else # only term_result1 is used
                if !(term_result1 == "0" || term_result1 == "")
                    if first_term1
                        result = term_result1;
                        first_term1 = false;
                    else
                        result *= " .+ " * term_result1;
                    end
                end
            end
            
        end # term loop
        
        if lorr == LHS && vors == "surface" # LHS surface returns two pieces per dof
            code *= "result = [" * result1 * "    " * result2 * "];\n";
            
        elseif vors == "surface" # RHS surface
            code *= "result = [" * result * "];\n";
            code *= "
if els[2][1] == eid && els[1][1] != els[2][1]
    result = -result; # Since this flux is applied to element eid, make sure it's going in the right direction.
end
return result;
"
        else # volume (LHS and RHS)
            code *= "result = [" * result * "];\n";
        end
    end
    
    return code;
end

function generate_term_calculation_fv_julia(term, var, lorr, vors)
    # This will return 3 things: for volume terms, the first is the result.
    # For surface terms, there is one for each side.
    # The third is the var_ind
    result1 = "";
    result2 = "";
    # Note: separate_factors return test and trial info, but FV will not have any test functions.
    # var_part refers to the unknown variable part.
    # example:
    # 0.1*D1__u_1*_FACENORMAL1_1  ->  0.1 * normal[1]                    on LHS
    #                             ->  0.1 * (coef_D1xu_1 .* normal[1])   on RHS
    if lorr == LHS
        (test_part, var_part, coef_part, test_ind, var_ind) = separate_factors(term, var);
        which_side = get_face_side_info(var_part); # 0=none, 1/2=side 1/2, 3=average, 4=neighborhood
        
        if coef_part === nothing # Only an unknown variable
            # handle negative sign
            if typeof(var_part) == Expr && (var_part.args[1] === :.- || var_part.args[1] === :-)
                if vors == "surface"
                    if which_side == 0 || which_side == 3 # average them
                        result1 = "-0.5";
                        result2 = "-0.5";
                    elseif which_side == 1
                        result1 = "-1";
                    elseif which_side == 2
                        result2 = "-1";
                    else
                        # neighborhoods are not supported for implicit steppers
                    end
                else
                    result1 = "-1";
                end
                
            else
                if vors == "surface"
                    if which_side == 0 || which_side == 3 # average them
                        result1 = "0.5";
                        result2 = "0.5";
                    elseif which_side == 1
                        result1 = "1";
                    elseif which_side == 2
                        result2 = "1";
                    else
                        # neighborhoods are not supported for implicit steppers
                    end
                else
                    result1 = "1";
                end
            end
            
        else
            if vors == "surface"
                # This is tricky. First make two copies of the coef_part.
                # Since the coef_part could actually contain a variable hidden in an expression like a conditional,
                # Set the correct side to 1 and the other side to 0, or if no side, set both to 0.5
                # Also make side 2 negative... because I need to make upwinding work.
                # This is a bad way to rig this up, but I'm out of ideas here.
                coef_side1 = copy(coef_part);
                coef_side2 = copy(coef_part);
                coef_side1 = replace_lhs_surface_var_entities(coef_side1, var, 1);
                coef_side2 = replace_lhs_surface_var_entities(coef_side2, var, 2);
                if which_side == 0 || which_side == 3 # average them
                    result1 = string(replace_entities_with_symbols(coef_side1));
                    result2 = string(replace_entities_with_symbols(coef_side2));
                elseif which_side == 1
                    result1 = string(replace_entities_with_symbols(coef_side1));
                elseif which_side == 2
                    result2 = string(replace_entities_with_symbols(coef_side2));
                else
                    # neighborhoods are not supported for implicit steppers
                end
            else
                result1 = string(replace_entities_with_symbols(coef_part));
            end
        end
        
    else
        (test_part, var_part, coef_part, test_ind, var_ind) = separate_factors(term);
        # RHS: coef_part
        if !(coef_part === nothing)
            result1 = string(replace_entities_with_symbols(coef_part));
        else
            result1 = "0";
        end
    end
    
    return (result1, result2, var_ind);
end

# A special function for handling LHS surface parts.
# Swaps variable entities with a number depending on side.
# negate will make side 2 negative.
function replace_lhs_surface_var_entities(ex, var, side, negate=true)
    if typeof(ex) == Expr
        for i=1:length(ex.args)
            ex.args[i] = replace_lhs_surface_var_entities(ex.args[i], var, side, negate);
        end
    elseif typeof(ex) == SymEntity
        if is_unknown_var(ex, var)
            which_side = get_face_side_info(ex);
            if which_side == 0 || which_side == 3
                ex = 0.5;
            elseif which_side == side
                ex = 1;
            else
                ex = 0;
            end
            if side == 2 && negate
                ex = -ex;
            end
        end
    else
        # number?
    end
    return ex;
end

# Generate the assembly loop structures and insert the content
function generate_assembly_loop_fv_julia(var, indices, parallel_type)
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
        parallel_type = ["none", parallel_type];
    end
    
    code = 
"# Label things that were allocated externally
if is_explicit
    sourcevec = allocated_vecs[1];
    fluxvec = allocated_vecs[2];
    facefluxvec = allocated_vecs[3];
    face_done = allocated_vecs[4];
else
    lhsmatV = allocated_vecs[3];
    sourcevec = allocated_vecs[4]; # In the implicit case, these vectors are for the RHS
    fluxvec = allocated_vecs[5];
    facefluxvec = allocated_vecs[6];
    face_done = allocated_vecs[7];
    
    for i=1:length(lhsmatV)
        lhsmatV[i] = 0;
    end
end

nel = fv_grid.nel_owned;
dofs_squared = dofs_per_node*dofs_per_node;
"
    # generate the loop structures
    code *= "# Loops\n"
    loop_start = ""
    loop_end = ""
    for i=1:length(indices)
        if indices[i] == "elements" || indices[i] == "cells"
            loop_start *= "for eid in elemental_order\n";
            loop_start *= "
for i=1:length(face_done)
    face_done[i] = 0; # Reset face_done in elemental loop
end\n";
            loop_end *= "end # loop for elements\n";
        else
            indexer_index = 0;
            for j=1:length(indexers)
                if indices[i].symbol === indexers[j].symbol
                    indexer_index += j;
                    break;
                end
            end
            loop_start *= "for indexing_variable_"*string(indices[i].symbol)*" in Finch.indexers["*string(indexer_index)*"].range\n";
            loop_start *= "    Finch.indexers["*string(indexer_index)*"].value = indexing_variable_"*string(indices[i].symbol)*";\n";
            loop_end *= "end # loop for "*string(indices[i].symbol)*"\n";
        end
    end
    
    # The index offset depends on the order in which indexers were used during variable creation.
    # They must all be the same for now.
    ind_offset = ""
    prev_ind = ""
    if typeof(var) <: Array
        var_indices = var[1].indexer;
    else
        var_indices = var.indexer;
    end
    if !(typeof(var_indices) <: Array)
        var_indices = [var_indices];
    end
    for i=1:length(var_indices)
        if length(prev_ind) == 0
            ind_offset *= "(indexing_variable_"*string(var_indices[i].symbol);
            prev_ind = string(length(var_indices[i].range));
        else
            ind_offset *= " + "*prev_ind*" * (indexing_variable_"*string(var_indices[i].symbol)*" - 1";
            prev_ind = string(length(var_indices[i].range));
        end
    end
    for i=1:length(var_indices)
        ind_offset *= ")";
    end
    
    # index_offset = "dofs_per_loop * ("*ind_offset*" - 1 + dofs_per_node * (eid - 1)) + 1"
    # face_index_offset = "dofs_per_loop * ("*ind_offset*" - 1 + dofs_per_node * (fid - 1)) + 1"
    index_offset = "dofs_per_node * (eid - 1) + ("*ind_offset*" - 1) * dofs_per_loop + 1";
    face_index_offset = "dofs_per_node * (fid - 1) + ("*ind_offset*" - 1) * dofs_per_loop + 1";
    
    code *= loop_start;
    
    # insert the content
    code *= "
    # determine index in global vector
    index_offset = "*index_offset*";
    first_ind = index_offset; # First global index for this loop
    last_ind = index_offset + dofs_per_loop-1; # last global index
    
    # Zero the result vectors for this element
    sourcevec[first_ind:last_ind] .= 0;
    fluxvec[first_ind:last_ind] .= 0;
    
    ##### Source integrated over the cell #####
    # Compute RHS volume integral
    if !(source_rhs === nothing)
        
        sourceargs = (var, eid, 0, grid_data, geo_factors, fv_info, refel, t, dt);
        source = source_rhs.func(sourceargs; "*index_args*");
        
        # Add to global source vector
        sourcevec[first_ind:last_ind] .= source;
    end
    
    # Compute LHS volume integral = (block)diagonal components
    if !(source_lhs === nothing) && !is_explicit
        sourceargs = (var, eid, 0, fv_grid, fv_geo_factors, fv_info, refel, t, dt);
        source = source_lhs.func(sourceargs; "*index_args*");
        # Add to global matrix
        for di = 1:dofs_per_loop
            first = (eid-1)*dofs_squared + "*ind_offset*"; # This needs to change for multi dof TODO
            last = first + dofs_per_loop - 1;
            lhsmatV[first:last] = source[((di-1)*dofs_per_loop + 1):(di*dofs_per_loop)];
        end
    end

    ##### Flux integrated over the faces #####
    # Loop over this element's faces.
    for i=1:refel.Nfaces
        fid = grid_data.element2face[i, eid];
        do_face_here = (face_done[fid] < dofs_per_node);
        if do_face_here
            face_done[fid] += dofs_per_loop;
        end
        # Only one element on either side is available here. For more use parent/child version.
        (leftel, rightel) = grid_data.face2element[:,fid];
        if rightel == 0
            neighborhood = [[leftel],[]];
        else
            neighborhood = [[leftel],[rightel]];
        end
        
        # determine index in global face vector
        face_index_offset = "*face_index_offset*";
        
        if !(flux_rhs === nothing)
            if do_face_here
                
                fluxargs = (var, eid, fid, neighborhood, grid_data, geo_factors, fv_info, refel, t, dt);
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
        
        # LHS surface integral
        if !(flux_lhs === nothing) && !is_explicit
            # These update the (block)diagonal for eid and one (block of)off diagonal for neighborID
            if do_face_here
                face_done[fid] = 1; # Possible race condition.
                neighborID = (rightel==eid ? leftel : rightel);
                eid_flux_index = (leftel==eid ? 1 : 2);
                neighbor_flux_index = (leftel==eid ? 2 : 1);
                nfirst_ind = (neighborID-1)*dofs_per_node + "*ind_offset*";
                nlast_ind = neighborID*dofs_per_node + "*ind_offset*" - 1;
                
                fluxargs = (var, eid, fid, neighborhood, fv_grid, fv_geo_factors, fv_info, refel, t, dt);
                flux = flux_lhs.func(fluxargs; "*index_args*"); println(flux);
                flux = flux_lhs.func(fluxargs; "*index_args*") .* fv_geo_factors.area[fid];
                # insert the matrix elements
                for di = 1:dofs_per_loop
                    # The eid components
                    first = nel * dofs_squared + (fid-1)*dofs_squared*4 + "*ind_offset*";
                    last = first + dofs_per_loop - 1;
                    lhsmatV[first:last] = flux[((di-1)*dofs_per_loop + 1):(di*dofs_per_loop), eid_flux_index] ./ fv_geo_factors.volume[eid];
                    
                    # The neighborID components
                    if neighborID > 0
                        # first = nel * dofs_squared + (fid-1)*dofs_squared*4 + dofs_squared + "*ind_offset*";
                        # last = first + dofs_per_loop - 1;
                        first += dofs_squared;
                        last += dofs_squared;
                        lhsmatV[first:last] = flux[((di-1)*dofs_per_loop + 1):(di*dofs_per_loop), neighbor_flux_index] ./ fv_geo_factors.volume[eid];
                    end
                end
                
                # While we're here, set the components for the neighborID rows as well.
                if neighborID > 0
                    for di = 1:dofs_per_loop
                        # The eid components
                        first = nel * dofs_squared + (fid-1)*dofs_squared*4 + dofs_squared*2 + "*ind_offset*";
                        last = first + dofs_per_loop - 1;
                        lhsmatV[first:last] = -flux[((di-1)*dofs_per_loop + 1):(di*dofs_per_loop), eid_flux_index] ./ fv_geo_factors.volume[neighborID];
                        
                        # The neighborID components
                        first += dofs_squared;
                        last += dofs_squared;
                        
                        lhsmatV[first:last] = -flux[((di-1)*dofs_per_loop + 1):(di*dofs_per_loop), neighbor_flux_index] ./ fv_geo_factors.volume[neighborID];
                    end
                end
                
            else
                # This flux has either been computed or is being computed by another thread.
                # Everything was set there.
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
                        
                        bflux = FVSolver.FV_flux_bc_rhs_only_simple(prob.bc_func[var[vi].index, fbid][compo], fid, t) .* geo_factors.area[fid];
                        
                        fluxvec[index_offset + dofind-1] += (bflux - facefluxvec[face_index_offset + dofind-1]) ./ geo_factors.volume[eid];
                        facefluxvec[face_index_offset + dofind-1] = bflux;
                    elseif prob.bc_type[var[vi].index, fbid] == DIRICHLET
                        # Set variable array and handle after the face loop
                        var[vi].values[compo,eid] = FVSolver.evaluate_bc(prob.bc_func[var[vi].index, fbid][compo], eid, fid, t);
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
                    
                    bflux = FVSolver.FV_flux_bc_rhs_only_simple(prob.bc_func[var.index, fbid][d], fid, t) .* geo_factors.area[fid];
                    
                    
                    fluxvec[index_offset] += (bflux - facefluxvec[face_index_offset]) ./ geo_factors.volume[eid];
                    facefluxvec[face_index_offset] = bflux;
                elseif prob.bc_type[var.index, fbid] == DIRICHLET
                    # Set variable array and handle after the face loop
                    var.values[d,eid] = FVSolver.evaluate_bc(prob.bc_func[var.index, fbid][d], eid, fid, t);
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