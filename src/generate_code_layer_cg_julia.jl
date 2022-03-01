#=
Code generation functions for Julia.
=#

# Extract the input args. This must match the arguments passed by the solver.
function handle_input_args_cg_julia(lorr, vors)
    code =
"var =       args[1];
eid =       args[2];
fid =       args[3];
grid =      args[4];
geo_facs =  args[5];
refel =     args[6]
time =      args[7];
dt =        args[8]
"
    # A trick for uniform grids to avoid repeated work
    if vors == "volume" && config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
        code *= "stiffness = args[9]; # set of stiffness matrices for each dimension\n";
        code *= "mass =      args[10]; # mass matrix\n";
    end
    
    return code;
end

# Allocate, compute, or fetch all needed values
function prepare_needed_values_cg_julia(entities, var, lorr, vors)
    # Only gather the pieces that are used. See the end of this function.
    if vors == "volume"
        # el, loc2glb, nodex, detj, J
        piece_needed = zeros(Bool, 5);
    else # surface
        # (els, loc2glb, nodex, frefelind, face2glb, facex, normal, fdetJ, vol_J)
        piece_needed = zeros(Bool, 9);
    end
    need_deriv_matrix = false;
    
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
        if is_test_function(entities[i])
            # Assign it a transpose quadrature matrix
            if vors == "volume"
                if length(entities[i].derivs) > 0
                    xyzchar = ["x","y","z"];
                    for di=1:length(entities[i].derivs)
                        code *= cname * " = TRQ"*string(entities[i].derivs[di])*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of test function\n";
                    end
                    need_deriv_matrix = true;
                else
                    code *= cname * " = refel.Q'; # test function.\n";
                end
            else
                if length(entities[i].derivs) > 0
                    xyzchar = ["x","y","z"];
                    for di=1:length(entities[i].derivs)
                        code *= cname * " = TRQ"*string(entities[i].derivs[di])*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of test function\n";
                        need_deriv_matrix = true;
                    end
                else
                    code *= cname * " = (refel.surf_Q[frefelind[1]]); # test function.\n";
                end
            end
            
        elseif is_unknown_var(entities[i], var) && lorr == LHS
            if vors == "volume"
                if length(entities[i].derivs) > 0
                    xyzchar = ["x","y","z"];
                    for di=1:length(entities[i].derivs)
                        code *= cname * " = RQ"*string(entities[i].derivs[di])*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of trial function\n";
                    end
                    need_deriv_matrix = true;
                else
                    code *= cname * " = refel.Q; # trial function.\n";
                end
            else
                if length(entities[i].derivs) > 0
                    xyzchar = ["x","y","z"];
                    for di=1:length(entities[i].derivs)
                        code *= cname * " = RQ"*string(entities[i].derivs[di])*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of trial function\n";
                        need_deriv_matrix = true;
                    end
                else
                    code *= cname * " = refel.surf_Q[frefelind[1]]; # trial function.\n";
                end
            end
            
        else
            # Is coefficient(number or function) or variable(array)?
            (ctype, cval) = get_coef_val(entities[i]);
            if ctype == -1
                # It was a special symbol like dt
                if vors == "surface"
                    if entities[i].name == "FACENORMAL1"
                        code *= cname * " = normal["*string(entities[i].index)*"]; # normal vector component\n"
                        piece_needed[7] = true;
                    else entities[i].name == "FACENORMAL2"
                        code *= cname * " = -normal["*string(entities[i].index)*"]; # reverse normal vector component\n"
                        piece_needed[7] = true;
                    end
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
                # coef_n_i = zeros(refel.Np);
                # for coefi = 1:refel.Np
                #     coef_k_i[coefi] = (Finch.genfunctions[cval]).func(x[1,coefi], x[2,coefi],x[3,coefi],time);
                # end
                ######################################
                # Note: for LHS, this is only Nfp values. For RHS, this is Np values, but only evaluated at face
                # cargs = "(nodex[coefi], 0, 0, time)";
                # if config.dimension == 2
                #     cargs = "(nodex[1, coefi], nodex[2, coefi], 0, time)";
                # elseif config.dimension == 3
                #     cargs = "(nodex[1, coefi], nodex[2, coefi], nodex[3, coefi], time)";
                # end
                # cargs1 = "(nodex[1][coefi], 0, 0, time)";
                # if config.dimension == 2
                #     cargs1 = "(nodex[1][1, coefi], nodex[1][2, coefi], 0, time)";
                # elseif config.dimension == 3
                #     cargs1 = "(nodex[1][1, coefi], nodex[1][2, coefi], nodex[1][3, coefi], time)";
                # end
                # cargs2 = "(nodex[2][coefi], 0, 0, time)";
                # if config.dimension == 2
                #     cargs2 = "(nodex[2][1, coefi], nodex[2][2, coefi], 0, time)";
                # elseif config.dimension == 3
                #     cargs2 = "(nodex[2][1, coefi], nodex[2][2, coefi], nodex[2][3, coefi], time)";
                # end
                # cargsf = "(facex[coefi], 0, 0, time)";
                # if config.dimension == 2
                #     cargsf = "(facex[1, coefi], facex[2, coefi], 0, time)";
                # elseif config.dimension == 3
                #     cargsf = "(facex[1, coefi], facex[2, coefi], facex[3, coefi], time)";
                # end
                
                nodesymbol = "nodex";
                nodesymbol1 = "nodex[1]";
                nodesymbolf = "facex";
                coef_index = get_coef_index(entities[i]);
                
                if vors == "volume"
                    code *= cname * " = zeros(refel.Np);\n";
                    # code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargs * " end\n";
                    code *= "for coefi = 1:refel.Np " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", "*nodesymbol*"[:,coefi], time, eid, fid) end\n";
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*" and interpolate at quadrature points.\n";
                        end
                        need_deriv_matrix = true;
                    else
                        code *= cname * " = refel.Q * " * cname * "; # Interpolate at quadrature points.\n";
                    end
                    piece_needed[3] = true;
                    
                else # surface
                    # If derivatives are needed, must evaluate at all volume nodes.
                    if length(entities[i].derivs) > 0
                        code *= cname * " = zeros(refel.Np);\n";
                        # code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargs1 * " end\n";
                        code *= "for coefi = 1:refel.Np " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", "*nodesymbol1*"[:,coefi], time, eid, fid) end\n";
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RD"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                        code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]]; # extract face values only.";
                        need_deriv_matrix = true;
                        piece_needed[3] = true;
                        
                    else # no derivatives, only need surface nodes
                        code *= cname * " = zeros(Nfp);\n";
                        # code *= "for coefi = 1:Nfp " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargsf * " end\n";
                        code *= "for coefi = 1:Nfp " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", "*nodesymbolf*"[:,coefi], time, eid, fid) end\n";
                        piece_needed[6] = true;
                    end
                    
                    # Interpolate at surface quadrature points
                    code *= cname * " = refel.surf_Q[frefelind[1]][:,refel.face2local[frefelind[1]]] * " * cname * "; # Interpolate at quadrature points.\n";
                end
                
            elseif ctype == 3 # a known variable value
                # This generates something like: coef_u_1 = (Finch.variables[1]).values[1, loc2glb]
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
                    if variables[cval].discretization == FV
                        code *= cname * " = fill((Finch.variables["*string(cval)*"]).values["*indstr*", eid], length(loc2glb));\n";
                    else
                        code *= cname * " = (Finch.variables["*string(cval)*"]).values["*indstr*", loc2glb];\n";
                    end
                    
                    # Apply any needed derivative operators.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*"\n";
                        end
                        need_deriv_matrix = true;
                    else
                        code *= cname * " = refel.Q * " * cname * "; # Interpolate at quadrature points.\n";
                    end
                    piece_needed[2] = true;
                    
                else # surface
                    # If derivatives are needed, must evaluate at all volume nodes.
                    if length(entities[i].derivs) > 0
                        if variables[cval].discretization == FV
                            code *= cname * " = fill((Finch.variables["*string(cval)*"]).values["*indstr*", eid], length(loc2glb[1]));\n";
                        else
                            code *= cname * " = (Finch.variables["*string(cval)*"]).values["*indstr*", loc2glb[1]];\n";
                        end
                        
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RD"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                        #     code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]; # extract face values only.";
                        need_deriv_matrix = true;
                        piece_needed[2] = true;
                        
                    else # no derivatives
                        if variables[cval].discretization == FV
                            code *= cname * " = fill((Finch.variables["*string(cval)*"]).values["*indstr*", eid], length(loc2glb[1]));\n";
                        else
                            code *= cname * " = (Finch.variables["*string(cval)*"]).values["*indstr*", loc2glb[1]];\n";
                        end
                        
                        # piece_needed[5] = true;
                        piece_needed[2] = true;
                    end
                    
                    #     code *= cname * " = refel.surf_Q[frefelind[1]][:,refel.face2local[frefelind[1]]] * " * cname * "; # Interpolate at quadrature points.\n";
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
                    piece_needed[6] = true;
                else
                    piece_needed[3] = true;
                end
                # cargs = "("*nodesymbol*"[coefi], 0, 0, time)";
                # if config.dimension == 2
                #     cargs = "("*nodesymbol*"[1, coefi], "*nodesymbol*"[2, coefi], 0, time)";
                # elseif config.dimension == 3
                #     cargs = "("*nodesymbol*"[1, coefi], "*nodesymbol*"[2, coefi], "*nodesymbol*"[3, coefi], time)";
                # end
                
                nodesymbol = "nodex";
                nodesymbolf = "facex"
                
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
                            code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*" and interpolate at quadrature points.\n";
                        end
                        need_deriv_matrix = true;
                    else
                        code *= cname * " = refel.Q * " * cname * "; # Interpolate at quadrature points.\n";
                    end
                    
                else # surface
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        code *= cname * " = zeros(refel.Np);\n";
                        # code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.coefficients["*string(cval)*"]).value["*indstr*"].func" * cargs * " end\n";
                        code *= "for coefi = 1:refel.Np " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(cval)*"], ["*indstr*"], "*nodesymbol*"[:,coefi], time, eid, fid) end\n";
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                        code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]]; # extract face values only.";
                        need_deriv_matrix = true;
                        
                    else # no derivatives, only need surface nodes
                        code *= cname * " = zeros(refel.Nfp[frefelind[1]]);\n";
                        # code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = (Finch.coefficients["*string(cval)*"]).value["*indstr*"].func" * cargs * " end\n";
                        code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(cval)*"], ["*indstr*"], "*nodesymbolf*"[:,coefi], time, eid, fid) end\n";
                    end
                    # Interpolate at surface quadrature points
                    code *= cname * " = refel.surf_Q[frefelind[1]][:,refel.face2local[frefelind[1]]] * " * cname * "; # Interpolate at quadrature points.\n";
                end
                
            end
        end # if coefficient
        #code *= "println(\"did \"*string("*string(i)*"))\n";
    end # entity loop
    
    # Add the needed pieces
    if vors == "volume"
        # el, loc2glb, nodex, detj, J
        piece_needed[4] = true; #I'm basically always going to need this.
        needed_pieces = "";
        if piece_needed[1]
            needed_pieces *= "el = eid;\n"
        end
        if piece_needed[2] || piece_needed[3]
            needed_pieces *= "loc2glb = grid.loc2glb[:,eid];    # local to global map\n"
        end
        if piece_needed[3]
            needed_pieces *= "nodex = grid.allnodes[:,loc2glb]; # node coordinates\n"
        end
        if piece_needed[4]
            needed_pieces *= "detj = geo_factors.detJ[eid];     # geometric factors\n"
            needed_pieces *= "wdetj = refel.wg .* detj;         # quadrature weights\n"
        end
        if piece_needed[5] || need_deriv_matrix
            needed_pieces *= "J = geo_factors.J[eid];           # geometric factors\n"
        end
        
        if need_deriv_matrix
            needed_pieces *= 
"
# Note on derivative matrices:
# RQn are quadrature matrices for the derivatives of the basis functions
# with Jacobian factors. They are made like this.
# |RQ1|   | rx sx tx || Qx |
# |RQ2| = | ry sy ty || Qy |
# |RQ3|   | rz sz tz || Qz |

"
            use_full_deriv_mat = lorr==RHS;
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
        
    else # surface
        piece_needed[4] = true; #I'm basically always going to need this.
        piece_needed[8] = true; #I'm basically always going to need this.
        needed_pieces = 
"if grid.face2element[1, fid] == eid # The normal points out of e"*(piece_needed[7] ? "\n\tnormal = grid.facenormals[:, fid];" : "")*"
    neighbor = grid.face2element[2, fid];
    frefelind = [grid.faceRefelInd[1,fid], grid.faceRefelInd[2,fid]]; # refel based index of face in both elements
else # The normal points into e"*(piece_needed[7] ? "\n\tnormal = -grid.facenormals[:, fid];" : "")*"
    neighbor = grid.face2element[1, fid];
    frefelind = [grid.faceRefelInd[2,fid], grid.faceRefelInd[1,fid]]; # refel based index of face in both elements
end
if neighbor == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
    neighbor = eid;
    frefelind[2] = frefelind[1];
end
Nfp = refel.Nfp[frefelind[1]]; # number of face nodes\n
";
        # (els, loc2glb, nodex, frefelind, face2glb, facex, normal, fdetJ, vol_J)
        if piece_needed[1]
            needed_pieces *= "els = (eid, neighbor); # indices of elements on both sides\n"
        end
        if piece_needed[2] || piece_needed[3]
            needed_pieces *= "loc2glb = (grid.loc2glb[:,eid], grid.loc2glb[:, neighbor]); # volume local to global\n"
        end
        if piece_needed[3]
            needed_pieces *= "nodex = (grid.allnodes[:,loc2glb[1][:]], grid.allnodes[:,loc2glb[2][:]]); # volume node coordinates\n"
        end
        if piece_needed[5] || piece_needed[6]
            needed_pieces *= "face2glb = grid.face2glb[:,:,fid];         # global index for face nodes for each side of each face\n"
        end
        if piece_needed[6]
            needed_pieces *= "facex = grid.allnodes[:, face2glb[:, 1]];  # face node coordinates\n"
        end
        if piece_needed[8]
            needed_pieces *= "face_detJ = geo_facs.face_detJ[fid];       # detJ on face\n"
            needed_pieces *= "wdetj = refel.surf_wg[1] .* face_detJ;     # quadrature weights\n"
        end
        if piece_needed[9] || need_deriv_matrix
            needed_pieces *= "vol_J = (geo_facs.J[eid], geo_facs.J[neighbor]);\n"
        end
        
        if need_deriv_matrix
            needed_pieces *= 
"
# Note on derivative matrices:
# RQn are quadrature matrices for the derivatives of the basis functions
# with Jacobian factors. They are made like this.
# |RQ1|   | rx sx tx || Qx |
# |RQ2| = | ry sy ty || Qy |
# |RQ3|   | rz sz tz || Qz |

"
            use_full_deriv_mat = lorr==RHS;
            if config.dimension == 1
                needed_pieces *= "RQ1 = refel.surf_Qr[frefelind[1]] .* vol_J[1].rx[1]";
                needed_pieces *= "TRQ1 = RQ1'";
                
            elseif config.dimension == 2
                needed_pieces *= "(RQ1,RQ2,RD1,RD2) = build_face_deriv_matrix(refel, frefelind[1], vol_J[1], $use_full_deriv_mat)";
                needed_pieces *= "(TRQ1,TRQ1) = (RQ1',RQ1')";
                
            elseif config.dimension == 3
                needed_pieces *= "(RQ1,RQ2,RQ3,RD1,RD2,RD3) = build_face_deriv_matrix(refel, frefelind[1], vol_J[1], $use_full_deriv_mat)";
                needed_pieces *= "(TRQ1,TRQ2,TRQ3) = (RQ1',RQ2',RQ3')";
            end
        end
    end
    
    code = needed_pieces * "\n" * code;
    
    return code;
end

function make_elemental_computation_cg_julia(terms, var, dofsper, offset_ind, lorr, vors)
    # Here is where I make some assumption about the form of the expression.
    # Since it was expanded by the parser it should look like a series of terms: t1 + t2 + t3...
    # Where each term is multiplied by one test function component, and if LHS, involves one unknown component.
    # The submatrix modified by a term is determined by these, so go through the terms and divide them
    # into their submatrix expressions. 
    # Each term will look something like 
    # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
    # RHS: test_part * (weight_part .* coef_part)
    code = "";
    
    # Allocate the vector or matrix to be returned if needed
    if dofsper > 1
        if lorr == RHS
            code *= "element_vector = zeros(refel.Np * "*string(dofsper)*"); # Allocate the returned vector.\n"
        else
            code *= "element_matrix = zeros(refel.Np * "*string(dofsper)*", refel.Np * "*string(dofsper)*"); # Allocate the returned matrix.\n"
        end
    end
    
    # Separate the factors of each term into test, trial, coef and form the calculation
    if dofsper > 1
        # Submatrices or subvectors for each component
        if lorr == LHS
            submatrices = Array{String, 2}(undef, dofsper, dofsper);
        else # RHS
            submatrices = Array{String, 1}(undef, dofsper);
        end
        for smi=1:length(submatrices)
            submatrices[smi] = "";
        end
        
        if typeof(var) <: Array
            for vi=1:length(var) # variables
                # Process the terms for this variable
                for ci=1:length(terms[vi]) # components
                    for i=1:length(terms[vi][ci])
                        (term_result, test_ind, trial_ind) = generate_term_calculation_cg_julia(terms[vi][ci][i], var, lorr, vors);
                        
                        # println(terms)
                        # println(terms[vi])
                        # println(terms[vi][ci])
                        # println(terms[vi][ci][i])
                        # println(term_result * " : "*string(test_ind)*", "*string(trial_ind))
                        
                        # Find the appropriate submatrix for this term
                        submati = offset_ind[vi] + test_ind;
                        submatj = trial_ind;
                        if lorr == LHS
                            submat_ind = submati + dofsper * (submatj-1);
                        else
                            submat_ind = submati;
                        end
                        
                        
                        if length(submatrices[submat_ind]) > 1
                            submatrices[submat_ind] *= " .+ " * term_result;
                        else
                            submatrices[submat_ind] = term_result;
                        end
                    end
                end
                
            end # vi
            
        else # only one variable
            # Process the terms for this variable
            for ci=1:length(terms) # components
                for i=1:length(terms[ci])
                    (term_result, test_ind, trial_ind) = generate_term_calculation_cg_julia(terms[ci][i], var, lorr, vors);
                    
                    # Find the appropriate submatrix for this term
                    if lorr == LHS
                        submat_ind = test_ind + dofsper * (trial_ind-1);
                    else
                        submat_ind = test_ind;
                    end
                    
                    if length(submatrices[submat_ind]) > 1
                        submatrices[submat_ind] *= " .+ " * term_result;
                    else
                        submatrices[submat_ind] = term_result;
                    end
                end
            end
            
        end
        
        # Put the submatrices together into element_matrix or element_vector
        if lorr == LHS
            for emi=1:dofsper
                for emj=1:dofsper
                    if length(submatrices[emi, emj]) > 1
                        rangei = "("*string(emi-1)*"*refel.Np + 1):("*string(emi)*"*refel.Np)";
                        rangej = "("*string(emj-1)*"*refel.Np + 1):("*string(emj)*"*refel.Np)";
                        code *= "element_matrix["*rangei*", "*rangej*"] = " * submatrices[emi,emj] * "\n";
                    end
                end
            end
            code *= "return element_matrix;\n"
            
        else # RHS
            for emi=1:dofsper
                if length(submatrices[emi]) > 1
                    rangei = "("*string(emi-1)*"*refel.Np + 1):("*string(emi)*"*refel.Np)";
                    code *= "element_vector["*rangei*"] = " * submatrices[emi] * "\n";
                end
            end
            code *= "return element_vector;\n"
        end
        
    else # one dof
        terms = terms[1];
        if lorr == LHS
            result = "zeros(refel.Np, refel.Np)";
        else
            result = "zeros(refel.Np)";
        end
        
        #process each term
        for i=1:length(terms)
            (term_result, test_ind, trial_ind) = generate_term_calculation_cg_julia(terms[i], var, lorr, vors);
            
            if i > 1
                result *= " .+ " * term_result;
            else
                result = term_result;
            end
        end
        code *= "return " * result * ";\n";
    end
    
    return code;
end

function generate_term_calculation_cg_julia(term, var, lorr, vors)
    precomputed_mass_stiffness = config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE;
    function extract_one_entity(ex)
        if typeof(ex) == Expr
            if ex.head == :call && (ex.args[1] == :- || ex.args[1] == :.-) && length(ex.args) == 2 && typeof(ex.args[2]) == SymEntity
                return (ex.args[2], true);
            else
                return nothing;
            end
        elseif typeof(ex) == SymEntity
            return (ex, false);
        end
        # return nothing
    end
    
    result = "";
    
    if lorr == LHS
        (test_part, trial_part, coef_part, test_ind, trial_ind) = separate_factors(term, var);
        # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
        if !(coef_part === nothing)
            result = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
                    string(replace_entities_with_symbols(coef_part)) * ") * " * 
                    string(replace_entities_with_symbols(trial_part));
        else # no coef_part
            using_precomputed = false;
            if precomputed_mass_stiffness
                test_ent= extract_one_entity(test_part); # defined above
                trial_ent = extract_one_entity(trial_part);
                if !(test_ent === nothing || trial_ent === nothing)
                    if length(test_ent[1].derivs) == 0 && length(trial_ent[1].derivs) == 0
                        if test_ent[2] || trial_ent[2]
                            result = "-mass";
                        else
                            result = "mass";
                        end
                        using_precomputed = true;
                    elseif length(test_ent[1].derivs) == 1 && length(trial_ent[1].derivs) == 1 && test_ent[1].derivs[1] == trial_ent[1].derivs[1]
                        if test_ent[2] || trial_ent[2]
                            result = "-stiffness["*string(test_ent[1].derivs[1])*"]";
                        else
                            result = "stiffness["*string(test_ent[1].derivs[1])*"]";
                        end
                        using_precomputed = true;
                    end
                end
            end
            if !using_precomputed
                result = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
                    string(replace_entities_with_symbols(trial_part));
            end
        end
    else
        (test_part, trial_part, coef_part, test_ind, trial_ind) = separate_factors(term);
        # RHS: test_part * (weight_part .* coef_part)
        if !(coef_part === nothing)
            result = string(replace_entities_with_symbols(test_part)) * " * (wdetj .* " * 
                    string(replace_entities_with_symbols(coef_part)) * ")";
        else
            result = string(replace_entities_with_symbols(test_part)) * " * (wdetj)";
        end
    end
    
    return (result, test_ind, trial_ind);
end

# Generate the assembly loop structures and insert the content
function generate_assembly_loop_cg_julia(var, indices)
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

"
    # generate the loop structures
    code *= "# Loops\n"
    loop_start = ""
    loop_end = ""
    for i=1:length(indices)
        if indices[i] == "elements" || indices[i] == "cells"
            loop_start *= "for eid in elemental_order\n";
            # loop_start *= "face_done .= 0; # Reset face_done in elemental loop\n";
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
    
    code *= loop_start;
    
    # insert the content
    code *= "
    loc2glb = grid_data.loc2glb[:,eid]; # global indices of this element's nodes
    
    # determine index in global vector
    dof_index = (loc2glb .- 1) .* dofs_per_node .+ "*ind_offset*";
    #println(\"dof index \"*string(dof_index))
    
    if !rhs_only
        Astart = (eid-1)*Np*Np*dofs_per_node + 1; # The segment of AI, AJ, AV for this element
        #println(\"Astart \"*string(Astart))
    end
    
    volargs = (var, eid, 0, grid_data, geo_factors, refel, t, dt, stiffness, mass);
    linchunk = linear.func(volargs; "*index_args*");  # get the elemental linear part
    if !rhs_only
        bilinchunk = bilinear.func(volargs; "*index_args*"); # the elemental bilinear part
    end
    
    if dofs_per_node == 1
        b[loc2glb] += linchunk;
        if !rhs_only
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
        b[dof_index] += linchunk;
        if !rhs_only
            #print(\"offsets \")
            for jj=1:Np
                offset = Astart + Np*Np*("*ind_offset*"-1) + (jj-1)*Np - 1;
                #print(\", \"*string(offset))
                for ii=1:Np
                    AI[offset + ii] = dof_index[ii];
                    AJ[offset + ii] = dof_index[jj];
                    AV[offset + ii] = bilinchunk[ii,jj];
                end
            end
            #println()
        end
    end
"
    # close loops
    code *= loop_end;
    
    # finish
    code *=
"
loop_time = Base.Libc.time() - loop_time;

if !rhs_only
    # Build the sparse A. Uses default + to combine overlaps
    #println(\"AI \"*string(AI))
    #println(\"AJ \"*string(AJ))
    #println(\"AV \"*string(AV))
    A = sparse(AI[1:nel*Np*Np*dofs_per_node], AJ[1:nel*Np*Np*dofs_per_node], AV[1:nel*Np*Np*dofs_per_node]);
end

# Boundary conditions
bc_time = Base.Libc.time();
if rhs_only
    b = CGSolver.apply_boundary_conditions_rhs_only(var, b, t);
else
    (A, b) = CGSolver.apply_boundary_conditions_lhs_rhs(var, A, b, t);
end
bc_time = Base.Libc.time() - bc_time;

if rhs_only
    return b;
    
else
    log_entry(\"Elemental loop time:     \"*string(loop_time));
    log_entry(\"Boundary condition time: \"*string(bc_time));
    #println(size(A))
    #println(size(b))
    #display(Array(A))
    #println(b)
    return (A, b);
end
"
    # println(code)
    return code;
end