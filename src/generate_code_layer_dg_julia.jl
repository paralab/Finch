#=
Code generation functions for Julia.
=#

# Extract the input args. This must match the arguments passed by the solver.
function handle_input_args_dg_julia(lorr, vors)
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

# If needed, build derivative matrices
function build_derivative_matrices_dg_julia(lorr, vors)
    return "";
    code = 
"
# Note on derivative matrices:
# RQn are vandermond matrices for the derivatives of the basis functions
# with Jacobian factors. They are made like this.
# |RQ1|   | rx sx tx || Qx |
# |RQ2| = | ry sy ty || Qy |
# |RQ3|   | rz sz tz || Qz |

"
    use_full_deriv_mat = lorr==RHS;
    if config.dimension == 1
        if vors == "volume"
            code *= "(RQ1, RD1) = build_deriv_matrix(refel, J);\n";
            code *= "TRQ1 = RQ1';\n"
        else
            code *= "RQ1_1 = refel.surf_Qr[frefelind[1]] .* vol_J1.rx[1]";
            code *= "RQ1_2 = refel.surf_Qr[frefelind[2]] .* vol_J2.rx[1]";
            code *= "TRQ1_1 = RQ1_1'";
            code *= "TRQ1_2 = RQ1_2'";
        end
        
    elseif config.dimension == 2
        if vors == "volume"
            code *= "(RQ1, RQ2, RD1, RD2) = build_deriv_matrix(refel, J);\n";
            code *= "(TRQ1, TRQ2) = (RQ1', RQ2');\n"
        else
            code *= "(RQ1_1,RQ2_1,RD1_1,RD2_1) = build_face_deriv_matrix(refel, frefelind[1], vol_J1, $use_full_deriv_mat)";
            code *= "(RQ1_2,RQ2_2,RD1_2,RD2_2) = build_face_deriv_matrix(refel, frefelind[2], vol_J2, $use_full_deriv_mat)";
            code *= "(TRQ1_1,TRQ1_2) = (RQ1_1',RQ1_2')";
            code *= "(TRQ2_1,TRQ2_2) = (RQ2_1',RQ2_2')";
        end
        
    elseif config.dimension == 3
        if vors == "volume"
            code *= "(RQ1, RQ2, RQ3, RD1, RD2, RD3) = build_deriv_matrix(refel, J);\n";
            code *= "(TRQ1, TRQ2, TRQ3) = (RQ1', RQ2', RQ3');\n"
        else
            code *= "(RQ1_1,RQ2_1,RQ3_1,RD1_1,RD2_1,RD3_1) = build_face_deriv_matrix(refel, frefelind[1], vol_J1, $use_full_deriv_mat)";
            code *= "(RQ1_2,RQ2_2,RQ3_2,RD1_2,RD2_2,RD3_2) = build_face_deriv_matrix(refel, frefelind[2], vol_J2, $use_full_deriv_mat)";
            code *= "(TRQ1_1,TRQ2_1,TRQ3_1) = (RQ1_1',RQ2_1',RQ3_1')";
            code *= "(TRQ1_2,TRQ2_2,TRQ3_2) = (RQ1_2',RQ2_2',RQ3_2')";
        end
    end
    
    return code;
end

# Allocate, compute, or fetch all needed values
function prepare_needed_values_dg_julia(entities, var, lorr, vors)
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
    for i=1:length(entities)
        cname = make_entity_name(entities[i]);
        side1_entity = copy(entities[i]);
        side2_entity = copy(entities[i]);
        push!(side1_entity.flags, "DGSIDE1");
        push!(side2_entity.flags, "DGSIDE2");
        cnameside1 = make_entity_name(side1_entity);
        cnameside2 = make_entity_name(side2_entity);
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
                dgside = 0;
                if "DGSIDE1" in entities[i].flags
                    dgside = 1;
                elseif "DGSIDE2" in entities[i].flags
                    dgside = 2;
                end
                if length(entities[i].derivs) > 0
                    xyzchar = ["x","y","z"];
                    for di=1:length(entities[i].derivs)
                        if dgside == 0
                            code *= cnameside1 * " = TRQ"*string(entities[i].derivs[di])*"_1; # d/d"*xyzchar[entities[i].derivs[di]]*" of test function on side 1\n";
                            code *= cnameside2 * " = TRQ"*string(entities[i].derivs[di])*"_2; # d/d"*xyzchar[entities[i].derivs[di]]*" of test function on side 2\n";
                        elseif dgside == 1
                            code *= cname * " = TRQ"*string(entities[i].derivs[di])*"_1; # d/d"*xyzchar[entities[i].derivs[di]]*" of test function on side 1\n";
                        elseif dgside == 2
                            code *= cname * " = TRQ"*string(entities[i].derivs[di])*"_2; # d/d"*xyzchar[entities[i].derivs[di]]*" of test function on side 2\n";
                        end
                        need_deriv_matrix = true;
                    end
                else
                    if dgside == 0
                        code *= cnameside1 * " = (refel.surf_Q[frefelind[1]])'; # test function side 1.\n";
                        code *= cnameside2 * " = (refel.surf_Q[frefelind[2]])'; # test function side 2.\n";
                    elseif dgside == 1
                        code *= cname * " = (refel.surf_Q[frefelind[1]])'; # test function side 1.\n";
                    elseif dgside == 2
                        code *= cname * " = (refel.surf_Q[frefelind[2]])'; # test function side 2.\n";
                    end
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
                dgside = 0;
                if "DGSIDE1" in entities[i].flags
                    dgside = 1;
                elseif "DGSIDE2" in entities[i].flags
                    dgside = 2;
                end
                if length(entities[i].derivs) > 0
                    xyzchar = ["x","y","z"];
                    for di=1:length(entities[i].derivs)
                        if dgside == 0
                            code *= cnameside1 * " = RQ"*string(entities[i].derivs[di])*"_1; # d/d"*xyzchar[entities[i].derivs[di]]*" of trial function on side 1\n";
                            code *= cnameside2 * " = RQ"*string(entities[i].derivs[di])*"_2; # d/d"*xyzchar[entities[i].derivs[di]]*" of trial function on side 2\n";
                        elseif dgside == 1
                            code *= cname * " = RQ"*string(entities[i].derivs[di])*"_1; # d/d"*xyzchar[entities[i].derivs[di]]*" of trial function on side 1\n";
                        elseif dgside == 2
                            code *= cname * " = RQ"*string(entities[i].derivs[di])*"_2; # d/d"*xyzchar[entities[i].derivs[di]]*" of trial function on side 2\n";
                        end
                        need_deriv_matrix = true;
                    end
                else
                    if dgside == 0
                        code *= cnameside1 * " = refel.surf_Q[frefelind[1]]; # trial function side 1.\n";
                        code *= cnameside2 * " = refel.surf_Q[frefelind[2]]; # trial function side 2.\n";
                    elseif dgside == 1
                        code *= cname * " = refel.surf_Q[frefelind[1]]; # trial function side 1.\n";
                    elseif dgside == 2
                        code *= cname * " = refel.surf_Q[frefelind[2]]; # trial function side 2.\n";
                    end
                end
            end
            
        else
            # Is coefficient(number or function) or variable(array)?
            (ctype, cval) = get_coef_val(entities[i]);
            if ctype == -1
                # It was a special symbol like dt
                if entities[i].name == "FACENORMAL1"
                    code *= cname * " = normal["*string(entities[i].index)*"]; # normal vector component\n"
                    piece_needed[7] = true;
                else entities[i].name == "FACENORMAL2"
                    code *= cname * " = -normal["*string(entities[i].index)*"]; # reverse normal vector component\n"
                    piece_needed[7] = true;
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
                nodesymbol2 = "nodex[2]";
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
                    dgside = 0;
                    if "DGSIDE1" in entities[i].flags
                        dgside = 1;
                    elseif "DGSIDE2" in entities[i].flags
                        dgside = 2;
                    end
                    # If derivatives are needed, must evaluate at all volume nodes.
                    if length(entities[i].derivs) > 0
                        code *= cname * " = zeros(refel.Np);\n";
                        if dgside == 2
                            # code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargs2 * " end\n";
                            code *= "for coefi = 1:refel.Np " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", "*nodesymbol2*"[:,coefi], time, eid, fid) end\n";
                        else
                            # code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargs1 * " end\n";
                            code *= "for coefi = 1:refel.Np " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", "*nodesymbol1*"[:,coefi], time, eid, fid) end\n";
                        end
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            if dgside == 2
                                code *= cname * " = RD"*string(entities[i].derivs[di])*"_2 * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                            else
                                code *= cname * " = RD"*string(entities[i].derivs[di])*"_1 * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                            end
                        end
                        if dgside == 2
                            code *= cname * " = " * cname * "[refel.face2local[frefelind[2]]; # extract face values only.";
                        else
                            code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]; # extract face values only.";
                        end
                        need_deriv_matrix = true;
                        piece_needed[3] = true;
                        
                    else # no derivatives, only need surface nodes
                        code *= cname * " = zeros(Nfp);\n";
                        # code *= "for coefi = 1:Nfp " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargsf * " end\n";
                        code *= "for coefi = 1:refel.Nfp " * cname * "[coefi] = evaluate_coefficient(Finch.coefficients["*string(coef_index)*"], "*string(entities[i].index)*", "*nodesymbolf*"[:,coefi], time, eid, fid) end\n";
                        piece_needed[6] = true;
                    end
                    
                    # Interpolate at surface quadrature points
                    if dgside == 2
                        code *= cname * " = refel.surf_Q[frefelind[2]][:,refel.face2local[frefelind[2]]] * " * cname * "; # Interpolate at quadrature points.\n";
                    else
                        code *= cname * " = refel.surf_Q[frefelind[1]][:,refel.face2local[frefelind[1]]] * " * cname * "; # Interpolate at quadrature points.\n";
                    end
                end
                
            elseif ctype == 3 # a known variable value
                # This generates something like: coef_u_1 = copy((Finch.variables[1]).values[1, loc2glb])
                if vors == "volume"
                    code *= cname * " = copy((Finch.variables["*string(cval)*"]).values["*string(entities[i].index)*", loc2glb]);\n";
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
                    dgside = 0;
                    if "DGSIDE1" in entities[i].flags
                        dgside = 1;
                    elseif "DGSIDE2" in entities[i].flags
                        dgside = 2;
                    end
                    # If derivatives are needed, must evaluate at all volume nodes.
                    if length(entities[i].derivs) > 0
                        if dgside == 0
                            code *= cnameside1 * " = copy((Finch.variables["*string(cval)*"]).values["*string(entities[i].index)*", loc2glb[1]]);\n";
                            code *= cnameside2 * " = copy((Finch.variables["*string(cval)*"]).values["*string(entities[i].index)*", loc2glb[2]]);\n";
                        elseif dgside == 1
                            code *= cname * " = copy((Finch.variables["*string(cval)*"]).values["*string(entities[i].index)*", loc2glb[1]]);\n";
                        else
                            code *= cname * " = copy((Finch.variables["*string(cval)*"]).values["*string(entities[i].index)*", loc2glb[2]]);\n";
                        end
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            if dgside == 0
                                code *= cnameside1 * " = RD"*string(entities[i].derivs[di])*"_1 * " * cnameside1 * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                                code *= cnameside2 * " = RD"*string(entities[i].derivs[di])*"_2 * " * cnameside2 * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                            elseif dgside == 1
                                code *= cname * " = RD"*string(entities[i].derivs[di])*"_1 * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                            elseif dgside == 2
                                code *= cname * " = RD"*string(entities[i].derivs[di])*"_2 * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                            end
                        end
                        # if dgside == 0
                        #     code *= cnameside1 * " = " * cnameside1 * "[refel.face2local[frefelind[1]]; # extract face values only.";
                        #     code *= cnameside2 * " = " * cnameside2 * "[refel.face2local[frefelind[2]]; # extract face values only.";
                        # elseif dgside == 1
                        #     code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]; # extract face values only.";
                        # elseif dgside == 2
                        #     code *= cname * " = " * cname * "[refel.face2local[frefelind[2]]; # extract face values only.";
                        # end
                        need_deriv_matrix = true;
                        piece_needed[2] = true;
                        
                    else # no derivatives
                        if dgside == 0
                            code *= cnameside1 * " = copy((Finch.variables["*string(cval)*"]).values["*string(entities[i].index)*", loc2glb[1]]);\n";
                            code *= cnameside2 * " = copy((Finch.variables["*string(cval)*"]).values["*string(entities[i].index)*", loc2glb[2]]);\n";
                        elseif dgside == 1
                            code *= cname * " = copy((Finch.variables["*string(cval)*"]).values["*string(entities[i].index)*", loc2glb[1]]);\n";
                        else
                            code *= cname * " = copy((Finch.variables["*string(cval)*"]).values["*string(entities[i].index)*", loc2glb[2]]);\n";
                        end
                        # piece_needed[5] = true;
                        piece_needed[2] = true;
                    end
                    
                    # # Interpolate at surface quadrature points
                    # if dgside == 0
                    #     code *= cnameside1 * " = refel.surf_Q[frefelind[1]][:,refel.face2local[frefelind[1]]] * " * cnameside1 * "; # Interpolate at quadrature points.\n";
                    #     code *= cnameside2 * " = refel.surf_Q[frefelind[2]][:,refel.face2local[frefelind[2]]] * " * cnameside2 * "; # Interpolate at quadrature points.\n";
                    # elseif dgside == 1
                    #     code *= cname * " = refel.surf_Q[frefelind[1]][:,refel.face2local[frefelind[1]]] * " * cname * "; # Interpolate at quadrature points.\n";
                    # else
                    #     code *= cname * " = refel.surf_Q[frefelind[2]][:,refel.face2local[frefelind[2]]] * " * cname * "; # Interpolate at quadrature points.\n";
                    # end
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
# RQn are vandermond matrices for the derivatives of the basis functions
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
            needed_pieces *= 
"if neighbor == eid
    face2glb[:,2] = face2glb[:,1];
end\n"
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
# RQn are vandermond matrices for the derivatives of the basis functions
# with Jacobian factors. They are made like this.
# |RQ1|   | rx sx tx || Qx |
# |RQ2| = | ry sy ty || Qy |
# |RQ3|   | rz sz tz || Qz |

"
            use_full_deriv_mat = lorr==RHS;
            if config.dimension == 1
                needed_pieces *= "RQ1_1 = refel.surf_Qr[frefelind[1]] .* vol_J[1].rx[1]\n";
                needed_pieces *= "RQ1_2 = refel.surf_Qr[frefelind[2]] .* vol_J[2].rx[1]\n";
                needed_pieces *= "TRQ1_1 = RQ1_1'\n";
                needed_pieces *= "TRQ1_2 = RQ1_2'\n";
                
            elseif config.dimension == 2
                needed_pieces *= "(RQ1_1,RQ2_1,RD1_1,RD2_1) = build_face_deriv_matrix(refel, frefelind[1], vol_J[1], $use_full_deriv_mat)\n";
                needed_pieces *= "(RQ1_2,RQ2_2,RD1_2,RD2_2) = build_face_deriv_matrix(refel, frefelind[2], vol_J[2], $use_full_deriv_mat)\n";
                needed_pieces *= "(TRQ1_1,TRQ1_2) = (RQ1_1',RQ1_2')\n";
                needed_pieces *= "(TRQ2_1,TRQ2_2) = (RQ2_1',RQ2_2')\n";
                
            elseif config.dimension == 3
                needed_pieces *= "(RQ1_1,RQ2_1,RQ3_1,RD1_1,RD2_1,RD3_1) = build_face_deriv_matrix(refel, frefelind[1], vol_J[1], $use_full_deriv_mat)\n";
                needed_pieces *= "(RQ1_2,RQ2_2,RQ3_2,RD1_2,RD2_2,RD3_2) = build_face_deriv_matrix(refel, frefelind[2], vol_J[2], $use_full_deriv_mat)\n";
                needed_pieces *= "(TRQ1_1,TRQ2_1,TRQ3_1) = (RQ1_1',RQ2_1',RQ3_1')\n";
                needed_pieces *= "(TRQ1_2,TRQ2_2,TRQ3_2) = (RQ1_2',RQ2_2',RQ3_2')\n";
            end
        end
    end
    
    code = needed_pieces * "\n" * code;
    
    return code;
end

function make_elemental_computation_dg_julia(terms, var, dofsper, offset_ind, lorr, vors)
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
    # There will be four pieces for LHS and two for RHS.
    if dofsper > 1
        if vors == "volume"
            if lorr == RHS
                code *= "element_vector = zeros(refel.Np * "*string(dofsper)*"); # Allocate the returned vector.\n"
            else
                code *= "element_matrix = zeros(refel.Np * "*string(dofsper)*", refel.Np * "*string(dofsper)*"); # Allocate the returned matrix.\n"
            end
        else # surface
            if lorr == RHS
                code *= "element_vectorL = zeros(Nfp * "*string(dofsper)*"); # Allocate return vector for element 1.\n"
                code *= "element_vectorR = zeros(Nfp * "*string(dofsper)*"); # Allocate return vector for element 2.\n"
            else
                code *= "element_matrixLL = zeros(Nfp * "*string(dofsper)*", Nfp * "*string(dofsper)*"); # Allocate return matrix for LL.\n"
                code *= "element_matrixLR = zeros(Nfp * "*string(dofsper)*", Nfp * "*string(dofsper)*"); # Allocate return matrix for LR.\n"
                code *= "element_matrixRL = zeros(Nfp * "*string(dofsper)*", Nfp * "*string(dofsper)*"); # Allocate return matrix for RL.\n"
                code *= "element_matrixRR = zeros(Nfp * "*string(dofsper)*", Nfp * "*string(dofsper)*"); # Allocate return matrix for RR.\n"
            end
        end
        
    else
        if vors == "volume"
            if lorr == RHS
                code *= "element_vector = zeros(refel.Np); # Allocate the returned vector.\n"
            else
                code *= "element_matrix = zeros(refel.Np, refel.Np); # Allocate the returned matrix.\n"
            end
        else # surface
            if lorr == RHS
                code *= "element_vectorL = zeros(Nfp); # Allocate return vector for element 1.\n"
                code *= "element_vectorR = zeros(Nfp); # Allocate return vector for element 2.\n"
            else
                code *= "element_matrixLL = zeros(Nfp, Nfp); # Allocate return matrix for LL.\n"
                code *= "element_matrixLR = zeros(Nfp, Nfp); # Allocate return matrix for LR.\n"
                code *= "element_matrixRL = zeros(Nfp, Nfp); # Allocate return matrix for RL.\n"
                code *= "element_matrixRR = zeros(Nfp, Nfp); # Allocate return matrix for RR.\n"
            end
        end
        
    end
    
    # Separate the factors of each term into test, trial, coef and form the calculation
    if vors == "volume"
        if lorr == LHS
            submatrices = Array{String, 2}(undef, dofsper, dofsper);
        else # RHS
            submatrices = Array{String, 1}(undef, dofsper);
        end
        for smi=1:length(submatrices)
            submatrices[smi] = "";
        end
    else # surface
        if lorr == LHS
            submatricesll = Array{String, 2}(undef, dofsper, dofsper);
            submatriceslr = Array{String, 2}(undef, dofsper, dofsper);
            submatricesrl = Array{String, 2}(undef, dofsper, dofsper);
            submatricesrr = Array{String, 2}(undef, dofsper, dofsper);
            submatrices = [submatricesll, submatriceslr, submatricesrl, submatricesrr];
        else # RHS
            submatricesl = Array{String, 1}(undef, dofsper);
            submatricesr = Array{String, 1}(undef, dofsper);
            submatrices = [submatricesl, submatricesr];
        end
        for smi=1:length(submatrices)
            for smj=1:length(submatrices[smi])
                submatrices[smi][smj] = "";
            end
        end
    end
    
    if typeof(var) <: Array # multiple variables
        for vi=1:length(var) # variables
            # Process the terms for this variable
            for ci=1:length(terms[vi]) # components
                for i=1:length(terms[vi][ci])
                    (term_result, test_ind, trial_ind) = generate_term_calculation_dg_julia(terms[vi][ci][i], var, lorr, vors);
                    
                    # println(terms[vi][ci][i])
                    # println(term_result * " : "*string(test_ind)*", "*string(trial_ind))
                    
                    # Find the appropriate submatrix for this term
                    if vors == "volume"
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
                        
                    else # surface
                        submati = offset_ind[vi] + test_ind;
                        submatj = trial_ind;
                        if lorr == LHS
                            submat_ind = submati + dofsper * (submatj-1);
                            smparts = 4;
                        else
                            submat_ind = submati;
                            smparts = 2;
                        end
                        
                        for smi=1:smparts
                            if length(submatrices[smi][submat_ind]) > 1
                                if length(term_result[smi]) > 0
                                    submatrices[smi][submat_ind] *= " .+ " * term_result[smi];
                                end
                            else
                                submatrices[smi][submat_ind] = term_result[smi];
                            end
                        end
                    end
                end
            end
        end # vi
        
    else # only one variable
        # Process the terms for this variable
        for ci=1:length(terms) # components
            for i=1:length(terms[ci])
                (term_result, test_ind, trial_ind) = generate_term_calculation_dg_julia(terms[ci][i], var, lorr, vors);
                
                # Find the appropriate submatrix for this term
                if vors == "volume"
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
                    
                else # surface
                    if lorr == LHS
                        submat_ind = test_ind + dofsper * (trial_ind-1);
                        smparts = 4;
                    else
                        submat_ind = test_ind;
                        smparts = 2;
                    end
                    
                    for smi=1:smparts
                        if length(submatrices[smi][submat_ind]) > 1
                            if length(term_result[smi]) > 0
                                submatrices[smi][submat_ind] *= " .+ " * term_result[smi];
                            end
                        else
                            submatrices[smi][submat_ind] = term_result[smi];
                        end
                    end
                end
            end
        end
    end
    
    # Put the submatrices together into element_matrix or element_vector
    if lorr == LHS
        if vors == "volume"
            for emi=1:dofsper
                for emj=1:dofsper
                    if length(submatrices[emi, emj]) > 1
                        rangei = "(("*string(emi-1)*")*refel.Np + 1):("*string(emi)*"*refel.Np)";
                        rangej = "(("*string(emj-1)*")*refel.Np + 1):("*string(emj)*"*refel.Np)";
                        code *= "element_matrix["*rangei*", "*rangej*"] = " * submatrices[emi,emj] * "\n";
                    end
                end
            end
            code *= "return element_matrix;\n"
             
        else # surface
            mat_tags = ["LL", "LR", "RL", "RR"];
            frefs = [["1","1"],["1","2"],["2","1"],["2","2"]];
            for emi=1:dofsper
                for emj=1:dofsper
                    for smpart=1:4
                        if length(submatrices[smpart][emi, emj]) > 1
                            rangei = "(("*string(emi-1)*")*Nfp + 1):("*string(emi)*"*Nfp)";
                            rangej = "(("*string(emj-1)*")*Nfp + 1):("*string(emj)*"*Nfp)";
                            code *= "element_matrix"*mat_tags[smpart]*"["*rangei*", "*rangej*"] = (" * submatrices[smpart][emi,emj] * ")[refel.face2local[frefelind["*frefs[smpart][1]*"]], refel.face2local[frefelind["*frefs[smpart][2]*"]]]\n";
                        end
                    end
                end
            end
            code *= "return [element_matrixLL, element_matrixLR, element_matrixRL, element_matrixRR];\n"
        end
        
    else # RHS
        if vors == "volume"
            for emi=1:dofsper
                if length(submatrices[emi]) > 1
                    rangei = "(("*string(emi-1)*")*refel.Np + 1):("*string(emi)*"*refel.Np)";
                    code *= "element_vector["*rangei*"] = " * submatrices[emi] * "\n";
                end
            end
            code *= "return element_vector;\n"
            
        else # surface
            mat_tags = ["L", "R"];
            frefs = ["1","2"];
            for emi=1:dofsper
                for smpart=1:2
                    if length(submatrices[smpart][emi]) > 1
                        rangei = "(("*string(emi-1)*")*Nfp + 1):("*string(emi)*"*Nfp)";
                        code *= "element_vector"*mat_tags[smpart]*"["*rangei*"] = (" * submatrices[smpart][emi] * ")[refel.face2local[frefelind["*frefs[smpart]*"]]]\n";
                    end
                end
            end
            code *= "return [element_vectorL, element_vectorR];\n"
        end
    end
    
    return code;
end

function generate_term_calculation_dg_julia(term, var, lorr, vors)
    result = "";
    
    if vors == "volume"
        if lorr == LHS
            (test_part, trial_part, coef_part, test_ind, trial_ind) = separate_factors(term, var);
            # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
            if !(coef_part === nothing)
                result = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
                        string(replace_entities_with_symbols(coef_part)) * ") * " * 
                        string(replace_entities_with_symbols(trial_part));
            else # no coef_part
                result = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
                        string(replace_entities_with_symbols(trial_part));
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
        
    else # surface
        if lorr == LHS
            result = ["","","",""];
            (test_part, trial_part, coef_part, test_ind, trial_ind) = separate_factors(term, var);
            
            # Determine which result, LL, LR, RL, RR to add to based on test and trial parts.
            # For example, SIDE1 test -> L and SIDE2 trial -> R will contribute to LR.
            # If there's no SIDE flag, diagonals are possible.
            test_side = 0;
            trial_side = 0;
            if "DGSIDE1" in test_part.flags
                test_side = 1;
            elseif "DGSIDE2" in test_part.flags
                test_side = 2;
            else
                test_part_side1 = copy(test_part);
                test_part_side2 = copy(test_part);
                push!(test_part_side1.flags, "DGSIDE1");
                push!(test_part_side2.flags, "DGSIDE2");
            end
            if "DGSIDE1" in trial_part.flags
                trial_side = 1;
                # also apply side1 to known variables
                coef_part_side1 = add_flag_to_var_entities(copy(coef_part), variables, "DGSIDE1", nevermind="DGSIDE", copy_ent=true);
            elseif "DGSIDE2" in trial_part.flags
                trial_side = 2;
                # also apply side2 to known variables
                coef_part_side2 = add_flag_to_var_entities(copy(coef_part), variables, "DGSIDE2", nevermind="DGSIDE", copy_ent=true);
            else
                trial_part_side1 = copy(trial_part);
                trial_part_side2 = copy(trial_part);
                push!(trial_part_side1.flags, "DGSIDE1");
                push!(trial_part_side2.flags, "DGSIDE2");
                # also apply side to known variables
                coef_part_side1 = add_flag_to_var_entities(copy(coef_part), variables, "DGSIDE1", nevermind="DGSIDE", copy_ent=true);
                coef_part_side2 = add_flag_to_var_entities(copy(coef_part), variables, "DGSIDE2", nevermind="DGSIDE", copy_ent=true);
            end
            
            # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
            if test_side == 1
                if trial_side == 1
                    # LL only
                    if !(coef_part === nothing)
                        result[1] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
                                string(replace_entities_with_symbols(coef_part_side1)) * ") * " * 
                                string(replace_entities_with_symbols(trial_part));
                    else # no coef_part
                        result[1] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
                                string(replace_entities_with_symbols(trial_part));
                    end
                    
                elseif trial_side == 2
                    # LR only
                    if !(coef_part === nothing)
                        result[2] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
                                string(replace_entities_with_symbols(coef_part_side2)) * ") * " * 
                                string(replace_entities_with_symbols(trial_part));
                    else # no coef_part
                        result[2] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
                                string(replace_entities_with_symbols(trial_part));
                    end
                    
                else
                    # LL only and must change trial symbols to side1
                    if !(coef_part === nothing)
                        result[1] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
                                string(replace_entities_with_symbols(coef_part_side1)) * ") * " * 
                                string(replace_entities_with_symbols(trial_part_side1));
                    else # no coef_part
                        result[1] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
                                string(replace_entities_with_symbols(trial_part_side1));
                    end
                end
                
                
            elseif test_side == 2
                if trial_side == 1
                    # RL only
                    if !(coef_part === nothing)
                        result[3] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
                                string(replace_entities_with_symbols(coef_part_side1)) * ") * " * 
                                string(replace_entities_with_symbols(trial_part));
                    else # no coef_part
                        result[3] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
                                string(replace_entities_with_symbols(trial_part));
                    end
                    
                elseif trial_side == 2
                    # RR only
                    if !(coef_part === nothing)
                        result[4] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
                                string(replace_entities_with_symbols(coef_part_side2)) * ") * " * 
                                string(replace_entities_with_symbols(trial_part));
                    else # no coef_part
                        result[4] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
                                string(replace_entities_with_symbols(trial_part));
                    end
                    
                else
                    # RR only and must change trial symbols to side2
                    if !(coef_part === nothing)
                        result[4] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
                                string(replace_entities_with_symbols(coef_part_side2)) * ") * " * 
                                string(replace_entities_with_symbols(trial_part_side2));
                    else # no coef_part
                        result[4] = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
                                string(replace_entities_with_symbols(trial_part_side2));
                    end
                end
                
            else # no test side specified. Diagonals are possible
                if trial_side == 1
                    # LL only
                    if !(coef_part === nothing)
                        result[1] = string(replace_entities_with_symbols(test_part_side1)) * " * diagm(wdetj .* " * 
                                string(replace_entities_with_symbols(coef_part_side1)) * ") * " * 
                                string(replace_entities_with_symbols(trial_part));
                    else # no coef_part
                        result[1] = string(replace_entities_with_symbols(test_part_side1)) * " * diagm(wdetj) * " * 
                                string(replace_entities_with_symbols(trial_part));
                    end
                    
                elseif trial_side == 2
                    # RR only
                    if !(coef_part === nothing)
                        result[4] = string(replace_entities_with_symbols(test_part_side2)) * " * diagm(wdetj .* " * 
                                string(replace_entities_with_symbols(coef_part_side2)) * ") * " * 
                                string(replace_entities_with_symbols(trial_part));
                    else # no coef_part
                        result[4] = string(replace_entities_with_symbols(test_part_side2)) * " * diagm(wdetj) * " * 
                                string(replace_entities_with_symbols(trial_part));
                    end
                    
                else
                    # LL and RR and must change both symbols
                    if !(coef_part === nothing)
                        result[1] = string(replace_entities_with_symbols(test_part_side1)) * " * diagm(wdetj .*" * 
                                string(replace_entities_with_symbols(coef_part_side1)) * ") * " * 
                                string(replace_entities_with_symbols(trial_part_side1));
                        result[4] = string(replace_entities_with_symbols(test_part_side2)) * " * diagm(wdetj .* " * 
                                string(replace_entities_with_symbols(coef_part_side2)) * ") * " * 
                                string(replace_entities_with_symbols(trial_part_side2));
                    else # no coef_part
                        result[1] = string(replace_entities_with_symbols(test_part_side1)) * " * diagm(wdetj) * " * 
                                string(replace_entities_with_symbols(trial_part_side1));
                        result[4] = string(replace_entities_with_symbols(test_part_side2)) * " * diagm(wdetj) * " * 
                                string(replace_entities_with_symbols(trial_part_side2));
                    end
                end
            end
            
        else # RHS
            result = ["",""];
            (test_part, trial_part, coef_part, test_ind, trial_ind) = separate_factors(term);
            # figure out what side, if any, the coefficient part is on
            function get_coef_side(ex)
                side = 0;
                if typeof(ex) == Expr
                    for i=1:length(ex.args)
                        side = max(side, get_coef_side(ex.args[i]));
                    end
                elseif typeof(ex) == SymEntity
                    if "DGSIDE1" in ex.flags
                        side = 1;
                    elseif "DGSIDE2" in ex.flags
                        side = 2;
                    end
                end
                return side;
            end
            coef_side = get_coef_side(coef_part);
            if coef_side == 0
                coef_part_side1 = add_flag_to_var_entities(copy(coef_part), variables, "DGSIDE1", nevermind="DGSIDE", copy_ent=true);
                coef_part_side2 = add_flag_to_var_entities(copy(coef_part), variables, "DGSIDE2", nevermind="DGSIDE", copy_ent=true);
            end
            test_side = 0;
            if "DGSIDE1" in test_part.flags
                test_side = 1;
            elseif "DGSIDE2" in test_part.flags
                test_side = 2;
            else
                test_part_side1 = copy(test_part);
                test_part_side2 = copy(test_part);
                push!(test_part_side1.flags, "DGSIDE1");
                push!(test_part_side2.flags, "DGSIDE2");
            end
            
            # RHS: test_part * (weight_part .* coef_part)
            if test_side == 1
                if !(coef_part === nothing)
                    if coef_side == 1
                        result[1] = string(replace_entities_with_symbols(test_part)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] * " * 
                            string(replace_entities_with_symbols(coef_part)) * "))";
                    elseif coef_side == 2
                        result[1] = string(replace_entities_with_symbols(test_part)) * " * (wdetj .* (refel.surf_Q[frefelind[2]] * " * 
                            string(replace_entities_with_symbols(coef_part)) * "))";
                    else
                        result[1] = string(replace_entities_with_symbols(test_part)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] * " * 
                            string(replace_entities_with_symbols(coef_part_side1)) * " .+ refel.surf_Q[frefelind[2]] * " *
                            string(replace_entities_with_symbols(coef_part_side2)) * "))";
                    end
                    
                else
                    result[1] = string(replace_entities_with_symbols(test_part)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] .+ refel.surf_Q[frefelind[2]]))";
                end
            elseif test_side == 2
                if !(coef_part === nothing)
                    if coef_side == 1
                        result[2] = string(replace_entities_with_symbols(test_part)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] * " * 
                            string(replace_entities_with_symbols(coef_part)) * "))";
                    elseif coef_side == 2
                        result[2] = string(replace_entities_with_symbols(test_part)) * " * (wdetj .* (refel.surf_Q[frefelind[2]] * " * 
                            string(replace_entities_with_symbols(coef_part)) * "))";
                    else
                        result[2] = string(replace_entities_with_symbols(test_part)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] * " * 
                            string(replace_entities_with_symbols(coef_part_side1)) * " .+ refel.surf_Q[frefelind[2]] * " *
                            string(replace_entities_with_symbols(coef_part_side2)) * "))";
                    end
                else
                    result[2] = string(replace_entities_with_symbols(test_part)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] .+ refel.surf_Q[frefelind[2]]))";
                end
            else # no test side specified
                if !(coef_part === nothing)
                    if coef_side == 1
                        result[1] = string(replace_entities_with_symbols(test_part_side1)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] * " * 
                            string(replace_entities_with_symbols(coef_part)) * "))";
                        result[2] = string(replace_entities_with_symbols(test_part_side2)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] * " * 
                            string(replace_entities_with_symbols(coef_part)) * "))";
                    elseif coef_side == 2
                        result[1] = string(replace_entities_with_symbols(test_part_side1)) * " * (wdetj .* (refel.surf_Q[frefelind[2]] * " * 
                            string(replace_entities_with_symbols(coef_part)) * "))";
                        result[2] = string(replace_entities_with_symbols(test_part_side2)) * " * (wdetj .* (refel.surf_Q[frefelind[2]] * " * 
                            string(replace_entities_with_symbols(coef_part)) * "))";
                    else
                        result[1] = string(replace_entities_with_symbols(test_part_side1)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] * " * 
                            string(replace_entities_with_symbols(coef_part_side1)) * "))";
                        result[2] = string(replace_entities_with_symbols(test_part_side2)) * " * (wdetj .* (refel.surf_Q[frefelind[2]] * " *
                            string(replace_entities_with_symbols(coef_part_side2)) * "))";
                    end
                else
                    result[1] = string(replace_entities_with_symbols(test_part_side1)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] .+ refel.surf_Q[frefelind[2]]))";
                    result[2] = string(replace_entities_with_symbols(test_part_side2)) * " * (wdetj .* (refel.surf_Q[frefelind[1]] .+ refel.surf_Q[frefelind[2]]))";
                end
            end
            
            
        end
    end
    
    
    return (result, test_ind, trial_ind);
end
