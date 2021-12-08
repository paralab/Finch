###############################################################################################################
# generate for cachesim
###############################################################################################################
#=
Want

l 2342
s 512 8
l 512 8

to be turned into

cs.load(2342)  # Loads one byte from address 2342
cs.store(512, length=8)  # Stores 8 bytes to addresses 512-519
cs.load(512, length=8)  # Loads from address 512 until (exclusive) 520 (eight bytes)
...
cs.force_write_back()
cs.print_stats()

The arrays are indexed as such:
1 A
2 b
3 u
4 Ak
5 bk
6 Q
7 RQ1
8 RQ2
9 RQ3
10 TRQ1
11 TRQ2
12 TRQ3
13 RD1
14 RD2
15 RD3
17 wdetJ
18+ any other needed arrays
=#

function generate_code_layer_cachesim(var, entities, terms, lorr, vors)
    multivar = typeof(var) <:Array;
    varcount = 1;
    dofsper = 0;
    offset_ind = [0];
    if multivar
        varcount = length(var);
        offset_ind = zeros(Int, varcount);
        dofsper = length(var[1].symvar);
        for i=2:length(var)
            offset_ind[i] = dofsper;
            dofsper = dofsper + length(var[i].symvar);
        end
    else
        dofsper = length(var.symvar);
    end
    
    code = handle_input_args_cachesim(lorr, vors);
    code *= "\n";
    code *= prepare_needed_values_cachesim(entities, var, lorr, vors);
    code *= "\n";
    code *= make_elemental_computation_cachesim(terms, var, dofsper, offset_ind, lorr, vors);
    
    return (code, code_string_to_expr(code));
end

# Extract the input args. This must match the arguments passed by the solver.
function handle_input_args_cachesim(lorr, vors)
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
    return code;
end

# Allocate, compute, or fetch all needed values
function prepare_needed_values_cachesim(entities, var, lorr, vors)
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
                        # code *= cname * " = TRQ"*string(entities[i].derivs[di])*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of test function\n";
                        code *= "Finch.cachesim_load_range(7)"*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of test function\n";
                    end
                    need_deriv_matrix = true;
                else
                    code *= "Finch.cachesim_load_range(7)"*"; # test function.\n";
                end
            else
                if length(entities[i].derivs) > 0
                    xyzchar = ["x","y","z"];
                    for di=1:length(entities[i].derivs)
                        code *= "Finch.cachesim_load_range(7)"*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of test function\n";
                        need_deriv_matrix = true;
                    end
                else
                    code *= "Finch.cachesim_load_range(7)"*"; # test function.\n";
                end
            end
            
        elseif is_unknown_var(entities[i], var) && lorr == LHS
            if vors == "volume"
                if length(entities[i].derivs) > 0
                    xyzchar = ["x","y","z"];
                    for di=1:length(entities[i].derivs)
                        code *= "Finch.cachesim_load_range(7)"*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of trial function\n";
                    end
                    need_deriv_matrix = true;
                else
                    code *= "Finch.cachesim_load_range(7)"*"; # trial function.\n";
                end
            else
                if length(entities[i].derivs) > 0
                    xyzchar = ["x","y","z"];
                    for di=1:length(entities[i].derivs)
                        code *= "Finch.cachesim_load_range(7)"*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of trial function\n";
                        need_deriv_matrix = true;
                    end
                else
                    code *= "Finch.cachesim_load_range(7)"*"; # trial function.\n";
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
                cargs = "(nodex[coefi], 0, 0, time)";
                if config.dimension == 2
                    cargs = "(nodex[1, coefi], nodex[2, coefi], 0, time)";
                elseif config.dimension == 3
                    cargs = "(nodex[1, coefi], nodex[2, coefi], nodex[3, coefi], time)";
                end
                cargs1 = "(nodex[1][coefi], 0, 0, time)";
                if config.dimension == 2
                    cargs1 = "(nodex[1][1, coefi], nodex[1][2, coefi], 0, time)";
                elseif config.dimension == 3
                    cargs1 = "(nodex[1][1, coefi], nodex[1][2, coefi], nodex[1][3, coefi], time)";
                end
                cargs2 = "(nodex[2][coefi], 0, 0, time)";
                if config.dimension == 2
                    cargs2 = "(nodex[2][1, coefi], nodex[2][2, coefi], 0, time)";
                elseif config.dimension == 3
                    cargs2 = "(nodex[2][1, coefi], nodex[2][2, coefi], nodex[2][3, coefi], time)";
                end
                cargsf = "(facex[coefi], 0, 0, time)";
                if config.dimension == 2
                    cargsf = "(facex[1, coefi], facex[2, coefi], 0, time)";
                elseif config.dimension == 3
                    cargsf = "(facex[1, coefi], facex[2, coefi], facex[3, coefi], time)";
                end
                if vors == "volume"
                    code *= "# "*cname * " = zeros(refel.Np);\n";
                    code *= "# "*"for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargs * " end\n";
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= "# "*cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*" and interpolate at quadrature points.\n";
                        end
                        need_deriv_matrix = true;
                    else
                        code *= "# "*cname * " = refel.Q * " * cname * "; # Interpolate at quadrature points.\n";
                    end
                    piece_needed[3] = true;
                    
                else # surface
                    # If derivatives are needed, must evaluate at all volume nodes.
                    if length(entities[i].derivs) > 0
                        code *= "# "*cname * " = zeros(refel.Np);\n";
                        code *= "# "*"for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargs1 * " end\n";
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= "# "*cname * " = RD"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                        code *= "# "*cname * " = " * cname * "[refel.face2local[frefelind[1]]]; # extract face values only.";
                        need_deriv_matrix = true;
                        piece_needed[3] = true;
                        
                    else # no derivatives, only need surface nodes
                        code *= "# "*cname * " = zeros(Nfp);\n";
                        code *= "# "*"for coefi = 1:Nfp " * cname * "[coefi] = (Finch.genfunctions["*string(cval)*"]).func" * cargsf * " end\n";
                        piece_needed[6] = true;
                    end
                    
                    # Interpolate at surface quadrature points
                    code *= "# "*cname * " = refel.surf_Q[frefelind[1]][:,refel.face2local[frefelind[1]]] * " * cname * "; # Interpolate at quadrature points.\n";
                end
                
            elseif ctype == 3 # a known variable value
                # This generates something like: coef_u_1 = copy((Finch.variables[1]).values[1, loc2glb])
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
                    # code *= cname * " = copy((Finch.variables["*string(cval)*"]).values["*indstr*", loc2glb]);\n";
                    csarrayind = 17+cval;
                    code *= "Finch.cachesim_load_range($csarrayind, loc2glb)";
                    # Apply any needed derivative operators.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= "# "*cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*"\n";
                        end
                        need_deriv_matrix = true;
                    else
                        code *= "# "*cname * " = refel.Q * " * cname * "; # Interpolate at quadrature points.\n";
                    end
                    piece_needed[2] = true;
                    
                else # surface
                    # If derivatives are needed, must evaluate at all volume nodes.
                    csarrayind = 17+cval;
                    code *= "Finch.cachesim_load_range($csarrayind, loc2glb[1])";
                    if length(entities[i].derivs) > 0
                        # code *= cname * " = copy((Finch.variables["*string(cval)*"]).values["*indstr*", loc2glb[1]]);\n";
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= "# "*cname * " = RD"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                        #     code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]; # extract face values only.";
                        need_deriv_matrix = true;
                        piece_needed[2] = true;
                        
                    else # no derivatives
                        code *= "# "*cname * " = copy((Finch.variables["*string(cval)*"]).values["*indstr*", loc2glb[1]]);\n";
                        # piece_needed[5] = true;
                        piece_needed[2] = true;
                    end
                    
                    #     code *= cname * " = refel.surf_Q[frefelind[1]][:,refel.face2local[frefelind[1]]] * " * cname * "; # Interpolate at quadrature points.\n";
                end
                
            elseif ctype == 4 # an indexed coefficient
                # # This generates something like:
                # ######################################
                # # coef_n_1 = zeros(refel.Np); # allocate
                # # for coefi = 1:refel.Np
                # #     coef_k_1[coefi] = (Finch.coefficients[cval]).value[INDEX_VAL_i].func(x[1,coefi], x[2,coefi],x[3,coefi],time); # evaluate at nodes
                # # end
                # ######################################
                # if vors == "surface" && length(entities[i].derivs) == 0
                #     nodesymbol = "facex"
                #     piece_needed[6] = true;
                # else
                #     nodesymbol = "nodex"
                #     piece_needed[3] = true;
                # end
                # cargs = "("*nodesymbol*"[coefi], 0, 0, time)";
                # if config.dimension == 2
                #     cargs = "("*nodesymbol*"[1, coefi], "*nodesymbol*"[2, coefi], 0, time)";
                # elseif config.dimension == 3
                #     cargs = "("*nodesymbol*"[1, coefi], "*nodesymbol*"[2, coefi], "*nodesymbol*"[3, coefi], time)";
                # end
                
                # indstr = "";
                # for indi=1:length(entities[i].index)
                #     if indi>1
                #         indstr *= ",";
                #     end
                #     indstr *= "INDEX_VAL_"*entities[i].index[indi];
                # end
                
                # if vors == "volume"
                #     code *= cname * " = zeros(refel.Np);\n";
                #     code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.coefficients["*string(cval)*"]).value["*indstr*"].func" * cargs * " end\n";
                #     # Apply any needed derivative operators. Interpolate at quadrature points.
                #     if length(entities[i].derivs) > 0
                #         xyzchar = ["x","y","z"];
                #         for di=1:length(entities[i].derivs)
                #             code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                #                     "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*" and interpolate at quadrature points.\n";
                #         end
                #         need_deriv_matrix = true;
                #     else
                #         code *= cname * " = refel.Q * " * cname * "; # Interpolate at quadrature points.\n";
                #     end
                    
                # else # surface
                #     # Apply any needed derivative operators. Interpolate at quadrature points.
                #     if length(entities[i].derivs) > 0
                #         code *= cname * " = zeros(refel.Np);\n";
                #         code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Finch.coefficients["*string(cval)*"]).value["*indstr*"].func" * cargs * " end\n";
                #         xyzchar = ["x","y","z"];
                #         for di=1:length(entities[i].derivs)
                #             code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                #                     "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                #         end
                #         code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]]; # extract face values only.";
                #         need_deriv_matrix = true;
                        
                #     else # no derivatives, only need surface nodes
                #         code *= cname * " = zeros(refel.Nfp[frefelind[1]]);\n";
                #         code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = (Finch.coefficients["*string(cval)*"]).value["*indstr*"].func" * cargs * " end\n";
                #     end
                #     # Interpolate at surface quadrature points
                #     code *= cname * " = refel.surf_Q[frefelind[1]][:,refel.face2local[frefelind[1]]] * " * cname * "; # Interpolate at quadrature points.\n";
                # end
                
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
            needed_pieces *= "# Build derivative matrices"
            
            needed_pieces *= "Finch.cachesim_load_range(13)\n";
            needed_pieces *= "Finch.cachesim_load_range(14)\n";
            needed_pieces *= "Finch.cachesim_load_range(15)\n";
            needed_pieces *= "Finch.cachesim_load_range(16)\n";
            if config.dimension == 1
                needed_pieces *= "Finch.cachesim_load_range(7)\n";
                needed_pieces *= "Finch.cachesim_load_range(10)\n";
            elseif config.dimension == 2
                needed_pieces *= "Finch.cachesim_load_range(7)\n";
                needed_pieces *= "Finch.cachesim_load_range(8)\n";
                needed_pieces *= "Finch.cachesim_load_range(10)\n";
                needed_pieces *= "Finch.cachesim_load_range(11)\n";
            elseif config.dimension == 3
                needed_pieces *= "Finch.cachesim_load_range(7)\n";
                needed_pieces *= "Finch.cachesim_load_range(8)\n";
                needed_pieces *= "Finch.cachesim_load_range(9)\n";
                needed_pieces *= "Finch.cachesim_load_range(10)\n";
                needed_pieces *= "Finch.cachesim_load_range(11)\n";
                needed_pieces *= "Finch.cachesim_load_range(12)\n";
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
            needed_pieces *= "Finch.cachesim_load_range(13)";
            needed_pieces *= "Finch.cachesim_load_range(14)";
            needed_pieces *= "Finch.cachesim_load_range(15)";
            needed_pieces *= "Finch.cachesim_load_range(16)";
            if config.dimension == 1
                needed_pieces *= "Finch.cachesim_load_range(7)";
                needed_pieces *= "Finch.cachesim_load_range(10)";
            elseif config.dimension == 2
                needed_pieces *= "Finch.cachesim_load_range(7)";
                needed_pieces *= "Finch.cachesim_load_range(8)";
                needed_pieces *= "Finch.cachesim_load_range(10)";
                needed_pieces *= "Finch.cachesim_load_range(11)";
            elseif config.dimension == 3
                needed_pieces *= "Finch.cachesim_load_range(7)";
                needed_pieces *= "Finch.cachesim_load_range(8)";
                needed_pieces *= "Finch.cachesim_load_range(9)";
                needed_pieces *= "Finch.cachesim_load_range(10)";
                needed_pieces *= "Finch.cachesim_load_range(11)";
                needed_pieces *= "Finch.cachesim_load_range(12)";
            end
        end
    end
    
    code = needed_pieces * "\n" * code;
    
    return code;
end

function make_elemental_computation_cachesim(terms, var, dofsper, offset_ind, lorr, vors)
    # don't actually do anything here. maybe later
    return "";
end

# Generate the assembly loop structures and insert the content
function generate_assembly_loop_cachesim(indices)
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
    ind_offset = ""
    prev_ind = ""
    for i=1:length(indices)
        if indices[i] == "elements" || indices[i] == "cells"
            loop_start *= "for eid in elemental_order\n";
            # loop_start *= "face_done .= 0; # Reset face_done in elemental loop\n";
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