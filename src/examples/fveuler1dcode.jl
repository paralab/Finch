#=
Generated functions for FVeuler1d
=#

# No lhs volume set for q1

# function lhs_surface_function_for_q1(args; kwargs...)
# var =       args[1];
# eid =       args[2];
# fid =       args[3];
# grid =      args[4];
# geo_facs =  args[5];
# fv_data =   args[6];
# refel =     args[7];
# time =      args[8];
# dt =        args[9];

# if grid.face2element[1, fid] == eid # The normal points out of e
#     neighbor = grid.face2element[2, fid];
# else # The normal points into e
#     neighbor = grid.face2element[1, fid];
# end
# if neighbor == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
#     neighbor = eid;
# end



# cell_matrix = zeros(refel.Nfp[frefelind[1]] * 3, 3); # Allocate for returned matrix.
# return cell_matrix;

# end #lhs_surface_function_for_q1

# No lhs volume set for q2

# function lhs_surface_function_for_q2(args; kwargs...)
# var =       args[1];
# eid =       args[2];
# fid =       args[3];
# grid =      args[4];
# geo_facs =  args[5];
# fv_data =   args[6];
# refel =     args[7];
# time =      args[8];
# dt =        args[9];

# if grid.face2element[1, fid] == eid # The normal points out of e
#     neighbor = grid.face2element[2, fid];
# else # The normal points into e
#     neighbor = grid.face2element[1, fid];
# end
# if neighbor == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
#     neighbor = eid;
# end



# cell_matrix = zeros(refel.Nfp[frefelind[1]] * 3, 3); # Allocate for returned matrix.
# return cell_matrix;

# end #lhs_surface_function_for_q2

# No lhs volume set for q3

# function lhs_surface_function_for_q3(args; kwargs...)
# var =       args[1];
# eid =       args[2];
# fid =       args[3];
# grid =      args[4];
# geo_facs =  args[5];
# fv_data =   args[6];
# refel =     args[7];
# time =      args[8];
# dt =        args[9];

# if grid.face2element[1, fid] == eid # The normal points out of e
#     neighbor = grid.face2element[2, fid];
# else # The normal points into e
#     neighbor = grid.face2element[1, fid];
# end
# if neighbor == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
#     neighbor = eid;
# end



# cell_matrix = zeros(refel.Nfp[frefelind[1]] * 3, 3); # Allocate for returned matrix.
# return cell_matrix;

# end #lhs_surface_function_for_q3

# No rhs volume set for q1

function rhs_surface_function_for_q1(args; kwargs...)
var =       args[1];
eid =       args[2];
fid =       args[3];
grid =      args[4];
geo_facs =  args[5];
fv_data =   args[6];
refel =     args[7];
time =      args[8];
dt =        args[9];

left_cells = kwargs[1]; # for higher order
right_cells = kwargs[2];

normal = grid.facenormals[:,fid];
if grid.face2element[1, fid] == eid # The normal points out of e
    neighbor = grid.face2element[2, fid];
    if normal[1] < 0
        # the left side is neighbor
        left_cells = kwargs[2];
        right_cells = kwargs[1];
        els = [neighbor, eid];
    else
        # left is eid
        els = [eid, neighbor];
    end
else # The normal points into e
    neighbor = grid.face2element[1, fid];
    if normal[1] > 0
        # the left side is neighbor
        left_cells = kwargs[2];
        right_cells = kwargs[1];
        els = [neighbor, eid];
    else
        # left is eid
        els = [eid, neighbor];
    end
end
if neighbor == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
    neighbor = eid;
    els = [eid, neighbor]; # indices of elements on both sides\
end

# Reconstruct values on both sides
# left_cells = [els[1]];
# right_cells = [els[2]];
cellx_left = fv_info.cellCenters[:, left_cells];
cellx_right = fv_info.cellCenters[:, right_cells];
facex = fv_info.faceCenters[:, fid];

cellval_r_1_left = Finch.variables[1].values[1, left_cells];
cellval_r_1_right = Finch.variables[1].values[1, right_cells];
(coef__r_1_left, coef__r_1_right) = Finch.FV_reconstruct_value_left_right(cellx_left, cellx_right, cellval_r_1_left, cellval_r_1_right, facex, limiter="vanleer");

cellval_u_1_left = Finch.variables[2].values[1, left_cells];
cellval_u_1_right = Finch.variables[2].values[1, right_cells];
(coef__u_1_left, coef__u_1_right) = Finch.FV_reconstruct_value_left_right(cellx_left, cellx_right, cellval_u_1_left, cellval_u_1_right, facex, limiter="vanleer");

cellval_p_1_left = Finch.variables[3].values[1, left_cells];
cellval_p_1_right = Finch.variables[3].values[1, right_cells];
(coef__p_1_left, coef__p_1_right) = Finch.FV_reconstruct_value_left_right(cellx_left, cellx_right, cellval_p_1_left, cellval_p_1_right, facex, limiter="vanleer");

coef__gamma_1 = 1.4;

# compute left and right moving flux on both sides
U_left = [coef__r_1_left, coef__u_1_left, coef__p_1_left];
U_right = [coef__r_1_right, coef__u_1_right, coef__p_1_right];
(Fm_left, Fp_left) = Finch.callback_functions[1].func(U_left);
(Fm_right, Fp_right) = Finch.callback_functions[1].func(U_right);

# Combine
if els[1] == eid
    result = -(Fp_left .+ Fm_right);
else # switch sign
    result = (Fp_left .+ Fm_right);
end

if els[1] == els[2] && normal[1] < 0
    result = -result;
end

return result;

end #rhs_surface_function_for_q1


function rhs_surface_function_for_q2(args; kwargs...)
    var =       args[1];
eid =       args[2];
fid =       args[3];
grid =      args[4];
geo_facs =  args[5];
fv_data =   args[6];
refel =     args[7];
time =      args[8];
dt =        args[9];

left_cells = kwargs[1]; # for higher order
right_cells = kwargs[2];

normal = grid.facenormals[:,fid];
if grid.face2element[1, fid] == eid # The normal points out of e
    neighbor = grid.face2element[2, fid];
    if normal[1] < 0
        # the left side is neighbor
        left_cells = kwargs[2];
        right_cells = kwargs[1];
        els = [neighbor, eid];
    else
        # left is eid
        els = [eid, neighbor];
    end
else # The normal points into e
    neighbor = grid.face2element[1, fid];
    if normal[1] > 0
        # the left side is neighbor
        left_cells = kwargs[2];
        right_cells = kwargs[1];
        els = [neighbor, eid];
    else
        # left is eid
        els = [eid, neighbor];
    end
end
if neighbor == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
    neighbor = eid;
    els = [eid, neighbor]; # indices of elements on both sides\
end

# Reconstruct values on both sides
# left_cells = [els[1]];
# right_cells = [els[2]];
cellx_left = fv_info.cellCenters[:, left_cells];
cellx_right = fv_info.cellCenters[:, right_cells];
facex = fv_info.faceCenters[:, fid];

cellval_r_1_left = Finch.variables[1].values[1, left_cells];
cellval_r_1_right = Finch.variables[1].values[1, right_cells];
(coef__r_1_left, coef__r_1_right) = Finch.FV_reconstruct_value_left_right(cellx_left, cellx_right, cellval_r_1_left, cellval_r_1_right, facex, limiter="vanleer");

cellval_u_1_left = Finch.variables[2].values[1, left_cells];
cellval_u_1_right = Finch.variables[2].values[1, right_cells];
(coef__u_1_left, coef__u_1_right) = Finch.FV_reconstruct_value_left_right(cellx_left, cellx_right, cellval_u_1_left, cellval_u_1_right, facex, limiter="vanleer");

cellval_p_1_left = Finch.variables[3].values[1, left_cells];
cellval_p_1_right = Finch.variables[3].values[1, right_cells];
(coef__p_1_left, coef__p_1_right) = Finch.FV_reconstruct_value_left_right(cellx_left, cellx_right, cellval_p_1_left, cellval_p_1_right, facex, limiter="vanleer");

coef__gamma_1 = 1.4;

# compute left and right moving flux on both sides
U_left = [coef__r_1_left, coef__u_1_left, coef__p_1_left];
U_right = [coef__r_1_right, coef__u_1_right, coef__p_1_right];
(Fm_left, Fp_left) = Finch.callback_functions[1].func(U_left);
(Fm_right, Fp_right) = Finch.callback_functions[1].func(U_right);

# Combine
if els[1] == eid
    result = -(Fp_left .+ Fm_right);
else # switch sign
    result = (Fp_left .+ Fm_right);
end

if els[1] == els[2] && normal[1] < 0
    result = -result;
end

return result;
    
end #rhs_surface_function_for_q2

# No assembly function set for q2

# No rhs volume set for q3

function rhs_surface_function_for_q3(args; kwargs...)
    var =       args[1];
    eid =       args[2];
    fid =       args[3];
    grid =      args[4];
    geo_facs =  args[5];
    fv_data =   args[6];
    refel =     args[7];
    time =      args[8];
    dt =        args[9];
    
    left_cells = kwargs[1]; # for higher order
    right_cells = kwargs[2];
    
    normal = grid.facenormals[:,fid];
    if grid.face2element[1, fid] == eid # The normal points out of e
        neighbor = grid.face2element[2, fid];
        if normal[1] < 0
            # the left side is neighbor
            left_cells = kwargs[2];
            right_cells = kwargs[1];
            els = [neighbor, eid];
        else
            # left is eid
            els = [eid, neighbor];
        end
    else # The normal points into e
        neighbor = grid.face2element[1, fid];
        if normal[1] > 0
            # the left side is neighbor
            left_cells = kwargs[2];
            right_cells = kwargs[1];
            els = [neighbor, eid];
        else
            # left is eid
            els = [eid, neighbor];
        end
    end
    if neighbor == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
        neighbor = eid;
        els = [eid, neighbor]; # indices of elements on both sides\
    end
    
    # Reconstruct values on both sides
    # left_cells = [els[1]];
    # right_cells = [els[2]];
    cellx_left = fv_info.cellCenters[:, left_cells];
    cellx_right = fv_info.cellCenters[:, right_cells];
    facex = fv_info.faceCenters[:, fid];
    
    cellval_r_1_left = Finch.variables[1].values[1, left_cells];
    cellval_r_1_right = Finch.variables[1].values[1, right_cells];
    (coef__r_1_left, coef__r_1_right) = Finch.FV_reconstruct_value_left_right(cellx_left, cellx_right, cellval_r_1_left, cellval_r_1_right, facex, limiter="vanleer");
    
    cellval_u_1_left = Finch.variables[2].values[1, left_cells];
    cellval_u_1_right = Finch.variables[2].values[1, right_cells];
    (coef__u_1_left, coef__u_1_right) = Finch.FV_reconstruct_value_left_right(cellx_left, cellx_right, cellval_u_1_left, cellval_u_1_right, facex, limiter="vanleer");
    
    cellval_p_1_left = Finch.variables[3].values[1, left_cells];
    cellval_p_1_right = Finch.variables[3].values[1, right_cells];
    (coef__p_1_left, coef__p_1_right) = Finch.FV_reconstruct_value_left_right(cellx_left, cellx_right, cellval_p_1_left, cellval_p_1_right, facex, limiter="vanleer");
    
    coef__gamma_1 = 1.4;
    
    # compute left and right moving flux on both sides
    U_left = [coef__r_1_left, coef__u_1_left, coef__p_1_left];
    U_right = [coef__r_1_right, coef__u_1_right, coef__p_1_right];
    (Fm_left, Fp_left) = Finch.callback_functions[1].func(U_left);
    (Fm_right, Fp_right) = Finch.callback_functions[1].func(U_right);
    
    # Combine
    if els[1] == eid
        result = -(Fp_left .+ Fm_right);
    else # switch sign
        result = (Fp_left .+ Fm_right);
    end
    
    if els[1] == els[2] && normal[1] < 0
        result = -result;
    end
    
    return result;
end #rhs_surface_function_for_q3

# No assembly function set for q3

