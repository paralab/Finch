# Geometric factors
export build_geometric_factors, geometric_factors, 
        build_derivative_matrix, build_deriv_matrix, build_face_deriv_matrix, get_quadrature_point_coords

include("tensor_ops.jl");

# #=
# # Stores the elemental jacobian
# # Used by the "geometric_factors" function
# =#
# struct Jacobian
#     rx::Vector
#     ry::Vector
#     rz::Vector
#     sx::Vector
#     sy::Vector
#     sz::Vector
#     tx::Vector
#     ty::Vector
#     tz::Vector
# end

# #=
# Stores all of the geometric factors for the grid.
# =#
# struct GeometricFactors
#     J::Vector{Jacobian}        # Jacobian for each element
#     detJ::Vector      # Determinant of Jacobian for each element
    
#     # These below are only computed if needed, otherwise empty arrays
#     volume::Vector    # Volume of each element (used by FV)
#     face_detJ::Array  # Determinant of Jacobian for each face (used by DG and FV)
#     area::Vector      # Area of each face (used by FV)
# end

import Base.copy
function copy(j::Jacobian)
    return Jacobian(copy(j.rx),copy(j.ry),copy(j.rz),copy(j.sx),copy(j.sy),copy(j.sz),copy(j.tx),copy(j.ty),copy(j.tz));
end

# Construct the geometric factors from a grid+refel
function build_geometric_factors(refel::Refel, grid::Grid; do_face_detj::Bool=false, 
                                do_vol_area::Bool=false, constant_jacobian::Bool=false)
    ftype = finch_state.config.float_type;
    nel = size(grid.loc2glb, 2);
    totalfaces = size(grid.facenormals, 2);
    dim = size(grid.allnodes, 1);
    nodes_per_element = refel.Np;
    qnodes_per_element = refel.Nqp;
    
    if constant_jacobian
        detJ = zeros(ftype, 1, nel);
    else
        detJ = zeros(ftype, qnodes_per_element, nel);
    end
    
    J = Vector{Jacobian{ftype}}(undef, nel);
    
    if do_vol_area
        volume = zeros(ftype, nel);
        area = zeros(ftype, totalfaces);
    else
        volume = zeros(0);
        area = zeros(0);
    end
    
    for e=1:nel
        glb = grid.loc2glb[:,e];
        xe = grid.allnodes[:,glb[:]];
        
        (e_detJ, e_J) = geometric_factors(refel, xe, constantJ=constant_jacobian);
        
        if constant_jacobian
            detJ[e] = e_detJ;
        else
            detJ[:,e] = e_detJ;
        end
        J[e] = e_J;
        
        if do_vol_area
            volume[e] = (2 ^ dim) * detJ[e];
        end
    end
    
    if do_face_detj
        if dim == 1
            # Faces in 1D are just points
            # detj = 1, area = 1
            face_detj = zeros(ftype, totalfaces);
            face_detj .= 1;
            if do_vol_area
                area .= 1;
            end
            
        else
            # Build a refel for the face
            face_faces = 2; # all 2d faces are lines
            face_const_J = true;
            if dim == 3
                if refel.Nfaces == 6 # hex
                    face_faces = 4; # faces are quads
                    face_const_J = constant_jacobian;
                elseif refel.Nfaces == 4 # tet
                    face_faces = 3; # faces are triangles
                end
            end
            face_refel = build_refel(dim-1, refel.N, face_faces, finch_state.config.elemental_nodes);
            face_detj = zeros(ftype, totalfaces);
            for fi=1:totalfaces
                xf = grid.allnodes[:,grid.face2glb[:,1,fi]];
                normal = grid.facenormals[:,fi];
                
                fdj = geometric_factors_face(face_refel, xf, normal, face_const_J);
                
                if face_const_J
                    face_detj[fi] = fdj[1];
                else
                    # I need to figure out what exactly to do in this situation
                    # the face has a non-constant jacobian
                    # For now I'll average it, but this is not correct.
                    face_detj[fi] = sum(fdj)/length(fdj);
                end
                
                if do_vol_area
                    area[fi] = (2 ^ (dim-1)) * face_detj[fi];
                end
            end
        end
        
    else
        face_detj = zeros(ftype,0);
    end
    
    return GeometricFactors(J, detJ, volume, face_detj, area);
end

function geometric_factors(refel::Refel, pts::Matrix{Float64}; constantJ::Bool=false, do_J::Bool=true)
    # pts = element node global coords
    # detJ = determinant(J)
    # J = Jacobian
    if refel.dim == 0
        detJ = [1];
        if do_J
            J = Jacobian([1.0],zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0));
        end
        
    elseif refel.dim == 1
        if length(pts) == 1
            # 0D face refels can only have 1 point
            detJ = [1];
            if do_J
                J = Jacobian([1],zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0));
            end
        else
            xr  = refel.Dg*pts[:];
            if constantJ
                detJ = xr[1];
                rx = [1 / detJ];
            else
                detJ = xr[:];
                rx = 1 ./ detJ;
            end
            if do_J
                J = Jacobian(rx,zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0));
            end
        end
        
    elseif refel.dim == 2
        if refel.Nfaces == 3 # triangle
            xr = refel.Qr*pts[1,:];
            xs = refel.Qs*pts[1,:];
            yr = refel.Qr*pts[2,:];
            ys = refel.Qs*pts[2,:];
            
        else # quad
            (xr, xs) = tensor_grad2(refel.Dg, pts[1,:][:]);
            (yr, ys) = tensor_grad2(refel.Dg, pts[2,:][:]);
        end
        
        if constantJ
            detJ = -xs[1]*yr[1] + xr[1]*ys[1];
            xr = [xr[1]];
            xs = [xs[1]];
            yr = [yr[1]];
            ys = [ys[1]];
        else
            detJ = -xs.*yr + xr.*ys;
        end
        
        if do_J
            rx =  ys./detJ;
            sx = -yr./detJ;
            ry = -xs./detJ;
            sy =  xr./detJ;
            J = Jacobian(rx,ry,zeros(0),sx,sy,zeros(0),zeros(0),zeros(0),zeros(0));
        end
        
    else
        if refel.Nfaces == 4 # tetrahedron
            xr = refel.Qr*pts[1,:];
            xs = refel.Qs*pts[1,:];
            xt = refel.Qt*pts[1,:];
            yr = refel.Qr*pts[2,:];
            ys = refel.Qs*pts[2,:];
            yt = refel.Qt*pts[2,:];
            zr = refel.Qr*pts[3,:];
            zs = refel.Qs*pts[3,:];
            zt = refel.Qt*pts[3,:];
            
        else # hexahedron
            (xr, xs, xt) = tensor_grad3(refel.Dg, pts[1,:][:]);
            (yr, ys, yt) = tensor_grad3(refel.Dg, pts[2,:][:]);
            (zr, zs, zt) = tensor_grad3(refel.Dg, pts[3,:][:]);
        end
        
        if constantJ
            detJ = xr[1]*(ys[1]*zt[1]-zs[1]*yt[1]) - yr[1]*(xs[1]*zt[1]-zs[1]*xt[1]) + zr[1]*(xs[1]*yt[1]-ys[1]*xt[1]);
            xr = [xr[1]];
            xs = [xs[1]];
            xt = [xt[1]];
            yr = [yr[1]];
            ys = [ys[1]];
            yt = [yt[1]];
            zr = [zr[1]];
            zs = [zs[1]];
            zt = [zt[1]];
        else
            detJ = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
        end
        
        if do_J
            rx =  (ys.*zt - zs.*yt)./detJ;
            ry = -(xs.*zt - zs.*xt)./detJ;
            rz =  (xs.*yt - ys.*xt)./detJ;
            
            sx = -(yr.*zt - zr.*yt)./detJ;
            sy =  (xr.*zt - zr.*xt)./detJ;
            sz = -(xr.*yt - yr.*xt)./detJ;
            
            tx =  (yr.*zs - zr.*ys)./detJ;
            ty = -(xr.*zs - zr.*xs)./detJ;
            tz =  (xr.*ys - yr.*xs)./detJ;
            J = Jacobian(rx,ry,rz,sx,sy,sz,tx,ty,tz);
        end
    end
    
    if do_J
        return (detJ,J);
    else
        return detJ;
    end
end

# NOTE: Only returns the Jacobian determinant. Use the volume Jacobian for derivatives.
# Also, this refel is a lower dimensional refel
function geometric_factors_face(refel::Refel, pts::Matrix{Float64}, normal::Vector{Float64}, constantJ::Bool=true)
    # pts = face node global coords
    if refel.dim == 0
        # A 1D face is just a point
        detJ = [1];
        
    elseif refel.dim == 1
        # A line. detJ = length/2
        if constantJ || size(pts,2)==2
            dx = pts[1,end] - pts[1,1];
            dy = pts[2,end] - pts[2,1];
            detJ = [0.5 * sqrt(dx*dx+dy*dy)];
        else
            # distance from point 1
            np = size(pts,2)
            distance = zeros(np);
            for i=2:np
                distance[i] = sqrt((pts[1,i]-pts[1,1]) * (pts[1,i]-pts[1,1]) + (pts[2,i]-pts[2,1]) * (pts[2,i]-pts[2,1]));
            end
            detJ = refel.Dg*distance;
        end
        
    elseif refel.dim == 2
        # rotate face ito x-y plane
        newpts = flatten_face(normal, pts);
        
        if refel.Nfaces == 3 # triangle
            np = size(pts,2);
            if constantJ
                xr = 0.0
                xs = 0.0
                yr = 0.0
                ys = 0.0
                for i=1:np
                    xr += refel.Qr[1,i]*newpts[1,i];
                    xs += refel.Qs[1,i]*newpts[1,i];
                    yr += refel.Qr[1,i]*newpts[2,i];
                    ys += refel.Qs[1,i]*newpts[2,i];
                end
                
                detJ = -xs * yr + xr * ys;
            else
                xr = refel.Qr*newpts[1,:];
                xs = refel.Qs*newpts[1,:];
                yr = refel.Qr*newpts[2,:];
                ys = refel.Qs*newpts[2,:];
                
                detJ = -xs.*yr + xr.*ys;
            end
            
        else # quad
            (xr, xs) = tensor_grad2(refel.Dg, newpts[1,:][:]);
            (yr, ys) = tensor_grad2(refel.Dg, newpts[2,:][:]);
            
            if constantJ
                detJ = [-xs[1] * yr[1] + xr[1] * ys[1]];
            else
                detJ = -xs.*yr + xr.*ys;
            end
        end
        
        
    end
    
    return detJ;
end

# Transforms a 3D surface into a 2D space with the first vertex at the origin
function flatten_face(normal::Vector{Float64}, pts::Matrix{Float64})
    newpts = copy(pts);
    if abs(normal[1]) + abs(normal[2]) < 1e-15
        # The face is already in the x-y plane
        if normal[3] < 0
             newpts[1,:] .*= 1;
        end
        
    else # need to transform
        np = size(pts,2);
        # move first point to origin
        for i=2:np
            newpts[1,i] -= newpts[1,1];
            newpts[2,i] -= newpts[2,1];
            newpts[3,i] -= newpts[3,1];
        end
        # find angles to rotate so that normal is along z axis
        # This is safe because we know normal is not along z axis now
        dxy = sqrt(normal[2]*normal[2] + normal[1]*normal[1]);
        theta = asin(normal[2] / dxy);
        phi = asin(1 / dxy);
        
        # Rotation matrix
        Rz = [cos(-theta) -sin(-theta) 0.0; sin(-theta) cos(-theta) 0.0 ; 0.0 0.0 1.0];
        Ry = [cos(phi) 0.0 sin(phi); 0.0 1.0 0.0; -sin(phi) 0.0 cos(phi)];
        Ryz = Ry * Rz;
        
        # Apply to all points
        for i=2:np
            newpts[:,i] = Ryz * newpts[:,i];
        end
    end
    
    return newpts[1:2,:];
end

#############################################################################################################

# builds one derivative matrix in place
function build_derivative_matrix(refel::Refel, geofacs::GeometricFactors, direction::Int, eid::Int, type::Int, mat::Matrix{Float64})
    N = size(mat,2);
    M = size(mat,1);
    J = geofacs.J[eid];
    const_j = (length(J.rx) == 1);
    if type == 0
        refel_dr = refel.Qr;
        refel_ds = refel.Qs;
        refel_dt = refel.Qt;
    else
        refel_dr = refel.Ddr;
        refel_ds = refel.Dds;
        refel_dt = refel.Ddt;
    end
    @inbounds begin
    if refel.dim == 1
        # Multiply rows of dr by J.rx
        if const_j
            for i=1:N # loop over columns
                for j=1:M # loop over rows
                    mat[j,i] = J.rx[1] * refel_dr[j,i];
                end
            end
        else
            for i=1:N # loop over columns
                for j=1:M # loop over rows
                    mat[j,i] = J.rx[j] * refel_dr[j,i];
                end
            end
        end
        
    elseif refel.dim == 2
        if direction == 1
            if const_j
                for i=1:N # loop over columns
                    for j=1:M # loop over rows
                        mat[j,i] = J.rx[1] * refel_dr[j,i] + J.sx[1] * refel_ds[j,i];
                    end
                end
            else
                for i=1:N # loop over columns
                    for j=1:M # loop over rows
                        mat[j,i] = J.rx[j] * refel_dr[j,i] + J.sx[j] * refel_ds[j,i];
                    end
                end
            end
            
        else
            if const_j
                for i=1:N # loop over columns
                    for j=1:M # loop over rows
                        mat[j,i] = J.ry[1] * refel_dr[j,i] + J.sy[1] * refel_ds[j,i];
                    end
                end
            else
                for i=1:N # loop over columns
                    for j=1:M # loop over rows
                        mat[j,i] = J.ry[j] * refel_dr[j,i] + J.sy[j] * refel_ds[j,i];
                    end
                end
            end
            
        end
        
    elseif refel.dim == 3
        if direction == 1
            if const_j
                for i=1:N # loop over columns
                    for j=1:M # loop over rows
                        mat[j,i] = J.rx[1] * refel_dr[j,i] + J.sx[1] * refel_ds[j,i] + J.tx[1] * refel_dt[j,i];
                    end
                end
            else
                for i=1:N # loop over columns
                    for j=1:M # loop over rows
                        mat[j,i] = J.rx[j] * refel_dr[j,i] + J.sx[j] * refel_ds[j,i] + J.tx[j] * refel_dt[j,i];
                    end
                end
            end
            
        elseif direction == 2
            if const_j
                for i=1:N # loop over columns
                    for j=1:M # loop over rows
                        mat[j,i] = J.ry[1] * refel_dr[j,i] + J.sy[1] * refel_ds[j,i] + J.ty[1] * refel_dt[j,i];
                    end
                end
            else
                for i=1:N # loop over columns
                    for j=1:M # loop over rows
                        mat[j,i] = J.ry[j] * refel_dr[j,i] + J.sy[j] * refel_ds[j,i] + J.ty[j] * refel_dt[j,i];
                    end
                end
            end
            
        else
            if const_j
                for i=1:N # loop over columns
                    for j=1:M # loop over rows
                        mat[j,i] = J.rz[1] * refel_dr[j,i] + J.sz[1] * refel_ds[j,i] + J.tz[1] * refel_dt[j,i];
                    end
                end
            else
                for i=1:N # loop over columns
                    for j=1:M # loop over rows
                        mat[j,i] = J.rz[j] * refel_dr[j,i] + J.sz[j] * refel_ds[j,i] + J.tz[j] * refel_dt[j,i];
                    end
                end
            end
        end
    end
    
    end # inbounds
end

# function build_deriv_matrix(refel, J)
#     if refel.dim == 1
#         RQ1 = zeros(size(refel.Q));
#         RD1 = zeros(size(refel.Q));
#         # Multiply rows of Qr and Ddr by J.rx
#         for i=1:size(RQ1,2) # loop over columns
#             RQ1[:,i] = J.rx .* refel.Qr[:,i];
#             RD1[:,i] = J.rx .* refel.Ddr[:,i];
#         end
#         return (RQ1,RD1);
        
#     elseif refel.dim == 2
#         RQ1 = zeros(size(refel.Q));
#         RQ2 = zeros(size(refel.Q));
#         RD1 = zeros(size(refel.Q));
#         RD2 = zeros(size(refel.Q));
#         for i=1:size(RQ1,2)
#             RQ1[:,i] = J.rx .* refel.Qr[:,i] + J.sx .* refel.Qs[:,i];
#             RQ2[:,i] = J.ry .* refel.Qr[:,i] + J.sy .* refel.Qs[:,i];
#             RD1[:,i] = J.rx .* refel.Ddr[:,i] + J.sx .* refel.Dds[:,i];
#             RD2[:,i] = J.ry .* refel.Ddr[:,i] + J.sy .* refel.Dds[:,i];
#         end
#         return (RQ1, RQ2, RD1, RD2);
        
#     elseif refel.dim == 3
#         RQ1 = zeros(size(refel.Q));
#         RQ2 = zeros(size(refel.Q));
#         RQ3 = zeros(size(refel.Q));
#         RD1 = zeros(size(refel.Q));
#         RD2 = zeros(size(refel.Q));
#         RD3 = zeros(size(refel.Q));
#         for i=1:size(RQ1,2)
#             RQ1[:,i] = J.rx .* refel.Qr[:,i] + J.sx .* refel.Qs[:,i] + J.tx .* refel.Qt[:,i];
#             RQ2[:,i] = J.ry .* refel.Qr[:,i] + J.sy .* refel.Qs[:,i] + J.ty .* refel.Qt[:,i];
#             RQ3[:,i] = J.rz .* refel.Qr[:,i] + J.sz .* refel.Qs[:,i] + J.tz .* refel.Qt[:,i];
#             RD1[:,i] = J.rx .* refel.Ddr[:,i] + J.sx .* refel.Dds[:,i] + J.tx .* refel.Ddt[:,i];
#             RD2[:,i] = J.ry .* refel.Ddr[:,i] + J.sy .* refel.Dds[:,i] + J.ty .* refel.Ddt[:,i];
#             RD3[:,i] = J.rz .* refel.Ddr[:,i] + J.sz .* refel.Dds[:,i] + J.tz .* refel.Ddt[:,i];
#         end
#         return (RQ1, RQ2, RQ3, RD1, RD2, RD3);
#     end
# end

# Build the regular deriv matrices, then extract the relevant face parts
function build_face_deriv_matrix(refel, face, J, full = false)
    if refel.dim == 1
        RQ1 = J.rx[1]*refel.surf_Qr[face];
        if full
            RD1 = J.rx[1]*refel.Ddr;
        else
            RD1 = J.rx[1]*refel.surf_Ddr[face];
        end
        
        return (RQ1,RD1);
        
    elseif refel.dim == 2
        RQ1 = J.rx[1]*refel.surf_Qr[face] + J.sx[1]*refel.surf_Qs[face];
        RQ2 = J.ry[1]*refel.surf_Qr[face] + J.sy[1]*refel.surf_Qs[face];
        if full
            RD1 = J.rx[1]*refel.Ddr + J.sx[1]*refel.Dds;
            RD2 = J.ry[1]*refel.Ddr + J.sy[1]*refel.Dds;
        else
            RD1 = J.rx[1]*refel.surf_Ddr[face] + J.sx[1]*refel.surf_Dds[face];
            RD2 = J.ry[1]*refel.surf_Ddr[face] + J.sy[1]*refel.surf_Dds[face];
        end
        
        return (RQ1, RQ2, RD1, RD2);
        
    elseif refel.dim == 3
        RQ1 = J.rx[1]*refel.surf_Qr[face] + J.sx[1]*refel.surf_Qs[face] + J.tx[1]*refel.surf_Qt[face];
        RQ2 = J.ry[1]*refel.surf_Qr[face] + J.sy[1]*refel.surf_Qs[face] + J.ty[1]*refel.surf_Qt[face];
        RQ3 = J.ry[1]*refel.surf_Qr[face] + J.sy[1]*refel.surf_Qs[face] + J.tz[1]*refel.surf_Qt[face];
        if full
            RD1 = J.rx[1]*refel.Ddr + J.sx[1]*refel.Dds + J.tx[1]*refel.Ddt;
            RD2 = J.ry[1]*refel.Ddr + J.sy[1]*refel.Dds + J.ty[1]*refel.Ddt;
            RD3 = J.ry[1]*refel.Ddr + J.sy[1]*refel.Dds + J.tz[1]*refel.Ddt;
        else
            RD1 = J.rx[1]*refel.surf_Ddr[face] + J.sx[1]*refel.surf_Dds[face] + J.tx[1]*refel.surf_Ddt[face];
            RD2 = J.ry[1]*refel.surf_Ddr[face] + J.sy[1]*refel.surf_Dds[face] + J.ty[1]*refel.surf_Ddt[face];
            RD3 = J.ry[1]*refel.surf_Ddr[face] + J.sy[1]*refel.surf_Dds[face] + J.tz[1]*refel.surf_Ddt[face];
        end
        
        return (RQ1, RQ2, RQ3, RD1, RD2, RD3);
    end
end

###################################################################################################################################
function geometric_factors_cachesim(refel, pts)
    # pts = element node global coords
    # J = detJ
    # D = Jacobian
    if refel.dim == 1
        xr  = refel.Dg*pts[:];
        J = xr[:];
        rx = 1 ./ J;
        D = Jacobian(rx,zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0));
        cachesim_load_range(10);
        cachesim_store_range(13);
        cachesim_store_range(16);
        
    elseif refel.dim == 2
        (xr, xs) = tensor_grad2(refel.Dg, pts[1,:][:]);
        (yr, ys) = tensor_grad2(refel.Dg, pts[2,:][:]);
        J = -xs.*yr + xr.*ys;
        
        rx =  ys./J;
        sx = -yr./J;
        ry = -xs./J;
        sy =  xr./J;
        D = Jacobian(rx,ry,zeros(0),sx,sy,zeros(0),zeros(0),zeros(0),zeros(0));
        
        cachesim_load_range(10);
        cachesim_load_range(11);
        cachesim_store_range(13);
        cachesim_store_range(14);
        cachesim_store_range(16);
        
    else
        (xr, xs, xt) = tensor_grad3(refel.Dg, pts[1,:][:]);
        (yr, ys, yt) = tensor_grad3(refel.Dg, pts[2,:][:]);
        (zr, zs, zt) = tensor_grad3(refel.Dg, pts[3,:][:]);
        J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
        
        rx =  (ys.*zt - zs.*yt)./J;
        ry = -(xs.*zt - zs.*xt)./J;
        rz =  (xs.*yt - ys.*xt)./J;
        
        sx = -(yr.*zt - zr.*yt)./J;
        sy =  (xr.*zt - zr.*xt)./J;
        sz = -(xr.*yt - yr.*xt)./J;
        
        tx =  (yr.*zs - zr.*ys)./J;
        ty = -(xr.*zs - zr.*xs)./J;
        tz =  (xr.*ys - yr.*xs)./J;
        D = Jacobian(rx,ry,rz,sx,sy,sz,tx,ty,tz);
        
        cachesim_load_range(10);
        cachesim_load_range(11);
        cachesim_load_range(12);
        cachesim_store_range(13);
        cachesim_store_range(14);
        cachesim_store_range(15);
        cachesim_store_range(16);
    end
    
    return (J,D);
end

function build_deriv_matrix_cachesim(refel, J)
    if refel.dim == 1
        RQ1 = zeros(size(refel.Q));
        RD1 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j];
                RD1[i,j] = J.rx[i]*refel.Ddr[i,j];
            end
        end
        return (RQ1,RD1);
        
        cachesim_load_range(13);
        cachesim_load_range(7);
        cachesim_load_range(10);
        
    elseif refel.dim == 2
        RQ1 = zeros(size(refel.Q));
        RQ2 = zeros(size(refel.Q));
        RD1 = zeros(size(refel.Q));
        RD2 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j] + J.sx[i]*refel.Qs[i,j];
                RQ2[i,j] = J.ry[i]*refel.Qr[i,j] + J.sy[i]*refel.Qs[i,j];
                RD1[i,j] = J.rx[i]*refel.Ddr[i,j] + J.sx[i]*refel.Dds[i,j];
                RD2[i,j] = J.ry[i]*refel.Ddr[i,j] + J.sy[i]*refel.Dds[i,j];
            end
        end
        return (RQ1, RQ2, RD1, RD2);
        
        cachesim_load_range(7);
        cachesim_load_range(8);
        cachesim_load_range(9);
        cachesim_load_range(10);
        cachesim_load_range(11);
        cachesim_load_range(12);
        cachesim_load_range(13);
        cachesim_load_range(14);
        cachesim_load_range(15);
        
    elseif refel.dim == 3
        RQ1 = zeros(size(refel.Q));
        RQ2 = zeros(size(refel.Q));
        RQ3 = zeros(size(refel.Q));
        RD1 = zeros(size(refel.Q));
        RD2 = zeros(size(refel.Q));
        RD3 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j] + J.sx[i]*refel.Qs[i,j] + J.tx[i]*refel.Qt[i,j];
                RQ2[i,j] = J.ry[i]*refel.Qr[i,j] + J.sy[i]*refel.Qs[i,j] + J.ty[i]*refel.Qt[i,j];
                RQ3[i,j] = J.rz[i]*refel.Qr[i,j] + J.sz[i]*refel.Qs[i,j] + J.tz[i]*refel.Qt[i,j];
                RD1[i,j] = J.rx[i]*refel.Ddr[i,j] + J.sx[i]*refel.Dds[i,j] + J.tx[i]*refel.Ddt[i,j];
                RD2[i,j] = J.ry[i]*refel.Ddr[i,j] + J.sy[i]*refel.Dds[i,j] + J.ty[i]*refel.Ddt[i,j];
                RD3[i,j] = J.rz[i]*refel.Ddr[i,j] + J.sz[i]*refel.Dds[i,j] + J.tz[i]*refel.Ddt[i,j];
            end
        end
        return (RQ1, RQ2, RQ3, RD1, RD2, RD3);
    end
end
