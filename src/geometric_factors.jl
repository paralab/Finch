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
function build_geometric_factors(refel, grid::Grid; do_face_detj::Bool=true, 
                                do_vol_area::Bool=false, constant_jacobian::Bool=false)
    ftype = finch_state.config.float_type;
    nel = size(grid.loc2glb, 2);
    totalfaces = size(grid.facenormals, 2);
    dim = size(grid.allnodes, 1);
    
    # Find maximal values for sizing arrays
    if typeof(refel) <: Array
        nodes_per_element = refel[1].Np;
        qnodes_per_element = refel[1].Nqp;
        for i=2:length(refel)
            nodes_per_element = max(nodes_per_element, refel[i].Np);
            qnodes_per_element = max(qnodes_per_element, refel[i].Nqp);
        end
        order = refel[1].N;
        
    else # one element type
        nodes_per_element = refel.Np;
        qnodes_per_element = refel.Nqp;
        order = refel.N;
        
        nvertex = size(grid.glbvertex, 1)
        ve = zeros(dim, nvertex);
        etype = 1;
        if dim == 1
            etype = 2;
        elseif dim == 2
            if refel.Nfaces == 3
                etype = 3;
            else
                etype = 4;
            end
        elseif dim == 3
            if refel.Nfaces == 4
                etype = 5;
            else
                etype = 6;
            end
        end
    end
    
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
        volume = zeros(ftype, 0);
        area = zeros(ftype, 0);
    end
    
    #
    xe = zeros(dim, nodes_per_element);
    ve = zeros(dim, size(grid.glbvertex, 1));
    
    # loop over elements to build their geo facs
    for e=1:nel
        # it's possible to have different element types with different refels
        if typeof(refel) <: Array
            refeli = refel[grid.refel_ind[e]];
            nodes_per_element = refeli.Np;
            qnodes_per_element = refeli.Nqp;
            nvertex = size(grid.glbvertex, 1)
            for i=1:nvertex
                if grid.glbvertex[i,e]==0
                    nvertex = i-1;
                    break;
                end
            end
            etype = 1;
            if dim == 1
                etype = 2;
            elseif dim == 2
                if refeli.Nfaces == 3
                    etype = 3;
                else
                    etype = 4;
                end
            elseif dim == 3
                if refeli.Nfaces == 4
                    etype = 5;
                else
                    etype = 6;
                end
            end
        else # only one element type
            refeli = refel;
        end
        
        # get node coords
        for i=1:nodes_per_element
            for j=1:dim
                xe[j,i] = grid.allnodes[j,grid.loc2glb[i,e]];
            end
        end
        for i=1:nvertex
            for j=1:dim
                ve[j,i] = grid.allnodes[j,grid.glbvertex[i,e]];
            end
        end
        
        (e_detJ, e_J) = geometric_factors(refeli, xe, constantJ=constant_jacobian);
        
        if constant_jacobian
            detJ[e] = e_detJ;
        else
            detJ[1:nodes_per_element,e] = e_detJ;
        end
        
        J[e] = e_J;
        
        if do_vol_area
            if constant_jacobian
                if etype == 4
                    volume[e] = 4 * e_detJ;
                elseif etype == 6
                    volume[e] = 8 * e_detJ;
                else
                    volume[e] = element_volume(etype, ve);
                end
            else
                volume[e] = element_volume(etype, ve);
            end
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
            face_type = 2;
            if dim == 3
                if refeli.Nfaces == 6 # hex
                    face_faces = 4; # faces are quads
                    face_const_J = constant_jacobian;
                    face_type = 4;
                elseif refeli.Nfaces == 4 # tet
                    face_faces = 3; # faces are triangles
                    face_type = 3;
                end
            end
            face_refel = build_refel(dim-1, order, face_faces, finch_state.config.elemental_nodes);
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
                    # area[fi] = (2 ^ (dim-1)) * face_detj[fi];
                    area[fi] = face_area(dim, face_type, xf);
                end
            end
        end
        
    else
        face_detj = zeros(ftype,0);
    end
    
    return GeometricFactors(J, detJ, volume, face_detj, area);
end

function geometric_factors(refel::Refel, pts::Matrix; constantJ::Bool=false, do_J::Bool=true)
    # pts = element node global coords
    # detJ = determinant(J)
    # J = Jacobian
    FT = finch_state.config.float_type;
    np = refel.Np;
    if refel.dim == 0
        detJ = [FT(1.0)];
        if do_J
            J = Jacobian([FT(1.0)],zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0));
        end
        
    elseif refel.dim == 1
        if np == 1
            # 0D face refels can only have 1 point
            detJ = [FT(1.0)];
            if do_J
                J = Jacobian([FT(1.0)],zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0));
            end
        else
            xr  = refel.Dg*pts[:];
            if constantJ
                detJ = FT(0.0);
                for i=1:np
                    detJ += refel.Dg[1,i] * pts[1,i];
                end
                if do_J
                    rx = [FT(1 / detJ)];
                end
                
            else
                detJ = zeros(FT, np);
                for j=1:np
                    for i=1:np
                        detJ[j] += refel.Dg[j,i] * pts[1,i];
                    end
                end
                if do_J
                    rx = FT(1 ./ detJ);
                end
            end
            if do_J
                J = Jacobian(rx,zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0));
            end
        end
        
    elseif refel.dim == 2
        # if refel.Nfaces == 3 # triangle
            if constantJ
                xr = FT(0.0); xs = FT(0.0); yr = FT(0.0); ys = FT(0.0);
                for i=1:np
                    xr += refel.Qr[1,i] * pts[1,i];
                    xs += refel.Qs[1,i] * pts[1,i];
                    yr += refel.Qr[1,i] * pts[2,i];
                    ys += refel.Qs[1,i] * pts[2,i];
                end
            else
                xr = zeros(FT, np);
                xs = zeros(FT, np);
                yr = zeros(FT, np);
                ys = zeros(FT, np);
                for j=1:np
                    for i=1:np
                        xr[j] += refel.Qr[j,i] * pts[1,i];
                        xs[j] += refel.Qs[j,i] * pts[1,i];
                        yr[j] += refel.Qr[j,i] * pts[2,i];
                        ys[j] += refel.Qs[j,i] * pts[2,i];
                    end
                end
            end
            # xr = refel.Qr*pts[1,:];
            # xs = refel.Qs*pts[1,:];
            # yr = refel.Qr*pts[2,:];
            # ys = refel.Qs*pts[2,:];
            
        # else # quad
        #     (xr, xs) = tensor_grad2(refel.Dg, pts[1,:][:]);
        #     (yr, ys) = tensor_grad2(refel.Dg, pts[2,:][:]);
        # end
        
        if constantJ
            detJ = FT(abs(-xs*yr + xr*ys));
            if do_J
                xr = [FT(xr)];
                xs = [FT(xs)];
                yr = [FT(yr)];
                ys = [FT(ys)];
            end
        else
            detJ = abs.(-xs.*yr + xr.*ys);
        end
        
        if do_J
            rx =  ys./detJ;
            sx = -yr./detJ;
            ry = -xs./detJ;
            sy =  xr./detJ;
            J = Jacobian(rx,ry,zeros(FT,0),sx,sy,zeros(FT,0),zeros(FT,0),zeros(FT,0),zeros(FT,0));
        end
        
    else
        # if refel.Nfaces == 4 # tetrahedron
            if constantJ
                xr = FT(0.0); xs = FT(0.0); xt = FT(0.0);
                yr = FT(0.0); ys = FT(0.0); yt = FT(0.0);
                zr = FT(0.0); zs = FT(0.0); zt = FT(0.0);
                for i=1:np
                    xr += refel.Qr[1,i] * pts[1,i];
                    xs += refel.Qs[1,i] * pts[1,i];
                    xt += refel.Qt[1,i] * pts[1,i];
                    yr += refel.Qr[1,i] * pts[2,i];
                    ys += refel.Qs[1,i] * pts[2,i];
                    yt += refel.Qt[1,i] * pts[2,i];
                    zr += refel.Qr[1,i] * pts[3,i];
                    zs += refel.Qs[1,i] * pts[3,i];
                    zt += refel.Qt[1,i] * pts[3,i];
                end
            else
                xr = zeros(FT, np);
                xs = zeros(FT, np);
                xt = zeros(FT, np);
                yr = zeros(FT, np);
                ys = zeros(FT, np);
                yt = zeros(FT, np);
                zr = zeros(FT, np);
                zs = zeros(FT, np);
                zt = zeros(FT, np);
                for j=1:np
                    for i=1:np
                        xr[j] += refel.Qr[j,i] * pts[1,i];
                        xs[j] += refel.Qs[j,i] * pts[1,i];
                        xt[j] += refel.Qt[j,i] * pts[1,i];
                        yr[j] += refel.Qr[j,i] * pts[2,i];
                        ys[j] += refel.Qs[j,i] * pts[2,i];
                        yt[j] += refel.Qt[j,i] * pts[2,i];
                        zr[j] += refel.Qr[j,i] * pts[3,i];
                        zs[j] += refel.Qs[j,i] * pts[3,i];
                        zt[j] += refel.Qt[j,i] * pts[3,i];
                    end
                end
            end
            # xr = refel.Qr*pts[1,:];
            # xs = refel.Qs*pts[1,:];
            # xt = refel.Qt*pts[1,:];
            # yr = refel.Qr*pts[2,:];
            # ys = refel.Qs*pts[2,:];
            # yt = refel.Qt*pts[2,:];
            # zr = refel.Qr*pts[3,:];
            # zs = refel.Qs*pts[3,:];
            # zt = refel.Qt*pts[3,:];
            
        # else # hexahedron
        #     (xr, xs, xt) = tensor_grad3(refel.Dg, pts[1,:][:]);
        #     (yr, ys, yt) = tensor_grad3(refel.Dg, pts[2,:][:]);
        #     (zr, zs, zt) = tensor_grad3(refel.Dg, pts[3,:][:]);
        # end
        
        if constantJ
            detJ = FT(abs(xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt)));
            if do_J
                xr = [FT(xr)];
                xs = [FT(xs)];
                xt = [FT(xt)];
                yr = [FT(yr)];
                ys = [FT(ys)];
                yt = [FT(yt)];
                zr = [FT(zr)];
                zs = [FT(zs)];
                zt = [FT(zt)];
            end
        else
            detJ = abs.(xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt));
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
function geometric_factors_face(refel::Refel, pts::Matrix{FT}, normal::Vector{FT}, constantJ::Bool=true) where FT<:AbstractFloat
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
        
        # if refel.Nfaces == 3 # triangle
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
            
        # else # quad
        #     (xr, xs) = tensor_grad2(refel.Dg, newpts[1,:][:]);
        #     (yr, ys) = tensor_grad2(refel.Dg, newpts[2,:][:]);
            
        #     if constantJ
        #         detJ = [-xs[1] * yr[1] + xr[1] * ys[1]];
        #     else
        #         detJ = -xs.*yr + xr.*ys;
        #     end
        # end
        
        
    end
    
    return detJ;
end

# Transforms a 3D surface into a 2D space with the first vertex at the origin
function flatten_face(normal::Vector{FT}, pts::Matrix{FT}) where FT<:AbstractFloat
    newpts = copy(pts);
    if abs(normal[1]) + abs(normal[2]) < 1e-15
        # The face is already in the x-y plane
        if normal[3] < 0
             newpts[1,:] .*= -1;
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
        # dxy = sqrt(normal[2]*normal[2] + normal[1]*normal[1]);
        dxy = 1.0;
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

# computes the volume of an element with the given vertex coords
# etype:
# 1=point, 2=line, 3=triangle, 4=quad, 5=tet, 6=hex, 7=prism
function element_volume(etype::Int, pts::Matrix{FT}) where FT<:AbstractFloat
    dim = size(pts,1);
    if etype == 1
        return 1.0; # this is meaningless
    elseif etype == 2
        if dim == 1
            return abs(pts[1,2]-pts[1,1]);
        else
            d1 = (pts[1,2] - pts[1,1]);
            d2 = (pts[2,2] - pts[2,1]);
            v = d1*d1 + d2*d2;
            if dim == 3
                d1 = (pts[3,2] - pts[3,1]);
                v += d1*d1;
            end
            return sqrt(v)
        end
        
    elseif etype == 3 # triangle
        if dim == 2
            # abs(Ax(By - Cy) + Bx(Cy - Ay) + Cx(Ay - By) )/2
            return 0.5 * abs(pts[1,1]*(pts[2,2]-pts[2,3]) + pts[1,2]*(pts[2,3]-pts[2,1]) + pts[1,3]*(pts[2,1]-pts[2,2]));
        else
            # |v1 X v2| / 2
            v1 = [pts[i,2]-pts[i,1] for i=1:3];
            v2 = [pts[i,3]-pts[i,1] for i=1:3];
            c1 = v1[2]*v2[3] - v1[3]*v2[2];
            c2 = v1[3]*v2[1] - v1[1]*v2[3];
            c3 = v1[1]*v2[2] - v1[2]*v2[1];
            return 0.5 * sqrt(c1*c1 + c2*c2 + c3*c3);
        end
        
    elseif etype == 4 # quad
        # treat as two triangles. Assume points 1 and 3 are diagonally oposite
        # Find diagonally opposite vertex from 1
        dx2 = pts[1,2]-pts[1,1];
        dy2 = pts[2,2]-pts[2,1];
        dx3 = pts[1,3]-pts[1,1];
        dy3 = pts[2,3]-pts[2,1];
        dx4 = pts[1,4]-pts[1,1];
        dy4 = pts[2,4]-pts[2,1];
        
        # Guess 3
        s = sqrt(dx3*dx3 + dy3*dy3);
        top2 = (-dx2*dy3/s + dy2*dx3/s) > 0;
        top4 = (-dx4*dy3/s + dy4*dx3/s) > 0;
        if top2 != top4
            # It was 3
            return element_volume(3, pts[:,[1,2,3]]) + element_volume(3, pts[:,[1,3,4]]);
        end
        
        # Guess 4
        s = sqrt(dx4*dx4 + dy4*dy4);
        top2 = (-dx2*dy4/s + dy2*dx4/s) > 0;
        top3 = (-dx3*dy4/s + dy3*dx4/s) > 0;
        if top2 != top3
            # It was 4
            return element_volume(3, pts[:,[1,2,4]]) + element_volume(3, pts[:,[1,3,4]]);
        end
        
        # It was 2
        return element_volume(3, pts[:,[1,2,3]]) + element_volume(3, pts[:,[1,2,4]]);
        
    elseif etype == 5 # tet
        # (a X b . c)/6
        a = [pts[i,2]-pts[i,1] for i=1:3];
        b = [pts[i,3]-pts[i,1] for i=1:3];
        c = [pts[i,4]-pts[i,1] for i=1:3];
        axb = [a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1]];
        axbdc = axb[1] * c[1] + axb[2] * c[2] + axb[3] * c[3];
        return abs(axbdc)/6;
        
    elseif etype == 6 # hex
        # [7-1, 2-1, 3-6] + [7-1, 5-1, 6-8] + [7-1, 4-1, 8-3]
        a = [pts[i,7]-pts[i,1] for i=1:3];
        b = [pts[i,2]-pts[i,1] for i=1:3];
        c = [pts[i,3]-pts[i,6] for i=1:3];
        bxc = [b[2]*c[3] - b[3]*c[2], b[3]*c[1] - b[1]*c[3], b[1]*c[2] - b[2]*c[1]];
        v = a[1] * bxc[1] + a[2] * bxc[2] + a[3] * bxc[3];
        
        b = [pts[i,5]-pts[i,1] for i=1:3];
        c = [pts[i,6]-pts[i,8] for i=1:3];
        bxc = [b[2]*c[3] - b[3]*c[2], b[3]*c[1] - b[1]*c[3], b[1]*c[2] - b[2]*c[1]];
        v += a[1] * bxc[1] + a[2] * bxc[2] + a[3] * bxc[3];
        
        b = [pts[i,4]-pts[i,1] for i=1:3];
        c = [pts[i,8]-pts[i,3] for i=1:3];
        bxc = [b[2]*c[3] - b[3]*c[2], b[3]*c[1] - b[1]*c[3], b[1]*c[2] - b[2]*c[1]];
        v += a[1] * bxc[1] + a[2] * bxc[2] + a[3] * bxc[3];
        
        return abs(v);
    end
    printerr("Can't compute volume for unknown element type: "*string(etype), fatal=true);
    return 1.0; 
end

# computes the area of a face with the given vertex coords
# face type:
# 1=point, 2=line, 3=triangle, 4=quad
function face_area(dim::Int, ftype::Int, pts::Matrix{FT}) where FT<:AbstractFloat
    if ftype == 1
        return 1.0; # this is meaningless
    elseif ftype == 2
        if dim == 2 # dim=1 doesn't make sense
            d1 = (pts[1,2] - pts[1,1]);
            d2 = (pts[2,2] - pts[2,1]);
            v = d1*d1 + d2*d2;
            return sqrt(v)
        elseif dim == 3 # This is an edge, not a face, but let's do it anyway
            d1 = (pts[1,2] - pts[1,1]);
            d2 = (pts[2,2] - pts[2,1]);
            d3 = (pts[3,2] - pts[3,1]);
            v = d1*d1 + d2*d2 + d3*d3;
            return sqrt(v)
        end
        
    elseif ftype == 3 # triangle
        if dim == 3
            # |v1 X v2| / 2
            v1 = [pts[i,2]-pts[i,1] for i=1:3];
            v2 = [pts[i,3]-pts[i,1] for i=1:3];
            c1 = v1[2]*v2[3] - v1[3]*v2[2];
            c2 = v1[3]*v2[1] - v1[1]*v2[3];
            c3 = v1[1]*v2[2] - v1[2]*v2[1];
            return 0.5 * sqrt(c1*c1 + c2*c2 + c3*c3);
        end
        
    elseif ftype == 4 # quad
        # treat as two triangles. Assume points 1 and 3 are diagonally oposite
        return face_area(dim, 3, pts[:,1:3]) + face_area(dim, 3, pts[:,[1,3,4]]);
    end
    printerr("Can't compute area for unknown face type: "*string(ftype), fatal=true);
    return 1.0; 
end

#############################################################################################################

# builds one derivative matrix in place
function build_derivative_matrix(refel::Refel, geofacs::GeometricFactors, direction::Int, eid::Int, type::Int, mat::Union{SubArray, Matrix{FT}}) where FT<:AbstractFloat
    N = size(mat,2);
    M = size(mat,1);
    J = geofacs.J[eid];
    const_j = (length(J.rx) == 1);
    if direction == 1
        Jr = J.rx;
        Js = J.sx;
        Jt = J.tx;
    elseif direction == 2
        Jr = J.ry;
        Js = J.sy;
        Jt = J.ty;
    else
        Jr = J.rz;
        Js = J.sz;
        Jt = J.tz;
    end
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
            mat .= J.rx[1] .* refel_dr;
        else
            for i=1:N
                for j=1:M
                    mat[j,i] = J.rx[j] * refel_dr[j,i];
                end
            end
        end
        
    elseif refel.dim == 2
        if const_j
            mat .= Jr[1] .* refel_dr .+ Js[1] .* refel_ds;
        else
            for i=1:N
                for j=1:M
                    mat[j,i] = Jr[j] * refel_dr[j,i] + Js[j] * refel_ds[j,i];
                end
            end
        end
        
    elseif refel.dim == 3
        if const_j
            mat .= Jr[1] .* refel_dr .+ Js[1] .* refel_ds .+ Jt[1] .* refel_dt;
        else
            for i=1:N # loop over columns
                for j=1:M # loop over rows
                    mat[j,i] = Jr[j] * refel_dr[j,i] + Js[j] * refel_ds[j,i] + Jt[j] * refel_dt[j,i];
                end
            end
        end
    end
    
    end # inbounds
end

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
