# Geometric factors
export GeometricFactors, Jacobian
export build_geometric_factors, geometric_factors, geometric_factors_face, 
        build_derivative_matrix, build_deriv_matrix, build_face_deriv_matrix, get_quadrature_point_coords

include("tensor_ops.jl");

#=
# Stores the elemental jacobian
# Used by the "geometric_factors" function
=#
struct Jacobian
    rx::Vector
    ry::Vector
    rz::Vector
    sx::Vector
    sy::Vector
    sz::Vector
    tx::Vector
    ty::Vector
    tz::Vector
end

#=
Stores all of the geometric factors for the grid.
=#
struct GeometricFactors
    J::Array{Jacobian,1}        # Jacobian for each element
    detJ::Array      # Determinant of Jacobian for each element
    
    # These below are only computed if needed, otherwise empty arrays
    volume::Vector    # Volume of each element (used by FV)
    face_detJ::Array  # Determinant of Jacobian for each face (used by DG and FV)
    area::Vector      # Area of each face (used by FV)
end

import Base.copy
function copy(j::Jacobian)
    return Jacobian(copy(j.rx),copy(j.ry),copy(j.rz),copy(j.sx),copy(j.sy),copy(j.sz),copy(j.tx),copy(j.ty),copy(j.tz));
end

# Construct the geometric factors from a grid+refel
function build_geometric_factors(refel, grid; do_face_detj=false, do_vol_area=false, constant_jacobian=false)
    nel = size(grid.loc2glb, 2);
    totalfaces = size(grid.facenormals, 2);
    dim = size(grid.allnodes, 1);
    nodes_per_element = refel.Np;
    qnodes_per_element = refel.Nqp;
    
    if constant_jacobian
        detJ = zeros(config.float_type, 1, nel);
    else
        detJ = zeros(config.float_type, qnodes_per_element, nel);
    end
    
    J = Array{Jacobian,1}(undef, nel);
    
    if do_vol_area
        volume = zeros(config.float_type, nel);
        area = zeros(config.float_type, totalfaces);
    else
        volume = [];
        area = [];
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
        face_detj = zeros(config.float_type, totalfaces);
        for fi=1:totalfaces
            xf = grid.allnodes[:,grid.face2glb[:,1,fi]];
            
            (fdj, fj) = geometric_factors_face(refel, grid.faceRefelInd[1,fi], xf);
            
            face_detj[fi] = fdj[1];
            
            if do_vol_area
                area[fi] = (2 ^ (dim-1)) * fdj[1];
            end
        end
    else
        face_detj = [];
    end
    
    return GeometricFactors(J, detJ, volume, face_detj, area);
end

function geometric_factors(refel, pts; constantJ = false)
    # pts = element node global coords
    # detJ = determinant(J)
    # J = Jacobian
    if refel.dim == 0
        detJ = [1];
        J = Jacobian([1],[],[],[],[],[],[],[],[]);
        
    elseif refel.dim == 1
        if length(pts) == 1
            # 0D face refels can only have 1 point
            detJ = [1];
            J = Jacobian([1],[],[],[],[],[],[],[],[]);
        else
            xr  = refel.Dg*pts[:];
            if constantJ
                detJ = xr[1];
                rx = [1 / detJ];
            else
                detJ = xr[:];
                rx = 1 ./ detJ;
            end
            J = Jacobian(rx,[],[],[],[],[],[],[],[]);
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
        
        rx =  ys./detJ;
        sx = -yr./detJ;
        ry = -xs./detJ;
        sy =  xr./detJ;
        J = Jacobian(rx,ry,[],sx,sy,[],[],[],[]);
        
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
    
    return (detJ,J);
end

# NOTE: Only returns the Jacobian determinant. Use the volume Jacobian for derivatives.
function geometric_factors_face(refel, face, pts)
    # pts = face node global coords
    if refel.dim == 0
        detJ = [1];
        
    elseif refel.dim == 1
        # 1D face can only have 1 point
        detJ = [1];
        
    elseif refel.dim == 2
        dx = pts[1,end] - pts[1,1];
        dy = pts[2,end] - pts[2,1];
        detJ = 0.5 * sqrt(dx*dx+dy*dy);
        
    else
        # TODO this assumes rectangular faces
        println("Warning: face geometric factors assume regular hex mesh. See geometric_factors_face() in geometric_factors.jl")
        dx = abs(pts[1,end] - pts[1,1]);
        dy = abs(pts[2,end] - pts[2,1]);
        dz = abs(pts[3,end] - pts[3,1]);
        # TODO this assumes hex aligned with axes
        d1 = max(dx, dy);
        d2 = max(min(dx, dy), dz);
        
        detJ = 0.25 * (d1*d2);
    end
    
    return (detJ, 0);
end

# builds one derivative matrix in place
function build_derivative_matrix(refel::Refel, geofacs, direction::Int, eid::Int, type::Int, mat::Matrix)
    @timeit timer_output "deriv-mat" begin
    N = size(mat,2);
    M = size(mat,1);
    J = geofacs.J[eid];
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
        for i=1:N # loop over columns
            mat[1:M,i] .= J.rx .* refel_dr[1:M,i];
        end
        
    elseif refel.dim == 2
        if direction == 1
            for i=1:N # loop over columns
                mat[1:M,i] .= J.rx .* refel_dr[1:M,i] .+ J.sx .* refel_ds[1:M,i];
            end
        else
            for i=1:N # loop over columns
                mat[1:M,i] .= J.ry .* refel_dr[1:M,i] .+ J.sy .* refel_ds[1:M,i];
            end
        end
        
    elseif refel.dim == 3
        if direction == 1
            for i=1:N # loop over columns
                mat[1:M,i] .= J.rx .* refel_dr[1:M,i] .+ J.sx .* refel_ds[1:M,i] .+ J.tx .* refel_dt[1:M,i];
            end
        elseif direction == 2
            for i=1:N # loop over columns
                mat[1:M,i] .= J.ry .* refel_dr[1:M,i] .+ J.sy .* refel_ds[1:M,i] .+ J.ty .* refel_dt[1:M,i];
            end
        else
            for i=1:N # loop over columns
                mat[1:M,i] .= J.rz .* refel_dr[1:M,i] .+ J.sz .* refel_ds[1:M,i] .+ J.tz .* refel_dt[1:M,i];
            end
        end
    end
    
    end # inbounds
    end # timer
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
        D = Jacobian(rx,[],[],[],[],[],[],[],[]);
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
        D = Jacobian(rx,ry,[],sx,sy,[],[],[],[]);
        
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
