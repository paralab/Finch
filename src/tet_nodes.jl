#=
Find nodes for tet elements.
=#
include("tet_quadrature_table.jl");

# Build global nodes for a triangular element with refel and vertices v
function tetrahedron_element_nodes(refel, v)
    return  tetrahedron_refel_to_xyz(refel.r[:,1], refel.r[:,2], refel.r[:,3], v);
end

# Set up refel nodal array
function tetrahedron_refel_nodes!(refel)
    (eqx, eqy, eqz) = tetrahedron_equilateral_nodes(refel.N);
    (r, s, t) = tetrahedron_equilateral_to_rst(eqx,eqy,eqz);
    
    # Adjust things within 1e-15 of the -1 boundaries to be exactly -1 and close to 0 to be 0
    tet_snap_to_bdry!(r, s, t, 1e-14);
    
    refel.r = zeros(refel.Np,3);
    refel.wr = zeros(refel.Np);
    
    refel.r[:,1] = r;
    refel.r[:,2] = s;
    refel.r[:,3] = t;
    
    # quadrature nodes/weights from a table
    xyzw = tetrahedron_quadrature_nodes_weights(refel.N+1);
    refel.g = xyzw[:,1:3];
    refel.wg = xyzw[:,4];
    
    # face node maps
    tol = 1e-12;
    tf1(x) = abs(x[1] + 1) < tol;
    tf2(x) = abs(x[2] + 1) < tol;
    tf3(x) = abs(x[3] + 1) < tol;
    tf4(x) = abs(x[1] + x[2] + x[3] + 1) < tol;
    refel.face2local = [get_face2local_map(refel.r, tf1),
                        get_face2local_map(refel.r, tf2),
                        get_face2local_map(refel.r, tf3),
                        get_face2local_map(refel.r, tf4)];
    
    # Surface quadrature nodes/weights are not ready. TODO
    if finch_state.config.solver_type == DG
        printerr("Surface quadrature for tets is not ready. Sorry.", fatal=true);
    end
end

# Purpose  : Compute (x,y,z) nodes in equilateral tetrahedron for polynomial of order N  
function tetrahedron_equilateral_nodes(N)
    # optimal alpha values for N up to 16 (from book)
    alpopt = [0, 0, 0, 0.1002, 1.1332, 1.5608, 1.3413, 1.2577, 1.1603, 1.10153, 0.6080, 0.4523, 0.8856, 0.8717, 0.9655];
            
    # Set optimized parameter, alpha, depending on order N
    if (N<16)
        alpha = alpopt[N];
    else
        alpha = 1;
    end

    # total number of nodes
    Np = Int64((N+1)*(N+2)*(N+3)/6);

    # Create equidistributed nodes on equilateral tetrahedron
    L1 = zeros(Np);
    L2 = zeros(Np);
    L3 = zeros(Np);
    L4 = zeros(Np);
    sk = 1;
    for n=1:(N+1)
        for m=1:(N+2-n)
            for q=1:(N+3-n-m)
                L1[sk] = (n-1)/N; # (1 + (-1 + (n-1)*2/N))/2
                L2[sk] = (m-1)/N; # (1 + (-1 + (m-1)*2/N))/2
                L3[sk] = 1 - (q+m+n-3)/N; # -(1 + -1 + (q-1)*2/N + -1 + (m-1)*2/N + -1 + (n-1)*2/N)/2
                L4[sk] = (q-1)/N; # (1 + (-1 + (q-1)*2/N))/2
                sk = sk+1;
            end
        end
    end
    
    # set vertices of tetrahedron
    invsqrt3 = 0.5773502691896258; # 1/sqrt(3)=0.5773502691896258
    invsqrt6 = 0.4082482904638631; # 1/sqrt(6)=0.4082482904638631
    v1 = [-1, -invsqrt3,   -invsqrt6]; 
    v2 = [ 1, -invsqrt3,   -invsqrt6]; 
    v3 = [ 0,  2*invsqrt3, -invsqrt6]; 
    v4 = [ 0,  0,          3*invsqrt6];

    # orthogonal axis tangents on faces 1-4
    t1 = zeros(4,3);
    t2 = zeros(4,3)
    t1[1,:] = v2.-v1;
    t1[2,:] = v2.-v1;
    t1[3,:] = v3.-v2;
    t1[4,:] = v3.-v1;
    t2[1,:] = v3.-0.5 .* (v1.+v2);
    t2[2,:] = v4.-0.5 .* (v1.+v2);
    t2[3,:] = v4.-0.5 .* (v2.+v3);
    t2[4,:] = v4.-0.5 .* (v1.+v3);
    
    for n=1:4 # normalize tangents 
        t1[n,:] = t1[n,:] ./ norm(t1[n,:]); 
        t2[n,:] = t2[n,:] ./ norm(t2[n,:]); 
    end 
    
    # Warp and blend for each face (accumulated in shiftXYZ)
    XYZ = zeros(Np,3);
    shift = zeros(Np,3);
    for i=1:Np
        for j=1:3
            XYZ[i,j] = L3[i]*v1[j] + L4[i]*v2[j] + L2[i]*v3[j] + L1[i]*v4[j]; # form undeformed coordinates
        end
    end
    
    for face=1:4 
        if face==1     La = L1; Lb = L2; Lc = L3; Ld = L4;
        elseif face==2 La = L2; Lb = L1; Lc = L3; Ld = L4;
        elseif face==3 La = L3; Lb = L1; Lc = L4; Ld = L2;
        elseif face==4 La = L4; Lb = L1; Lc = L3; Ld = L2;
        end
      
        # compute warp tangential to face
        (warp1, warp2) = tetrahedron_warpfactor(N, Np, alpha, Lb, Lc, Ld); 
        
        blend = Lb.*Lc.*Ld;   # compute volume blending
      
        denom = (Lb .+ .5 .* La).*(Lc .+ .5 .* La).*(Ld .+ .5 .* La);   # modify linear blend
        tol = 1e-10
        for ids=1:Np
            if denom[ids] > tol
                blend[ids] = (1 .+ (alpha.*La[ids]).^2).*blend[ids]./denom[ids];
            end
        end
        
        # compute warp & blend
        shift = shift .+ (blend.*warp1) * t1[face,:]' .+ (blend.*warp2) * t2[face,:]';
        
        # fix face warp 
        for ids=1:Np
            numtrue = 0;
            if Lb[ids]>tol
                numtrue += 1;
            end
            if Lc[ids]>tol
                numtrue += 1;
            end
            if Ld[ids]>tol
                numtrue += 1;
            end
            if La[ids]<tol && (numtrue < 3)
                shift[ids,:] = warp1[ids]*t1[face,:] + warp2[ids]*t2[face,:];
            end
        end
    end
    
    XYZ = XYZ .+ shift;
    
    return (XYZ[:,1], XYZ[:,2], XYZ[:,3]);
end

# Compute scaled warp function at order N based on rout interpolation nodes
function tetrahedron_warpfactor(N, Np, alpha, L1, L2, L3)
    # Compute LGL and equidistant node distribution
    (LGLr, w) = jacobi_LGL_quad(N);
    
    # 2) compute blending function at each node for each edge
    blend1 = L2.*L3; 
    blend2 = L1.*L3; 
    blend3 = L1.*L2;
    
    # 3) amount of warp for each node, for each edge
    warp1 = zeros(Np);
    warp2 = zeros(Np);
    warp3 = zeros(Np);
    
    dx = zeros(Np);
    dy = zeros(Np);
    
    xeq = zeros(N+1);
    for i=1:N+1
        xeq[i] = -1 + 2*(N+1-i)/N;
    end
    
    for wi=1:Np
        for i=1:N+1
            d1 = -LGLr[i] - xeq[i];
            d2 = -LGLr[i] - xeq[i];
            d3 = -LGLr[i] - xeq[i];
            for j=2:N
                if i!=j 
                    d1 = d1.*(L3[wi]-L2[wi]-xeq[j])/(xeq[i]-xeq[j]);
                    d2 = d2.*(L1[wi]-L3[wi]-xeq[j])/(xeq[i]-xeq[j]);
                    d3 = d3.*(L2[wi]-L1[wi]-xeq[j])/(xeq[i]-xeq[j]);
                end
            end
            
            if i!=1
                d1 = -d1/(xeq[i]-xeq[1]);
                d2 = -d2/(xeq[i]-xeq[1]);
                d3 = -d3/(xeq[i]-xeq[1]);
            end
    
            if i!=(N+1)
                d1 = d1/(xeq[i]-xeq[N+1]);
                d2 = d2/(xeq[i]-xeq[N+1]);
                d3 = d3/(xeq[i]-xeq[N+1]);
            end
    
            warp1[wi] = warp1[wi]+d1;
            warp2[wi] = warp2[wi]+d3;
            warp3[wi] = warp3[wi]+d3;
        end
        
        # 4) combine blend & warp
        warp1[wi] = 4 * blend1[wi] * warp1[wi] * (1 + (alpha * L1[wi])^2);
        warp2[wi] = 4 * blend2[wi] * warp2[wi] * (1 + (alpha * L2[wi])^2);
        warp3[wi] = 4 * blend3[wi] * warp3[wi] * (1 + (alpha * L3[wi])^2);
        
        # 5) evaluate shift in equilateral triangle
        # dx[wi] = 1*warp1[wi] + cos(2*pi/3)*warp2[wi] + cos(4*pi/3)*warp3[wi];
        # dy[wi] = 0*warp1[wi] + sin(2*pi/3)*warp2[wi] + sin(4*pi/3)*warp3[wi];
        dx[wi] = warp1[wi] - 0.5*warp2[wi] - 0.5*warp3[wi];
        dy[wi] = 0.8660254037844387*warp2[wi] - 0.8660254037844387*warp3[wi];
    end
    
    return (dx, dy);
end

# From (x,y,z) in equilateral tetrahedron to (r,s,t) coordinates in standard tetrahedron
function tetrahedron_equilateral_to_rst(x,y,z)
    invsqrt3 = 0.5773502691896258; # 1/sqrt(3)=0.5773502691896258
    invsqrt6 = 0.4082482904638631; # 1/sqrt(6)=0.4082482904638631
    v1 = [-1, -invsqrt3,   -invsqrt6]; 
    v2 = [ 1, -invsqrt3,   -invsqrt6]; 
    v3 = [ 0,  2*invsqrt3, -invsqrt6]; 
    v4 = [ 0,  0,          3*invsqrt6];

    # back out right tet nodes
    rhs = zeros(3,length(x));
    rhs[1,:] = x .- (0.5 * (v2[1]+v3[1]+v4[1]-v1[1]));
    rhs[2,:] = y .- (0.5 * (v2[2]+v3[2]+v4[2]-v1[2]));
    rhs[3,:] = z .- (0.5 * (v2[3]+v3[3]+v4[3]-v1[3]));
    
    A = [0.5 .* (v2 .- v1) 0.5 .* (v3 .- v1) 0.5 .* (v4 .- v1)];
    RST = A \ rhs;
    
    return (RST[1,:], RST[2,:], RST[3,:]);
end

# # From (r,s,t) coordinates in reference tetrahedron to (x,y,z) in tetrahedron with vertices v
# # v is a 3x4 array [x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4]
# function tetrahedron_refel_to_xyz(r, s, t, v)
#     d = v[:,1];
#     A = [v[:,2].-d v[:,3].-d v[:,4].-d];
    
#     np = length(r);
#     mv = zeros(3,np);
#     for i=1:np
#         tmp = [(r[i]+1)/2, (s[i]+1)/2, (t[i]+1)/2];
#         mv[:,i] = A*tmp + d;
#     end
    
#     x = mv[1,:]
#     y = mv[2,:]
#     z = mv[3,:]
    
#     return (x, y, z);
# end

# function get_face2local_map(r, compare)
#     n = size(r,1);
#     map = zeros(Int,n);
#     nf = 0;
#     for i=1:n
#         if compare(r[i,:])
#             nf = nf+1;
#             map[nf] = i;
#         end
#     end
    
#     return map[1:nf];
# end

function tet_snap_to_bdry!(r, s, t, tol)
    n = length(r);
    for i=1:n
        allowance = [1,1,1];
        if abs(r[i]+1) < tol
            r[i] = -1;
            allowance[1] = 0;
        end
        if abs(s[i]+1) < tol
            s[i] = -1;
            allowance[2] = 0;
        end
        if abs(t[i]+1) < tol
            t[i] = -1;
            allowance[3] = 0;
        end
        
        if abs(r[i] + s[i] + t[i] + 1) < tol
            offset = r[i] + s[i] + t[i] + 1;
            r[i] -= offset * allowance[1]/sum(allowance);
            s[i] -= offset * allowance[2]/sum(allowance);
            t[i] -= offset * allowance[3]/sum(allowance);
        end
    end
end