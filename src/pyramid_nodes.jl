#=
Find nodes for tet elements.
=#
include("pyramid_quadrature_table.jl");

# Set up refel nodal array
function pyramid_refel_nodes!(refel)

    beta = 0
    # target nodes on surface - default to Nodes3D optimized blend for
    # conformity with tet nodes
    (mapr, maps, mapt) = pyramidSurfaceNodes3D( refel.N );

    # find coefficients mapping from equispaced to warped surface
    (rbc, sbc, tbc) = pyramidSurfaceEquiNodes3D( refel.N );
    (Vbc, v_ids, etri_ids, equad_ids, ftri_ids, fquad_ids) = JVandermonde3D( refel.N, rbc, sbc, tbc, beta );

    ids = vcat( v_ids, etri_ids, equad_ids, ftri_ids, fquad_ids);
    maprst = [mapr maps mapt];

    println( typeof( Vbc[ :, ids ] ) )
    println( size( Vbc[:, ids] ) )
    println( size( maprst ) )
    mapcrst = Vbc[ :, ids] \ maprst;

    # evaluate map at equispaced volumed nodes
    (req, seq, teq) = pyramidEquiNodes3D( refel.N );
    (Veq, v_ids, etri_ids, equad_ids, ftri_ids, fquad_ids) = JVandermonde3D(refel.N, req, seq, teq, beta);

    ids = vcat(v_ids, etri_ids, equad_ids, ftri_ids, fquad_ids);
    rst = Veq[ :, ids ] * mapcrst; # final coordinates of all nodes in triangle

    r = rst[ :, 1 ]; s = rst[ :, 2 ]; t = rst[ :, 3 ];
    
    refel.r = zeros( refel.Np, 3 );
    refel.wr = zeros( refel.Np );
    
    refel.r[:,1] = r;
    refel.r[:,2] = s;
    refel.r[:,3] = t;
    
    # quadrature nodes/weights from a table
    xyzw = pyramid_quadrature_nodes_weights( refel.N + 1 );
    refel.g = xyzw[ :, 1:3 ];
    refel.wg = xyzw[ :, 4 ];
    
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
        printerr("Surface quadrature for pyramids is not ready. Sorry.", fatal=true);
    end
end


# returns pyramid surface nodes for conformity with a tet/hex element.
# Warp and Blend nodes on the triangular faces, GLL tensor product nodes on
# the quadrilateral base. 

function pyramidSurfaceNodes3D( N )

    alphastore = [0, 0, 0, 0.1002, 1.1332, 1.560, 1.3413, 1.2577, 1.1603, 1.10153, 0.6080, 0.4523, 0.8856, 0.8717, 0.9655];
    alphaS = alphastore[ N ];  # surface blending 

    # build surface nodes for a degree 5 pyramid trace
    # [-1,1]x[-1,1]x[0,1] with top vertex at (0,0,1)

    (x, y, z) = Nodes3D(N, alphaS);
    (r, s, t) = xyztorst(x, y, z);
    inds = findall( tVal -> ( abs( tVal + 1 ) < 1e-8 ), t );
    r = r[ inds ];
    s = s[ inds ];

    Np = length(r);

    (gr, vals) = jacobi_LGL_quad( N );

    # (r2d, s2d) = meshgrid(gr);
    println( size(gr) )
    println(typeof(gr))
    r2d = gr' .* ones(N + 1);
    s2d = ones(1, N + 1) .* gr;
    t2d = zeros(N + 1, N + 1);

    r2d = reshape( r2d, (N + 1) * (N + 1), 1);
    s2d = reshape( s2d, (N + 1) * (N + 1), 1);
    t2d = reshape( t2d, (N + 1) * (N + 1), 1);

    I = ones(Np, 1);

    # Face 2
    x1 = -1; y1 = -1; z1 = +0;
    x2 =  1; y2 = -1; z2 = +0;
    x3 =  0; y3 =  0; z3 = +1;
    x3d = [ r2d; 0.5*(-(r + s)*x1 + (1 .+ r)*x2 + (1 .+ s)*x3) ];
    y3d = [ s2d; 0.5*(-(r + s)*y1 + (1 .+ r)*y2 + (1 .+ s)*y3) ];
    z3d = [ t2d; 0.5*(-(r + s)*z1 + (1 .+ r)*z2 + (1 .+ s)*z3) ];

    # Face 3
    x1 =  1; y1 = -1; z1 = +0;
    x2 =  1; y2 = +1; z2 = +0;
    x3 =  0; y3 =  0; z3 = +1;
    x3d = [ x3d; 0.5*(-(r + s)*x1 + (1 .+ r)*x2 + (1 .+ s)*x3) ];
    y3d = [ y3d; 0.5*(-(r + s)*y1 + (1 .+ r)*y2 + (1 .+ s)*y3) ];
    z3d = [ z3d; 0.5*(-(r + s)*z1 + (1 .+ r)*z2 + (1 .+ s)*z3) ];

    # Face 4
    x1 =  1; y1 = +1; z1 = +0;
    x2 = -1; y2 = +1; z2 = +0;
    x3 =  0; y3 =  0; z3 = +1;
    x3d = [ x3d; 0.5*(-( r + s )*x1 + ( 1 .+ r )*x2 + ( 1 .+ s )*x3) ];
    y3d = [ y3d; 0.5*(-( r + s )*y1 + ( 1 .+ r )*y2 + ( 1 .+ s )*y3) ];
    z3d = [ z3d; 0.5*(-( r + s )*z1 + ( 1 .+ r )*z2 + ( 1 .+ s )*z3) ];

    # Face 5
    x1 =  -1; y1 = -1; z1 = +0;
    x2 =  -1; y2 = +1; z2 = +0;
    x3 =  0; y3 =  0; z3 = +1;
    x3d = [ x3d; 0.5*(-( r + s )*x1 + ( 1 .+ r )*x2 + ( 1 .+ s )*x3) ];
    y3d = [ y3d; 0.5*(-( r + s )*y1 + ( 1 .+ r )*y2 + ( 1 .+ s )*y3) ];
    z3d = [ z3d; 0.5*(-( r + s )*z1 + ( 1 .+ r )*z2 + ( 1 .+ s )*z3) ];

    # remove duplicates
    tol = 1e-9;
    cnt = 1;
    x = x3d[1]; y = y3d[1]; z = z3d[1];
    for n = 2:length(x3d)
        d = (x3d[n] .- x).^2 + (y3d[n] .- y).^2 + (z3d[n] .- z).^2;
        
        if( minimum( d ) > tol )
            x = [ x; x3d[ n ] ];
            y = [ y; y3d[ n ] ];
            z = [ z; z3d[ n ] ];
        end
    end

    return (x, y, z)
end

function Nodes3D( p, alphain = nothing )

    # function [X,Y,Z] = Nodes3D(p)
    # Purpose: compute Warp & Blend nodes
    #  input:    p=polynomial order of interpolant
    #  output: X,Y,Z vectors of node coordinates in equilateral tetrahedron
    
    # choose optimized blending parameter
    alphastore = [0;0;0;0.1002; 1.1332;1.5608;1.3413;1.2577;1.1603;1.10153;0.6080;0.4523;0.8856;0.8717;0.9655];
    if alphain == nothing
        if p <= 15 
            alpha = alphastore[p];
        else
            alpha = 1.;
        end
    else
        alpha = alphain;
    end
    
    # total number of nodes and tolerance
    N = (Int)( (p + 1)*(p + 2)*(p + 3)/6 ); tol = 1e-10;
    
    ( r, s, t ) = EquiNodes3D(p); # create equidistributed nodes
    L1 = ( 1 .+ t )/2; L2 = ( 1 .+ s )/2; L3 = -( 1 .+ r .+ s .+ t )/2; L4 = ( 1 .+ r )/2;
    
    # set vertices of tetrahedron
    v1 = [-1 -1/sqrt(3) -1/sqrt(6)]; v2 = [ 1 -1/sqrt(3) -1/sqrt(6)];
    v3 = [ 0  2/sqrt(3) -1/sqrt(6)]; v4 = [ 0  0         3/sqrt(6)];
    
    # orthogonal axis tangents on faces 1-4
    t1 = zeros( 4, 3 );
    t2 = zeros( 4, 3 );
    t3 = zeros( 4, 3 );
    t4 = zeros( 4, 3 );

    t1[ 1, : ] = v2 - v1;              t1[ 2, : ] = v2 - v1;
    t1[ 3, : ] = v3 - v2;              t1[ 4, : ] = v3 - v1;
    t2[ 1, : ] = v3 - 0.5 * (v1 + v2); t2[ 2, : ] = v4 - 0.5 * (v1 + v2);
    t2[ 3, : ] = v4 - 0.5 * (v2 + v3); t2[ 4, : ] = v4 - 0.5 * (v1 + v3);
    
    for n = 1:4 # normalize tangents
        t1[ n, : ] = t1[ n, : ] / norm( t1[ n, : ] ); t2[ n, : ] = t2[ n, : ] / norm( t2[ n, : ] );
    end
    
    # Warp and blend for each face (accumulated in shiftXYZ)
    XYZ = L3 * v1 + L4 * v2 + L2 * v3 + L1 * v4; # form undeformed coordinates
    shift = zeros( size(XYZ) );
    for face = 1:4
        if face == 1 
            La = L1; Lb = L2; Lc = L3; Ld = L4; 
        end
        if face == 2 
            La = L2; Lb = L1; Lc = L3; Ld = L4; 
        end
        if face == 3 
            La = L3; Lb = L1; Lc = L4; Ld = L2; 
        end
        if face == 4 
            La = L4; Lb = L1; Lc = L3; Ld = L2; 
        end;
        
        # println( size(La) )
        # println( size(Lb) )

        # compute warp tangential to face
        (warp1, warp2) = WarpShiftFace3D(p, alpha, alpha, La, Lb, Lc, Ld);
        
        blend = Lb.*Lc.*Ld;   # compute volume blending
        
        denom = ( Lb .+ .5*La ).*( Lc .+ .5*La ).*( Ld .+ .5*La );   # modify linear blend
        ids = findall( denomVal -> ( denomVal > tol ), denom );
        blend[ids] = (1 .+ sign(alpha) * ( alpha .* La[ids] ).^2 ).*blend[ids] ./ denom[ids];
        
        # compute warp & blend
        # println( size( blend ) );
        # println( size( warp1 ) );
        # println( size(t1) );

        shift = shift .+ (blend.*warp1)*t1[ face:face, : ] .+ (blend.*warp2)*t2[ face:face, : ];
        
        # fix face warp
        ids = []
        for idx = 1:length(La)
            if La[idx] < tol && ( (Lb[idx] > tol) + (Lc[idx] > tol) + (Ld[idx] > tol) < 3 )
                ids = [ ids; idx ];
            end
        end

        shift[ ids, : ] = warp1[ids] * t1[ face:face, : ] + warp2[ids] * t2[ face:face, : ];
    end

    XYZ = XYZ .+ shift;
    X = XYZ[ :, 1 ]; Y = XYZ[ :, 2 ]; Z = XYZ[ :, 3 ];
    return (X, Y, Z);
end

function EquiNodes3D(N)

    # function [X,Y,Z] = EquiNodes3D(N)
    # Purpose: compute the equidistributed nodes on the reference tetrahedron
    
    # total number of nodes
    Np = (Int)( (N + 1) * (N + 2) * (N + 3)/6 );
    
    # 2) create equidistributed nodes on equilateral triangle
    X = zeros( Np, 1 ); Y = zeros( Np, 1 ); Z = zeros( Np, 1 ); 
    
    sk = 1;
    for n = 1:(N + 1)
      for m = 1:(N + 2 - n)
        for q = 1:(N + 3 - n - m)
            
            X[sk] = -1 + ( q - 1 ) * 2 / N; 
            Y[sk] = -1 + ( m - 1 ) * 2 / N; 
            Z[sk] = -1 + ( n - 1 ) * 2 / N;
            sk = sk + 1;

        end
      end
    end

    return (X, Y, Z)
end

function WarpShiftFace3D(p,pval, pval2, L1,L2,L3,L4)

    # function [warpx, warpy] = WarpShiftFace3D(p,pval, pval2, L1,L2,L3,L4)     
    # Purpose: compute warp factor used in creating 3D Warp & Blend nodes
    
    (dtan1, dtan2) = evalshift(p, pval, L2, L3, L4);
    warpx = dtan1; 
    warpy = dtan2;
    
    return (warpx, warpy)

end

function evalshift(p, pval, L1, L2, L3)  

    # function [dx, dy] = evalshift(p, pval, L1, L2, L3)  
    # Purpose: compute two-dimensional Warp & Blend transform
    
    # 1) compute Gauss-Lobatto-Legendre node distribution
    (gaussX, val) = jacobi_LGL_quad( p );
    gaussX = -gaussX;
     
    # 2) compute blending function at each node for each edge
    blend1 = L2.*L3; 
    blend2 = L1.*L3; 
    blend3 = L1.*L2;
    
    # 3) amount of warp for each node, for each edge
    warpfactor1 = 4*evalwarp(p, gaussX, L3-L2); 
    warpfactor2 = 4*evalwarp(p, gaussX, L1-L3); 
    warpfactor3 = 4*evalwarp(p, gaussX, L2-L1); 
    
    # 4) combine blend & warp
    warp1 = blend1.*warpfactor1.*(1 .+ (pval*L1).^2);
    warp2 = blend2.*warpfactor2.*(1 .+ (pval*L2).^2);
    warp3 = blend3.*warpfactor3.*(1 .+ (pval*L3).^2);
    
    # 5) evaluate shift in equilateral triangle
    dx = 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
    dy = 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;
    return (dx, dy)
end    


function evalwarp(p, xnodes, xout)

    # function warp = evalwarp(p, xnodes, xout)
    # Purpose: compute one-dimensional edge warping function
    
    warp = zeros( size(xout) );
    xeq = zeros( p + 1 );

    for i = 1:(p + 1)
      xeq[i] = -1 + 2*( p + 1 - i ) / p;
    end
    
    for i = 1:(p + 1)
      d = ( xnodes[i] - xeq[i] );

      for j = 2:p
        if i != j
            d = d.*( xout .- xeq[j] ) / ( xeq[i] - xeq[j] );
        end
      end
      
      if i != 1
        d = -d / ( xeq[i] - xeq[1] );
      end
    
      if i != ( p + 1 )
        d = d/( xeq[i] - xeq[ p + 1 ] );
      end
    
      warp = warp .+ d;
    end

    return warp
end    

function xyztorst(X, Y, Z)

    # function [r,s,t] = xyztorst(x, y, z)
    # Purpose : Transfer from (x,y,z) in equilateral tetrahedron
    #           to (r,s,t) coordinates in standard tetrahedron
    
    v1 = [-1 -1/sqrt(3) -1/sqrt(6)];
    v2 = [ 1 -1/sqrt(3) -1/sqrt(6)];
    v3 = [ 0  2/sqrt(3) -1/sqrt(6)]; 
    v4 = [ 0  0/sqrt(3)  3/sqrt(6)];
    
    # back out right tet nodes
    rhs = [X'; Y'; Z'] - 0.5 * (v2' + v3' + v4' - v1') * ones( 1, length(X) );
    A = [0.5*(v2 - v1)' 0.5*(v3 - v1)' 0.5*(v4 - v1)'];
    RST = A \ rhs;
    r = RST[1:1, :]'; s = RST[2:2, :]'; t = RST[3:3, :]';
    return (r, s, t);

end

# constructs equispaced nodes on the pyramid surface. 

function pyramidSurfaceEquiNodes3D(N)

    # build surface nodes for a degree 5 pyramid trace
    # [-1,1]x[-1,1]x[0,1] with top vertex at (0,0,1)
    # start with W&B nodes for triangle
    #[r,s] = Nodes2D(N);
    #[r,s] = xytors(r,s);

    (r, s, t) = EquiNodes3D(N);

    inds = findall( tVal -> ( abs( tVal + 1 ) < 1e-8 ), t );
    r = r[ inds ];
    s = s[ inds ];

    Np = length(r);

    gr = LinRange(-1, 1, N + 1);

    r2d = gr' .* ones( N + 1 );
    s2d = ones( 1, N + 1 ) .* gr;
    t2d = zeros( N + 1, N + 1 );

    r2d = reshape( r2d, (N + 1) * (N + 1), 1);
    s2d = reshape( s2d, (N + 1) * (N + 1), 1);
    t2d = reshape( t2d, (N + 1) * (N + 1), 1);

    I = ones( Np, 1 );

    # Face 2
    x1 = -1; y1 = -1; z1 = +0;
    x2 =  1; y2 = -1; z2 = +0;
    x3 =  0; y3 =  0; z3 = +1;
    x3d = [r2d;0.5*(-(r .+ s)*x1 + (1 .+ r)*x2 + (1 .+ s)*x3)];
    y3d = [s2d;0.5*(-(r .+ s)*y1 + (1 .+ r)*y2 + (1 .+ s)*y3)];
    z3d = [t2d;0.5*(-(r .+ s)*z1 + (1 .+ r)*z2 + (1 .+ s)*z3)];


    # Face 3
    x1 =  1; y1 = -1; z1 = +0;
    x2 =  1; y2 = +1; z2 = +0;
    x3 =  0; y3 =  0; z3 = +1;
    x3d = [x3d;0.5*(-(r .+ s)*x1 + (1 .+ r)*x2 + (1 .+ s)*x3)];
    y3d = [y3d;0.5*(-(r .+ s)*y1 + (1 .+ r)*y2 + (1 .+ s)*y3)];
    z3d = [z3d;0.5*(-(r .+ s)*z1 + (1 .+ r)*z2 + (1 .+ s)*z3)];

    # Face 4
    x1 =  1; y1 = +1; z1 = +0;
    x2 = -1; y2 = +1; z2 = +0;
    x3 =  0; y3 =  0; z3 = +1;
    x3d = [x3d;0.5*(-(r .+ s)*x1 + (1 .+ r)*x2 + (1 .+ s)*x3)];
    y3d = [y3d;0.5*(-(r .+ s)*y1 + (1 .+ r)*y2 + (1 .+ s)*y3)];
    z3d = [z3d;0.5*(-(r .+ s)*z1 + (1 .+ r)*z2 + (1 .+ s)*z3)];
    
    # Face 5
    x1 =  -1; y1 = -1; z1 = +0;
    x2 =  -1; y2 = +1; z2 = +0;
    x3 =  0; y3 =  0; z3 = +1;
    x3d = [x3d;0.5*(-(r .+ s)*x1 + (1 .+ r)*x2 + (1 .+ s)*x3)];
    y3d = [y3d;0.5*(-(r .+ s)*y1 + (1 .+ r)*y2 + (1 .+ s)*y3)];
    z3d = [z3d;0.5*(-(r .+ s)*z1 + (1 .+ r)*z2 + (1 .+ s)*z3)];

    # remove duplicates
    tol = 1e-10;
    cnt = 1;
    x = x3d[1]; 
    y = y3d[1]; 
    z = z3d[1];

    for n = 2:length(x3d)
        d = ( x3d[n] .- x ).^2 + ( y3d[n] .- y ).^2 + ( z3d[n] .- z ).^2;
        
        if( minimum(d) > tol )
            x = [ x; x3d[n] ];
            y = [ y; y3d[n] ];
            z = [ z; z3d[n] ];
        end
    end

    return (x, y, z);

end

function JVandermonde3D(N, r, s, t, beta)

    # pyramid
    #
    #    4------3
    #    |\    /|
    #    | \  / |
    #    |  5   |
    #    | /  \ |
    #    |/    \|
    #    1------2
    
    # alphastore = [0;0;0;0.1002; 1.1332;1.5608;1.3413;1.2577;1.1603;...
    #     1.10153;0.6080;0.4523;0.8856;0.8717;0.9655];
    # alphaT = alphastore(N);
    
    tol = 1e-10;
    V = [];
    
    # vertex functions
    V = zeros( 5, length(r) );
    V[ :, 1 ] = .25 * ( 1 .- r .- s .- t .+ r.*s./(1 .- t .+ tol) );
    V[ :, 2 ] = .25 * ( 1 .+ r .- s .- t .- r.*s./(1 .- t .+ tol) );
    V[ :, 3 ] = .25 * ( 1 .+ r .+ s .- t .+ r.*s./(1 .- t .+ tol) );
    V[ :, 4 ] = .25 * ( 1 .- r .+ s .- t .- r.*s./(1 .- t .+ tol) );
    V[ :, 5 ] = t;
    
    v_ids = 1:5;
    
    sk = 6;
    
    # edge function for triangular edges
    edges = [1 5; 2 5; 3 5; 4 5];
    etri_ids = [];
    for e = 1:4
        i1 = edges[ e, 1 ];
        i2 = edges[ e, 2 ];
        
        if e == 1
            blend1 = V[ :, 2 ];
            blend2 = V[ :, 4 ];
        elseif e == 2
            blend1 = V[ :, 1 ];
            blend2 = V[ :, 3 ];
        elseif e == 3
            blend1 = V[ :, 4 ];
            blend2 = V[ :, 2 ];
        elseif e == 4
            blend1 = V[ :, 1 ];
            blend2 = V[ :, 3 ];
        end

        blend2 = (1 .+ beta*(blend1).^2).*(1 .+ beta*(blend2).^2); # blend towards opposite vertex on faces
        
        for i = 0:(N - 2)
            xi = V[ :, i1 ] - V[ :, i2 ];
            V[ :, sk ] = blend2.*V[ :, i1 ].*V[ :, i2 ].*jacobi_polynomial(xi, 1, 1, i);
            etri_ids = [etri_ids sk];
            sk = sk + 1;
        end
    end
    
    # edge functions for base
    edges = [1 2; 2 3; 3 4; 4 1];
    equad_ids = [];
    for e = 1:4
        i1 = edges[ e, 1 ];
        i2 = edges[ e, 2 ];
        
        #     op_vert = setdiff(1:5,edges(e,:));
        #V(:,op_vert(1)).*V(:,op_vert(2)).*V(:,op_vert(3));
        blend2 = (1 .+ beta*( V[ :, 5 ] ).^2 ); # blend towards top vertex
        
        for i=0:( N - 2 )
            xi = V[ :, i1 ] - V[ :, i2 ];
            V[ :, sk ] = blend2.*V[ :, i1 ] .* V[ :, i2 ].*JacobiP(xi, 1, 1, i);
            equad_ids = [equad_ids sk];
            sk = sk + 1;
        end
    end
    
    #triangular faces
    faces = [1 2 5; 2 3 5; 3 4 5; 4 1 5];
    ftri_ids = [];
    for f = 1:4
        i1 = faces[ f, 1 ];
        i2 = faces[ f, 2 ];
        i3 = faces[ f, 3 ];
        
        # bubble edge blend
        #         op_vert = setdiff(1:5,faces(f,:));  blend2 = V(:,op_vert(1)).*V(:,op_vert(2));
        
        # plane edge blend
        if f == 1
            blend2 = ( ( 1 .+ s ) - t )/2;
        elseif f == 2
            blend2 = ( ( 1 .- r ) - t )/2;
        elseif f == 3
            blend2 = ( ( 1 .- s ) - t )/2;
        elseif f == 4
            blend2 = ( ( 1 .+ r ) - t ) / 2;
        end
        #blend2 = 1-V(:,5); % blend to base
        
        blend2 = 1 .+ beta*( blend2 ).^2;
        #blend2 = beta*(V(:,5)-.5); % blend to base - switches negative at midway pt.
        
        
        L1 = V[ :, i1 ];    L2 = V[ :, i2 ];    L3 = V[ :, i3 ];
        (x, y) = eqbarytoxy(L1, L2, L3);
        (rr, ss) = xytors( x, y );
        Vf = Vandermonde2D( N - 3, rr, ss );
        for i = 1:size(Vf,2)
            V[ :, sk ] = blend2.*V[ :, i1 ].*V[ :, i2 ].*V[ :, i3 ].*Vf[ :, i ];
            ftri_ids = [ftri_ids sk];
            sk = sk + 1;
        end
    end
    
    # square face on the bottom
    blend2 = 1 .+ beta * ( V[ :, 5 ] ).^2;
    fquad_ids = [];
    for i = 0:(N - 2)
        for j = 0:(N - 2)
            V[ :, sk ] = blend2 .* V[ :, 1 ].*V[ :, 2 ].*V[ :, 3 ].*V[ :, 4 ].*jacobi_polynomial(r, 1, 1, i) .* jacobi_polynomial(s, 1, 1, j);
            fquad_ids = [fquad_ids sk];
            sk = sk + 1;
        end
    end
    
    return (V, v_ids, etri_ids, equad_ids, ftri_ids, fquad_ids)
end

# convert from reference triangle coordinates to equilateral triangle
# coordinates 
#
# function [x, y] = rstoxy(r,s)
#
function eqbarytoxy(L1, L2, L3)

    # equilateral vertices
    v1 = 2*[-.5 -sqrt(3)/6]';
    v2 = 2*[.5 -sqrt(3)/6]';
    v3 = 2*[0 sqrt(3)/3]';

    x = zeros( length( L1 ), 1 );
    y = zeros( length( L1 ), 1 );

    for i = 1:length(L1)
        XY = v1 * L1[i] + v2 * L2[i] + v3 * L3[i];
        x[i] = XY[1];
        y[i] = XY[2];
    end

    return (x, y)

end

function xytors(x, y)

    # function [r,s] = xytors(x, y)
    # Purpose : From (x,y) in equilateral triangle to (r,s) coordinates in standard triangle
    
    L1 = ( sqrt(3.0)*y .+ 1.0 )/3.0;
    L2 = (-3.0*x - sqrt(3.0)*y .+ 2.0)/6.0;
    L3 = ( 3.0*x - sqrt(3.0)*y .+ 2.0)/6.0;
    
    r = -L2 + L3 - L1; s = -L2 - L3 + L1;

    return (r, s)

end
    
function Vandermonde2D(N, r, s);

    # function [V2D] = Vandermonde2D(N, r, s);
    # Purpose : Initialize the 2D Vandermonde Matrix,  V_{ij} = phi_j(r_i, s_i);
    
    V2D = zeros( length(r), Int( (N + 1) * (N + 2)/2 ) );
    
    # Transfer to (a,b) coordinates
    (a, b) = rstoab(r, s);
    
    # build the Vandermonde matrix

    sk = 1;
    for i = 0:N
      for j = 0:(N - i)
        V2D[:, sk] = Simplex2DP(a,b,i,j);
        sk = sk + 1;
      end
    end

    return V2D

end
    

function rstoab(r, s)

    # function [a,b] = rstoab(r,s)
    # Purpose : Transfer from (r,s) -> (a,b) coordinates in triangle
    
    Np = length(r); 
    a = zeros(Np, 1);

    for n = 1:Np

        if( s[n] != 1 )
            a[n] = 2 * ( 1 + r[n] )/( 1 - s[n] ) - 1;
        else 
            a[n] = -1; 
        end

    end
    
    b = s;
    return (a, b)

end    

function Simplex2DP(a, b, i, j);

    # function [P] = Simplex2DP(a,b,i,j);
    # Purpose : Evaluate 2D orthonormal polynomial
    #           on simplex at (a,b) of order (i,j).
    
    h1 = jacobi_polynomial(a, 0, 0, i); 
    h2 = jacobi_polynomial(b, 2*i + 1, 0, j);
    P = sqrt(2.0) * h1.*h2.*(1 - b).^i;

    return P;

end

# Conical construction of equispaced nodes on the pyramid
# optional arg: t, specifies levels at which to place nodes

function pyramidEquiNodes3D(N)

    t = LinRange(0, 1, N + 1);

    x = []; 
    y = []; 
    z= [];

    for level = 0:N
        a = ( 1 - t[level + 1] );
        if level < N        
            r1D = LinRange( -a, a, N + 1 - level );
        else
            r1D = 0;
        end

        r = r1D' .* ones( N + 1 - level );
        s = ones( 1, N + 1 - level ) .* r1D;
        
        r = reshape( r, (N + 1 - level) * (N + 1 - level), 1);
        s = reshape( s, (N + 1 - level) * (N + 1 - level), 1);

        x = [x; r];
        y = [y; s];    
        z = [z; t[ level + 1 ] * ones( size( r ) ) ];

    end

    return (x, y, z)
end
