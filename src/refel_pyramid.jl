#=
# A tetrehedral reference element.
# Notation follows that used in Nodal Discontinuous Galerkin Methods
# by Hesthaven and Warburton.
#  https://link.springer.com/book/10.1007/978-0-387-72067-8
#
=#

function build_pyramid_refel(refel)
    # refel has already been created, but needs quadrature matrices
    (refel.V, refel.Dr, refel.Ds, refel.Dt) = pyramid_vandermonds(refel, refel.r);
    refel.invV = inv(refel.V);
    
    (refel.Vg, DVgr, DVgs, DVgt) = pyramid_vandermonds(refel, refel.g);
    #refel.invVg = inv(refel.Vg);
    
    refel.Q = refel.Vg*refel.invV;
    refel.Qr = DVgr*refel.invV;
    refel.Qs = DVgs*refel.invV;
    refel.Qt = DVgt*refel.invV;
    refel.Ddr = refel.Dr*refel.invV;
    refel.Dds = refel.Ds*refel.invV;
    refel.Ddt = refel.Dt*refel.invV;

    println("Start")
    println("Q")
    println(refel.Q)
    println("Qr")
    println(refel.Qr)
    println("Qs")
    println(refel.Qs)
    println("Qt")
    println(refel.Qt)
    println("End")

    # println(refel.r)
    
    return refel;
end

function pyramid_vandermonds(refel, r)
    Np = refel.Np;      # number of nodal points
    Nrp = size(r,1);    # number of quadrature points may be different
    V = zeros(Nrp, Np);
    gradVr = zeros(Nrp, Np);
    gradVs = zeros(Nrp, Np);
    gradVt = zeros(Nrp, Np);
    
    # Transfer (r,s,t) to (a,b,c) coordinates
    a = zeros(Nrp); 
    b = zeros(Nrp);
    for ni=1:Nrp  
        if abs( r[ni,3] - 1 ) > 1e-6
            a[ni] = ( r[ni,1] )/( 1 - r[ ni, 3 ] );
        else
            a[ni] = 0;
        end
        if abs( r[ni,3] - 1 ) > 1e-6
            b[ni] = ( r[ni,2] )/( 1 - r[ni,3] );
        else
            b[ni] = 0;
        end
    end
    c = 2 * r[:,3] .- 1;
    
    # println( "a is ")
    # println( a )

    # println( "b is ")
    # println( b )

    # println( "c is ")
    # println( c )
    
    # build the Vandermonde and gradVandermond matrix
    sk = 1;
    for i = 0:refel.N

        h1 = jacobi_polynomial(a, 0 , 0 , i);
        dfa = grad_jacobi_polynomial(a, 0, 0, i);

        for j = 0:(refel.N)

            h2 = jacobi_polynomial(b, 0, 0, j);
            dgb = grad_jacobi_polynomial(b, 0, 0, j);
            muij = max( i, j );
            
            for k = 0:(refel.N - muij)

                h3 = jacobi_polynomial( c, 2*muij + 2, 0, k );
                dhc = grad_jacobi_polynomial( c, 2*muij + 2, 0, k);
                
                V[:,sk] = 2^( ( 2 * muij + 3 ) / 2 ) .* h1 .* h2 .* h3 .* ( 0.5 * (1 .- c) ).^( muij );
                
                # r-derivative
                dmodedr = dfa.*h2.*h3;
                if muij > 0
                    dmodedr = dmodedr.*( ( 0.5*(1 .- c) ).^( muij - 1 ) );
                end
                
                # s-derivative
                dmodeds = h1 .* dgb .* h3;
                if muij > 0
                    dmodeds = dmodeds.*( ( 0.5 * (1 .- c) ).^( muij - 1 ) );
                end
                
                # t-derivative
                tmp = 0

                if muij > 0
                    tmp = tmp .+ dfa .* h2 .* h3 .* a
                    tmp = tmp .+ h1 .* dgb .* h3 .* b
                    tmp = tmp .* ( 0.5 * (1 .- c ) ) .^ (muij - 1)
                end

                tmp = tmp .+ 2 * h1 .* h2 .* dhc .* ( (0.5 * (1 .- c)) .^ muij )

                if muij > 0
                    tmp = tmp .- muij .* h1 .* h2 .* h3 .* ( (0.5 * (1 .- c)) .^ (muij - 1) )
                end

                dmodedt = tmp;
                
                # Normalize
                dmodedr = 2^( ( 2 * muij + 3 )/2 ) .* dmodedr; 
                dmodeds = 2^( ( 2 * muij + 3 )/2 ) .* dmodeds;
                dmodedt = 2^( ( 2 * muij + 3 )/2 ) .* dmodedt;
                
                gradVr[ :,sk ] = dmodedr;
                gradVs[ :,sk ] = dmodeds;
                gradVt[ :,sk ] = dmodedt;
                
                sk = sk + 1;
            end
        end
    end
   
    # println("GRAD VR")
    # println( gradVr )

    # println("GRAD VS")
    # println( gradVs )

    # println("GRAD VT")
    # println( gradVt )
    # println("END")

    return (V, gradVr, gradVs, gradVt);
end