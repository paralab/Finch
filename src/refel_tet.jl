#=
# A tetrehedral reference element.
# Notation follows that used in Nodal Discontinuous Galerkin Methods
# by Hesthaven and Warburton.
#  https://link.springer.com/book/10.1007/978-0-387-72067-8
#
=#

function build_tetrahedron_refel(refel)
    # refel has already been created, but needs quadrature matrices
    (refel.V, refel.Dr, refel.Ds, refel.Dt) = tetrahedron_vandermonds(refel, refel.r);
    refel.invV = inv(refel.V);
    
    (refel.Vg, DVgr, DVgs, DVgt) = tetrahedron_vandermonds(refel, refel.g);
    #refel.invVg = inv(refel.Vg);
    
    refel.Q = refel.Vg*refel.invV;
    refel.Qr = DVgr*refel.invV;
    refel.Qs = DVgs*refel.invV;
    refel.Qt = DVgt*refel.invV;
    refel.Ddr = refel.Dr*refel.invV;
    refel.Dds = refel.Ds*refel.invV;
    refel.Ddt = refel.Dt*refel.invV;
    
    return refel;
end

function tetrahedron_vandermonds(refel, r)
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
        if r[ni,2]+r[ni,3] != 0
            a[ni] = 2*(1+r[ni,1])/(-r[ni,2]-r[ni,3])-1;
        else
            a[ni] = -1;
        end
        if r[ni,3] != 1
            b[ni] = 2*(1+r[ni,2])/(1-r[ni,3])-1;
        else
            b[ni] = -1;
        end
    end
    c = r[:,3];
    
    # build the Vandermonde and gradVandermond matrix
    sk = 1;
    for i=0:refel.N
        h1 = jacobi_polynomial(a,0,0,i);
        dfa = grad_jacobi_polynomial(a, 0, 0, i);
        for j=0:(refel.N - i)
            h2 = jacobi_polynomial(b,2*i+1,0,j);
            dgb = grad_jacobi_polynomial(b, 2*i+1,0, j);
            for k=0:(refel.N - i - j)
                h3 = jacobi_polynomial(c,2*(i+j)+2,0,k);
                dhc = grad_jacobi_polynomial(b, 2*(i+j)+2,0, k);
                
                V[:,sk] = 2.8284271247461903 .* h1 .* h2 .* h3 .* (1 .- b).^i .* (1 .- c).^(i+j);
                
                # r-derivative
                dmodedr = dfa.*h2.*h3;
                if i>0
                    dmodedr = dmodedr.*((0.5 .* (1 .- b)).^(i-1));
                end
                if (i+j)>0
                    dmodedr = dmodedr.*((0.5 .* (1 .- c)).^(i+j-1));
                end
                
                # s-derivative
                dmodeds = dmodedr .* (0.5 .* (1 .+ a));
                tmp = dgb .* (0.5 .* (1 .- b)).^i;
                if i>0
                    tmp = tmp .+ (-0.5*i) .* (h2 .* (0.5 .* (1 .- b)).^(i-1));
                end
                if i+j>0
                    tmp = tmp .* (0.5 .* (1 .- c)).^(i+j-1); 
                end
                tmp = tmp .* h3 .* h1;
                dmodeds = dmodeds .+ tmp;
                
                # t-derivative
                dmodedt = 0.5 .* (1 .+ a) .* dmodedr .+ 0.5 .* (1 .+ b) .* tmp;
                tmp = dhc .* (0.5 .* (1 .- c)).^(i+j);
                if i+j>0
                    tmp = tmp .- (0.5*(i+j)) .* h3 .* (0.5 .* (1 .- c)).^(i+j-1); 
                end
                tmp = h1 .* (h2 .* tmp) .* (0.5 .* (1 .- b)).^i;
                dmodedt = dmodedt .+ tmp;
                
                # Normalize
                dmodedr = 2^(2*i+j+1.5) .* dmodedr; 
                dmodeds = 2^(2*i+j+1.5) .* dmodeds;
                dmodedt = 2^(2*i+j+1.5) .* dmodedt;
                
                gradVr[:,sk] = dmodedr;
                gradVs[:,sk] = dmodeds;
                gradVt[:,sk] = dmodedt;
                
                sk = sk+1;
            end
        end
    end
    
    return (V, gradVr, gradVs, gradVt);
end