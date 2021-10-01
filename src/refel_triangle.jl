#=
# A triangular reference element.
# Notation follows that used in Nodal Discontinuous Galerkin Methods
# by Hesthaven and Warburton.
#  https://link.springer.com/book/10.1007/978-0-387-72067-8
#
=#

function build_triangle_refel(refel)
    # refel has already been created, but needs quadrature matrices
    (refel.V, refel.Dr, refel.Ds) = triangle_vandermonds(refel, refel.r);
    refel.invV = inv(refel.V);
    
    (refel.Vg, DVgr, DVgs) = triangle_vandermonds(refel, refel.g);
    #refel.invVg = inv(refel.Vg);
    
    refel.Q = refel.Vg*refel.invV;
    refel.Qr = DVgr*refel.invV;
    refel.Qs = DVgs*refel.invV;
    refel.Ddr = refel.Dr*refel.invV;
    refel.Dds = refel.Ds*refel.invV;
    
    ## Surfaces
    nfaces = 3;
    refel.surf_V = fill!(Array{Array{Float64}}(undef, nfaces), []);
    #refel.surf_gradV = fill!(Array{Array{Float64}}(undef, nfaces), []);
    refel.surf_Vg = fill!(Array{Array{Float64}}(undef, nfaces), []);
    #refel.surf_gradVg = fill!(Array{Array{Float64}}(undef, nfaces), []);
    
    refel.surf_Q = fill!(Array{Array{Float64}}(undef, nfaces), []);
    refel.surf_Qr = fill!(Array{Array{Float64}}(undef, nfaces), []);
    refel.surf_Qs = fill!(Array{Array{Float64}}(undef, nfaces), []);
    refel.surf_Ddr = fill!(Array{Array{Float64}}(undef, nfaces), []);
    refel.surf_Dds = fill!(Array{Array{Float64}}(undef, nfaces), []);
    
    for fi=1:nfaces
        (refel.surf_V[fi], surf_gradVr, surf_gradVs) = triangle_vandermonds(refel, refel.surf_r[fi]);
        (refel.surf_Vg[fi], surf_gradVgr, surf_gradVgs) = triangle_vandermonds(refel, refel.surf_g[fi]);
        
        refel.surf_Q[fi] = refel.surf_V[fi] * refel.invV;
        refel.surf_Qr[fi] = surf_gradVgr * refel.invV;
        refel.surf_Qs[fi] = surf_gradVgs * refel.invV;
        refel.surf_Ddr[fi] = surf_gradVr * refel.invV;
        refel.surf_Dds[fi] = surf_gradVs * refel.invV;
    end
    
    return refel;
end

function triangle_vandermonds(refel, r)
    Np = refel.Np;      # number of nodal points
    Nrp = size(r,1);    # number of quadrature points may be different
    V = zeros(Nrp, Np);
    gradVr = zeros(Nrp, Np);
    gradVs = zeros(Nrp, Np);
    
    # Transfer (r,s) to (a,b) coordinates
    a = zeros(Nrp,1);
    for ni=1:Nrp
        if r[ni,2] != 1
            a[ni] = 2*(1+r[ni,1])/(1-r[ni,2])-1;
        else 
            a[ni] = -1; 
        end
    end
    b = r[:,2];
    
    # build the Vandermonde and gradVandermond matrix
    sk = 1;
    for i=0:refel.N
        h1 = jacobi_polynomial(a,0,0,i);
        dfa = grad_jacobi_polynomial(a, 0, 0, i);
        for j=0:refel.N - i
            h2 = jacobi_polynomial(b,2*i+1,0,j);
            dgb = grad_jacobi_polynomial(b, 2*i+1,0, j);
            V[:,sk] = 1.4142135623730951 .* h1 .* h2 .* (1 .- b).^i;
            
            # r-derivative
            # d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
            dmodedr = dfa.*h2;
            if i>0
                dmodedr = dmodedr.*((0.5 .* (1 .- b)).^(i-1));
            end
            # s-derivative
            # d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
            dmodeds = dfa.*(h2.*(0.5 .* (1 .+ a)));
            if i>0
                dmodeds = dmodeds.*((0.5 .* (1 .- b)).^(i-1));
            end
            
            tmp = dgb.*((0.5 .* (1 .- b)).^i);
            if i>0
                tmp = tmp .- 0.5*i .* h2 .* ((0.5 .* (1 .- b)).^(i-1));
            end
            dmodeds = dmodeds + h1.*tmp;
            
            # Normalize
            dmodedr = 2^(i+0.5) .* dmodedr; 
            dmodeds = 2^(i+0.5) .* dmodeds;
            
            gradVr[:,sk] = dmodedr;
            gradVs[:,sk] = dmodeds;
            
            sk = sk+1;
        end
    end
    
    return (V, gradVr, gradVs);
end