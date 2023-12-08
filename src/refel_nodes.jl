#=
# Compute the node locations for the reference element
=#
include("jacobi_gauss_quad.jl");
include("triangle_nodes.jl");
include("tet_nodes.jl");

function refel_nodes!(refel, nodetype)
    if refel.dim == 0
        # 0D is a point so some of this doesn't make sense
        refel.r1d = [0];
        refel.wr1d = [1];
        refel.g1d = [0];
        refel.wg1d = [1];
        refel.r = zeros(1,1);
        refel.wr = [1];
        refel.g = zeros(1,1);
        refel.wg = [1];
        
    elseif refel.dim == 1
        # 1D has line segments
        if nodetype == UNIFORM
            refel.r1d = Array(-1:(2/(refel.Np-1)):1);
            refel.wr1d = ones(length(refel.r1d)) ./ length(refel.r1d);
        elseif nodetype == GAUSS
            (r,w) = jacobi_gauss_quad(0,0,refel.N);
            refel.r1d = r;
            refel.wr1d = w;
        elseif nodetype == LOBATTO
            if refel.N == 1
                refel.r1d = [-1; 1];
                refel.wr1d = [1; 1];
            else
                # (r,w) = jacobi_gauss_quad(1,1,refel.N-2);
                # refel.r1d = [-1; r ; 1];
                
                # # compute the weights
                # w = jacobi_polynomial(refel.r1d, 0, 0, refel.N)
                # adgammaN = (2*refel.N + 1) / (refel.N * (refel.N + 1))
                # w = w.*w
                # w = adgammaN./w
                
                # refel.wr1d = w;
                (refel.r1d, refel.wr1d) = jacobi_LGL_quad(refel.N);
            end
        end
        # Then find Gauss points
        (g,w) = jacobi_gauss_quad(0,0,refel.N);
        refel.g1d = g;
        refel.wg1d = w;
        
        # 1D is easy
        refel.g = reshape(refel.g1d,:,1);
        refel.wg = refel.wg1d;
        refel.r = reshape(refel.r1d,:,1);
        refel.wr = refel.wr1d;
        
        # surface is just points
        mone = ones(1,1);
        mnone = -mone;
        refel.face2local = [[1], [refel.Np]];
        refel.surf_r = Vector{Matrix{Float64}}(undef,2);
        refel.surf_wr = Vector{Vector{Float64}}(undef,2);
        refel.surf_g = Vector{Matrix{Float64}}(undef,2);
        refel.surf_wg = Vector{Vector{Float64}}(undef,2);
        
        refel.surf_r[1] = mnone; refel.surf_r[2] = mone;
        refel.surf_wr[1] = [1]; refel.surf_wr[2] = [1];
        refel.surf_g[1] = mnone; refel.surf_g[2] = mone;
        refel.surf_wg[1] = [1]; refel.surf_wg[2] = [1];
        
    elseif refel.dim == 2
        # 2D has triangles and quads
        if refel.Nfaces == 3 # triangles
            triangle_refel_nodes!(refel);
        else # quads
            if nodetype == UNIFORM
                refel.r1d = Array(-1 : (2 / ( refel.N ) ) : 1);
                refel.wr1d = ones(length(refel.r1d)) ./ length(refel.r1d);
            elseif nodetype == GAUSS
                (r,w) = jacobi_gauss_quad(0,0,refel.N);
                refel.r1d = r;
                refel.wr1d = w;
            elseif nodetype == LOBATTO
                if refel.N == 1
                    refel.r1d = [-1; 1];
                    refel.wr1d = [1; 1];
                else
                    (r,w) = jacobi_gauss_quad(1,1,refel.N-2);
                    refel.r1d = [-1; r ; 1];
                    
                    # compute the weights
                    w = jacobi_polynomial(refel.r1d, 0, 0, refel.N);
                    adgammaN = (2*refel.N + 1) / (refel.N * (refel.N + 1));
                    w = w.*w;
                    w = adgammaN./w;
                    
                    refel.wr1d = w;
                end
            end
            # Then find Gauss quadrature points
            (g,w) = jacobi_gauss_quad(0,0,refel.N);
            refel.g1d = g;
            refel.wg1d = w;
            
            # r and w
            refel.r = zeros(refel.Np,2);
            refel.wr = zeros(refel.Np);
            refel.g = zeros(refel.Np,2);
            refel.wg = zeros(refel.Np);
            n1d = length(refel.r1d);
            for j=1:n1d
                for i=1:n1d
                    k = (j-1)*n1d + i;
                    refel.r[k,:] = [refel.r1d[i]; refel.r1d[j]];
                    refel.wr[k] = refel.wr1d[i] * refel.wr1d[j];
                    refel.g[k,:] = [refel.g1d[i]; refel.g1d[j]];
                    refel.wg[k] = refel.wg1d[i] * refel.wg1d[j];
                end
            end
            
            ### surface lines ################
            # First get the face2local maps
            tol = 1e-8;
            qf1(x) = abs(x[1] + 1) < tol;
            qf2(x) = abs(x[2] + 1) < tol;
            qf3(x) = abs(x[1] - 1) < tol;
            qf4(x) = abs(x[2] - 1) < tol;
            refel.face2local = [get_face2local_map(refel.r, qf1),
                                get_face2local_map(refel.r, qf2),
                                get_face2local_map(refel.r, qf3),
                                get_face2local_map(refel.r, qf4)];
            
            tmp = zeros(2,0);
            refel.surf_r = [tmp, tmp, tmp, tmp];
            refel.surf_g = [tmp, tmp, tmp, tmp];
            
            for fi=1:4
                refel.surf_r[fi] = refel.r[refel.face2local[fi], :];
                refel.surf_g[fi] = refel.g[refel.face2local[fi], :];
            end
            refel.surf_wr = [refel.wr1d, refel.wr1d, refel.wr1d, refel.wr1d];
            refel.surf_wg = [refel.wg1d, refel.wg1d, refel.wg1d, refel.wg1d];
        end
        
        
    elseif refel.dim == 3
        # 3D has tets, hexs and prisms
        if refel.Nfaces == 4 # tets
            tetrahedron_refel_nodes!(refel);
        else # hexs
            if nodetype == UNIFORM
                refel.r1d = Array( -1 : (2 / ( refel.N ) ) : 1 );
                refel.wr1d = ones(length(refel.r1d)) ./ length(refel.r1d);
            elseif nodetype == GAUSS
                (r,w) = jacobi_gauss_quad(0,0,refel.N);
                refel.r1d = r;
                refel.wr1d = w;
            elseif nodetype == LOBATTO
                if refel.N == 1
                    refel.r1d = [-1; 1];
                    refel.wr1d = [1; 1];
                else
                    (r,w) = jacobi_gauss_quad(1,1,refel.N-2);
                    refel.r1d = [-1; r ; 1];
                    
                    # compute the weights
                    w = jacobi_polynomial(refel.r1d, 0, 0, refel.N);
                    adgammaN = (2*refel.N + 1) / (refel.N * (refel.N + 1));
                    w = w.*w;
                    w = adgammaN./w;
                    
                    refel.wr1d = w;
                end
            end
            # Then find Gauss points
            (g,w) = jacobi_gauss_quad(0,0,refel.N);
            refel.g1d = g;
            refel.wg1d = w;
            
            # r and w
            refel.r = zeros(refel.Np,3);
            refel.wr = zeros(refel.Np);
            refel.g = zeros(refel.Np,3);
            refel.wg = zeros(refel.Np);
            n1d = length(refel.r1d);
            for k=1:n1d
                for j=1:n1d
                    for i=1:n1d
                        ind = (k-1)*n1d*n1d + (j-1)*n1d + i;
                        refel.r[ind,:] = [refel.r1d[i]; refel.r1d[j]; refel.r1d[k]];
                        refel.wr[ind] = refel.wr1d[i] * refel.wr1d[j] * refel.wr1d[k];
                        refel.g[ind,:] = [refel.g1d[i]; refel.g1d[j]; refel.g1d[k]];
                        refel.wg[ind] = refel.wg1d[i] * refel.wg1d[j] * refel.wg1d[k];
                    end
                end
            end
            
            ### surface quads ################
            # First get the face2local maps
            tol = 1e-8;
            hf1(x) = abs(x[1] + 1) < tol;
            hf2(x) = abs(x[2] + 1) < tol;
            hf3(x) = abs(x[3] + 1) < tol;
            hf4(x) = abs(x[1] - 1) < tol;
            hf5(x) = abs(x[2] - 1) < tol;
            hf6(x) = abs(x[3] - 1) < tol;
            refel.face2local = [get_face2local_map(refel.r, hf1),
                                get_face2local_map(refel.r, hf2),
                                get_face2local_map(refel.r, hf3),
                                get_face2local_map(refel.r, hf4),
                                get_face2local_map(refel.r, hf5),
                                get_face2local_map(refel.r, hf6)];
            
            tmp = zeros(3,0);
            refel.surf_r = [tmp, tmp, tmp, tmp, tmp, tmp];
            refel.surf_g = [tmp, tmp, tmp, tmp, tmp, tmp];
            
            for fi=1:6
                refel.surf_r[fi] = refel.r[refel.face2local[fi], :];
                refel.surf_g[fi] = refel.g[refel.face2local[fi], :];
            end
            
            # The weights are a little more tricky
            fwr = zeros(refel.Nfp[1]);
            fwg = zeros(refel.Nfp[1]);
            n1d = length(refel.r1d);
            for j=1:n1d
                for i=1:n1d
                    k = (j-1)*n1d + i;
                    fwr[k] = refel.wr1d[i] * refel.wr1d[j];
                    fwg[k] = refel.wg1d[i] * refel.wg1d[j];
                end
            end
            refel.surf_wr = [fwr, fwr, fwr, fwr, fwr, fwr];
            refel.surf_wg = [fwg, fwg, fwg, fwg, fwg, fwg];
        end
        
    else
        # Not ready for 4D
    end
end

# Determines surface nodes using a supplied function.
# Returns a face to local map. The order of which will be omnipotent.
function get_face2local_map(r, compare)
    n = size(r,1);
    map = zeros(Int,n);
    nf = 0;
    for i=1:n
        if compare(r[i,:])
            nf = nf+1;
            map[nf] = i;
        end
    end
    
    return map[1:nf];
end

# Returns the global node locations for elemental nodes
# Has size Np,dim
function get_node_coords(vx, elnodes)
    x = [];
    if refel.dim == 1
        hx = vx[2] - vx[1];
        x = vx[1] .+ (elnodes .+ 1) .* (hx*0.5);
    elseif refel.dim == 2
        # for now assume rectangular quads
        x = zeros(size(elnodes));
        hx = max(vx[2,1]-vx[1,1], vx[3,1]-vx[1,1]);
        hy = max(vx[2,2]-vx[1,2], vx[3,2]-vx[1,2]);
        x[:,1] = vx[1,1] .+ (elnodes[:,1] .+ 1) .*(hx*0.5);
        x[:,2] = vx[1,2] .+ (elnodes[:,2] .+ 1) .*(hy*0.5);
    elseif refel.dim == 3
        # for now assume regular hexs
        x = zeros(size(elnodes));
        hx = max(vx[2,1]-vx[1,1], vx[3,1]-vx[1,1]);
        hy = max(vx[2,2]-vx[1,2], vx[3,2]-vx[1,2]);
        x[:,1] = vx[1,1] .+ (elnodes[:,1] .+ 1) .*(hx*0.5);
        x[:,2] = vx[1,2] .+ (elnodes[:,2] .+ 1) .*(hy*0.5);
        x[:,3] = vx[1,3] .+ (elnodes[:,3] .+ 1) .*(hy*0.5);
    end
    
    return x;
end