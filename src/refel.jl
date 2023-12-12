#=
# A reference element.
# Notation follows that used in Nodal Discontinuous Galerkin Methods
# by Hesthaven and Warburton.
#  https://link.springer.com/book/10.1007/978-0-387-72067-8
#
# Can be used for 1, 2, 3, 4 dimensions
=#
include("jacobi_polynomial.jl");
include("refel_nodes.jl");
include("refel_triangle.jl");
include("refel_tet.jl");

# mutable struct Refel
#     dim::Int                # Dimension
#     N::Int                  # Order of polynomials
#     Np::Int                 # Number of nodes
#     Nqp::Int                # Number of quadrature points
#     Nfaces::Int             # Number of faces
#     Nfp::Array{Int,1}       # Number of nodes for each face

#     ######################################
#     # Volume nodes and quadrature matrices
#     ######################################
#     r1d::Array{Float64}     # Node coordinates in 1D
#     r::Array{Float64}       # dim-dim Node coordinates

#     wr1d::Array{Float64}    # r1d gll Quadrature weights
#     wr::Array{Float64}      # r gll Quadrature weights

#     g1d::Array{Float64}     # 1D Gauss points
#     wg1d::Array{Float64}    # 1D Gauss weights

#     g::Array{Float64}       # dim-dim Gauss points
#     wg::Array{Float64}      # dim-dim Gauss weights

#     V::Array{Float64}       # basis at r
#     gradV::Array{Float64}   # grad of basis at r
#     invV::Array{Float64}    # Inverse V

#     Vg::Array{Float64}      # basis at Gauss
#     gradVg::Array{Float64}  # grad of basis at g
#     invVg::Array{Float64}   # Inverse Vg

#     Dr::Array{Float64}      # Differentiation matrix for r
#     Ds::Array{Float64}      # Differentiation matrix for s
#     Dt::Array{Float64}      # Differentiation matrix for t
#     Dg::Array{Float64}      # Differentiation matrix for g

#     # Useful quadrature matrices for the volume integrals
#     Q1d::Array{Float64}     # 1D quadrature matrix: like Vg*invV
#     Q::Array{Float64}       # dim-dim quadrature matrix
#     Qr::Array{Float64}      # quad of derivative matrix: like gradVg*invV
#     Qs::Array{Float64}      # 
#     Qt::Array{Float64}      # 

#     Ddr::Array{Float64}      # Derivatives at the elemental nodes, not quadrature nodes
#     Dds::Array{Float64}      # 
#     Ddt::Array{Float64}      #

#     #######################################
#     # Surface nodes and quadrature matrices
#     #######################################
#     face2local::Array{Array{Int}}       # maps face nodes to local indices

#     surf_r::Array{Array{Float64}}       # surface node coordinates
#     surf_wr::Array{Array{Float64}}      # surface gll weights

#     surf_g::Array{Array{Float64}}       # surface Gauss points
#     surf_wg::Array{Array{Float64}}      # surface Gauss weights

#     surf_V::Array{Array{Float64}}       # basis at surf_r
#     surf_gradV::Array{Array{Float64}}   # grad of basis at surf_r
#     surf_invV::Array{Array{Float64}}    # Inverse surf_V

#     surf_Vg::Array{Array{Float64}}      # basis at surf_g
#     surf_gradVg::Array{Array{Float64}}  # grad of basis at surf_g
#     surf_invVg::Array{Array{Float64}}   # Inverse surf_Vg

#     surf_Dr::Array{Array{Float64}}      # Differentiation matrix for surf_r
#     surf_Ds::Array{Array{Float64}}      # Differentiation matrix for surf_s
#     surf_Dt::Array{Array{Float64}}      # Differentiation matrix for surf_t
#     surf_Dg::Array{Array{Float64}}      # Differentiation matrix for surf_g

#     surf_Q::Array{Array{Float64}}       # quadrature matrix
#     surf_Qr::Array{Array{Float64}}      # derivative quadrature matrix
#     surf_Qs::Array{Array{Float64}}      # 
#     surf_Qt::Array{Array{Float64}}      # 

#     surf_Ddr::Array{Array{Float64}}     # Derivatives at the elemental nodes, not quadrature nodes
#     surf_Dds::Array{Array{Float64}}     # 
#     surf_Ddt::Array{Array{Float64}}     #

#     # Constructor needs at least this information
#     Refel(dim, order, nnodes, nfaces, nfp) = new(
#         dim,
#         order,
#         nnodes,
#         -1,
#         nfaces,
#         nfp,
#         [],[],
#         [],[],
#         [],[],
#         [],[],
#         [],[],[],
#         [],[],[],
#         [],[],[],[],
#         [],[],[],[],[],
#         [],[],[],
#         [[]],
#         [[]],[[]],
#         [[]],[[]],
#         [[]],[[]],[[]],
#         [[]],[[]],[[]],
#         [[]],[[]],[[]],[[]],
#         [[]],[[]],[[]],[[]],
#         [[]],[[]],[[]]
#     )
# end

import Base.copy
function copy(ref::Refel)
    T = finch_state.config.float_type;
    newref = Refel(T, ref.dim, ref.N, ref.Np, ref.Nfaces, ref.Nfp);
    newref.Nqp = ref.Nqp;
    newref.r1d = copy(ref.r1d);
    newref.r = copy(ref.r);
    newref.wr1d = copy(ref.wr1d);
    newref.wr = copy(ref.wr);
    newref.g1d = copy(ref.g1d);
    newref.wg1d = copy(ref.wg1d);
    newref.g = copy(ref.g);
    newref.wg = copy(ref.wg);
    newref.V = copy(ref.V);
    newref.gradV = copy(ref.gradV);
    newref.invV = copy(ref.invV);
    newref.Vg = copy(ref.Vg);
    newref.gradVg = copy(ref.gradVg);
    newref.invVg = copy(ref.invVg);
    newref.Dr = copy(ref.Dr);
    newref.Ds = copy(ref.Ds);
    newref.Dt = copy(ref.Dt);
    newref.Dg = copy(ref.Dg);
    newref.Q1d = copy(ref.Q1d);
    newref.Q = copy(ref.Q);
    newref.Qr = copy(ref.Qr);
    newref.Qs = copy(ref.Qs);
    newref.Qt = copy(ref.Qt);
    newref.Ddr = copy(ref.Ddr);
    newref.Dds = copy(ref.Dds);
    newref.Ddt = copy(ref.Ddt);
    
    return newref;
end

function build_refel(dimension, order, nfaces, nodetype)
    # Check for errors
    if (dimension == 1 && nfaces != 2) || (dimension > 1 && nfaces < dimension+1)
        println("Error: buildRefel(dimension, order, nfaces, nodetype), check for valid parameters.");
        return nothing;
    end
    
    # Number of points determined by order element type
    if (dimension == 0) 
        Np = 1; #point
        Nfp = []; #no faces 
        
    elseif (dimension == 1) # line segment
        Np = order+1; 
        Nfp = [1, 1]; # face points
        
    elseif (dimension == 2 && nfaces == 3) # triangle
        Np = (Int)( ( order + 1 )*( order + 2 )/2);
        Nfp = [order + 1, order + 1, order + 1]; # face lines
        
    elseif (dimension == 2 && nfaces == 4) # quad
        Np = (Int)( (order + 1) * (order + 1) );
        Nfp = [order + 1, order + 1, order + 1, order + 1]; # face lines
        
    elseif (dimension == 3 && nfaces == 4)  # tet
        Np = (Int)((order + 1)*(order + 2)*(order + 3)/6);
        M = (Int)((order + 1)*(order + 2)/2);
        Nfp = [M, M, M, M]; # face triangles
        
    elseif (dimension == 3 && nfaces == 6)  # hex
        Np = (Int)((order + 1)*(order + 1)*(order + 1)); # hex
        M = (Int)((order + 1)*(order + 1));
        Nfp = [M, M, M, M, M, M]; # face quads
        
    elseif (dimension == 4)  # ??
        Np = (Int)((order + 1)*(order + 2)*(order + 3)*(order + 4)/24); # ??
        # TODO
    end
    
    refel = Refel( finch_state.config.float_type, dimension, order, Np, nfaces, Nfp);
    
    # Get nodes on the reference element
    refel_nodes!(refel, nodetype);
    refel.Nqp = length(refel.wg);
    
    #  0D refels
    if dimension == 0
        refel.V = ones(1,1);
        refel.gradV = zeros(1,1);
        refel.invV = ones(1,1);
        refel.Vg = ones(1,1);
        refel.gradVg = zeros(1,1);
        refel.invVg = ones(1,1);
        refel.Dr = zeros(1,1);
        refel.Dg = zeros(1,1);
        refel.Q1d = zeros(1,1);
        refel.Q = ones(1,1);
        refel.Qr = zeros(1,1);
        refel.Qs = zeros(1,1);
        refel.Qt = zeros(1,1);
        refel.Ddr = zeros(1,1);
        refel.Dds = zeros(1,1);
        refel.Ddt = zeros(1,1);
        
        return refel;
    end
    
    if (dimension == 2 && nfaces == 3) # triangle
        refel = build_triangle_refel(refel);
        
    elseif (dimension == 3 && nfaces == 4) # tet
        refel = build_tetrahedron_refel(refel);
        
    else # line, quad, hex
        # Vandermonde matrix and grad,inv
        # Values of basis functions and derivs at points
        refel.V = zeros(order + 1, order + 1);
        refel.gradV = zeros(order+1, order+1);
        # Gauss versions
        refel.Vg = zeros(order+1, order+1);
        refel.gradVg = zeros(order+1, order+1);
        # surface versions
        refel.surf_V = Array{Array{Float64}}(undef, nfaces);
        refel.surf_gradV = Array{Array{Float64}}(undef, nfaces);
        refel.surf_Vg = Array{Array{Float64}}(undef, nfaces);
        refel.surf_gradVg = Array{Array{Float64}}(undef, nfaces);
        for fi = 1:nfaces
            refel.surf_V[fi] = zeros(refel.Nfp[fi], Np);
            refel.surf_gradV[fi] = zeros(refel.Nfp[fi], Np);
            refel.surf_Vg[fi] = zeros(refel.Nfp[fi], Np);
            refel.surf_gradVg[fi] = zeros(refel.Nfp[fi], Np);
        end
        
        # nodal versions
        refel.V[:, :] = jacobi_polynomial( refel.r1d, 0, 0, refel.N, polynomialStartIdx = 1 )'

        # Previous Evaluation
        # for i = 1:refel.N + 1
        #     refel.V[:,i] = jacobi_polynomial(refel.r1d, 0, 0, i-1);
        # end

        polynomialOrders = collect( 1:refel.N )
        derivativePrefix = sqrt.( polynomialOrders .* ( polynomialOrders .+ 1 ) )

        refel.gradV[ :, 2:end ] = ( derivativePrefix .* jacobi_polynomial( refel.r1d, 1, 1, refel.N - 1, polynomialStartIdx = 1 ) )'

        # Previous Evaluation
        # for i=1:refel.N
        #     refel.gradV[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(refel.r1d, 1, 1, i-1);
        # end

        refel.invV = inv(refel.V);
        
        # Gauss versions
        refel.Vg = jacobi_polynomial(refel.g1d, 0, 0, refel.N, polynomialStartIdx = 1)';

        # Previous Evaluation
        # for i=1:refel.N+1
        #     refel.Vg[:,i] = jacobi_polynomial(refel.g1d, 0, 0, i-1);
        # end

        refel.gradVg[ :, 2:end ] = ( derivativePrefix .* jacobi_polynomial( refel.g1d, 1, 1, refel.N - 1, polynomialStartIdx = 1 ) )'

        # Previous Evaluation
        # for i=1:refel.N
        #     refel.gradVg[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(refel.g1d, 1, 1, i-1);
        # end
        refel.invVg = inv(refel.Vg);
        
        # Differentiation matrices
        refel.Dr = refel.gradV*refel.invV;
        refel.Dg = refel.gradVg*refel.invV;
        
        refel.Q1d = refel.Vg*refel.invV;
        
        if dimension == 1
            # volume
            refel.Q = refel.Q1d;
            refel.Qr = refel.Dg;
            refel.Ddr = refel.Dr;
            #surface
            # r and g are the same for this case because there's only one point
            refel.surf_Q = [refel.V[[1],:] * refel.invV, refel.V[[Np],:] * refel.invV];
            refel.surf_Qr = [refel.gradV[[1],:] * refel.invV, refel.gradV[[Np],:] * refel.invV];
            refel.surf_Ddr = refel.surf_Qr; 
            
        elseif dimension == 2
            # volume
            ident = Matrix( 1.0*I, order + 1, order + 1 );
            refel.Q = kron(refel.Q1d, refel.Q1d);
            refel.Qr = kron(refel.Q1d, refel.Dg);
            refel.Qs = kron(refel.Dg, refel.Q1d);
            refel.Ddr = kron(ident, refel.Dr);
            refel.Dds = kron(refel.Dr, ident);
            # surface
            refel.surf_Q = Array{Array{Float64}}(undef, nfaces);
            refel.surf_Qr = Array{Array{Float64}}(undef, nfaces);
            refel.surf_Qs = Array{Array{Float64}}(undef, nfaces);
            refel.surf_Ddr = Array{Array{Float64}}(undef, nfaces);
            refel.surf_Dds = Array{Array{Float64}}(undef, nfaces);
            ## Surfaces will use the full matrices rather than using tensor products
            fullinvV = kron(refel.invV, refel.invV);

            for fi = 1:nfaces
                surf_gradVr = zeros( refel.Nfp[fi], Np );
                surf_gradVs = zeros( refel.Nfp[fi], Np );
                surf_gradVgr = zeros( refel.Nfp[fi], Np );
                surf_gradVgs = zeros( refel.Nfp[fi], Np );
                for ni = 1 : size( refel.surf_r[fi], 1 )
                    # nodal versions

                    # Faster Evaluation
                    xEvals = jacobi_polynomial(refel.surf_r[fi][ ni, 1 ], 0, 0, refel.N, polynomialStartIdx = 1 )'
                    yEvals = jacobi_polynomial(refel.surf_r[fi][ ni, 2 ], 0, 0, refel.N, polynomialStartIdx = 1 )'

                    # X varies fastest in linearized kron operation
                    refel.surf_V[fi][ni, :] = kron( yEvals, xEvals )                      

                    # Previous Evaluation, X varies fastest in inner loop
                    # for i = 1 : refel.N + 1
                    #     for j = 1 : refel.N + 1
                    #         ind = (j - 1) * (refel.N + 1) + i;
                    #         refel.surf_V[fi][ni, ind] = (jacobi_polynomial(refel.surf_r[fi][ni,1], 0, 0, i-1) .* jacobi_polynomial(refel.surf_r[fi][ni,2], 0, 0, j-1))[1];
                    #     end
                    # end

                    # New Evaluation for surface derivatives
                    xderivEvals = ( derivativePrefix .* jacobi_polynomial(refel.surf_r[fi][ni, 1], 1, 1, refel.N - 1, polynomialStartIdx = 1) )'
                    yderivEvals = ( derivativePrefix .* jacobi_polynomial(refel.surf_r[fi][ni, 2], 1, 1, refel.N - 1, polynomialStartIdx = 1) )'
   
                    indicesToUpdate = [ ( j - 1 ) * (refel.N + 1) + i for j in 1:(refel.N + 1) for i in 2:(refel.N + 1) ]
                    surf_gradVr[ ni, indicesToUpdate ] = kron( yEvals, xderivEvals )

                    indicesToUpdate = [ ( j - 1 ) * (refel.N + 1) + i for j in 2:(refel.N + 1) for i in 1:(refel.N + 1) ]
                    surf_gradVs[ ni, indicesToUpdate ] = kron( yderivEvals, xEvals ) 

                    # Bug Fix for surface derivatives.
                    # Consider polynomial orders 0 -> N in both X and Y directions

                    # For multiplied X Y polynomials, if X direction varies fastest, 
                    # order in X direction varies from 0 -> N fastest with index

                    # Hence, if i == 1, polynomial derivative with X is 0 and needs to be skipped.
                    # For all other indexes, derivative value needs to be calculated regardless of j direction
                    
                    # Vice versa observation in Y direction
                    # for i = 1:refel.N + 1
                    #     for j = 1:refel.N + 1

                    #         ind = ( j - 1 ) * (refel.N + 1) + i;
                    #         ivalDeriv = i - 1
                    #         jvalDeriv = j - 1

                    #         if( i > 1 )  

                    #             surf_gradVr[ ni, ind ] = sqrt( ivalDeriv * ( ivalDeriv + 1 ) ) * 
                    #                 (jacobi_polynomial( refel.surf_r[fi][ ni, 1 ], 1, 1, ivalDeriv - 1 ) .*
                    #                 jacobi_polynomial( refel.surf_r[fi][ ni, 2 ], 0, 0, jvalDeriv ) )[1];

                    #         end

                    #         if( j > 1 )

                    #             surf_gradVs[ ni, ind ] = sqrt( jvalDeriv * ( jvalDeriv + 1 ) ) *
                    #             (jacobi_polynomial( refel.surf_r[fi][ ni, 1 ], 0, 0, ivalDeriv ) .*
                    #             jacobi_polynomial( refel.surf_r[fi][ ni, 2 ], 1, 1, jvalDeriv - 1) )[1];

                    #         end
                    #     end
                    # end
                    
                    # Previous Evaluation of surface derivatives (with bug)
                    # for i=1:refel.N
                    #     for j=1:refel.N
                    #         ind = (j)*(refel.N) + i + 1;
                    #         surf_gradVr[ni,ind] = sqrt(i*(i+1)) * (jacobi_polynomial(refel.surf_r[fi][ni,1], 1, 1, i-1) .* jacobi_polynomial(refel.surf_r[fi][ni,2], 0, 0, j-1))[1];
                    #         surf_gradVs[ni,ind] = sqrt(j*(j+1)) * (jacobi_polynomial(refel.surf_r[fi][ni,1], 0, 0, i-1) .* jacobi_polynomial(refel.surf_r[fi][ni,2], 1, 1, j-1))[1];
                    #     end
                    # end

                    # Faster Evaluation
                    xEvals = jacobi_polynomial(refel.surf_g[fi][ ni, 1 ], 0, 0, refel.N, polynomialStartIdx = 1 )'
                    yEvals = jacobi_polynomial(refel.surf_g[fi][ ni, 2 ], 0, 0, refel.N, polynomialStartIdx = 1 )'

                    # X varies fastest in linearized kron operation
                    refel.surf_Vg[fi][ni, :] = kron( yEvals, xEvals )

                    # Gauss versions, Previous Evaluation, X varies fastest in inner loop
                    # for i=1:refel.N+1
                    #     for j=1:refel.N+1
                    #         ind = (j-1)*(refel.N+1) + i;
                    #         refel.surf_Vg[fi][ni,ind] = (jacobi_polynomial(refel.surf_g[fi][ni,1], 0, 0, i-1) .* jacobi_polynomial(refel.surf_g[fi][ni,2], 0, 0, j-1))[1];
                    #     end
                    # end

                    # New Evaluation for surface derivatives
                    xderivEvals = ( derivativePrefix .* jacobi_polynomial(refel.surf_g[fi][ni, 1], 1, 1, refel.N - 1, polynomialStartIdx = 1) )'
                    yderivEvals = ( derivativePrefix .* jacobi_polynomial(refel.surf_g[fi][ni, 2], 1, 1, refel.N - 1, polynomialStartIdx = 1) )'
   
                    indicesToUpdate = [ ( j - 1 ) * (refel.N + 1) + i for j in 1:(refel.N + 1) for i in 2:(refel.N + 1) ]
                    surf_gradVgr[ ni, indicesToUpdate ] = kron( yEvals, xderivEvals )

                    indicesToUpdate = [ ( j - 1 ) * (refel.N + 1) + i for j in 2:(refel.N + 1) for i in 1:(refel.N + 1) ]
                    surf_gradVgs[ ni, indicesToUpdate ] = kron( yderivEvals, xEvals ) 

                    # Bug Fix for surface derivatives (Gauss Version).
                    # Consider polynomial orders 0 -> N in both X and Y directions

                    # For multiplied X Y polynomials, if X direction varies fastest, 
                    # order in X direction varies from 0 -> N fastest with index

                    # Hence, if i == 1, polynomial derivative with X is 0 and needs to be skipped.
                    # For all other indexes, derivative value needs to be calculated regardless of j direction
                    
                    # Vice versa observation in Y direction
                    # for i = 1:refel.N + 1
                    #     for j = 1:refel.N + 1

                    #         ind = ( j - 1 ) * (refel.N + 1) + i;
                    #         ivalDeriv = i - 1
                    #         jvalDeriv = j - 1

                    #         if( i > 1 )  

                    #             surf_gradVgr[ ni, ind ] = sqrt( ivalDeriv * ( ivalDeriv + 1 ) ) * 
                    #                 (jacobi_polynomial( refel.surf_g[fi][ ni, 1 ], 1, 1, ivalDeriv - 1 ) .*
                    #                 jacobi_polynomial( refel.surf_g[fi][ ni, 2 ], 0, 0, jvalDeriv ) )[1];

                    #         end

                    #         if( j > 1 )

                    #             surf_gradVgs[ ni, ind ] = sqrt( jvalDeriv * ( jvalDeriv + 1 ) ) *
                    #             (jacobi_polynomial( refel.surf_g[fi][ ni, 1 ], 0, 0, ivalDeriv ) .*
                    #             jacobi_polynomial( refel.surf_g[fi][ ni, 2 ], 1, 1, jvalDeriv - 1) )[1];

                    #         end
                    #     end
                    # end

                    # Gauss Versions, Previous Evaluation (With Bug)
                    # for i=1:refel.N
                    #     for j=1:refel.N
                    #         ind = (j)*(refel.N) + i + 1;
                    #         surf_gradVgr[ni,ind] = sqrt(i*(i+1)) * (jacobi_polynomial(refel.surf_g[fi][ni,1], 1, 1, i-1) .* jacobi_polynomial(refel.surf_g[fi][ni,2], 0, 0, j-1))[1];
                    #         surf_gradVgs[ni,ind] = sqrt(j*(j+1)) * (jacobi_polynomial(refel.surf_g[fi][ni,1], 0, 0, i-1) .* jacobi_polynomial(refel.surf_g[fi][ni,2], 1, 1, j-1))[1];
                    #     end
                    # end
                end
                
                refel.surf_Q[fi] = refel.surf_V[fi] * fullinvV;
                refel.surf_Qr[fi] = surf_gradVgr * fullinvV;
                refel.surf_Qs[fi] = surf_gradVgs * fullinvV;
                refel.surf_Ddr[fi] = surf_gradVr * fullinvV;
                refel.surf_Dds[fi] = surf_gradVs * fullinvV;
            end
            
        elseif dimension == 3
            ident = Matrix(1.0*I,order+1,order+1);
            refel.Q = kron(kron(refel.Q1d, refel.Q1d), refel.Q1d);
            refel.Qr = kron(kron(refel.Q1d, refel.Q1d), refel.Dg);
            refel.Qs = kron(kron(refel.Q1d, refel.Dg), refel.Q1d);
            refel.Qt = kron(kron(refel.Dg, refel.Q1d), refel.Q1d);
            refel.Ddr = kron(kron(ident, ident), refel.Dg);
            refel.Dds = kron(kron(ident, refel.Dg), ident);
            refel.Ddt = kron(kron(refel.Dg, ident), ident);
            
            # surface
            refel.surf_Q = Array{Array{Float64}}(undef, nfaces);
            refel.surf_Qr = Array{Array{Float64}}(undef, nfaces);
            refel.surf_Qs = Array{Array{Float64}}(undef, nfaces);
            refel.surf_Qt = Array{Array{Float64}}(undef, nfaces);
            refel.surf_Ddr = Array{Array{Float64}}(undef, nfaces);
            refel.surf_Dds = Array{Array{Float64}}(undef, nfaces);
            refel.surf_Ddt = Array{Array{Float64}}(undef, nfaces);
            ## Surfaces will use the full matrices rather than using tensor products
            fullinvV = kron(refel.invV, kron(refel.invV, refel.invV));
            
            for fi=1:nfaces
                surf_gradVr = zeros(refel.Nfp[fi], Np);
                surf_gradVs = zeros(refel.Nfp[fi], Np);
                surf_gradVt = zeros(refel.Nfp[fi], Np);
                surf_gradVgr = zeros(refel.Nfp[fi], Np);
                surf_gradVgs = zeros(refel.Nfp[fi], Np);
                surf_gradVgt = zeros(refel.Nfp[fi], Np);
                for ni=1:size(refel.surf_r[fi],1)
                    # nodal versions
                    xEvals = jacobi_polynomial(refel.surf_r[fi][ ni, 1 ], 0, 0, refel.N, polynomialStartIdx = 1 )'
                    yEvals = jacobi_polynomial(refel.surf_r[fi][ ni, 2 ], 0, 0, refel.N, polynomialStartIdx = 1 )'
                    zEvals = jacobi_polynomial(refel.surf_r[fi][ ni, 3 ], 0, 0, refel.N, polynomialStartIdx = 1 )'

                    # Faster Evaluation 
                    refel.surf_V[fi][ ni, : ] = kron( zEvals, kron( yEvals, xEvals ) )

                    # Previous Evaluation
                    for i=1:refel.N+1
                        for j=1:refel.N+1
                            for k=1:refel.N+1
                                ind = (k-1)*(refel.N+1)*(refel.N+1) + (j-1)*(refel.N+1) + i;
                                refel.surf_V[fi][ni,ind] = (jacobi_polynomial(refel.surf_r[fi][ni,1], 0, 0, i-1) .* jacobi_polynomial(refel.surf_r[fi][ni,2], 0, 0, j-1) .* jacobi_polynomial(refel.surf_r[fi][ni,3], 0, 0, k-1))[1];
                            end
                        end
                    end

                    # New/Faster Evaluation for surface derivatives
                    xderivEvals = ( derivativePrefix .* jacobi_polynomial(refel.surf_r[fi][ni, 1], 1, 1, refel.N - 1, polynomialStartIdx = 1) )'
                    yderivEvals = ( derivativePrefix .* jacobi_polynomial(refel.surf_r[fi][ni, 2], 1, 1, refel.N - 1, polynomialStartIdx = 1) )'
                    zderivEvals = ( derivativePrefix .* jacobi_polynomial(refel.surf_r[fi][ni, 3], 1, 1, refel.N - 1, polynomialStartIdx = 1) )'

                    indicesToUpdate = [  ( k - 1 ) * (refel.N + 1) * (refel.N + 1) + ( j - 1 ) * (refel.N + 1) + i 
                                    for k in 1:(refel.N + 1) 
                                    for j in 1:(refel.N + 1) 
                                    for i in 2:(refel.N + 1) ]

                    surf_gradVr[ ni, indicesToUpdate ] = kron( zEvals, kron( yEvals, xderivEvals ) )

                    indicesToUpdate = [  ( k - 1 ) * (refel.N + 1) * (refel.N + 1) + ( j - 1 ) * (refel.N + 1) + i 
                                    for k in 1:(refel.N + 1) 
                                    for j in 2:(refel.N + 1) 
                                    for i in 1:(refel.N + 1) ]

                    surf_gradVs[ ni, indicesToUpdate ] = kron( zEvals, kron( yderivEvals, xEvals ) )
                    
                    indicesToUpdate = [  ( k - 1 ) * (refel.N + 1) * (refel.N + 1) + ( j - 1 ) * (refel.N + 1) + i 
                                    for k in 2:(refel.N + 1) 
                                    for j in 1:(refel.N + 1) 
                                    for i in 1:(refel.N + 1) ]

                    surf_gradVt[ ni, indicesToUpdate ] = kron( zderivEvals, kron( yEvals, xEvals ) )

                    # Correct Evaluation (Longer Version with bug fix for readability)
                    # for i = 1:refel.N + 1
                    #     for j = 1:refel.N + 1
                    #         for k = 1:refel.N + 1

                    #             ind =  ( k - 1 ) * (refel.N + 1) * (refel.N + 1) + ( j - 1 ) * (refel.N + 1) + i;
                    #             ivalDeriv = i - 1
                    #             jvalDeriv = j - 1
                    #             kvalDeriv = k - 1

                    #             if( i > 1 )

                    #                 surf_gradVr[ ni, ind ] = sqrt( ivalDeriv * ( ivalDeriv + 1 ) ) * 
                    #                 (jacobi_polynomial( refel.surf_r[fi][ ni, 1 ], 1, 1, ivalDeriv - 1 ) .*
                    #                 jacobi_polynomial( refel.surf_r[fi][ ni, 2 ], 0, 0, jvalDeriv ) .*
                    #                 jacobi_polynomial( refel.surf_r[fi][ ni, 3 ], 0, 0, kvalDeriv ))[1];

                    #             end

                    #             if( j > 1 )

                    #                 surf_gradVs[ ni, ind ] = sqrt( jvalDeriv * ( jvalDeriv + 1 ) ) *
                    #                 (jacobi_polynomial( refel.surf_r[fi][ ni, 1 ], 0, 0, ivalDeriv ) .*
                    #                 jacobi_polynomial( refel.surf_r[fi][ ni, 2 ], 1, 1, jvalDeriv - 1) .*
                    #                 jacobi_polynomial( refel.surf_r[fi][ ni, 3 ], 0, 0, kvalDeriv) )[1];

                    #             end

                    #             if( k > 1 )

                    #                 surf_gradVt[ ni, ind ] = sqrt( kvalDeriv * ( kvalDeriv + 1 ) ) *
                    #                 (jacobi_polynomial( refel.surf_r[fi][ ni, 1 ], 0, 0, ivalDeriv ) .*
                    #                 jacobi_polynomial( refel.surf_r[fi][ ni, 2 ], 0, 0, jvalDeriv ) .*
                    #                 jacobi_polynomial( refel.surf_r[fi][ ni, 3 ], 1, 1, kvalDeriv - 1 ) )[1];

                    #             end
                    #         end
                    #     end
                    # end

                    # Previous Evaluation (with bug)
                    # for i = 1:refel.N 
                    #     for j = 1:refel.N 
                    #         for k = 1:refel.N 
                    #              ind = (k)*(refel.N)*(refel.N) + (j)*(refel.N) + i + 1;
                    #             surf_gradVr[ni,ind] = sqrt(i*(i+1)) * (jacobi_polynomial(refel.surf_r[fi][ni,1], 1, 1, i-1) .* jacobi_polynomial(refel.surf_r[fi][ni,2], 0, 0, j-1) .* jacobi_polynomial(refel.surf_r[fi][ni,3], 0, 0, k-1))[1];
                    #             surf_gradVs[ni,ind] = sqrt(j*(j+1)) * (jacobi_polynomial(refel.surf_r[fi][ni,1], 0, 0, i-1) .* jacobi_polynomial(refel.surf_r[fi][ni,2], 1, 1, j-1) .* jacobi_polynomial(refel.surf_r[fi][ni,3], 0, 0, k-1))[1];
                    #             surf_gradVt[ni,ind] = sqrt(k*(k+1)) * (jacobi_polynomial(refel.surf_r[fi][ni,1], 0, 0, i-1) .* jacobi_polynomial(refel.surf_r[fi][ni,2], 0, 0, j-1) .* jacobi_polynomial(refel.surf_r[fi][ni,3], 1, 1, k-1))[1];     
                    #         end
                    #     end
                    # end

                    # Gauss versions
                    xEvals = jacobi_polynomial(refel.surf_g[fi][ ni, 1 ], 0, 0, refel.N, polynomialStartIdx = 1 )'
                    yEvals = jacobi_polynomial(refel.surf_g[fi][ ni, 2 ], 0, 0, refel.N, polynomialStartIdx = 1 )'
                    zEvals = jacobi_polynomial(refel.surf_g[fi][ ni, 3 ], 0, 0, refel.N, polynomialStartIdx = 1 )'

                    # Faster Evaluation 
                    refel.surf_Vg[fi][ ni, : ] = kron( zEvals, kron( yEvals, xEvals ) )

                    # Previous Evaluation (Gauss Version)
                    # for i=1:refel.N+1
                    #     for j=1:refel.N+1
                    #         for k=1:refel.N+1
                    #             ind = (k-1)*(refel.N+1)*(refel.N+1) + (j-1)*(refel.N+1) + i;
                    #             refel.surf_Vg[fi][ni,ind] = (jacobi_polynomial(refel.surf_g[fi][ni,1], 0, 0, i-1) .* jacobi_polynomial(refel.surf_g[fi][ni,2], 0, 0, j-1) .* jacobi_polynomial(refel.surf_g[fi][ni,3], 0, 0, k-1))[1];
                    #         end
                    #     end
                    # end

                    # New Evaluation for surface derivatives
                    xderivEvals = ( derivativePrefix .* jacobi_polynomial(refel.surf_g[fi][ni, 1], 1, 1, refel.N - 1, polynomialStartIdx = 1) )'
                    yderivEvals = ( derivativePrefix .* jacobi_polynomial(refel.surf_g[fi][ni, 2], 1, 1, refel.N - 1, polynomialStartIdx = 1) )'
                    zderivEvals = ( derivativePrefix .* jacobi_polynomial(refel.surf_g[fi][ni, 3], 1, 1, refel.N - 1, polynomialStartIdx = 1) )'

                    indicesToUpdate = [  ( k - 1 ) * (refel.N + 1) * (refel.N + 1) + ( j - 1 ) * (refel.N + 1) + i 
                                    for k in 1:(refel.N + 1) 
                                    for j in 1:(refel.N + 1) 
                                    for i in 2:(refel.N + 1) ]

                    surf_gradVgr[ ni, indicesToUpdate ] = kron( zEvals, kron( yEvals, xderivEvals ) )

                    indicesToUpdate = [  ( k - 1 ) * (refel.N + 1) * (refel.N + 1) + ( j - 1 ) * (refel.N + 1) + i 
                                    for k in 1:(refel.N + 1) 
                                    for j in 2:(refel.N + 1) 
                                    for i in 1:(refel.N + 1) ]

                    surf_gradVgs[ ni, indicesToUpdate ] = kron( zEvals, kron( yderivEvals, xEvals ) )
                    
                    indicesToUpdate = [  ( k - 1 ) * (refel.N + 1) * (refel.N + 1) + ( j - 1 ) * (refel.N + 1) + i 
                                    for k in 2:(refel.N + 1) 
                                    for j in 1:(refel.N + 1) 
                                    for i in 1:(refel.N + 1) ]

                    surf_gradVgt[ ni, indicesToUpdate ] = kron( zderivEvals, kron( yEvals, xEvals ) )

                    # Correct Evaluation (Longer Version with bug fix for readability)
                    # for i = 1:refel.N + 1
                    #     for j = 1:refel.N + 1
                    #         for k = 1:refel.N + 1

                    #             ind =  ( k - 1 ) * (refel.N + 1) * (refel.N + 1) + ( j - 1 ) * (refel.N + 1) + i;
                    #             ivalDeriv = i - 1
                    #             jvalDeriv = j - 1
                    #             kvalDeriv = k - 1

                    #             if( i > 1 )

                    #                 surf_gradVgr[ ni, ind ] = sqrt( ivalDeriv * ( ivalDeriv + 1 ) ) * 
                    #                 (jacobi_polynomial( refel.surf_g[fi][ ni, 1 ], 1, 1, ivalDeriv - 1 ) .*
                    #                 jacobi_polynomial( refel.surf_g[fi][ ni, 2 ], 0, 0, jvalDeriv ) .*
                    #                 jacobi_polynomial( refel.surf_g[fi][ ni, 3 ], 0, 0, kvalDeriv ))[1];

                    #             end

                    #             if( j > 1 )

                    #                 surf_gradVgs[ ni, ind ] = sqrt( jvalDeriv * ( jvalDeriv + 1 ) ) *
                    #                 (jacobi_polynomial( refel.surf_g[fi][ ni, 1 ], 0, 0, ivalDeriv ) .*
                    #                 jacobi_polynomial( refel.surf_g[fi][ ni, 2 ], 1, 1, jvalDeriv - 1) .*
                    #                 jacobi_polynomial( refel.surf_g[fi][ ni, 3 ], 0, 0, kvalDeriv) )[1];

                    #             end

                    #             if( k > 1 )

                    #                 surf_gradVgt[ ni, ind ] = sqrt( kvalDeriv * ( kvalDeriv + 1 ) ) *
                    #                 (jacobi_polynomial( refel.surf_g[fi][ ni, 1 ], 0, 0, ivalDeriv ) .*
                    #                 jacobi_polynomial( refel.surf_g[fi][ ni, 2 ], 0, 0, jvalDeriv ) .*
                    #                 jacobi_polynomial( refel.surf_g[fi][ ni, 3 ], 1, 1, kvalDeriv - 1 ) )[1];

                    #             end
                    #         end
                    #     end
                    # end

                    # Previous Evaluation (With bug, gauss version)
                    # for i=1:refel.N
                    #     for j=1:refel.N
                    #         for k=1:refel.N
                    #             ind = (k)*(refel.N)*(refel.N) + (j)*(refel.N) + i + 1;
                    #             surf_gradVgr[ni,ind] = sqrt(i*(i+1)) * (jacobi_polynomial(refel.surf_g[fi][ni,1], 1, 1, i-1) .* jacobi_polynomial(refel.surf_g[fi][ni,2], 0, 0, j-1) .* jacobi_polynomial(refel.surf_g[fi][ni,3], 0, 0, k-1))[1];
                    #             surf_gradVgs[ni,ind] = sqrt(j*(j+1)) * (jacobi_polynomial(refel.surf_g[fi][ni,1], 0, 0, i-1) .* jacobi_polynomial(refel.surf_g[fi][ni,2], 1, 1, j-1) .* jacobi_polynomial(refel.surf_g[fi][ni,3], 0, 0, k-1))[1];
                    #             surf_gradVgt[ni,ind] = sqrt(k*(k+1)) * (jacobi_polynomial(refel.surf_g[fi][ni,1], 0, 0, i-1) .* jacobi_polynomial(refel.surf_g[fi][ni,2], 0, 0, j-1) .* jacobi_polynomial(refel.surf_g[fi][ni,3], 1, 1, k-1))[1];
                    #         end
                    #     end
                    # end
                end
                
                refel.surf_Q[fi] = refel.surf_V[fi] * fullinvV;
                refel.surf_Qr[fi] = surf_gradVgr * fullinvV;
                refel.surf_Qs[fi] = surf_gradVgs * fullinvV;
                refel.surf_Qt[fi] = surf_gradVgt * fullinvV;
                refel.surf_Ddr[fi] = surf_gradVr * fullinvV;
                refel.surf_Dds[fi] = surf_gradVs * fullinvV;
                refel.surf_Ddt[fi] = surf_gradVt * fullinvV;
            end
        end
    end
    
    return refel;
end

function custom_quadrature_refel(oldrefel, nodes, weights)
    # The supplied nodes will replace the gauss quadrature nodes and weights r.g and r.wg
    # NOTE: Elemental node and weight arrays will not be changed, so quadrature must be set to GAUSS.
    # NOTE: Only quadrature matrices are changed. Nothing else.
    refel = copy(oldrefel);
    refel.g = copy(nodes);
    refel.wg = copy(weights);
    
    # Vandermonde matrix and grad,inv for r will not change
    
    # Gauss versions will change
    Nn = size(nodes,2);
    refel.Vg = zeros(Nn, refel.N+1);
    refel.gradVg = zeros(Nn, refel.N+1);
    if refel.dim == 1
        for i=1:refel.N+1
            refel.Vg[:,i] = jacobi_polynomial(nodes[1,:], 0, 0, i-1);
        end
        for i=1:refel.N
            refel.gradVg[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[1,:], 1, 1, i-1);
        end
        
        # Differentiation matrices
        refel.Dg = refel.gradVg*refel.invV;
        refel.Q1d = refel.Vg*refel.invV;
        
        # Quadrature matrices
        refel.Q = refel.Q1d;
        refel.Qr = refel.Dg;
        refel.Ddr = refel.Dr;
        
    elseif refel.dim == 2
        # Now it's a little trickier, so forget the efficiency and do it carefully.
        Nn = size(nodes,2);
        tmp1 = zeros(Nn, refel.N+1);
        tmp2 = zeros(Nn, refel.N+1);
        Dtmp1 = zeros(Nn, refel.N+1);
        Dtmp2 = zeros(Nn, refel.N+1);
        
        for i=1:refel.N+1
            tmp1[:,i] = jacobi_polynomial(nodes[1,:], 0, 0, i-1);
            tmp2[:,i] = jacobi_polynomial(nodes[2,:], 0, 0, i-1);
        end
        for i=1:refel.N
            Dtmp1[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[1,:], 1, 1, i-1);
            Dtmp2[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[2,:], 1, 1, i-1);
        end
        
        # Use kron(Vg*invV, Vg,invV) = kron(Vg,Vg)*kron(invV,invV)
        VgXVg = zeros(Nn, (refel.N+1)*(refel.N+1));
        VgXDg = zeros(Nn, (refel.N+1)*(refel.N+1));
        DgXVg = zeros(Nn, (refel.N+1)*(refel.N+1));
        for i=1:Nn
            for j=1:(refel.N+1)
                for k=1:(refel.N+1)
                    ind = (k-1)*(refel.N+1) + j;
                    VgXVg[i,ind] = tmp1[i,j]*tmp2[i,k];
                    VgXDg[i,ind] = Dtmp1[i,j]*tmp2[i,k];
                    DgXVg[i,ind] = tmp1[i,j]*Dtmp2[i,k];
                end
            end
        end
        viXvi = kron(refel.invV, refel.invV);
        
        refel.Q = VgXVg * viXvi;
        refel.Qr = VgXDg * viXvi;
        refel.Qs = DgXVg * viXvi;
        
    elseif refel.dim == 3
        Nn = size(nodes,2);
        tmp1 = zeros(Nn, refel.N+1);
        tmp2 = zeros(Nn, refel.N+1);
        tmp3 = zeros(Nn, refel.N+1);
        Dtmp1 = zeros(Nn, refel.N+1);
        Dtmp2 = zeros(Nn, refel.N+1);
        Dtmp3 = zeros(Nn, refel.N+1);
        
        for i=1:refel.N+1
            tmp1[:,i] = jacobi_polynomial(nodes[1,:], 0, 0, i-1);
            tmp2[:,i] = jacobi_polynomial(nodes[2,:], 0, 0, i-1);
            tmp3[:,i] = jacobi_polynomial(nodes[3,:], 0, 0, i-1);
        end
        for i=1:refel.N
            Dtmp1[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[1,:], 1, 1, i-1);
            Dtmp2[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[2,:], 1, 1, i-1);
            Dtmp3[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[3,:], 1, 1, i-1);
        end
        
        viXviXvi = kron(kron(refel.invV, refel.invV), refel.invV);
        
        VgXVgXVg = zeros(Nn, (refel.N+1)*(refel.N+1)*(refel.N+1));
        VgXVgXDg = zeros(Nn, (refel.N+1)*(refel.N+1)*(refel.N+1));
        VgXDgXVg = zeros(Nn, (refel.N+1)*(refel.N+1)*(refel.N+1));
        DgXVgXVg = zeros(Nn, (refel.N+1)*(refel.N+1)*(refel.N+1));
        for ni=1:Nn
            for i=1:(refel.N+1)
                for j=1:(refel.N+1)
                    for k=1:(refel.N+1)
                        ind = ((k-1)*(refel.N+1) + (j-1))*(refel.N+1) + i;
                        VgXVgXVg[ni,ind] = tmp1[ni,i]*tmp2[ni,j]*tmp2[ni,k];
                        VgXVgXDg[ni,ind] = Dtmp1[ni,i]*tmp2[ni,j]*tmp2[ni,k];
                        VgXDgXVg[ni,ind] = tmp1[ni,i]*Dtmp2[ni,j]*tmp2[ni,k];
                        DgXVgXVg[ni,ind] = tmp1[ni,i]*tmp2[ni,j]*Dtmp2[ni,k];
                    end
                end
            end
        end
        
        refel.Q = VgXVgXVg * viXviXvi;
        refel.Qr = VgXVgXDg * viXviXvi;
        refel.Qs = VgXDgXVg * viXviXvi;
        refel.Qt = DgXVgXVg * viXviXvi;
    end
    
    return refel;
    
end

# Map 1D nodes to 2D face
# Assumes interval is scaled the same 
function map_face_nodes_2d(g1d, v1, v2)
    fnodes = zeros(2,length(g1d));
    dx = v2[1] - v1[1];
    dy = v2[2] - v1[2];
    dist = sqrt(dx*dx+dy*dy);
    xmult = dx/dist;
    ymult = dy/dist;
    
    fnodes[1,:] = v1[1] .+ (g1d .- v1[1]) .*xmult;
    fnodes[2,:] = v1[2] .+ (g1d .- v1[2]) .*ymult;
    
    return fnodes;
end

# Map 2D nodes to 3D face
# Assumes g2d[1] aligns with v1 and interval is scaled the same 
function map_face_nodes_3d(g2d, v1, v2, v3, v4=[])
    fnodes = zeros(3,length(g2d));
    if length(v4) > 0 # quads
        dx = 0
        dy = 0
        dz = 0
        dist = sqrt(dx*dx+dy*dy+dz*dz);
        xmult = dx/dist;
        ymult = dy/dist;
        zmult = dz/dist;
        
        fnodes[1,:] = v1[1] .+ (g2d .- v1[1]) .*xmult;
        fnodes[2,:] = v1[2] .+ (g2d .- v1[2]) .*ymult;
        fnodes[3,:] = v1[3] .+ (g2d .- v1[3]) .*zmult;
    else #triangles
        
    end
    
    
    return fnodes;
end