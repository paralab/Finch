#=
# Based on the file JacobiGQ.m from the book:
# Nodal Discontinuous Galerkin Method - Hesthaven, Warburton
#  https://link.springer.com/book/10.1007/978-0-387-72067-8
=#

# Compute the Nth order Gauss quadrature points, x, and weights, w,
# associated with the Jacobi polynomial of type (alpha,beta).
function jacobi_gauss_quad(alpha::Int, beta::Int, N::Int)
    if N==0
        x = [-(alpha-beta)/(alpha+beta+2)];
        w = [2];
        return (x,w);
    end
    
    function sqrtvec(v::Array)
        w = zeros(size(v));
        for i=1:length(v)
            w[i] = sqrt(v[i]);
        end
        return w;
    end
    
    # Form symmetric matrix from recurrence.
    #J = zeros(N+1,N+1);
    h1 = 2*(0:N).+(alpha+beta);
    J = diagm(-1/2*(alpha^2-beta^2)./(h1.+2)./h1) +
        diagm(1 => 2 ./(h1[1:N].+2).*sqrtvec((1:N).*((1:N).+(alpha+beta)).*
        ((1:N).+alpha).*((1:N).+beta)./(h1[1:N].+1)./(h1[1:N].+3)));
    if (alpha+beta == 0)
        J[1,1]=0.0;
    end
    J = J + J';
    
    # Find eigenvalues and vectors
    E = eigen(J);
    #x = E.values;
    w = (E.vectors[1,:]').^2 .*2^(alpha+beta+1)/(alpha+beta+1)*factorial(alpha)*
        factorial(beta)/factorial(alpha+beta);
    w = w[:];
    return (E.values,w);
end

function jacobi_LGL_quad(N::Int)
    if N == 1
        lgl = [-1; 1];
        w = [1; 1];
    else
        # compute the points
        (r,w) = jacobi_gauss_quad(1,1,N-2);
        lgl = [-1; r ; 1];
        
        # compute the weights
        w = jacobi_polynomial(lgl, 0, 0, N);
        adgammaN = (2*N + 1) / (N * (N + 1));
        w = w.*w;
        w = adgammaN./w;
    end
    
    return (lgl, w);
end