# Polyharmonic functions used below.
function polyharmonic_phi1(r) return r; end
function polyharmonic_phi2(r) return r>1e-16 ? r^2 * log(r) : r; end
function polyharmonic_phi3(r) return r^3; end
function polyharmonic_phi4(r) return r>1e-16 ? r^4 * log(r) : r^3; end

# Uses polyharmonic interpolation to interpolate at x using values u at centers c.
# sizes are:
# x is d x nx
# c is d x nc
# u is nc
# (optional) ord determines the type of basis function. default is 3. 1,2,3,4 are available
function polyharmonic_interp(x, c, u, ord=3)
    dim = size(x,1);
    nx = size(x,2); # number of places to evaluate
    
    # If no centers are provided, this is meaningless, but return zero.
    if length(c)==0
        return zeros(nx);
    end
    
    nc = size(c,2); # number of centers
    
    # Constant if only one center is given
    if nc == 1
        return u[1] .* ones(nx);
    end
    
    if ord == 1
        phi = polyharmonic_phi1;
    elseif ord == 2
        phi = polyharmonic_phi2;
    elseif ord == 3
        phi = polyharmonic_phi3;
    else
        phi = polyharmonic_phi4;
    end
    
    # build A matrix
    A = zeros(nc + dim+1, nc + dim+1);
    for i=1:nc
        for j=1:nc
            A[i,j] = phi(norm(c[:,j] - c[:,i]));
        end
    end
    # build B matrix (within A matrix)
    A[1:nc, nc+1] .= 1;
    A[nc+1, 1:nc] .= 1;
    A[1:nc, (nc+2):(nc+dim+1)] = c';
    A[(nc+2):(nc+dim+1), 1:nc] = c;
    
    # build RHS vector
    b = vcat(u, zeros(dim+1));
    
    wv = A\b; # coefficients for the polyharmonic
    
    sol = zeros(nx);
    for i=1:nx
        # evaluate at x[i]
        sol[i] = wv[nc+1];
        for k=1:dim
            sol[i] += wv[nc+k+1] * x[k,i];
        end
        for j=1:nc
            sol[i] += wv[j]*phi(norm(c[:,j] - x[:,i]));
        end
    end
    
    return sol;
end