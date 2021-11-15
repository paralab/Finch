# Polyharmonic functions used below.
# function polyharmonic_phi1(r) return r; end
# function polyharmonic_phi2(r) return r>1e-16 ? r^2 * log(r) : r; end
# function polyharmonic_phi3(r) return r^3; end
# function polyharmonic_phi4(r) return r>1e-16 ? r^4 * log(r) : r^3; end

# Uses polyharmonic interpolation to interpolate at x using values u at centers c.
# sizes are:
# x is d x nx
# c is d x nc
# u is nc
# (optional) ord determines the type of basis function. default is 3. 1,2,3,4 are available
# Note: This now redirects to a dimension specific function
function polyharmonic_interp(x, c, u; ord=3, limiter=0, neighbor_vals=[])
    dim = size(x,1);
    if dim==1
        return polyharmonic_interp_1d(x, c, u, ord=ord, limiter=limiter, neighbor_vals=neighbor_vals);
    elseif dim==2
        return polyharmonic_interp_2d(x, c, u, ord=ord, limiter=limiter, neighbor_vals=neighbor_vals);
    else
        return polyharmonic_interp_3d(x, c, u, ord=ord, limiter=limiter, neighbor_vals=neighbor_vals);
    end
end

# Uses polyharmonic interpolation to interpolate at x using values u at centers c.
# See above for input details.
function polyharmonic_interp_1d(x, c, u; ord=3, limiter=0, neighbor_vals=[])
    dim = 1;
    nx = size(x,2); # number of places to evaluate
    
    # If no centers are provided, this is meaningless, but return zero.
    if length(c)==0
        return zeros(nx);
    end
    
    nc = size(c,2); # number of centers
    
    # Constant if only one center is given
    if nc == 1
        if nx==1
            return u;
        else
            return u[1] .* ones(nx);
        end
        
    elseif nc == 2 # Linear if two centers are given
        du = u[2] - u[1];
        dx1 = c[2] - c[1];
        dx2 = x .- c[1];
        if abs(dx1) > 1e-20
            sol = u[1] .+ (du/dx1 .* dx2);
        else # overlapping points
            sol = u[1] .* ones(nx);
        end
        
    else # Use polyharmonics if more than two centers are given
        # build A matrix
        A = zeros(nc + dim+1, nc + dim+1);
        # Handle each order separately for efficiency
        if ord == 1
            for i=1:nc
                for j=(i+1):nc
                    A[i,j] = norm(c[:,j] - c[:,i]); # phi(norm(c[:,j] - c[:,i]))
                    A[j,i] = A[i,j];
                end
            end
        elseif ord == 2
            for i=1:nc
                for j=(i+1):nc
                    # A[i,j] = norm(c[:,j] - c[:,i]); # phi(norm(c[:,j] - c[:,i]))
                    tmp = sum((c[:,j] - c[:,i]).^2);
                    if tmp > 1e-30
                        A[i,j] = tmp * log(tmp) * 0.5;
                    else
                        A[i,j] = sqrt(tmp);
                    end
                    A[j,i] = A[i,j];
                end
            end
        elseif ord == 3
            for i=1:nc
                for j=(i+1):nc
                    A[i,j] = (norm(c[:,j] - c[:,i]))^3; # phi(norm(c[:,j] - c[:,i]))
                    A[j,i] = A[i,j];
                end
            end
        else # max order is 4
            for i=1:nc
                for j=(i+1):nc
                    # A[i,j] = norm(c[:,j] - c[:,i]); # phi(norm(c[:,j] - c[:,i]))
                    tmp = sum((c[:,j] - c[:,i]).^2);
                    if tmp > 1e-30
                        A[i,j] = tmp*tmp * log(tmp) * 0.5;
                    else
                        A[i,j] = tmp*sqrt(tmp);
                    end
                    A[j,i] = A[i,j];
                end
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
            if ord == 1
                for j=1:nc
                    sol[i] += wv[j]*norm(c[:,j] - x[:,i]); # wv[j]*phi(norm(c[:,j] - x[:,i]));
                end
            elseif ord == 2
                for j=1:nc
                    # sol[i] += wv[j]*phi(norm(c[:,j] - x[:,i]));
                    tmp = sum((c[:,j] - x[:,i]).^2);
                    if tmp > 1e-30
                        sol[i] += wv[j] * tmp * log(tmp) * 0.5;
                    else
                        sol[i] += wv[j] * sqrt(tmp);
                    end
                end
            elseif ord == 3
                for j=1:nc
                    sol[i] += wv[j]*norm(c[:,j] - x[:,i])^3;
                end
            else # 4+
                for j=1:nc
                    # sol[i] += wv[j]*phi(norm(c[:,j] - x[:,i]));
                    tmp = sum((c[:,j] - x[:,i]).^2);
                    if tmp > 1e-30
                        sol[i] += wv[j] * tmp * tmp * log(tmp) * 0.5;
                    else
                        sol[i] += wv[j] * tmp*sqrt(tmp);
                    end
                end
            end
        end # building sol
    end
    # At this point sol contains the extrapolated value
    
    # If a limiting is desired, limit the extrapolated values to the range of the neighboring cell values
    if limiter > 0 && length(neighbor_vals) == 2
        maxv = max(neighbor_vals[1], neighbor_vals[2]);
        minv = min(neighbor_vals[1], neighbor_vals[2]);
        
        for i=1:nx
            #sol[i] = max(sol[i], minv);
            #sol[i] = min(sol[i], maxv);
            if sol[i] < minv
                sol[i] = limiter * minv + (1-limiter) * sol[i];
            elseif sol[i] > maxv
                sol[i] = limiter * maxv + (1-limiter) * sol[i];
            end
        end
    end
    
    return sol;
end

# Uses polyharmonic interpolation to interpolate at x using values u at centers c.
# See above polyharmonic_interp() for input details
function polyharmonic_interp_2d(x, c, u; ord=3, limiter=0, neighbor_vals=[])
    dim = 2;
    nx = size(x,2); # number of places to evaluate
    
    # If no centers are provided, this is meaningless, but return zero.
    if length(c)==0
        return zeros(nx);
    end
    
    nc = size(c,2); # number of centers
    
    # Constant if only one center is given
    if nc == 1
        if nx==1
            return u;
        else
            return u[1] .* ones(nx);
        end
        
    elseif nc == 2 # NOT A UNIQUE INTERPOLATION, so assume the gradient is parallel to the vector between the two points
        sol = linear_interp(x, c, u, 2);
        
    else # Use polyharmonics if more than two centers are given and at least three create a nondegenerate triangle
        # First check for nondegeneracy
        degenerate = true;
        v1 = c[:,2] - c[:,1];
        for i=3:nc
            v2 = c[:,i] - c[:,1];
            if (v1[1]*v2[2] - v1[2]*v2[1]) > 1e-20 # rotate v2 by pi/2 and dot
                degenerate = false;
            end
        end
        
        if degenerate
            sol = linear_interp(x, c[:,1:2], u[1:2], 2);
            
        else
            # build A matrix
            A = zeros(nc + dim+1, nc + dim+1);
            # Handle each order separately for efficiency
            if ord == 1
                for i=1:nc
                    for j=(i+1):nc
                        A[i,j] = norm(c[:,j] - c[:,i]); # phi(norm(c[:,j] - c[:,i]))
                        A[j,i] = A[i,j];
                    end
                end
            elseif ord == 2
                for i=1:nc
                    for j=(i+1):nc
                        # A[i,j] = norm(c[:,j] - c[:,i]); # phi(norm(c[:,j] - c[:,i]))
                        tmp = sum((c[:,j] - c[:,i]).^2);
                        if tmp > 1e-30
                            A[i,j] = tmp * log(tmp) * 0.5;
                        else
                            A[i,j] = sqrt(tmp);
                        end
                        A[j,i] = A[i,j];
                    end
                end
            elseif ord == 3
                for i=1:nc
                    for j=(i+1):nc
                        A[i,j] = (norm(c[:,j] - c[:,i]))^3; # phi(norm(c[:,j] - c[:,i]))
                        A[j,i] = A[i,j];
                    end
                end
            else # max order is 4
                for i=1:nc
                    for j=(i+1):nc
                        # A[i,j] = norm(c[:,j] - c[:,i]); # phi(norm(c[:,j] - c[:,i]))
                        tmp = sum((c[:,j] - c[:,i]).^2);
                        if tmp > 1e-30
                            A[i,j] = tmp*tmp * log(tmp) * 0.5;
                        else
                            A[i,j] = tmp*sqrt(tmp);
                        end
                        A[j,i] = A[i,j];
                    end
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
                if ord == 1
                    for j=1:nc
                        sol[i] += wv[j]*norm(c[:,j] - x[:,i]); # wv[j]*phi(norm(c[:,j] - x[:,i]));
                    end
                elseif ord == 2
                    for j=1:nc
                        # sol[i] += wv[j]*phi(norm(c[:,j] - x[:,i]));
                        tmp = sum((c[:,j] - x[:,i]).^2);
                        if tmp > 1e-30
                            sol[i] += wv[j] * tmp * log(tmp) * 0.5;
                        else
                            sol[i] += wv[j] * sqrt(tmp);
                        end
                    end
                elseif ord == 3
                    for j=1:nc
                        sol[i] += wv[j]*norm(c[:,j] - x[:,i])^3;
                    end
                else # 4+
                    for j=1:nc
                        # sol[i] += wv[j]*phi(norm(c[:,j] - x[:,i]));
                        tmp = sum((c[:,j] - x[:,i]).^2);
                        if tmp > 1e-30
                            sol[i] += wv[j] * tmp * tmp * log(tmp) * 0.5;
                        else
                            sol[i] += wv[j] * tmp*sqrt(tmp);
                        end
                    end
                end
            end # building sol
        end
        
    end
    # At this point sol contains the extrapolated value
    
    # If a limiting is desired, limit the extrapolated values to the range of the neighboring cell values
    if limiter > 0 && length(neighbor_vals) == 2
        maxv = max(neighbor_vals[1], neighbor_vals[2]);
        minv = min(neighbor_vals[1], neighbor_vals[2]);
        
        for i=1:nx
            #sol[i] = max(sol[i], minv);
            #sol[i] = min(sol[i], maxv);
            if sol[i] < minv
                sol[i] = limiter * minv + (1-limiter) * sol[i];
            elseif sol[i] > maxv
                sol[i] = limiter * maxv + (1-limiter) * sol[i];
            end
        end
    end
    
    return sol;
end

# Uses polyharmonic interpolation to interpolate at x using values u at centers c.
# See above polyharmonic_interp() for input details
function polyharmonic_interp_3d(x, c, u; ord=3, limiter=0, neighbor_vals=[])
    dim = 2;
    nx = size(x,2); # number of places to evaluate
    
    # If no centers are provided, this is meaningless, but return zero.
    if length(c)==0
        return zeros(nx);
    end
    
    nc = size(c,2); # number of centers
    
    # Constant if only one center is given
    if nc == 1
        if nx==1
            return u;
        else
            return u[1] .* ones(nx);
        end
        
    else
        #TODO
        sol = u[1] .* ones(nx);
    end
    
    return sol;
end

# Linear interpolation used when polyharmonic won't work.
# It is assumed that there are two centers
function linear_interp(x, c, u, dim)
    nx = size(x,2); # number of places to evaluate
    if dim==1
        du = u[2] - u[1];
        dx1 = c[2] - c[1];
        dx2 = x .- c[1];
        if abs(dx1) > 1e-20
            sol = u[1] .+ (du/dx1 .* dx2);
        else # overlapping points
            sol = u[1] .* ones(nx);
        end
        
    elseif dim==2 # NOT A UNIQUE INTERPOLATION, so assume the gradient is parallel to the vector between the two points
        du = u[2] - u[1];
        dx1 = [c[1,2] - c[1,1]; c[2,2] - c[2,1]];
        if nx == 1
            dx2 = [x[1,1] - c[1,1]; x[2,1] - c[2,1]];
        else
            dx2 = similar(x);
            for xi=1:nx
                dx2[:,xi] = [x[1,xi] - c[1,1]; x[2,xi] - c[2,1]];
            end
        end
        if abs(dx1[1])+abs(dx1[2]) > 1e-20
            # u + grad(u) . dx2
            if nx == 1
                sol = [u[1] .+ (du/(dx1^2 + dx2^2) .* (dx1[1]*dx2[1] + dx1[2]*dx2[2]))];
            else
                sol = zeros(nx);
                for xi=1:nx
                    sol[xi] = u[1] .+ (du/(dx1[1]^2 + dx1[2]^2) .* (dx1[1]*dx2[1,xi] + dx1[2]*dx2[2,xi]));
                end
            end
        else # overlapping points
            sol = u[1] .* ones(nx);
        end
        
    else # 3D
        #TODO
        sol = u[1] .* ones(nx);
    end
    
    return sol;
end