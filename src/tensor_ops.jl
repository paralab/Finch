#=
# Some things based on HOMG Tensor class
=#

# 2D routines
function tensor_IAX(A::Matrix, x::Vector)
    N = size(A)[1];
    y = A * reshape(x, N, N);
    return reshape(y, :);
end

function tensor_AIX(A, x)
    N = size(A)[1];
    y = A * reshape(x, N, N)';
    y = y';
    return reshape(y, :);
end

# 3D routines
function tensor_IIAX(A, x)
    N = size(A)[1];
    y = A * reshape(x, N, N*N);
    return reshape(y, :);
end

function tensor_IAIX(A, x)
    N = size(A)[1];
    q = reshape(x, N, N, N);
    y = zeros(typeof(x[1]), N,N,N);
    for i=1:N
        y[i,:,:] = A * q[i,:,:];
    end
    return reshape(y, :);
end

function tensor_AIIX(A, x)
    N = size(A)[1];
    y = reshape(x, N*N, N) * A';
    return reshape(y, :);
end

function tensor_grad2(A, x)
    dx = tensor_IAX(A, x);
    dy = tensor_AIX(A, x);
    return (dx, dy);
end

function tensor_grad3(A, x)
    dx = tensor_IIAX(A, x);
    dy = tensor_IAIX(A, x);
    dz = tensor_AIIX(A, x);
    return (dx, dy, dz);
end