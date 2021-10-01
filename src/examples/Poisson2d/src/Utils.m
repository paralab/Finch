
classdef Utils
    
    methods(Static)
        % 2D routines
        function y = tensor_IAX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N);
            y = y(:);
        end

        function y = tensor_AIX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N)';
            y = y'; 
            y = y(:);
        end

        % 3D routines
        function y = tensor_IIAX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N*N);
            y = y(:);
        end

        function y = tensor_IAIX (A, x)
            N = size (A, 1);
            q = reshape(x, N, N, N);
            y = zeros(N,N,N);
            for i=1:N
                y(i,:,:) = A * squeeze( q(i,:,:) );
            end
            y = y(:);
        end

        function y = tensor_AIIX (A, x)
            N = size (A, 1);
            y = reshape(x, N*N, N) * A';
            y = y(:);
        end

        function du = tensor_grad(refel, u)
            du = zeros(length(u), refel.dim);
            if (refel.dim == 2)
                du(:,1) = Utils.tensor_IAX (refel.Dr, u);
                du(:,2) = Utils.tensor_AIX (refel.Dr, u);
            else
                du(:,1) = Utils.tensor_IIAX (refel.Dr, u);
                du(:,2) = Utils.tensor_IAIX (refel.Dr, u);
                du(:,3) = Utils.tensor_AIIX (refel.Dr, u);
            end
        end

        function [dx, dy] = tensor_grad2(A, x)
            dx = Utils.tensor_IAX (A, x);
            dy = Utils.tensor_AIX (A, x);
        end

        function [dx, dy, dz] = tensor_grad3(A, x)
            dx = Utils.tensor_IIAX (A, x);
            dy = Utils.tensor_IAIX (A, x);
            dz = Utils.tensor_AIIX (A, x);
        end


        function [J, D] = geometric_factors(refel, pts)

            if (refel.dim == 0)
                xr  = [1];
                J = xr;
            elseif (refel.dim == 1)
                xr  = refel.Dr*pts;
                J = xr;
            elseif (refel.dim == 2)
                if refel.Nfaces == 3 % triangle
                    xr = refel.Ddr*pts(1,:)';
                    xs = refel.Dds*pts(1,:)';
                    yr = refel.Ddr*pts(2,:)';
                    ys = refel.Dds*pts(2,:)';
                    J = -xs.*yr + xr.*ys;
                    J = J(1); 
                else % quad
                    [xr, xs] = Utils.tensor_grad2 (refel.Dg, pts(1,:));
                    [yr, ys] = Utils.tensor_grad2 (refel.Dg, pts(2,:));

                    J = -xs.*yr + xr.*ys;
                end

            else
                [xr, xs, xt] = Utils.tensor_grad3 (refel.Dg, pts(1,:));
                [yr, ys, yt] = Utils.tensor_grad3 (refel.Dg, pts(2,:));
                [zr, zs, zt] = Utils.tensor_grad3 (refel.Dg, pts(3,:));

                J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
            end

            if (nargout > 1)
                if (refel.dim == 1)
                    D.rx = 1./J;
                elseif (refel.dim == 2)
                    if refel.Nfaces == 3 % triangle
                        D.rx =  ys./J;
                        D.sx = -yr./J;
                        D.ry = -xs./J;
                        D.sy =  xr./J;
                    else % quad
                        D.rx =  ys./J;
                        D.sx = -yr./J;
                        D.ry = -xs./J;
                        D.sy =  xr./J;
                    end

                else
                    D.rx =  (ys.*zt - zs.*yt)./J;
                    D.ry = -(xs.*zt - zs.*xt)./J;
                    D.rz =  (xs.*yt - ys.*xt)./J;

                    D.sx = -(yr.*zt - zr.*yt)./J;
                    D.sy =  (xr.*zt - zr.*xt)./J;
                    D.sz = -(xr.*yt - yr.*xt)./J;

                    D.tx =  (yr.*zs - zr.*ys)./J;
                    D.ty = -(xr.*zs - zr.*xs)./J;
                    D.tz =  (xr.*ys - yr.*xs)./J;
                end

            end
        end
        
        % pack and unpack variables into a global vector
        % packed vars are in a vector, unpacked are in a cell table
        function pk = pack_vars(upk, pk)
            for i=1:Nvars
                pk(i:Nvars:length(pk)) = upk{i};
            end
        end
        
        function upk = unpack_vars(pk)
            upk = cell(Nvars,1);
            for i=1:Nvars
                upk{i} = pk(i:Nvars:length(pk));
            end
        end
    end
end


