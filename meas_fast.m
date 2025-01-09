function [y,dy] = meas_fast(nw,x,div)

%Generate nonlinear map of states and then the estimated output
fx = x(nw.i1,:).*x(nw.i2,:);
y  = [x ; nw.y0 + nw.A*x + nw.B*fx];

if nargout > 1
    ny = size(y,1);
    N  = size(y,2);
    nx = size(x,1);
    if div==1 %Generate derivatives of y wrt x (for all times)
        dydx = zeros(ny,nx,N);
        on = ones(1,N);
        for i=1:nx
            dxdx = zeros(nx,N);
            dxdx(i,:) = 1;
            dx1 = nw.i1==i;
            dx2 = nw.i2==i;
            dfdx = dx1.*x(nw.i2,:) + x(nw.i1,:).*dx2;
            dydx(:,i,:) = [dxdx; nw.A*dxdx + nw.B*dfdx];
        end
        dy=dydx;
    else %Generate derivatives of y wrt theta
        nt = 2*nw.num_lines;
        dydt = zeros(ny*N,nt);
        d    = zeros(size(y));
        for i=1:nt
            dA = sparse(nw.dA(:,:,i));
            dB = sparse(nw.dB(:,:,i));
            d(nx+1:end,:) = dA*x + dB*fx;
            dydt(:,i) = d(:);
        end
        dy=dydt;
    end
end