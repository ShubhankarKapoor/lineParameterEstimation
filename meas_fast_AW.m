function [y,dydx,dydt,dydxdt] = meas_fast_AW(nw,x)

%Generate nonlinear map of states and then the estimated output
fx = x(nw.i1,:).*x(nw.i2,:);
y  = [x ; nw.y0 + nw.A*x + nw.B*fx];

%Generate derivatives of y wrt x (for all times)
%Derivatives of y wrt states?

ny = size(y,1);
N  = size(y,2);
nx = size(x,1);
nt = 2*nw.num_lines;
dydx = zeros(ny,nx,N);
dxdx = zeros(nx,N,nx);
dx1  = zeros(size(nw.i1,1),nx);
dx2  = zeros(size(nw.i2,1),nx);
dfdx = zeros(size(nw.i1,1),N,nx);
for i=1:nx
	dxdx(i,:,i) = 1;
	dx1(:,i)    = nw.i1==i;
	dx2(:,i)    = nw.i2==i;
	dfdx(:,:,i) = dx1(:,i).*x(nw.i2,:) + x(nw.i1,:).*dx2(:,i);
	dydx(:,i,:) = [dxdx(:,:,i); nw.A*dxdx(:,:,i) + nw.B*dfdx(:,:,i)];
end

if nargout > 2	
	%Generate derivatives of y wrt theta
	dydt = zeros(ny,N,nt);
	d    = zeros(size(y));
	for i=1:nt
		dA = sparse(nw.dA(:,:,i));
		dB = sparse(nw.dB(:,:,i));
		d(nx+1:end,:) = dA*x + dB*fx;
		dydt(:,:,i) = d;
	end

	%Generate derivatives of dydx wrt theta
	dydxdt = zeros(ny,nx,N,nt);

	for j=1:nx
		sdxdx = sparse(dxdx(:,:,j));
		sdfdx = sparse(dfdx(:,:,j));
		for i=1:nt
			dA = sparse(nw.dA(:,:,i));
			dB = sparse(nw.dB(:,:,i));
			% dxdx            = zeros(nx,N);
			% dxdx(j,:)       = 1;
			% dx1             = nw.i1==j;
			% dx2             = nw.i2==j;
			% dfdx            = dx1.*x(nw.i2,:) + x(nw.i1,:).*dx2;
			tmp1  = dA*sdxdx;
			tmp2  = dB*sdfdx;
			%dydxdt(nx+1:end,j,:,i) = dA*dxdx(:,:,j) + dB*dfdx(:,:,j);
			dydxdt(nx+1:end,j,:,i) = tmp1 + tmp2;
		end
	end
end