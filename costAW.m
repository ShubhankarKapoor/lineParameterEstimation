function [cost,grad,J] = costAW(Z,G,O,vars)
num_lines = G.nw.num_lines;
G.theta = vars(1:num_lines*2);

is_grad = nargout > 1;

if length(vars) == num_lines*2+1
    G.sigma_v = vars(num_lines*2+1);
	G.sigma_p = vars(num_lines*2+1);
elseif length(vars) == num_lines*2+2
    G.sigma_p = vars(num_lines*2+1);
    G.sigma_v = vars(num_lines*2+2);
end

if is_grad
	[G.nw, G.nw.dA, G.nw.dB] = updateNWwithTheta(G.nw,G.theta,G.u0);
else
	G.nw = updateNWwithTheta(G.nw,G.theta,G.u0);
end

nx = size(G.mu,1); % num states
% construct R from sigma_p and sigma_v
G.R  = [G.sigma_p^2*ones(nx+2, 1); G.sigma_v^2*ones(nx/2+1, 1)];
G.sR = [G.sigma_p*ones(nx+2, 1); G.sigma_v*ones(nx/2+1, 1)];

if ~is_grad
	[yh,dydx] = meas_fast_AW(G.nw,G.mu);
else
	[yh,dydx,dydt,dydxdt] = meas_fast_AW(G.nw,G.mu);
end

cost = 0;
% cost1 = 0;
if is_grad
	grad = zeros(length(vars),1);
	Ny = 0;
	for t=1:Z.N
		ny = length(Z.y{t});
		Ny = Ny + ny + ny*nx;
	end
	% J    = zeros(Ny,length(vars));
	% e    = zeros(Ny,1);
	J = [];
end

istart = 1;
istop  = 0;
for t=1:Z.N
	ny = length(Z.y{t});
	istop = istop + ny + ny*nx;
    % inedxing on C
    indices    = Z.avail_index(:, t) == 1;
    err        = Z.y{t} - yh(indices, t);
    C          = squeeze(dydx(indices,:,t));

    % index R
    Sig        = G.sR(indices);
    diagSig    = diag(Sig);
    
    % higher order term in the cost function
    CP         = C*G.P(:,:,t);
    % CPt        = CP.';
    Sigerr     = diagSig\err;
    SigCP      = diagSig\CP;
	err1       = [Sigerr;SigCP(:)];

    % Gammat     = err*err'+ (C*G.P(:,:,t))*(G.P(:,:,t)'*C');
    % SigGam     = diagSig\(diagSig\Gammat);
    cost       = cost + 0.5*sum(log(Sig.*Sig)) + 0.5*(err1.'*err1);%0.5*trace(SigGam);
	% cost       = cost + 0.5*sum(log(Sig.*Sig)) + 0.5*trace(SigGam);


	%update the gradient
	if is_grad
		% e(istart:istop) = err1;
		Jt = zeros(ny+ny*nx,length(vars));
		for i=1:length(vars)
			if i <= num_lines*2
				derr    = -dydt(indices,t,i); % dG
				Sigderr = diagSig\derr;
				dC      = dydxdt(indices,:,t,i);
                dCP     = dC*G.P(:,:,t);
				SigdCP  = diagSig\dCP;
				derr1   = [Sigderr;SigdCP(:)];
				% grad(i) = grad(i) + 0.5*trace(diagSig\(diagSig\(derr*err.'))) + 0.5*trace(diagSig\(diagSig\(dCP*CP.')));
				% grad(i) = grad(i) + 0.5*trace(diagSig\(diagSig\(err*derr.'))) + 0.5*trace(diagSig\(diagSig\(CP*dCP.')));
                grad(i) = grad(i) + derr1.'*err1;%Sigerr.'*derr + trace(SigCP*dCP.');
				% if abs(derr1.'*err1 - trace(diagSig\(diagSig\(derr*err.' + CP*dCP.')))) > 1e-8
				% 	keyboard;
				% end
				% J(istart:istop,i) = derr1;
				Jt(:,i) = derr1;
			elseif i==num_lines*2+1	&& length(vars)==num_lines*2+2
				dSig    = [ones(nx+2, 1); zeros(nx/2+1, 1)];
				dSig    = dSig(indices);
				dSigerr = -diagSig\(diag(dSig)*(Sigerr));
				dSigCP  = -diagSig\(diag(dSig)*(SigCP));
				dSigerr1 = [dSigerr;dSigCP(:)];
                % grad(i) = grad(i) + sum(dSig./Sig./Sig) - 0.5*trace(diagSig\(diag(dSig)*SigGam));
				grad(i) = grad(i) + sum(dSig./Sig) + dSigerr1.'*err1;
				% J(istart:istop,i) = dSigerr1;
				Jt(:,i) = dSigerr1;
			elseif i==num_lines*2+2 && length(vars)==num_lines*2+2
				dSig = [zeros(nx+2, 1); ones(nx/2+1, 1)];
				dSig = dSig(indices);
				dSigerr = -diagSig\(diag(dSig)*(Sigerr));
				dSigCP = -diagSig\(diag(dSig)*(SigCP));
				dSigerr1 = [dSigerr;dSigCP(:)];
                % grad(i) = grad(i) + sum(dSig./Sig./Sig) - 0.5*trace(diagSig\(diag(dSig)*SigGam));
				grad(i) = grad(i) + sum(dSig./Sig) + dSigerr1.'*err1;
				% J(istart:istop,i) = dSigerr1;
				Jt(:,i) = dSigerr1;
			elseif i==num_lines*2+1 && length(vars)==num_lines*2+1
				dSig = [ones(nx+2, 1); ones(nx/2+1, 1)];
				dSig = dSig(indices);
				dSigerr = -diagSig\(diag(dSig)*(Sigerr));
				dSigCP = -diagSig\(diag(dSig)*(SigCP));
				dSigerr1 = [dSigerr;dSigCP(:)];
                % grad(i) = grad(i) + sum(dSig./Sig./Sig) - 0.5*trace(diagSig\(diag(dSig)*SigGam));
				grad(i) = grad(i) + sum(dSig./Sig) + dSigerr1.'*err1;
				% J(istart:istop,i) = dSigerr1;
				Jt(:,i) = dSigerr1;
			end
		end
		J = triu(qr([J;Jt]));
		J = J(1:length(vars),:); % dont get this
		istart = istart + ny + ny*nx;
	end
end

% keyboard
