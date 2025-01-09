function [xOpt,R] = optGN(fcn,Z,G,O,x)

x = x(:);
n = length(x);

% ub = G.ub;

% mu = 1;

lam = 1e-3;

for i=1:1000
    [f,g,J] = feval(fcn,Z,G,O,x);

	% f = f/mu - sum(log(x));% + log(ub-x));
	% g = g/mu - 1./x;% + 1./(ub-x);
	
	R         = triu(qr([J;lam*eye(n)]));%/sqrt(mu);diag(1./x)]));%;diag(1./(ub-x))]));
    % R         = R(1:n,1:n);
	% R = J;

    % p is the search direction?
    p = -(R.'\g);
    p = R\p;

	% keyboard

	% sl  = -p./x;
	% sl  = max(sl(sl>0));
	% if isempty(sl)
	% 	sl = 0.0;
	% end
	% su  = p./(ub-x);
	% su  = max(su(su>0));
	% if isempty(su)
	% 	su = 0.0;
	% end
    % alp = min(1.0,0.9999/max(sl,su));
	% alp = min(2.0,0.9999/sl);
	% ff = @(a) feval(fcn,Z,G,O,x+a*p);%/mu - sum(log(x+a*p));
	% alp = fminbnd(ff,0,2.0);
	% x = x + alp*p;

	alp = 1.0;
    for k=1:52
        [fa] = feval(fcn,Z,G,O,x+alp*p);
		% fa = fa/mu - sum(log(x+alp*p));%+log(ub-x-alp*p));
		if ~isfinite(fa)
			fa = inf;
		end
        if fa < f
            x = x+alp*p;
            break;
        end
        alp = alp/2;
    end
       
	if alp == 1.0
		lam = max(1e-8,lam/2);
	elseif alp < 0.25
		lam = 2*lam;
	end
	
    fprintf('Iteration Number = %5i, Cost = %10.2e,  Newton Decrement = %10.2e,  Alpha = %10.2e, Lambda = %10.2e\n',i,f,p'*g,alp,lam);
    
    if abs(p.'*g) < 1e-1
        break;
	end

	% if abs(p.'*g) < 1e-1 && mu < 1e-6
	% 	break;
	% end
end
xOpt = x;
