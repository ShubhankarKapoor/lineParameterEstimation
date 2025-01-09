function [cost] = mle_cost(Z, G, vars)
num_lines = G.nw.num_lines;
G.theta = vars(1:num_lines*2);

if size(vars,2) == num_lines*2+1
    G.sigma_v = vars(num_lines*2+1);
elseif size(vars,2) == num_lines*2+2
    G.sigma_p = vars(num_lines*2+1);
    G.sigma_v = vars(num_lines*2+2);
end

G.nw = updateNWwithTheta(G.nw,G.theta,G.u0);
nx = size(G.mu,1); % num states
% construct R from sigma_p and sigma_v
G.R = diag([G.sigma_p^2*ones(nx+2, 1); G.sigma_v^2*ones(nx/2+1, 1)]);

yh = meas_fast(G.nw,G.mu); % meas estimates

% LL = 0;
% for k=1:Z.N
%     e  = Z.y(:,k) - yh(:,k);
%     Se = G.R.'\e;
%     LL = LL + 0.5*(Se.'*Se + log(det(G.R))); % log-likelihood
% end

% compact way LL
E = Z.y - yh;
Se = G.R.' \ E;
% Calculate the log-likelihood contributions for each iteration
LL_contributions = sum(Se.^2) + log(det(G.R)) * ones(1, Z.N);
% Sum up the log-likelihood contributions
cost = 0.5*sum(LL_contributions);