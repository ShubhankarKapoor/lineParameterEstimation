function G = e_step(Z,G,O)

G.LL     = 0;
G.nw     = updateNWwithTheta(G.nw,G.theta,G.u0);
[y,dydx] = meas_fast(G.nw,G.mu,1);
nx       = size(G.mu,1);
G.R      = diag([G.sigma_p*ones(nx+2, 1); G.sigma_v*ones(nx/2+1, 1)]);

for k=1:Z.N
    ny       = size(Z.y{k},1);
    % extract the indices that are required
    indices     = Z.avail_index(:, k) == 1;
    e           = Z.y{k} - y(indices, k);
    C           = squeeze(dydx(indices,:,k));
    P           = G.P(:,:,k).';
    % index R
    R           = G.R(indices,indices);
    R1          = triu(qr([R zeros(ny,nx);P*C.' P])); % R for each time step
    % G.Ri(:,:,k) = R1(1:ny,1:ny);
    Se          = R1(1:ny,1:ny).'\e;
    G.mu(:,k)   = G.mu(:,k) + R1(1:ny,ny+1:end).'*Se; % state mean
    G.P(:,:,k)  = R1(ny+1:end,ny+1:end).'; % state vars covariance, be mindful this is P^(1/2) and not P
    G.LL        = G.LL + Se.'*Se + sum(log(diag(R).*diag(R))); % log-likelihood
end
G.dydx = dydx;

% Gmu = abs(G.mu);

%keyboard