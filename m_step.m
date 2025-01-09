function G = m_step(Z,G,O)
num_lines = G.nw.num_lines;
vars_init = [G.theta(:); G.sigma_p; G.sigma_v];

result = optGN('costAW',Z,G,O,vars_init(:));
result = abs(result);

G.theta   = result(1:num_lines*2);

if length(vars_init) == num_lines*2+1
    G.sigma_v = result(num_lines*2+1);
	G.sigma_p = result(num_lines*2+1);
elseif length(vars_init) == num_lines*2+2
    G.sigma_p = result(num_lines*2+1);
    G.sigma_v = result(num_lines*2+2);
end