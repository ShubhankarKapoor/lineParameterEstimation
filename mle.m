function [G_ml] = mle(Z, G_ml)
    num_lines = G_ml.nw.num_lines;
    vars_init = [G_ml.theta(:); G_ml.sigma_p; G_ml.sigma_v];
    lb = zeros(1,size(G_ml.theta,2)+2);
    ub = [ones(1, num_lines*2)*100, 1, 1];
    options      = optimoptions('fmincon', 'Display', 'iter');
    result       = fmincon(@(vars) mle_cost(Z,G_ml,vars),vars_init,[],[],[],[],lb,ub,[],options);
    G_ml.theta   = result(1:num_lines*2);
    G_ml.sigma_p = result(num_lines*2+1);
    G_ml.sigma_v = result(num_lines*2+2);
end