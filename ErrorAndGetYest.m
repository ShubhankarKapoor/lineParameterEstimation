function [G] = ErrorAndGetYest(Z,G,yt0,yt)
    % Estimates for all y including missng meas
    % Estimates for error
    num_lines    = G.nw.num_lines;
    G.nw         = updateNWwithTheta(G.nw,G.theta,G.u0);
    G.yEst       = meas_fast(G.nw,G.mu);
    G.errorTheta = abs((G.true_theta(:) - G.theta(:)) ./ G.true_theta(:)) * 100;
    errV         = abs(G.yEst(num_lines*2+3:end,:) - yt0(num_lines*2+3:end,:)); % =3 is correct
    G.errV       = errV(:);
    % G.errEst  = mean((G.yEst(75:end,:) - yt0(75:end,:)).^2,2);
    % G.errTrue = mean((yt(75:end,:) - yt0(75:end,:)).^2,2);
    % errTrue = mean((Z.y(11:end,:) - yt0(11:end,:)).^2,2);
    % errEM   = mean((yEM(11:end,:) - yt0(11:end,:)).^2,2);
    % figure()
    % plot([errTrue(:) errEst(:)])
    % legend('True', 'Est')
    % pause(0.0001)
end