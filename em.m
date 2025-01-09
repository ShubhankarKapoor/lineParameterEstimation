function G = em(Z,M,O)

G = M;
G.thetaHistory(1,:) = G.theta;
for i=1:O.numit
    %E step
    G = e_step(Z,G,O);
    
    %M step
    G = m_step(Z,G,O);
    G.thetaHistory(i+1,:) = G.theta;
    
    %Update progress
    fprintf('It# %4d, LL: %13.4e\n',i,G.LL);
    G.LLHistory(i) = G.LL;
    
    % plot
    figure(50031)
    plot([G.true_theta(:) G.theta(:)])
    legend('True','Est.')
    pause(1)

	if i>2 
		if G.LLHistory(i) >= G.LLHistory(i-1)
			break;
		end
	end
end
