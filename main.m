clear
% close all

rng(1) %1

% load network data
load saved_data37.mat;
nodes = 36;
% for different 5 node network configurations
% load saved_data5_straight.mat;
% load saved_data5.mat;
% load saved_data5_ac.mat
% nodes = 4;

%setup the network structure
nw = pre;
nw.num_lines = num_lines;
nw.indices = 1:num_lines;
nw.indices_v = 1:num_lines+1;
nw.V_node_idx = 1:num_lines+1;
nw.non_zib_index_array = nw.indices+1;
nw.comb_idx1 = nw.comb_idx1 + 1;
nw.comb_idx2 = nw.comb_idx2 + 1;
nw.elems_comb = nw.elems_comb+1;
nw.downstream_matrix_p = nw.downstream_matrix(:,nw.non_zib_index_array);
nw.nodes = nodes;
N = 600;

%set initial theta
theta0 = 0.018 + 0*(true_theta + 0.5*true_theta.*(rand(size(true_theta))-0.5));
% theta0 = true_theta;
nw    = updateNWwithTheta(nw,true_theta,1.0);
yt0    = load('y_true.mat').data; 

% load real dataset, consumer end load profile
% yt0 = load('ynextgen2.mat').y_real;
% if you need real world data used here please email me or request access to Nextgen
% random dataset
% randomIndices = randperm(size(yt0, 2), N);
% yt0 = yt0(:, randomIndices);

x_true = yt0(1:2*nodes,1:N);

%Monte carlo to investigate sensitivity to noise realisation
num_runs = 1;
ths      = [];
for ii=1:num_runs

    %Generate measured data 
    yt = abs(yt0 + 1e-3*randn(size(yt0)));

    % Calculate the percentage error for each element
    percentage_error = (abs(yt - yt0) ./ yt0) * 100;
    mean_percentage_error = mean(mean(percentage_error));
    % disp(mean_percentage_error);
    
    %data stuff
    Z.y = yt;
    Z.N = size(yt,2);
    Z.xt = x_true;
    
    %Model stuff
    n       = 2*num_lines;
    M.theta = theta0;
    M.u0    = 1.0;
    M.nw    = nw;
    M.mu    = Z.y(1:size(x_true,1),:);
    M.true_theta = true_theta; % to plot

    %Now setup measurement covariance (R) and initial state covariance (P)
    for i=1:Z.N
        % M.R(:,:,i) = 1*eye(size(yt,1)); %measurement noise
        M.P(:,:,i) = 1e-1*eye(n,n); %state noise/ prior
    end

    % noise params
    M.sigma_p = 1e-1;
    M.sigma_v = 1e-1; %0.001
    
    % MLE estimates with full noisy states and meas for comparison
    % [G_ml] = mle(Z, M);
    % break
    [Z, M] = extract_random_data(Z, M, 0.8); % extracts random data
    % [Z, M] = extract_data_same_node(Z,M);  % extracts data from same node
    
    %Options    
    O.numit = 20;

    M = m_step(Z,M,[]);

    % %Now run EM
    G = em(Z,M,O);

    figure()
    % plot([true_theta(:) theta0(:) G.theta(:)])
    plot([true_theta(:) G.theta(:)])
    legend('True', 'Est.')
	pause(0.0001)

    % figure()
    % plot([true_theta(:) G_ml.theta(:) G.theta(:)])
    % plot([true_theta(:) G.theta(:)])
    % legend('True', 'ML', 'Est.')
    % pause(0.0001)

    ths = [ths, G.theta(:)];
       
    % plot likelihood
    figure()
    plot(G.LLHistory, '-o')

end

% ERRORS %
% mse_ml = mse(true_theta(:) , G_ml.theta(:)); % MSE for ML
mse_em = mse(true_theta(:) , G.theta(:)); % MSE for EM

% Measaurement estimates and voltage errors 
[G] = ErrorAndGetYest(Z,G,yt0,yt);
% [G_ml] = ErrorAndGetYest(Z,G_ml,yt0,yt);
