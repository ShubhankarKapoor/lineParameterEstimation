function [Z, M] = extract_data_same_node(Z,M)
% update mu, P and y
% chose all measurements at the same node

num_missing    = floor(0.2*M.nw.num_lines); % num_lines used as a proxy for nodes with loads
% missing_index  = sort(randperm(M.nw.num_lines, num_missing));
% disp(missing_index);    
missing_index  = [3]; % node number missing - chosing manually here
% missing_index  = [2, 5, 8, 10, 20, 27, 30, 35]; % node number missing - chosing manually here
% missing_index  = [5, 8, 10, 20, 27, 30, 35]; % node number missing - chosing manually here
active_p_idx   = missing_index;
reactive_p_idx = missing_index+M.nw.nodes;
volt_indx      = missing_index+2+2*M.nw.nodes+1; % adding 1 because 37 meas as opposed to 36 & 2 for lineflows
missing_index  = [active_p_idx reactive_p_idx volt_indx];
avail_index    = ones(size(Z.y));
avail_index(missing_index,:) = 0;

% update mu
extracted_subarray = Z.y .* avail_index;
num_states         = size(M.mu,1);
mu                 = extracted_subarray(1:num_states,:);
M.mu               = mu; % prior mean

% indexing for P
index_P            = avail_index(1:num_states,:);
% Modify the diagonal elements using linear indexing
diagonal_indices = 1:(num_states + 1):(num_states * num_states);

% Preallocate a cell array to store the measurements for each column
num_cols = size(Z.y, 2);
extracted_cells = cell(1, num_cols);

for i=1:Z.N
    % State Update
    % modify diagoal elements of matrix P
    diagonal_indices_to_modify = index_P(:,i) == 0;
    % update prior variance on states
    slice_P = M.P(:,:,i);
    slice_P(diagonal_indices(diagonal_indices_to_modify)) = 1e-1; % higher value of prior
    M.P(:,:,i) = slice_P;
    
	mean_mu = mean(M.mu(index_P(:,i)==1,i));
	M.mu(diagonal_indices_to_modify,i) = mean_mu;

    % Measurement update
    % Extract measurements for the current column where avail_index is 1
    extracted_column = Z.y(avail_index(:, i) == 1, i);
    % Store the extracted column in the cell array
    extracted_cells{i} = extracted_column;
end
Z.y           = extracted_cells;
Z.avail_index = avail_index; % put the binary matrix as an arguement
end
