function [Z, M] = extract_random_data(Z,M,percentage)
rng(10)
num_states         = size(M.mu,1);
% update mu, P and y
% case of chosing measurements randomly
% percentage     = 0.8; % percentage of data available
disp(percentage);
total_elements = numel(Z.y);
num_avail      = floor(percentage * total_elements); % percentage of measurements available
avail_index    = sort(randperm(total_elements, num_avail));
[row_indices, col_indices] = ind2sub(size(Z.y), avail_index);
% Create the matrix where measurements are available
avail_index = zeros(size(Z.y));

% Set the elements at row_indices and col_indices to 1
avail_index(sub2ind(size(Z.y), row_indices, col_indices)) = 1; % binary array
% run a loop to make sure there are pline and qline measurements
% while any(avail_index(num_states+1,:) == 0) && any(avail_index(num_states+2,:) == 0)
%     num_changed = sum(avail_index(num_states+1,:)==0) + sum(avail_index(num_states+2,:)==0);
%     disp(num_changed)
%     % convert pline and qline to 1
%     avail_index(num_states+1,:)   = 1;
%     avail_index(num_states+2,:)  = 1;
%     % convert 1s to 0s to make sure still it is 80% availability
%     % the same number that you changed to 1 for rows 9 and 10
% 
%     % Find the indices of the non-zero elements in the array
%     [non_zero_rows, non_zero_cols] = find(avail_index ~= 0);
%     non_zero_indices = sub2ind(size(avail_index), non_zero_rows, non_zero_cols);
% 
%     % Shuffle the indices randomly so that we can select m of them
%     shuffled_indices = non_zero_indices(randperm(length(non_zero_indices)));
% 
%     % Set the first m non-zero elements to 0
%     avail_index(shuffled_indices(1:num_changed)) = 0;
% end

% update mu
extracted_subarray = Z.y .* avail_index;
mu                 = extracted_subarray(1:num_states,:);
M.mu               = mu;
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
    slice_P(diagonal_indices(diagonal_indices_to_modify)) =2e-1; % 5e-1; % highrt value of prior
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
