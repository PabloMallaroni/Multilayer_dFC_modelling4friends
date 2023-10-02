function T = all_window_allegiance_mean(M)
    % Calculates the mean allegiance matrix T for per time window time windows from a 3D matrix M
    % Across all optimisation runs
    %
    % Parameters:
    %   M: 3D matrix of size (n_opt x n_nod x n_win), where
    %      n_win is the number of time windows,
    %      n_nod is the number of nodes.
    %      n_opt is the number of optimisations. 
    %
    % Returns:
    %   T: 3D matrix of size (n_win x n_nod x n_nod), containing the mean allegiance matrices.

    [~,n_nod, n_win] = size(M);
    T = zeros(n_win, n_nod, n_nod);

    for i = 1:n_win
        T(i, :, :) = single_window_allegiance_mean(squeeze(M(:,:,i)));
    end
end
