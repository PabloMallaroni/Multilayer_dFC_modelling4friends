function T_mean = single_window_allegiance_mean(M)
    % Calculates the mean of single-window allegiance matrices.
    %
    % Parameters
    % ------------
    % M: n_opt x n_nod module assignment matrix
    %
    % Returns
    % ------------
    % T_mean: n_nod x n_nod mean allegiance matrix

    [n_opt, n_nod] = size(M);
    T = zeros(n_opt, n_nod, n_nod);

    for i = 1:n_opt
        T(i, :, :) = single_window_allegiance(M(i, :));
    end

    T_mean = squeeze(mean(T, 1));
end