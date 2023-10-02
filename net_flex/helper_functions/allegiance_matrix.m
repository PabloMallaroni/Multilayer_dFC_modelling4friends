function P = allegiance_matrix(M)
    % Calculates n x n allegiance matrix P from o x n x t matrix M,
    % where n represents nodes, and t represents time windows.
    % Each value on allegiance matrix P represents the probability that
    % node i and node j have been assigned to the same community.
    %
    % Parameters
    % ------------
    % M: n x t module assignment matrix
    %
    % Returns
    % ------------
    % P: n x n allegiance matrix
    
    n_nodes = size(M, 1);
    n_slices = size(M, 2);
    T = zeros(n_nodes, n_nodes);
    
    for i = 1:n_nodes
        for j = 1:i-1
            if i == j
                continue;
            else
                t = sum(M(i, :) == M(j, :));
                T(i, j) = t;
            end
        end
    end
    
    P = (T + T') / n_slices;
    P(1:n_nodes+1:end) = 1; % Fill diagonal with 1
end