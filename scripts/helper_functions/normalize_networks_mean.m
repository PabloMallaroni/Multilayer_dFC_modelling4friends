function normalizedMatrix = normalize_networks_mean(matrix, labels, n_iter)
    % Normalizes recruitment and integration values using a permutation approach.
    % Null module allegiance matrices are created by randomly permuting the
    % correspondence between regions of interest (ROIs) and large-scale systems.
    % Functional cartography measures are then calculated for all permuted matrices.
    % To obtain normalized values of recruitment and integration, they are divided
    % by the mean of the corresponding null distribution.
    % This procedure yields null distributions of recruitment and integration
    % coefficients resulting solely from the size of each system.
    %
    % Args:
    %   matrix: N x N matrix
    %   labels: N x 1 vector
    %   n_iter: int
    %

    n_networks = length(unique(labels));

    function nam = calculate_networks_mean(matrix, labels, n_networks)
        nam = zeros(n_networks, n_networks);

        for i = 1:n_networks
            for j = 1:n_networks
                indices_i = find(labels == i);
                indices_j = find(labels == j);
                nam(i, j) = mean(matrix(indices_i, indices_j), 'all');
            end
        end
    end

    nam = calculate_networks_mean(matrix, labels, n_networks);
    nam_null = zeros(n_networks, n_networks);
    labels_null = labels;

    for k = 1:n_iter
        labels_null = labels_null(randperm(length(labels_null)));
        nam_null = nam_null + calculate_networks_mean(matrix, labels_null, n_networks);
    end

    nam_null = nam_null / n_iter;

    normalizedMatrix = nam ./ nam_null;
end
