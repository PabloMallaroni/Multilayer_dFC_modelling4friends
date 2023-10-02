function PO = allegiance_matrix_opti(M)
    % Calculates n x n allegiance matrix P from o x n x t matrix M, where o represents module community detection 
    % optimizations, n represents nodes, and t represents time windows. Each value on allegiance matrix PO represents the 
    % probability that node i and node j have been assigned to the same 
    % community (Bassett et al., 2014; Mattar et al., 2015).
    
    % Parameters:
    %   M: o x n x t module assignment matrix

    % Returns:
    %   PO: n x n allegiance matrix P (mean across all optimizations)

    [n_opt, n_nodes, ~] = size(M);
    P = zeros(n_opt, n_nodes, n_nodes);
    
    for i = 1:n_opt
        P(i, :, :) = allegiance_matrix(squeeze(M(i, :, :)));
    end
    
    PO = squeeze(mean(P, 1));
end