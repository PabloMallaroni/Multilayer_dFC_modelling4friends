function T = single_window_allegiance(t)
    % Calculates the single-window allegiance matrix.
    %
    % Parameters
    % ------------
    % t: n_nod x 1 module assignment vector
    %
    % Returns
    % ------------
    % T: n_nod x n_nod allegiance matrix

    n_nod = length(t);
    T = zeros(n_nod, n_nod);

    for i = 1:n_nod
        for j = 1:n_nod
            if t(i) == t(j)
                T(i, j) = 1;
            else
                continue;
            end
        end
    end
end
