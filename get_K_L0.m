function K = get_K_L0(L0)
    persistent L0_vals K_vals

    if isempty(L0_vals)
        data = load('K_lookup_data.mat', 'L0_vals', 'K_vals');
        L0_vals = data.L0_vals;
        K_vals = data.K_vals;
    end

    % Bound the input (avoid extrapolation errors)
    L0 = max(min(L0, max(L0_vals)), min(L0_vals));

    % Find the closest L0 index
    for j = 1:2
        for k = 1:6
            K(j,k) = interp1(L0_vals, reshape(K_vals(j,k,:),1,[]), L0, 'linear', 'extrap');
        end
    end

    % Extract the corresponding K
    K = K_vals(:,:,idx);
end