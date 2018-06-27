function [ratio, weighted_se, a_est, group] = result_gen(group_est, a_out, group0, a0, N, N_frac, K, p)

a_est = zeros(K, p);
group = zeros(1, N);

for k = 1:K
    
	dif = sum( ( repmat(a_out(k,:), K, 1) - a0 ).^2, 2);
    gg = find(dif == min(dif));
    
    a_est(gg,:) = a_out(k, :);
    group(group_est(:, k) ~= 0) = gg;
    
end

if sum(group == 0)
    warning('Some individuals are not classfied!');
end

% Compute the correct ratio
ratio = sum(group == group0) / N;
% weighted RMSE
weighted_se = sum( (a_est(:, 1) - a0(:, 1)).^2 .* N_frac' );

end
