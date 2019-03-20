function [Result] = result_sum(correct_ratio, se_record, time_cvx, err)

[Rep, m] = size(correct_ratio);

err_num = sum(err);
time_cvx_sum = sum(time_cvx) ./ (Rep - err_num) .* Rep;


cr_result = ones(1,m);
rmse_result = ones(1,m);
for i = 1:m
    cr_result(i) = mean( correct_ratio(err(:, i) == 0,i) );
    rmse_result(i) = sqrt( mean( se_record(err(:, i) == 0,i) ) );
end

Result = [cr_result; rmse_result; time_cvx_sum];

end

