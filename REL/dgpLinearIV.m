function [y, X, Z] = dgpLinearIV(b, n, m)


% n: sample size
% m: number of IVs

%% instruments and disturbances

Z = random('Normal', 0, 1, n, m );
E = random('Normal', 0, 1, n, 3);
Sigma = [0.25, 0.15, 0.15; 0.15, 0.25, 0; 0.15, 0, 0.25];
E = E * Sigma;
e = E(:,1);
v1 = E(:,2);
v2 = E(:,3);

%% structural equation


x1 = Z(:, 1:2) * ones(2, 1) + v1;
x2 = Z(:, 3:4) * ones(2, 1) + v2;

X = [x1, x2];

y = X * b + e;

end
