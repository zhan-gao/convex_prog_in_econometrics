function [ H ] = MomentMatrix( y, X, Z, b )

% Prepare the H matrix in the numerical implementtation note
% return a m x n matrix

[n, m] = size(Z);

E = repmat(y - X*b, 1, m);
G = (Z .* E)';
sigma = repmat(std(G, 0, 2), 1, n);

H = G ./ sigma;


end

