function[L, weight] = REL_inner(b, y, X, Z, tau)

[n, m] = size(Z);
H = MomentMatrix(y, X, Z, b);

cvx_solver mosek
cvx_begin
    variable p(n,1);
    maximize( sum(log(p)) );
    subject to
        sum(p) == 1;
        p <= 1;
        p >= 0;
        H*p >= -tau*ones(m,1);
        H*p <= tau*ones(m,1);
cvx_end

L = sum(log(p));
weight = p;

end

