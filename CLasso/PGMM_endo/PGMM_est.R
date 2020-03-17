PGMM.mosek = function(N, t, yy, XX, Z, K, lambda, R, tol = 1e-4){
    
    # library("Rmosek")
    # library("SparseM")
    # library("Matrix")
    
    # PGMM estimation by the iterative algorithm
    
    # INPUT Arg:
    #   dimensions N, TT
    #   data y(TN * 1), X(TN * P), Z(TN*d)
    #   tuning parameter lambda
    #   maximum number of iteration R
    
    # OUTPUT:
    #   b_est: estimated beta (N*p)
    #   a_out: estimated alpha (K*p)
    #   group_est: estimated group identity (N*1)
    
    p = dim(X)[2]
    d = dim(Z)[2]
    
    y = matrix(0,d*N,1)
    X = matrix(0,d*N,p)
    
    for(i in 1:N){
        
        ind = ( (i-1)*t+1 ):(i*t)
        ind.d = ( (i-1)*d+1 ):(i*d)
        
        X.i = XX[ind, ]
        
        y[ind.d] = t(Z[ind, ]) %*% yy[ind]
        X[ind.d, ] = t(Z[ind,]) %*% XX[ind, ]
        
    }
    
    mosek.out = PLS.mosek(N, d, y, X, K, lambda, R, tol)
    
    return(mosek.out)
    
    
}