
DGP.static <- function(N, TT, p, N.cut, a0){
    
    # Static panel panel: DGP
    
    # INPUT Arg:
    # N: num of obs
    # TT: time horizon
    # N_cut: True group structure
    # a0: True parameter
    
    # Outputs:
    # X: TN * p, demeaned
    # y: TN * 1, demeaned
    # y_{i=1,t=1}, ..., y_{i=1, t=T}, y_{i=2, t=1}, ...
    
    
    e <- matrix( rnorm(N*TT*p), nrow = N*TT);
    u <- rnorm(N);
    uu <- rep( u, times = rep(TT,N) );
    XX <- 0.2 * matrix( rep(uu, times = p), nrow = N*TT) + e;
    
    epsl <- matrix( rnorm(N*TT), nrow = N*TT);
    y.temp <- XX %*% a0;
    yy <- matrix(0, N*TT, 1);
    for(kk in 1:K){
        cut <- c(0, N.cut);
        ind <- ( cut[kk]*N*TT+1) : (cut[kk+1]*N*TT);
        yy[ind, ] <- y.temp[ind, kk];
    }
    yy <- yy + epsl;
    
    y <- demean(yy, N, TT);
    X <- demean(XX, N, TT);
    
    result <- list(y = y, X = X);
    return(result)
    
}

demean <- function(yy, N, TT){
    
    # Output is the demeaned data with the same dimension as input
    # NT * 1 or NT * p
    
    if(dim(yy)[1] != N*TT) print("Error! Dimension of
                                 inputs in demean is wrong!")
    
    p <- dim(yy)[2];
    
    if( p == 1){
        y.temp <- matrix(yy, nrow = TT);
        m <- apply(y.temp, 2, mean);
        y.temp <- y.temp - matrix( rep(m, times = TT), nrow = TT, 
                                   ncol = N, byrow = TRUE);
        y <- matrix(y.temp, nrow = N*TT)
        return(y)
    }
    else{
        y <- matrix(0, N*TT, p);
        for(j in 1:p){
            y.temp <- matrix( yy[,j], nrow = TT);
            m <- apply(y.temp, 2, mean);
            y.temp <- y.temp - matrix( rep(m, times = TT), nrow = TT,
                                       ncol = N, byrow = TRUE);
            y[,j] <- matrix(y.temp, nrow = N*TT);
        }
        return(y)
    }
}























