PLS.optimx <- function(N, TT, y, X, K, lambda, R){

    # PLS estimation by the iterative algorithm
    # try OPTIMX package here
    
    # INPUT Arg:
    #   dimensions N, TT
    #   data y(TN * 1), X(TN * P)
    #   tuning parameter lambda
    #   maximum number of iteration R
    
    # OUTPUT:
    #   b_est: estimated beta (N*p)
    #   a_out: estimated alpha (K*p)
    #   group_est: estimated group identity (N*1)
    
    
    p <- dim(X)[2];
    # Set up the tolerence level
    tol <- 1e-4;
    
    # Use individual regression result as the initial value
    beta0 = numeric();
    for(i in 1:N){
        ind <- ( (i-1)*TT+1 ):(i*TT); 
        yy <- y[ind, ]
        XX <- X[ind, ]
        beta0 <- c( beta0, solve( t(XX) %*% XX ) %*% ( t(XX) %*% yy ) );
    }
    
    init <- c(beta0, c(1,1));
    
    b.out <- array( beta0, c(N,p,K) );
    a.out <- matrix(0,K,p);
    
    b.old <- matrix(1,N,p);
    a.old <- matrix(1,1,p); 
    
    for(r in 1:R){
        
        print(r)
        
        for(k in 1:K){
            
            # N * 1: consider it as gamma
            penalty.out <- pen.generate(b.out, a.out, N, p, K, k);
            
            # opt.out <- optimx(fn = obj, gr = obj.grad, par = init,
            #              penalty = penalty.out);
            
            opts = list(algorithm = "NLOPT_LD_LBFGS", xtol_rel = 1e-5, maxeval = 1000);
            nlopt.out <- nloptr(x0=init, eval_f = obj,eval_grad_f = obj.grad, 
                                opts = opts, penalty = penalty.out);
            
            coef.temp <- nlopt.out$solution;
            
            a.out[k, ] <- coef.temp[(p*N+1):(p*(N+1))];
            b.out[ , ,k] <- matrix(coef.temp[1:(N*p)], N, p, byrow = TRUE);
            
            
        }
        
        # Check the convergence criterion
        a.new <- a.out[K,];
        b.new <- b.out[ , ,K];
        
        if(criterion(a.old,a.new,b.old,b.new,tol)){
            break;
        }
        # Update
        a.old <- a.out[K,];
        b.old <- b.out[ , ,K];
    }
    
    # put b.out to nearest a.out and get the group estimation
    
    a.out.exp <- aperm( array(a.out, c(K,p,N)), c(3,2,1) );
    d.temp <- (b.out - a.out.exp)^2;
    dist <- sqrt( apply(d.temp, c(1,3), sum) );
    group.est <- apply(dist,1,which.min);
    
    # Post-Lasso estimation
    a.out <- post.lasso( group.est, y, X, K, p, N, TT );
    
    b.est <- matrix(999, N, p);
    for(i in 1:N){
        group <- group.est[i];
        b.est[i,] <- a.out[group, ];
    }
    
    result <- list( b.est = b.est, a.out = a.out, 
                    group.est = group.est );
    
    return(result)
}


##########################################################
obj <- function(coeff, penalty = NULL){
    B <- matrix( coeff[1:(p*N)], nrow = p);
    ee <- X %*% B;
    
    ii <- as.logical( c( rep( c( rep(1,TT),rep(0,TT*N) ), N-1), rep(1,TT) ) );
    ob <- ( sum((y - ee[ii])^2) )/( N*TT );
    
    ob <- ob + sum( sqrt( apply( matrix( ( coeff[1:(p*N)] - rep( coeff[(p*N+1):(p*(N+1))], N) )^2, 
                         nrow = N, byrow = TRUE ) , 1, sum) ) * penalty * (lambda/N) );
    return(ob)
    
}
##########################################################

obj.grad <- function(coeff, penalty = NULL){
    B <- matrix( coeff[1:(p*N)], nrow = p);
    ee <- X %*% B;
    ii <- as.logical( c( rep( c( rep(1,TT),rep(0,TT*N) ), N-1), rep(1,TT) ) );
    ob.grad.1 <- as.vector( t( apply( array( matrix( sum((y - ee[ii])^2), TT*N, p ) * X, c(TT, N, p) ),
           c(2,3), sum) ) ) * 2 / (N*TT);
    
    b.minus.a <- coeff[1:(p*N)] - rep( coeff[(p*N+1):(p*(N+1))], N);
    norm2 <- sqrt( apply( matrix( ( b.minus.a )^2, nrow = N, byrow = TRUE ) , 1, sum) );
    ob.grad.1 <- ob.grad.1 + ( b.minus.a / rep(norm2,rep(p,N)) ) * penalty * (lambda/N);
    
    ob.grad.2 <- apply( matrix( (-b.minus.a / rep(norm2,rep(p,N)) ) * penalty * (lambda/N),
            nrow = N, byrow = TRUE ), 2, sum );
    
    return( c(ob.grad.1,ob.grad.2) )
    
}

##########################################################
pen.generate <- function(b, a, N, p, K, kk){
    
    # generate the known part of the penalty term
    # output a N*1 vector
    
    a.out.exp <- aperm( array(a, c(K,p,N)), c(3,2,1) );
    p.temp <- ( b - a.out.exp )^2;
    p.norm <- sqrt( apply(p.temp, c(1,3), sum) );
    
    ind <- setdiff(1:K,kk);
    
    pen <- apply(p.norm[, ind], 1, prod);
    return(pen)
    
}

##########################################################
criterion <- function( a.old, a.new, b.old, b.new, tol ){
    
    d <- FALSE;
    
    a.cri <- sum( abs( a.old - a.new ) ) / (sum( abs(a.old) ) + 0.0001); 
    b.cri <- mean( abs( b.old - b.new ) ) / (mean( abs(b.old) ) + 0.0001);
    
    if(a.cri < tol & b.cri < tol){
        d <- TRUE;
    }
    
    return(d)
    
}

##########################################################
post.lasso <- function( group.est, y, X,  K, p, N, TT ){
    
    a.out <- matrix(0,K,p);
    
    for( k in 1:K ){
        Ind <- 1:N; group.ind <- Ind[group.est == k];
        data.ind <- as.numeric( sapply(group.ind, extend.ind) );
        yy <- y[data.ind,];
        XX <- X[data.ind,];
        
        a.out[k, ] <- solve( t(XX) %*% XX ) %*% ( t(XX) %*% yy );
    }
    return(a.out)
}

extend.ind <- function(i){
    return( ((i-1)*TT+1):(i*TT) );
}
