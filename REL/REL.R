
innerloop <- function(b, y = NULL, X = NULL, Z = NULL, tau = NULL){
    
    # Solve the inner loop optimization of Relaxed Empirical Likelihood by RMOSEK
    # max sum_i(log(p_i)) s.t. sum_i(p_i) = 1; |sum_i p_i*g_ij(b)|/ sigma_j <= tau
    
    # Give the data y, X, Z, this is a numerical function of b
    # tau is the tuning parameter
    
    # Packages: Rmosek, Matrix
    
    n <- nrow(Z); m <- ncol(Z);
    H <- MomentMatrix(y, X, Z, b);
    
    # Initialize the mosek problem
    Prob <- list(sense = "max");
    
    # Linear constraints
    Prob$A <- Matrix( rbind( rep(1,n), H ), sparse = TRUE );
    Prob$bc <- rbind( c(1, rep(-tau, m)), c(1, rep(tau, m)) );
    Prob$bx <- rbind( rep(0,n), rep(1,n) );
    
    # Linear coefficients of the objective
    Prob$c <- rep(0, n);
    # Logarithm part of the objective
    NUMOPRO <- n;
    opro <- matrix( list(), nrow = 5, ncol = NUMOPRO );
    rownames(opro) <- c("type", "j" , "f", "g", "h");
    for( i in 1:n ){
        opro[ , i] <- list("LOG", i, 1.0, 1.0, 0);
    }
    Prob$scopt = list( opro = opro )
    
    # Invoke Mosek
    mosek.out <- mosek(Prob, opts = list(verbose = 0, soldetail = 1));
    
    if(mosek.out$sol$itr$solsta == "OPTIMAL"){
        # Since the default of NLOPTR is to do minimization, need to set it as negative
        return( -mosek.out$sol$itr$pobjval );
    }else{
        print( "WARNING: Inner loop not optimized" );
        return( Inf );
    }
    
}



innerloop.cvxr <- function(b, y = NULL, X = NULL, Z = NULL, tau = NULL){
    
    # Solve the inner loop optimization of Relaxed Empirical Likelihood by cvxr
    # max sum_i(log(p_i)) s.t. sum_i(p_i) = 1; |sum_i p_i*g_ij(b)|/ sigma_j <= tau
    
    # Give the data y, X, Z, this is a numerical function of b
    # tau is the tuning parameter
    
    # Packages: cvxr
    
    n <- nrow(Z); m <- ncol(Z);
    H <- MomentMatrix(y, X, Z, b);
    
    p = Variable(n)
    
    constr = list( sum(p)==1, 
                   p>=0,
                   p<=1,
                   H%*%p >= -tau,
                   H%*%p <= tau )
    
    obj = sum( log(p) )
    obj = Maximize(obj)
    
    Prob = Problem(obj, constr)
    cvxr.out = solve(Prob)
    
    
    if(cvxr.out$status == "optimal"){
        # Since the default of NLOPTR is to do minimization, need to set it as negative
        return( -cvxr.out$value );
    }else{
        print( "WARNING: Inner loop not optimized" );
        return( Inf );
    }
    
}


REL <- function(y, X, Z, tau){
    
    # Package: nloptr
    
    # invoke nloptr
    opts = list(algorithm = "NLOPT_LN_NELDERMEAD", xtol_rel = 1e-5, maxeval = 5000);
    nlopt.out <- nloptr(x0 = c(0,0), eval_f = innerloop, opts = opts, y = y, X = X, Z = Z, tau = tau);
    
    if(nlopt.out$status != 0){
        print( "WARNING: Outer loop not optimized" )
    }
    
    return( nlopt.out$solution )

}


REL.cvxr = function(y, X, Z, tau){
    
    opts = list(algorithm = "NLOPT_LN_NELDERMEAD", xtol_rel = 1e-5, maxeval = 5000);
    nlopt.out <- nloptr(x0 = c(0,0), eval_f = innerloop.cvxr, opts = opts, y = y, X = X, Z = Z, tau = tau);
    
    if(nlopt.out$status != 0){
        print( "WARNING: Outer loop not optimized" )
    }
    
    return( nlopt.out$solution )
}

MomentMatrix <- function(y, X, Z, b){
    
    # Prepare the H matrix in the numerical implementation notes
    # return a m-by-n matrix
    
    n <- nrow(Z); m <- ncol(Z);
    
    E <- matrix( y - X %*% b, n, m ); G <- t( Z*E );
    sigma <- matrix( apply(G, 1, sd), m, n);
    H <- G / sigma; return(H);
}



###############################################################

# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
    Envir = as.environment(-1)
    
    if (length(r) > length(l))
        warning("RHS has more args than LHS. Only first", length(l), "used.")
    
    if (length(l) > length(r))  {
        warning("LHS has more args than RHS. RHS will be repeated.")
        r <- extendToMatch(r, l)
    }
    
    for (II in 1:length(l)) {
        do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
    }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
    s <- length(source)
    d <- length(destin)
    
    # Assume that destin is a length when it is a single number and source is not
    if(d==1 && s>1 && !is.null(as.numeric(destin)))
        d <- destin
    
    dif <- d - s
    if (dif > 0) {
        source <- rep(source, ceiling(d/s))[1:d]
    }
    return (source)
}

# Grouping the left hand side
g = function(...) {
    List = as.list(substitute(list(...)))[-1L]
    class(List) = 'lbunch'
    return(List)
}