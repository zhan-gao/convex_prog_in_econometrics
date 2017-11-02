group.coerce <- function(group.est, a.out, group0, a0, N, N.frac, K, p){
    
    a.est <- matrix(0,K,p);
    group <- rep(0,N);
    
    for(k in 1:K){
        dif <- apply( (matrix( a.out[k,], K, p, byrow = TRUE) - a0)^2, 1, sum);
        gg <- which.min(dif);
        
        a.est[gg, ] <- a.out[k, ];
        group[group.est == k] <- gg;
    }
    
    if( sum(group == 0) ){
        noquote("Error: some inidividual are not classfied!")
    }
    
    # correct ratio
    ratio <- sum( group == group0 ) / N;
    # rmse
    weighted.se <- sum( ( a.est[,1] - a0[,1] )^2 * N.frac );
    
    return( list(ratio = ratio, se = weighted.se,
                 a = a.est, group = group) );
}



###########################################################

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