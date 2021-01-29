lasso.mosek <- function(x, y, lambda, rtol = 1e-6, verb = 0, intercept = TRUE, scalex = TRUE){

    n <- nrow(x)
    p <- ncol(x)

    P = list( sense = "min" )
    P$bc <- rbind(c(y, -0.5, 0.5), c(y, -0.5, 0.5))

    P$cones <- matrix(list("QUAD",
                      c(n + 2 * p + 3, (2 * p + 1):(2 * p + n), n + 2 * p + 2)), 2, 1)
    rownames(P$cones) <- c("type", "sub")

    if(intercept){
        if(scalex){
            std = apply(x,2,sd)
            P$c = c(rep(lambda*std, 2), rep(0, n), 1/n, 0, 0, 0)
        }else{
            P$c = c(rep(lambda, 2*p), rep(0, n), 1/n, 0, 0, 0)
        }
        
        A <- as.matrix.csr(x)
        A <- cbind(A, -A, as(n, "matrix.diag.csr"), as.matrix.csr(0, n, 3), as.matrix.csr(rep(1,n),n,1) )
        A <- rbind(A, cbind(as.matrix.csr(0, 2, 2*p + n),
                           as.matrix.csr(c(-.5,-.5,1,0,0,1), 2, 3),
                           as.matrix.csr(0,2,1) ) )

        P$A <- as(A,"CsparseMatrix")

        P$bx <- rbind(c(rep(0, 2 * p), rep(-Inf, n), rep(0, 3), -Inf),
                      c(rep(Inf, 2 * p + n + 4)))

    }else{
        
        if(scalex){
            std = apply(x,2,sd)
            P$c = c(rep(lambda*std, 2), rep(0, n), 1/n, 0, 0)
        }else{
            P$c = c(rep(lambda, 2*p), rep(0, n), 1/n, 0, 0)
        }
        
        A <- as.matrix.csr(x)
        A <- cbind(A, -A, as(n, "matrix.diag.csr"), as.matrix.csr(0, n, 3))
        A <- rbind(A,cbind(as.matrix.csr(0, 2, 2*p + n),
                           as.matrix.csr(c(-.5,-.5,1,0,0,1), 2, 3)))
        P$A <- as(A,"CsparseMatrix")
        
        P$bx <- rbind(c(rep(0, 2 * p), rep(-Inf, n), rep(0, 3)),
                      c(rep(Inf, 2 * p + n + 3)))
    }

    P$dparam <- list(INTPNT_CO_TOL_REL_GAP=1e-5)
    z <- mosek(P, opts = list(verbose = verb))
    status <- z$sol$itr$solsta
    f <- z$sol$itr$xx
    bhat <- f[1:p] - f[(p + 1):(2 * p)]
    if(intercept) ahat  = f[length(f)]
    else ahat = NULL
    resid <- f[(2 * p + 1):(2 * p + n)]
    return( list(a = ahat, b=bhat) )
    
}




library("SparseM")
library("Rmosek")

n = 200
p = 20
x1 = rnorm(n*5) * 10;
x2 = runif(n*10, 0, 5)
x3 = rpois(n*5, 4)
X = matrix( c( x1, x2, x3 ), nrow = n, ncol = p )
beta = as.matrix( c(rep(1,4), rep(0,16)) )
e = as.matrix( rnorm(n) )
y = 0.25+ X %*% beta + e
lambda = 0.25


f <- lasso.mosek(X, y, lambda, intercept = TRUE)
