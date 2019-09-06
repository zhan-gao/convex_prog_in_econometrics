library("Rmosek")
library("SparseM")

# Main function
logistic_mosek <- function(y, x, rtol = 1e-6) {
    n <- nrow(x)
    p <- ncol(x)
    
    # maximize the standard logistic regression objective
    # variables:
    # t_1, ..., t_n
    # phi_1, ..., phi_n
    # u_1, ...,u_n
    # v_1, ..., v_n
    # beta_1, ..., beta_p
    # alpha_1, ...,alpha_2n
    # gamma_1, ..., gamma_n
    
    prob <- list(sense = "max")
    prob$c <- c(rep(1, n), y, rep(0, 5 * n + p))
    
    xx <- as.matrix.csr(-x)
    A <- rbind(cbind(
        as(matrix(0, n, n), "matrix.csr"),
        as(n, "matrix.diag.csr"),
        as(matrix(0, n, 2 * n), "matrix.csr"),
        xx,
        as(matrix(0, n, 3 * n), "matrix.csr")
    ),
    cbind(
        as(matrix(0, n, 2 * n), "matrix.csr"),
        as(n, "matrix.diag.csr"),
        as(n, "matrix.diag.csr"),
        as(matrix(0, n, 3 * n + p), "matrix.csr")
    ),
    cbind(
        as(n, "matrix.diag.csr"),
        as(n, "matrix.diag.csr"),
        as(matrix(0, n, 4 * n + p), "matrix.csr"),-as(n, "matrix.diag.csr")
    ))
    prob$A <- as(A, "CsparseMatrix")
    # Bounds for linear constraints
    prob$bc <- rbind(blc = c(rep(0, n), rep(-Inf, n), rep(0, n)),
                     buc = c(rep(0, n), rep(1, n), rep(0, n)))
    # Bounds for variables
    prob$bx <-
        rbind(
            blx = c(
                rep(-Inf, n),
                rep(-Inf, n),
                rep(0, 2 * n),
                rep(-Inf, p),
                rep(1, 2 * n),
                rep(-Inf, n)
            ),
            bux = c(
                rep(0, n),
                rep(Inf, n),
                rep(1, 2 * n),
                rep(Inf, p),
                rep(1, 2 * n),
                rep(Inf, n)
            )
        )
    
    # The cones
    NUMCONES <- 2 * n
    prob$cones <- matrix(list(), nrow = 2, ncol = NUMCONES)
    rownames(prob$cones) <- c("type", "sub")
    for (i in 1:n) {
        prob$cones[, i] <- list("PEXP", c(2 * n + i, 4 * n + p + i, 6 * n + p +
                                              i))
    }
    for (i in 1:n) {
        prob$cones[, n + i] <- list("PEXP", c(3 * n + i, 5 * n + p + i, i))
    }
    
    # Invoke the solver
    res <- mosek(prob, opts = list(verbose = 0))
    
    return(res$sol)
}

# --------------------------------
# Monte Carlo Experiments
# --------------------------------
dgp <- function(n, b0){
    
    p <- length(b0)
    x <- matrix(rnorm(n*p), n, p)
    z <- as.numeric(x %*% b0)
    
    pr <- exp(z)/(1+exp(z))        
    y <-  rbinom(n,1,pr)
    
    return(list(y = y, x = x))
}


set.seed(99)
b0 <- c(2,3)
n <- 100
p <- length(b0)

R <- 1000
bhat <- matrix(0, R, p)

for(r in 1:R){
    
    cat("Iteration ", r, " / ", R, "\n")
    
    d <- dgp(n, b0)
    x <- d$x
    y <- d$y
    
    res <- logistic_mosek(y, x)
    bhat[r, ] <- res$itr$xx[(4*n+1):(4*n+p)]    
}

plot(density(bhat[, 1]))
plot(density(bhat[, 2]))