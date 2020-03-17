# This chunk of code is intended to implement Pseudo Poisson Maximum 
# Likelihood (PPML) estimator by Rmosek solver.


PPML = function(X, y, sigma = 1, alpha = 0.05, c = 1.1, rtol = 1e-6){
  # Iuput arguments sigma, alpha and c are used to generate the 
  # tuning prameter lambda
  
  # obtain number of observation n and number of variables p
  n = nrow(X)
  p = ncol(X)
  # The way computing lambda is from Koenker's lasso code
  # lambda <- c * sigma * 2 * sqrt(n) * qnorm(1 - alpha/(2*p))
  lambda <- 0.01
  
  # Set up the problem
  # The optimization problem is minimizing the negative of log likelihood
  # plus the l1 penalty over v, beta+, beta-. The order of v, beta+, beta-
  # is followed in the following part of the function.
  
  P = list(sense = "min")
  # The linear part of the objective
  P$c = c(-y, rep(lambda, 2*p) )
  # Linear constraint
  A = as.matrix.csr(X)
  A = cbind( as(n,"matrix.diag.csr"), -A, A)
  P$A = as(A, "CsparseMatrix")
  P$bc = rbind( rep(0, n), rep(0, n)  )
  P$bx = rbind( c( rep(-Inf, n), rep(0, 2*p) ), c( rep(Inf, 2*p + n)  ) )
  
  # The Exponential part of objective
  NUMOPRO = n 
  opro = matrix(list(), nrow = 5, ncol = NUMOPRO)
  rownames(opro) = c("type", "j" , "f", "g", "h")
  for(i in 1:n){
    opro[,i] = list("EXP", i, 1/n, 1.0, 0)
  }
  P$scopt = list(opro=opro)
  
  #Tolerance
  P$dparam$intpnt_nl_tol_rel_gap <- rtol
  
  # Solve the optimization problem by mosek
  mosek_result = mosek(P)
  
  # Extract result from mosek_result
  status = mosek_result$sol$itr$solsta
  f = mosek_result$sol$itr$xx
  beta = f[(n+1):(n+p)] - f[(n+p+1):(n+2*p)]
  v = f[1:n]
  
  # Output result
  list(beta = beta, v = v, status = status)
}


