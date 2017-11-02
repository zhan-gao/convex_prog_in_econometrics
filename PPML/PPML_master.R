# Data generating process for PPML
library("SparseM")
library("Rmosek")

n = 200
p = 3


x1 = abs( rnorm(n*2) );
x2 = runif(n, 0, 5)


X = matrix( c( x1, x2 ), nrow = n, ncol = p )
beta = as.matrix( rep(1,3) );
#beta = as.matrix( c(rep(1,4), rep(0,16)) )
L = X %*% beta

# y = exp( L ) + as.matrix( rnorm(n, mean = 0, sd = 0.5) )

y = as.matrix( rpois(n, exp(L) ) ) + as.matrix( rnorm(n, mean = 0, sd = 0.5) )

source("PPML.R")

result <- PPML(X,y)