library("tidyverse")
library("Rmosek")
library("SparseM")
library("Matrix")
library("nloptr")
library("CVXR")

source("./CLasso/tools_func.R")
source("./CLasso/PLS_est.R")

tol <- 1e-4
MaxIter <- 500 # In Wei's code, R = 80: Need to check whether converge or not
K_max <- 4

data <- read.csv("./China GDP/Data_for_CLasso_0007.csv", header = FALSE)
# ==================================
# Variables included:
# ngdp : log term of norminal provincial GDP
# light : log term of night light
# ntax : log term of national tax
# expor : log term of export
# impor : log term of import
# ecos : log term of electricity consumption
# rail : log term of railway cargo volumn
# ==================================
colnames(data) <- c('year', 'code', 
                    'ngdp', 'light', 'ntax', 'expor', 'impor', 'ecos', 'rail')
N <- max(data$code)
TT <- max(data$year)

# CASE 1: (y, X1) all variables included
# CASE 2: (y, X2) drop light
# CASE 3: (y, X3) drop light and national tax
y <- as.matrix(data$ngdp)
X1 <- as.matrix(data[, -c(1:3)])
X2 <- as.matrix(data[, -c(1:4)])
X3 <- as.matrix(data[, -c(1:5)])

# Standardize data
yy <- data.normalization(y, N, TT)
XX1 <- data.normalization(X1, N, TT)
XX2 <- data.normalization(X2, N, TT)
XX3 <- data.normalization(X3, N, TT)

# Parameter candidates
c_min <- 0.001
c_max <- 1
num_c <- 100
ss <- (c_max / c_min) ^ (1 / (num_c-1) )
c_seq <- c_min * ss^(0:(num_c-1))

# Result container
IC_total = matrix(0, K_max, num_c)

# Apr 14:
# Need to add: 1/ Half panel jackknife to do bias correction
#              2/ variance estimation.


for(c in 1:num_c){
    for(K in 1:K_max){
        
        print( paste(as.character(ll), "th parameter; K = ", as.character(K)) )
        
        lambda <- as.numeric( c_seq[c] * var(y) * TT^(-1/3) )
        pls_out <- PLS.mosek(N, TT, y, X, K, lambda, MaxIter, tol)
        
    }
}
    
        