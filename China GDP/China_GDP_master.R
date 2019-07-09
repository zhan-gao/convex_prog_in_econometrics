library("classo")
library("CVXR")

tol <- 1e-4
MaxIter <- 80 # In Wei's code, R = 80: Need to check whether converge or not
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

case <- 1

y <- as.matrix(data$ngdp)
X1 <- as.matrix(data[, -c(1:3, 9)])
X2 <- as.matrix(data[, -c(1:4, 9)])
X3 <- as.matrix(data[, -c(1:5, 9)])

# Standardize data
yy_norm <- data.normalization(y, N, TT)
yy <- yy_norm$y
yy_raw <- yy_norm$y.raw

XX1_norm <- data.normalization(X1, N, TT)
XX1 <- XX1_norm$y
XX1_raw <- XX1_norm$y.raw

XX2_norm <- data.normalization(X2, N, TT)
XX2 <- XX2_norm$y
XX2_raw <- XX2_norm$y.raw

XX3_norm <- data.normalization(X3, N, TT)
XX3 <- XX3_norm$y
XX3_raw <- XX3_norm$y.raw


# Initial estimator
beta0_1 <- init_est(XX1, yy, TT)
beta0_2 <- init_est(XX2, yy, TT)
beta0_3 <- init_est(XX3, yy, TT)

if(case == 1){
    XX <- XX1
    XX_raw <- XX1_raw
    beta0 <- beta0_1
} else if (case == 2){
    XX <- XX2
    XX_raw <- XX2_raw
    beta0 <- beta0_2
} else if (case == 3){
    XX <- XX3
    XX_raw <- XX3_raw
    beta0 <- beta0_3
} else{
    stop("Wrong case number!")
}

p <- ncol(XX)

# Parameter candidates
c_min <- 0.001
c_max <- 0.01
num_c <- 10
ss <- (c_max / c_min) ^ (1 / (num_c-1) )
c_seq <- c_min * ss^(0:(num_c-1))

# Result container
IC_total <- matrix(0, K_max, num_c)
Time <- matrix(0, K_max, num_c)

for(c in 1:num_c){
    
    for(K in 1:K_max){
        
        print( paste(as.character(c), "/", as.character(num_c), "th parameter; K = ", as.character(K)) )
        
        if(K == 1){
            a <- lsfit(XX, yy, intercept = FALSE)$coefficients
            bias <- SPJ_PLS(TT, yy_raw, XX_raw)
            a_corr <- 2*a - bias
            
            IC_total[K, ] <- mean( (yy - XX %*% a_corr)^2 )
            
            next
        }
        
        lambda <- as.numeric( c_seq[c] * var(yy) * TT^(-1/3) )
        
        t0 <- Sys.time()
        pls_out <- PLS.mosek(N,
                             TT,
                             yy,
                             XX,
                             K,
                             lambda,
                             beta0,
                             MaxIter,
                             tol,
                             post_est = FALSE,
                             bias_corr = FALSE)
        Time[K, c] <- Sys.time() - t0
        
        Q <- rep(1e10, K)
        
        # Post-estimation
        for (k in 1:K){
            
            group_k <- (pls_out$group.est == k)
            
            if(sum(group_k) >= 2*p/TT){
                
                Ind <- 1:N
                group_ind <- Ind[group_k]
                data_ind <- as.numeric( sapply(group_ind, function(i){((i-1)*TT+1):(i*TT)}) )
                yy_k <- yy[data_ind]
                XX_k <- XX[data_ind, ]
                yy_raw_k <- yy_raw[data_ind]
                XX_raw_k <- XX_raw[data_ind, ]
                
                # bias correction
                bias_k <- SPJ_PLS(TT, yy_raw_k, XX_raw_k)
                a_k <- lsfit(XX_k, yy_k, intercept = FALSE)$coefficients
                a_corr_k <- 2*a_k - bias_k
            } else {
                a_corr_k <- pls_out$a.out[k, ]
            }
            
            Q[k] <- sum( (yy_k - XX_k %*% a_corr_k)^2 )
        }
        
        IC_total[K, c] <- sum(Q) / (N*TT)
    }
}

IC <- log(IC_total) + 2/3*(N*TT)^(-0.5) * p * matrix(rep(1:K_max, num_c), nrow = K_max)

min_ind <- which.min(IC)
c_opt <- c_seq[ceiling(min_ind / K_max)]
lambda_opt <- as.numeric( c_opt * var(yy) * TT^(-1/3) )
K_temp <- min_ind %% K_max
if (K_temp == 0) K_temp <- K_max
K_opt <- K_temp

result <- PLS.mosek(N, TT, yy, XX, K_opt, lambda_opt, beta0, MaxIter, tol)

group_result <- as.data.frame( readxl::read_xls("./China GDP/province.xls", col_names = FALSE) )
group_result[[2]] <- result$group_est
colnames(group_result) <- c("province","group")

write.csv(group_result, "./China GDP/group_result.csv")


