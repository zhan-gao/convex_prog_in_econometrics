# sessionInfo()
# R version 3.6.2 (2019-12-12)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 17763)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#     [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#     [1] Rmosek_9.0.96     classo_0.0.0.9100 CVXR_1.0          Matrix_1.2-18     SparseM_1.78     
# 
# loaded via a namespace (and not attached):
# [1] bit_1.1-15.2    compiler_3.6.2  R6_2.4.1        tools_3.6.2     gmp_0.5-13.6    Rcpp_1.0.3      bit64_0.9-7    
# [8] grid_3.6.2      ECOSolveR_0.5.3 Rmpfr_0.8-1     lattice_0.20-38

rm(list = ls())
# This is the master file for simulation.

library("Rmosek")
library("SparseM")
library("Matrix")
library("CVXR")

source("PLS_est.R")
source("tools_func.R")

p <- 2
N.cut <- c(0.3, 0.6, 1)
N.frac <- c(0.3, 0.3, 0.4)
K <- length(N.cut)
a0 <- matrix(c(0.4, 1.6, 1, 1, 1.6, 0.4), nrow = p)
Rep <- 5
MaxIter <- 500

correct.ratio <- se.record <- time.record <- array(0, c(6, Rep, 4))
case = 0

set.seed(200)

for (N in c(100, 200)) {
    
    group0 <- rep(1:K, N * N.frac)
    
    for (TT in c(15, 25, 50)) {
        
        case = case + 1
        
        for (r in 1:Rep) {
            
            print(paste("Iteration: Case", as.character(case), " - ", as.character(r), 
                        "/", as.character(Rep), sep = " "))
            
            filename = paste0("simu_data_small/pls_data_N_", N, "_T_", TT, 
                              "_r_", r, ".csv")
            D = read.csv(filename, header = TRUE)
            y = as.matrix(D[, 1])
            X = as.matrix(D[, -1])
            lambda <- as.numeric(0.5 * var(y)/(TT^(1/3)))
            
            # Rmosek
            t.begin <- Sys.time()
            pls <- PLS.mosek(N, TT, y, X, K, lambda, MaxIter)
            t.end <- Sys.time()
            
            time.record[case, r, 1] <- t.end - t.begin
            
            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), 
                                        N, N.frac, K, p)
            correct.ratio[case, r, 1] <- result.temp$ratio
            se.record[case, r, 1] <- result.temp$se
            
            # CVXR + ECOS
            t.begin <- Sys.time()
            pls <- PLS.cvxr(N, TT, y, X, K, lambda, MaxIter)
            t.end <- Sys.time()
            
            time.record[case, r, 2] <- t.end - t.begin
            
            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), 
                                        N, N.frac, K, p)
            correct.ratio[case, r, 2] <- result.temp$ratio
            se.record[case, r, 2] <- result.temp$se
            
            # CVXR + ECOS + nodcp
            t.begin <- Sys.time()
            pls <- PLS.cvxr.nodcp(N, TT, y, X, K, lambda, MaxIter)
            t.end <- Sys.time()
            
            time.record[case, r, 3] <- t.end - t.begin
            
            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), 
                                        N, N.frac, K, p)
            correct.ratio[case, r, 3] <- result.temp$ratio
            se.record[case, r, 3] <- result.temp$se

            # CVXR + Mosek
            t.begin <- proc.time()
            pls <- PLS.cvxr(N, TT, y, X, K, lambda, MaxIter, solver = "MOSEK")
            t.end <- proc.time()

            time.record[case, r, 4] <- (t.end - t.begin)[3]

            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0),
                                        N, N.frac, K, p)
            correct.ratio[case, r, 4] <- result.temp$ratio
            se.record[case, r, 4] <- result.temp$se

            print(time.record[case, r, ])
        }
    }
}


rmse = sqrt(apply(se.record, c(1, 3), mean))
cr = apply(correct.ratio, c(1, 3), mean)
time.sum = apply(time.record, c(1, 3), sum)

write.csv(rbind(rmse, cr, time.sum), "PLS_Result_comparison.csv")
save.image(file = "PLS_Result_comparison.RData")
