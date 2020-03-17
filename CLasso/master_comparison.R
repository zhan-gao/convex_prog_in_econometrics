rm(list = ls())
# This is the master file for simulation.

library("Rmosek")
library("SparseM")
library("Matrix")
library("CVXR")

source("./CLasso/PLS_est.R")
source("./CLasso/tools_func.R")

p <- 2
N.cut <- c(0.3, 0.6, 1)
N.frac <- c(0.3, 0.3, 0.4)
K <- length(N.cut)
a0 <- matrix(c(0.4, 1.6, 1, 1, 1.6, 0.4), nrow = p)
Rep <- 100
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
            
            filename = paste0("./CLasso/simu_data_full/pls_data_N_", N, "_T_", TT, 
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
            # 
            # # CVXR + Mosek
            # t.begin <- proc.time()
            # pls <- PLS.cvxr(N, TT, y, X, K, lambda, MaxIter, solver = "MOSEK")
            # t.end <- proc.time()
            # 
            # time.record[case, r, 4] <- (t.end - t.begin)[3]
            # 
            # result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), 
            #                             N, N.frac, K, p)
            # correct.ratio[case, r, 4] <- result.temp$ratio
            # se.record[case, r, 4] <- result.temp$se
            
            
            print(time.record[case, r, ])
        }
    }
}


rmse = sqrt(apply(se.record, c(1, 3), mean))
cr = apply(correct.ratio, c(1, 3), mean)
time.sum = apply(time.record, c(1, 3), sum)

write.csv(rbind(rmse, cr, time.sum), "PLS_Result_comparison.csv")
save.image(file = "PLS_Result_comparison.RData")
