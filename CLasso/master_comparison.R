rm(list = ls())
# This is the master file for simulation.

library("Rmosek")
library("SparseM")
library("Matrix")
library("nloptr")
library("CVXR")

source("./CLasso/DGP_static.R")
source("./CLasso/PLS_est.R")
source("./CLasso/tools_func.R")

p <- 2;
N.cut <- c(0.3,0.6,1);
N.frac <- c(0.3,0.3,0.4);
K <- length(N.cut);
a0 <- matrix(c(0.4, 1.6, 1,1, 1.6, 0.4), nrow = p);
Rep <- 30;
MaxIter <- 500;

correct.ratio <- se.record <- time.record <- array(0,c(6,Rep,3));
case = 0;

set.seed(200)

for( N in c(100,200) ){
    
    group0 <-  rep(1:K, N*N.frac);
    
    for( TT in c(10, 25, 50) ){
        
        case = case + 1;

        for(r in 1:Rep) {
            
            print( paste( "Iteration: Case",as.character(case)," - ", as.character(r), 
                          "/", as.character(Rep), sep = " ") );
            
            filename = paste0('./CLasso/simu_data/pls_data_N_', N, '_T_', TT, '_r_', r, '.csv')
            D = read.csv(filename, header = TRUE)
            y = as.matrix(D[ , 1])
            X = as.matrix(D[ , -1])
            lambda <- as.numeric( 0.5 * var(y) / (TT^(1/3)) )
            
            # Rmosek
            t.begin <- Sys.time()
            pls <- PLS.mosek(N, TT, y, X, K, lambda, MaxIter)
            t.end <- Sys.time()
            
            time.record[case, r, 1] <- t.end - t.begin
            
            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), N, N.frac, K, p )
            correct.ratio[case, r, 1] <- result.temp$ratio
            se.record[case, r, 1] <- result.temp$se
            
            # NLOPTR
            # WARNING: VERY SLOW
            t.begin <- Sys.time()
            pls <- PLS.nlopt(N, TT, y, X, K, lambda, MaxIter)
            t.end <- Sys.time()

            time.record[case, r, 2] <- t.end - t.begin

            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), N, N.frac, K, p )
            correct.ratio[case, r, 2] <- result.temp$ratio
            se.record[case, r, 2] <- result.temp$se
            
            # CVXR
            t.begin <- Sys.time()
            pls <- PLS.cvxr(N, TT, y, X, K, lambda, MaxIter)
            t.end <- Sys.time()
            
            time.record[case, r,3] <- t.end - t.begin
            
            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), N, N.frac, K, p )
            correct.ratio[case, r, 3] <- result.temp$ratio
            se.record[case, r, 3] <- result.temp$se
            
        }
    }
}


rmse = sqrt( apply(se.record, c(1,3), mean) )
cr = apply(correct.ratio, c(1,3), mean)
time.sum = apply(time.record, c(1,3), sum)

write.csv(rbind(rmse, cr, time.sum), "PLS_Result_comparison.csv")
save.image(file = "PLS_Result_comparison.RData")