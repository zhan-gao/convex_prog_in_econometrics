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
Rep <- 500;
MaxIter <- 500;

correct.ratio = se.record = time.record = matrix(0, 6, Rep)

case = 0

set.seed(200)

for( N in c(100,200) ){
    
    group0 <-  rep(1:K, N*N.frac);
    
    for( TT in c(15, 25, 50) ){
        
        case = case + 1;
        
        for( r in 1:Rep ){
            
            print( paste( "Iteration: Case",as.character(case)," - ", as.character(r), 
                          "/", as.character(Rep), sep = " ") );
            
            d <- DGP.static(N, TT, p, N.cut, a0)
            y <- d$y
            X <- d$X
            lambda <- as.numeric( 0.5 * var(y) / (TT^(1/3)) )
            
            t.begin <- Sys.time()
            pls <- PLS.mosek(N, TT, y, X, K, lambda, MaxIter)
            t.end <- Sys.time()

            time.record[case, r] <- t.end - t.begin

            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), N, N.frac, K, p )
            correct.ratio[case, r] <- result.temp$ratio
            se.record[case, r] <- result.temp$se
            
        }
    }
}

report <- function(ratio, se, Rep){
    ratio.mean <- apply(ratio, 1, mean);
    ratio.quantile <- apply(ratio, 1, quantile);
    RMSE <- sqrt( apply(se, 1, sum) / Rep )
    
    return( list(mean = ratio.mean, quantile = ratio.quantile, rmse = RMSE) )
}

Result.Rep = report(correct.ratio, se.record, Rep)
write.csv(rbind(Result.Rep$mean, Result.Rep$rmse, rowMeans(time.record)), "PLS_Result_rep.csv")
save.image(file = "PLS_Result_rep.RData")