    # This is the master file for simulation.

library("Rmosek")
library("SparseM")
library("Matrix")
library("nloptr")
library("CVXR")

source("DGP_static.R")
source("PLS_est.R")
source("tools_func.R")

p <- 2;
N.cut <- c(0.3,0.6,1);
N.frac <- c(0.3,0.3,0.4);
K <- length(N.cut);
a0 <- matrix(c(0.4, 1.6, 1,1, 1.6, 0.4), nrow = p);
Rep <- 30;
MaxIter <- 500;
    
correct.ratio <- se.record <- array(0,c(6,Rep,3));

case = 0;

time.mosek <- 0;
time.nlopt <- 0;
time.cvxr <- 0;


for( N in c(100,200) ){
    
    group0 <-  rep(1:K, N*N.frac);
    
    for( TT in c(10, 20, 40) ){
        
        case = case + 1;
        
        for( r in 1:Rep ){
            
            print( paste( "Iteration: Case",as.character(case)," - ", as.character(r), 
                          "/", as.character(Rep), sep = " ") );

            # To make sure the data used in the R environment and Matlab environment are the same.
            # We generate data first, and let programs in both environments read the data from csv files.
            
		    filename = paste0('./simu_data/pls_data_N_', N, '_T_', TT, '_r_', r, '.csv')
		    D = read.csv(filename, header = TRUE)
		    
		    y = as.matrix(D[ , 1])
		    X = as.matrix(D[ , -1])


            t.begin <- Sys.time();
            pls <- PLS.mosek(N, TT, y, X, K, lambda, MaxIter);
            t.end <- Sys.time(); time.record[case, r,1] <- t.end - t.begin; time.mosek <- time.mosek + ( t.end - t.begin );
            print(time.mosek)
            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), N, N.frac, K, p );
            correct.ratio[case, r, 1] <- result.temp$ratio;
            se.record[case, r, 1] <- result.temp$se;
            
            t.begin <- Sys.time();
            pls <- PLS.nlopt(N, TT, y, X, K, lambda, MaxIter);
            t.end <- Sys.time(); time.record[case, r,2] <- t.end - t.begin; time.nlopt <- time.nlopt + ( t.end - t.begin );
            print(time.nlopt)
            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), N, N.frac, K, p );
            correct.ratio[case, r, 2] <- result.temp$ratio;
            se.record[case, r, 2] <- result.temp$se;
            
            t.begin <- Sys.time();
            pls <- PLS.cvxr(N, TT, y, X, K, lambda, MaxIter);
            t.end <- Sys.time(); time.record[case, r,3] <- t.end - t.begin; time.cvxr <- time.cvxr + ( t.end - t.begin );
            print(time.cvxr)
            result.temp <- group.coerce(pls$group.est, pls$a.out, group0, t(a0), N, N.frac, K, p );
            correct.ratio[case, r, 3] <- result.temp$ratio;
            se.record[case, r, 3] <- result.temp$se;
        }
    }
}

report <- function(ratio, se, Rep){
    ratio.mean <- apply(ratio, 1, mean);
    ratio.quantile <- apply(ratio, 1, quantile);
    RMSE <- sqrt( apply(se, 1, sum) / Rep )
    
    return( list(mean = ratio.mean, quantile = ratio.quantile, rmse = RMSE) )
}

 
save.image(file = "Result_comparison.RData")