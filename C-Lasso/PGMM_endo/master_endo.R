library("Rmosek")
library("SparseM")
library("Matrix")
library("MASS")


source("DGP_endo.R")
source("PGMM_est.R")
source("PLS_est.R")
source("tools_func.R")



p = 2
d = 2
N.cut = c(0.3,0.6,1)
N.frac = c(0.3,0.3,0.4)
K = length(N.cut)
a0 = matrix(c(0.2, 1.8, 1,1, 1.8, 0.2), nrow = p)
Rep =500
MaxIter = 500


correct.ratio = matrix(0,6,Rep)
se.record = matrix(0,6,Rep)
time.record = matrix(0,6,Rep)


case = 0;


for( N in c(100,200) ){
    
    group0 =  rep(1:K, N*N.frac);
    
    for( t in c(15, 25, 50) ){
        
        case = case + 1
        
        for( r in 1:Rep ){
            
            print( paste( "Iteration: Case",as.character(case)," - ", as.character(r), 
                          "/", as.character(Rep), sep = " ") )
            
            data <- DGP.endo(N, t, p, d, N.cut, a0)
            
            g(y,X,Z) %=% data
            
            lambda = as.numeric( 0.5 * var(y) / (t^(1/3)) )
            
            t.begin = Sys.time()
            pls = PGMM.mosek(N, t, y, X, Z, K, lambda, MaxIter)
            t.end = Sys.time()
            time.record[case, r] = t.end - t.begin
            print(t.end - t.begin)
            
            result.temp = group.coerce(pls$group.est, pls$a.out, group0, 
                                        t(a0), N, N.frac, K, p )
            correct.ratio[case, r] <- result.temp$ratio;
            se.record[case, r] <- result.temp$se;
               
        }
    }
}

report <- function(ratio, se, Rep){
    ratio.mean <- apply(ratio, 1, mean);
    ratio.quantile <- apply(ratio, 1, quantile);
    RMSE <- sqrt( apply(se, 1, sum) / Rep )
    
    return( list(mean = ratio.mean, quantile = ratio.quantile, rmse = RMSE) )
}


save.image(file = "Result.RData")


print("mosek result:")
print( report(correct.ratio, se.record, Rep) )
print( sum(time.record))