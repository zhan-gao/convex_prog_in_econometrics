
source("./CLasso/DGP_static.R")

p <- 2
N.cut <- c(0.3,0.6,1)
N.frac <- c(0.3,0.3,0.4)
K <- length(N.cut)
a0 <- matrix(c(0.4, 1.6, 1,1, 1.6, 0.4), nrow = p)

set.seed(200)

for(N in c(100, 200)){
    for(TT in c(15, 25 ,50)){
 
        for(r in 1:500){
            
            d <- DGP.static(N, TT, p, N.cut, a0)
            y <- d$y
            X <- d$X
            
            file.name = paste0("./CLasso/simu_data/pls_data_N_", N, "_T_", TT, "_r_", r, ".csv")
            write.csv(cbind(y, X), file = file.name, row.names = FALSE, col.names = FALSE)
        }
        
    }
}




