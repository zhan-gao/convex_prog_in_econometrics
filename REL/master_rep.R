rm(list = ls())

library("Rmosek")
library("SparseM")
library("Matrix")
library("nloptr")
library("expm")
library("CVXR")

source("./REL/dgpLinearIV.R")
source("./REL/REL.R")

set.seed(99)

b0 = c(1,1) 
Rep = 500;

case = 0 

time.record = B =  matrix(0, Rep, 4)

for( n in c(120,240) ){

	for( m in c(80,160) ){
		
		case = case + 1

		# We follow the generic parameter tuning scheme in Section 4:
		# tau = c * sqrt( log(m) / n ), c = 0.5
		tau = 0.5 * sqrt( log(m) / n )

		for( r in 1:Rep ){

			print( paste( "Iteration: Case", as.character(case)," - ", as.character(r), 
                          "/", as.character(Rep), sep = " ") )
			# Generate data
			g(y, X, Z) %=% LinearIV(b0, n, m)
		    
			# Estimation (report only the estimation b_1)
			t.begin = Sys.time() 
			b_hat = REL(y, X, Z, tau)
			t.end = Sys.time();
			
			B[r, case] = b_hat[1]
			time.record[r, case] = t.end - t.begin
			print( t.end - t.begin )

		}
	}
}

bias.mosek = apply(B, 2, mean) - b0[1]
RMSE.mosek = sqrt( apply( (B - b0[1])^2, 2, mean) )

write.csv(cbind(bias.mosek, RMSE.mosek), "REL_Result_Rep.csv")

save.image(file = "REL_Result_Rep.RData")


