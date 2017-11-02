
library("Rmosek")
library("SparseM")
library("Matrix")
library("nloptr")
library("expm")
library("cvxr")

source("dgpLinearIV.R")
source("REL.R")

set.seed(99)

b0 <- c(1,1); Rep <- 500;

case <- 0; time.record <- array(0, c(Rep,4,2));

B <- array(0, c(Rep,4,2));

for( n in c(120,240) ){

	for( m in c(80,160) ){
		
		case <- case + 1;

		# We follow the generic parameter tuning scheme in Section 4:
		# tau = c * sqrt( log(m) / n ), c = 0.5
		tau <- 0.5 * sqrt( log(m) / n );

		for( r in 1:Rep ){

			print( paste( "Iteration: Case", as.character(case)," - ", as.character(r), 
                          "/", as.character(Rep), sep = " ") );
			# Generate data
			g(y, X, Z) %=% LinearIV(b0, n, m);
		    
			# Estimation (report only the estimation b_1)
			t.begin <- Sys.time(); b_hat <- REL(y, X, Z, tau); t.end <- Sys.time();
			B[r, case, 1] <- b_hat[1]; time.record[r, case, 1] <- t.end - t.begin;
			print( t.end - t.begin );
			
			t.begin <- Sys.time(); b_hat <- REL.cvxr(y, X, Z, tau); t.end <- Sys.time();
			B[r, case, 2] <- b_hat[1]; time.record[r, case, 2] <- t.end - t.begin;
			print( t.end - t.begin );

		}
	}
}

bias.mosek <- apply(B[,,1], 2, mean) - b0[1];
RMSE.mosek <- sqrt( apply( (B[,,1] - b0[1])^2, 2, mean) );
bias.cvxr <- apply(B[,,2], 2, mean) - b0[1];
RMSE.cvxr <- sqrt( apply( (B[,,2] - b0[1])^2, 2, mean) );

save.image(file = "Result_cvxr.RData")

print( "mosek result:" )
print( paste( "Bias: ", as.character(bias.mosek), sep  = " " ) )
print( paste( "RMSE: ", as.character(RMSE.mosek), sep  = " " ) )
print( paste( "Time: ", as.character(apply(time.record[,,1],2,sum)), sep = " " ) )


print( "CVXR result:" )
print( paste( "Bias: ", as.character(bias.cvxr), sep  = " " ) )
print( paste( "RMSE: ", as.character(RMSE.cvxr), sep  = " " ) )
print( paste( "Time: ", as.character(apply(time.record[,,2],2,sum)), sep = " " ) )
