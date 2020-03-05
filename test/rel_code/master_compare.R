# R version 3.6.2 (2019-12-12)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 17763)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United States.1252 
# [2] LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#     [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#     [1] CVXR_1.0      expm_0.999-4  nloptr_1.2.2  SparseM_1.78  Rmosek_9.0.96
# [6] Matrix_1.2-18
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.3      lattice_0.20-38 gurobi_9.0-0    gmp_0.5-13.6   
# [5] slam_0.1-47     grid_3.6.2      R6_2.4.1        Rmpfr_0.8-1    
# [9] tools_3.6.2     bit64_0.9-7     bit_1.1-15.2    compiler_3.6.2 
# [13] ECOSolveR_0.5.3


source("REL.R")
# Compare CVX + Mosek to Rmosek
# Read the contains genrated data: 100 replications
# Focus on innerloop only with beta = [.9; .9]

library("Rmosek")
library("SparseM")
library("Matrix")
library("nloptr")
library("expm")
library("CVXR")

b = c(0.9, 0.9)

T.Result = L.Result = array(0, c(100, 4, 4))

case = 0
for(n in c(120, 240)){
	for(m in c(80, 160)){

		case = case + 1

		filename = paste0('Data_', n, '_', m, '.csv')
		D = read.csv(filename, header = FALSE)

		tau = 0.5 * sqrt( log(m) / n )

		for(r in 1:100){

			cat("n = ", n , "; m = ", m, "; r = ", r, " / 100\n")

			id = ((r-1)*n+1):(r*n)
			y = D[id, 1]
			X = as.matrix(D[id, 2:3])
			Z = as.matrix(D[id, -(1:3)])

			t0 = Sys.time()
			L.mosek = -innerloop(b, y, X, Z, tau)
			t.mosek = Sys.time() - t0

			t0 = Sys.time()
			L.cvxr = -innerloop.cvxr(b, y, X, Z, tau)
			t.cvxr = Sys.time() - t0

			t0 = Sys.time()
			L.nloptr = -innerloop.nloptr(b, y, X, Z, tau)
			t.nloptr = Sys.time() - t0
			
			t0 = Sys.time()
			L.cvxr.nodcp = -innerloop.cvxr.nodcp(b, y, X, Z, tau)
			t.cvxr.nodcp = Sys.time() - t0
			

			T.Result[r, case, ] = c(t.mosek, t.cvxr, t.cvxr.nodcp, t.nloptr)
			L.Result[r, case, ] = c(L.mosek, L.cvxr, L.cvxr.nodcp, L.nloptr)

		}

	}
}

Result.mosek = cbind(L.Result[, , 1], T.Result[, , 1])
Result.cvxr = cbind(L.Result[, , 2], T.Result[, , 2])
Result.cvxr.nodcp = cbind(L.Result[, , 3], T.Result[, , 3])
Result.nloptr = cbind(L.Result[, , 4], T.Result[, , 4])

Time.Sum = apply(T.Result, c(2,3), sum)
colnames(Time.Sum) = c("Rmosek", "CVXR", "CVXR_nodcp", "nloptr")
write.csv(Time.Sum, "REL_time_compare_R.csv")
save.image("REL_Compare_Result.RData")
