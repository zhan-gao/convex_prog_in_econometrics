rm(list= ls())

source("./REL/REL.R")
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

T.Result = L.Result = array(0, c(100, 4, 3))

case = 0
for(n in c(120, 240)){
	for(m in c(80, 160)){

		case = case + 1

		filename = paste0('./REL/Data_', n, '_', m, '.csv')
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

			T.Result[r, case, ] = c(t.mosek, t.cvxr, t.nloptr)
			L.Result[r, case, ] = c(L.mosek, L.cvxr, L.nloptr)

		}

	}
}

Result.mosek = cbind(L.Result[, , 1], T.Result[, , 1])
Result.cvxr = cbind(L.Result[, , 2], T.Result[, , 2])
Result.nloptr = cbind(L.Result[, , 3], T.Result[, , 3])

Time.Sum = apply(T.Result, c(2,3), sum)
colnames(Time.Sum) = c("Rmosek", "CVXR", "nloptr")
write.csv(Time.Sum, "REL_time_compare_R.csv")
save.image("REL_Compare_Result.RData")
