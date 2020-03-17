DGP.endo = function(N, TT, p, d, N.cut, a0){

	# panel data with endogeneity: DGP
	# DGP 4 in appendix of SSP (2016)

	# require("MASS")
    
    # INPUT Arg:
    # N: num of obs
    # TT: time horizon
    # N_cut: True group structure
    # a0: True parameter
    
    # Outputs:
    # X: TN * p,
    # y: TN * 1
    # Z: TN * (d+1)
    # y_{i=1,t=1}, ..., y_{i=1, t=T}, y_{i=2, t=1}, ...


	K0 = length(N.cut)

	T1 = TT + 1


	E = mvrnorm( n = T1*N, 
				 mu = c(0,0), 
				 Sigma =  matrix(c(1,0.3,0.3,1),2,2) )
	g(epsl, e) %=% lapply( seq_len(ncol(E)), function(i) E[,1] )


	u = rnorm(N)
	uu = rep(u, times = rep(T1,N))

	ZZ = matrix(rnorm(T1*N*(d+1)), T1*N, d+1)

	g(X2, Z) %=% list( ZZ[,1], ZZ[,c(2,3)] )

	X1 = 0.2 * uu + 0.5 * rowSums(Z) + 0.5 * e
	XX = cbind(X1, X2)

	y.temp = XX %*% a0
	yy = matrix(0,T1*N,1)

	for( k in 1:K0 ){

		cut = c(0, N.cut)
		ind = ( cut[k]*N*T1+1) : (cut[k+1]*N*T1)
		yy[ind] = uu[ind] + y.temp[ind, k] + epsl[ind] 

	}

	id = rep( c(FALSE, rep(TRUE, TT)), N )
	id.minus = rep( c(rep(TRUE, TT), FALSE), N)

	Dy = yy[id, drop = FALSE] - yy[id.minus, drop = FALSE]
	DX = XX[id, ] - XX[id.minus, ]
	DZ = ZZ[id, ] - ZZ[id.minus, ]

	return( list(y = Dy, X = DX, Z = DZ) )

}