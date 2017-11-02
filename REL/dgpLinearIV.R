LinearIV <- function(b, n, m){
    
    # Input Args: 
    # b: ture parameter of interests
    # n: sample size, m: number of IVs
    
    # DGP check - July 31, 2017
    
	
	# IV and error terms
	Z <- matrix( rnorm(n*m), n, m );
	Sigma <- sqrtm( matrix( c(0.25, 0.15, 0.15, 0.15, 0.25, 0, 0.15, 0, 0.25), 3, 3 ) );
	E <- matrix( rnorm(n*3), n, 3 ) %*% Sigma;
	e0 <- E[ , 1]; e1 <- E[ , 2]; e2 <- E[ , 3];

	# Structural Equation

	x1 <- Z[ , 1:2] %*% matrix( c(0.5,0.5), 2, 1) + e1;
	x2 <- Z[ , 3:4] %*% matrix( c(0.5,0.5), 2, 1) + e2;
	X <- cbind(x1, x2);

	y <- X %*% b + e0;

	return( list(y = y, X = X, Z = Z) )
}



