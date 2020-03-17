library(Rmosek)
library(nloptr)

# this script apply the idea of empirical likelihood 
# through the R-MOSEk in the inner optimization

n = floor( T/blocksize )


################## functions 

inner = function(theta){
    prob$A = Matrix( rbind( 1,  t( gg.block ) ) ) # combine the sum of prob into the linear constraint
    
    prob$scopt <- list ( opro = opro  )
    
    r = mosek ( prob, opts = list(soldetail = 2, verbose = 0) )
    if (r$response$msg == "MSK_RES_ERR_USER_NLO_FUNC: The user-defined nonlinear function reported an error."){
      J = Inf}
    else{     J = -r$sol$itr$pobjval } # primal objective value 
}


############### execution 

mm = 8 # number of moments

prob = list(sense = "max")

prob$c = rep(0,n) 


prob$bc = rbind( blc = c( 1, rep(0,mm) ), buc = c( 1,  rep(0,mm) ) ) 
prob$bx = rbind( blx = rep(0,1), bux = rep(1,n) )

# the decision variables are the probablity 'p'
opro <- matrix ( list (), nrow =5, ncol = n )
rownames ( opro ) <- c(" type ","j","f","g","h")
for (i in 1:n){ opro[,i] = list("LOG", i, 1, 1, 0 ) }


#"Nelder-Mead" is also valid for a smooth problem,
# even if it takes a bit long time
# it is easier to use a brutal force method to implement the restriction
# on the parameter space than "L-BFGS-B"

opts = list("algorithm"="NLOPT_LN_NELDERMEAD","xtol_rel"=1.0e-7)
# opts = list("algorithm"="NLOPT_LN_SBPLX","xtol_rel"=1.0e-7)

res = nloptr( x0= theta0, 
                   eval_f=inner, probb = FALSE,
                   opts=opts)

theta.hat= res$solution
probs = inner(theta.hat, probb = TRUE)

print( res )


probs = probs/sum(probs)


LR = 2* (-sum(  log(probs) ) - ( n ) *log( n ) )
cat( "LR = ", LR, "p-value = ", 1-pchisq(LR, mm - 4), "\n" )
cat("\n sum probs", sum(probs), "\n")
