# > sessionInfo()
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
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] CVXR_1.0
# 
# loaded via a namespace (and not attached):
# [1] bit_1.1-15.2    compiler_3.6.2  R6_2.4.1        Matrix_1.2-18  
# [5] tools_3.6.2     gmp_0.5-13.6    Rmosek_9.0.96   Rcpp_1.0.3     
# [9] slam_0.1-47     bit64_0.9-7     grid_3.6.2      gurobi_9.0-0   
# [13] ECOSolveR_0.5.3 Rmpfr_0.8-1     lattice_0.20-38

# First time run
# user  system elapsed 
# 0.53    0.07    0.60 
# user  system elapsed 
# 0.30    0.04    0.33

# Second time run
# user  system elapsed 
# 0.29    0.00    0.29 
# user  system elapsed 
# 0.3     0.0     0.3 

suppressMessages(suppressWarnings(library(CVXR)))

data(cdiac)
y <- cdiac$annual
m <- length(y)
lambda <- 0.44
beta <- Variable(m)
obj <- 0.5 * sum((y - beta)^2) + lambda * sum(pos(diff(beta)))
prob <- Problem(Minimize(obj))

t0 <- proc.time()
soln <- solve(prob, solver = "ECOS")
betaHat <- soln$getValue(beta)
print(proc.time() - t0)

t0 <- proc.time()
prob_data <- get_problem_data(prob, solver = "ECOS")
ECOS_dims <- ECOS.dims_to_solver_dict(prob_data$data[["dims"]])
solver_output <- ECOSolveR::ECOS_csolve(c = prob_data$data[["c"]],
                                        G = prob_data$data[["G"]],
                                        h = prob_data$data[["h"]],
                                        dims = ECOS_dims,
                                        A = prob_data$data[["A"]],
                                        b = prob_data$data[["b"]])
direct_soln <- unpack_results(prob, solver_output, prob_data$chain, prob_data$inverse_data)
beta_hat <- direct_soln$getValue(beta)
print(proc.time() - t0)