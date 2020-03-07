suppressMessages(suppressWarnings(library(CVXR)))

data(cdiac)
y <- cdiac$annual
m <- length(y)
lambda <- 0.44
beta <- Variable(m)
obj <- 0.5 * sum((y - beta)^2) + lambda * sum(pos(diff(beta)))
prob <- Problem(Minimize(obj))

t0 <- Sys.time()
soln <- solve(prob, solver = "ECOS")
betaHat <- soln$getValue(beta)
print(Sys.time() - t0)

t0 <- Sys.time()
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
print(Sys.time() - t0)