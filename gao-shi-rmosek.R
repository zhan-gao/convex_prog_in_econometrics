# ==========================================================================
# The master file replicates the simulation results in
# "Implementing Convex Optimization in R: Two Econometric Examples" 
# by Zhan Gao and Zhentao Shi
# ==========================================================================

# ==========================================================================
# Replication of Table 1: CLasso Results
# ==========================================================================

# To Replicate the results in Table 1 with only Rmosek
# The results are in "PLS_Result_rep.csv"
source("./CLasso/master_rep.R") 

# Generate data explicitly since we need to compare across platforms 
# 	(R v.s Matlab v.s. Python)
# source("./CLasso/data_gen.R")

# For the comparison among all methods, we provide a small scale sample 
#	(30 replications)
# To generate the full 500 replication data, uncomment the data_gen.R line
# 	and generate data explicitly first, and then change the 
#	"./CLasso/master_comparison.R", "./CLasso/master_cvx.m" 
# 	"./CLasso/cvxpy_master.py" by change the variable Rep from 30 to 500

# The CVX results are generated in Matlab: Run "./CLasso/master_cvx.m" in Matlab
# The results are stored in "CVX_PLS_Result.csv"

# The CVXPY results are generated in Python: Run "./CLasso/cvxpy_master.py" in Python
# The results are stored in "python_result.csv"

# For CVXR and Rmosek
# The results are saved in "PLS_Result_comparison.csv"
source("./CLasso/master_comparison.R")

# ==========================================================================
# Replication of Table 2 and 3: Empirical application
# ==========================================================================

# The Rmosek replication
devtools::install_github("zhan-gao/classo", INSTALL_opts=c("--no-multiarch"))
source("./China GDP/China_GDP_master.R")

# The MATLAB implementation can be done by "./CLasso/China GDP/China gdp/main_ChinaGDP_PLS.m" in MATLAB

# ==========================================================================
# Replication of Table 4: REL Results
# ==========================================================================

# The bias and RMSE by Rmosek documented in Table 4 (left panel) are stored
# in "REL_Result_Rep.csv"
# The workplace is saved in "REL_Result_Rep.RData"
source("./REL/master_rep.R")

# ==========================================================================
# Replication of Table 3: REL inner-loop time comparison
# ==========================================================================

# Time documented in Table 3 is stored in "REL_time_compare_R.csv"
# The workplace is saved in "REL_Compare_Result.RData"
# The CVX results are generated in Matlab: Run "./REL/master_cvx.m" in Matlab
source("./REL/master_compare.R")

