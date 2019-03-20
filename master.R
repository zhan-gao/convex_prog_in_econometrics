# ==========================================================================
# The master file replicates the simulation results in
# "Two examples of convex-programming-based high-dimensional
#  econometric estimators" by Zhan Gao and Zhentao Shi
# ==========================================================================

# ==========================================================================
# Replication of Table 1: CLasso Results
# ==========================================================================

# CAVEAT: As mentioned in the paper, solvers other than Rmosek are very slow 
#		  running. It is recommended to try a small scale experiment and try
# 		  each case separately.

# To Replicate the results in Table 1 with only Rmosek
# The results are in "PLS_Result_rep.csv"
source("./CLasso/master_rep.R") 

# Generate data explicitly since we need to compare across platforms 
# 	(R v.s Matlab)

# source("./CLasso/data_gen.R")

# It will generate a full set of sample for 500 replications.
# Considering ECOS is sensitive to data input sometimes, for comparison reasons,
# we check whether it works with ECOS or not when generate data.

# For the comparison among all methods, we also provide a small scale sample 
#	(30 replications)
# To generate the full 500 replication data, uncomment the data_gen.R line
# 	and generate data explicitly first, and then change the 
#	"./CLasso/master_comparison.R" and "./CLasso/master_cvx.m" by change 
#	the variable Rep from 30 to 500

# The CVX results are generated in Matlab: Run "./CLasso/master_cvx.m" in Matlab
# The results are stored in "CVX_PLS_Result.csv"

# For comparison results
# The results are saved in "PLS_Result_comparison.csv"
source("./CLasso/master_comparison.R")

# ==========================================================================
# Replication of Table 2: REL Results
# ==========================================================================

# The bias and RMSE by Rmosek documented in Table 2 (left panel) are stored
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

