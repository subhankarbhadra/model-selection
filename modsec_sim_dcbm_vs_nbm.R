library(pbapply)

source("bkm_functions.R")
source("nbm.R")

# Simulates a network from NBM and performs model selection (DCBM vs. PABM) using Q_2
#
# Arguments:
#  n: Integer, number of nodes
#  k: Integer, number of communities
#  l: Integer, number of meta-communities
#  w: Numeric, homophily factor
#
# Returns: p-value of the test
#
exnbm <- function(n, k, l, w) {
  P <- Data_HBM(n, k, l, w)
  A <- rg_sample(P)
  
  maxiter <- 100
  nstart <- 25
  
  dcbm_vs_pabm(A, k, maxiter, nstart, 200)$pvalue
}

n <- 900
k <- 6

# set.seed(12345)
ff <- pbreplicate(100, exnbm(n, k, 2, 0.6)) 
mean(ff)
mean(ff < 0.05)
