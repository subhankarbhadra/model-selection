library(pbapply)

source("bkm_functions.R")

# Simulates a network from DCBM and performs model selection (DCBM vs. PABM) using Q_2
#
# Arguments:
#  n: Integer, number of nodes
#  k: Integer, number of communities
#  beta: Numeric, homophily factor
#  avgdeg: Integer, average degree
#
# Returns: p-value of the test
#
f <- function(n, k, beta, avgdeg) {
  wmat <- (1-beta)*diag(k) + beta*outer(rep(1, k), rep(1, k))
  t <- rplcon(n, xmin = 1, alpha = 5)
  comm <- rep(1:k, each = n/k)
  M <- wmat[comm, comm]
  diag(M) <- 0
  wmat <- wmat*(avgdeg/n)/mean(sapply(1:length(t), function(i) mean(t[i]*M[i,]*t)))
  
  maxiter <- 100
  nstart <- 25
  
  A <- DCBM.fast2(n, k, wmat, t, comm)
  
  dcbm_vs_pabm(A, k, maxiter, nstart, 200)$pvalue
}

#set.seed(12345)
n <- 600; k <- 5; beta <- 0.2; avgdeg <- 40
ff <- pbreplicate(100, f(n, k, beta, avgdeg))
mean(ff)
mean(ff < 0.05)
