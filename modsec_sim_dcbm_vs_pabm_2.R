library(pbapply)

source("bkm_functions.R")

# Simulates a network from a 2-community PABM and performs model selection (DCBM vs. PABM) using Q_2
#
# Arguments:
#  density: Numeric, network density
#
# Returns: p-value of the test
#
expabm1 <- function(density) {
  n <- 900
  k <- 2
  L <- cbind(c(rbeta(n/k, 2, 1), rbeta(n/k, 1, 2)),
             c(rbeta(n/k, 1, 2), rbeta(n/k, 2, 1)))
  P <- matrix(0, n, n)
  comm <- rep(1:2, each = n/k)
  P[comm == 1, comm == 1] <- tcrossprod(L[comm == 1, 1, drop = F])
  P[comm == 2, comm == 2] <- tcrossprod(L[comm == 2, 2, drop = F])
  P[comm == 1, comm == 2] <- L[comm == 1, 2, drop = F] %*% t(L[comm == 2, 1, drop = F])
  P[comm == 2, comm == 1] <- t(P[comm == 1, comm == 2])
  diag(P) <- 0
  
  P <- (density/mean(P))*P
  A <- rg_sample(P)
  
  maxiter <- 100
  nstart <- 25
  
  dcbm_vs_pabm(A, k, maxiter, nstart, 200)$pvalue
}

#set.seed(12345)
ff <- pbreplicate(100, expabm1(0.05))
mean(ff)
mean(ff < 0.05)
