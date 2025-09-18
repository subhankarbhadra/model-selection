library(pbapply)

source("bkm_functions.R")

# Simulates a network from a 3-community PABM and performs model selection (DCBM vs. PABM) using Q_2
#
# Arguments:
#  density: Numeric, network density
#
# Returns: p-value of the test
#
expabm2 <- function(density) {
  n <- 900
  k <- 3
  L <- cbind(c(rbeta(n/k, 2, 1), rbeta(n/k, 1, 2), rbeta(n/k, 1, 2)), 
             c(rbeta(n/k, 1, 2), rbeta(n/k, 2, 1), rbeta(n/k, 1, 2)),
             c(rbeta(n/k, 1, 2), rbeta(n/k, 1, 2), rbeta(n/k, 2, 1)))
  P <- matrix(0, n, n)
  comm <- rep(1:3, each = n/k)
  P[comm == 1, comm == 1] <- tcrossprod(L[comm == 1, 1, drop = F])
  P[comm == 2, comm == 2] <- tcrossprod(L[comm == 2, 2, drop = F])
  P[comm == 3, comm == 3] <- tcrossprod(L[comm == 3, 3, drop = F])
  P[comm == 1, comm == 2] <- L[comm == 1, 2, drop = F] %*% t(L[comm == 2, 1, drop = F])
  P[comm == 2, comm == 1] <- t(P[comm == 1, comm == 2])
  P[comm == 1, comm == 3] <- L[comm == 1, 3, drop = F] %*% t(L[comm == 3, 1, drop = F])
  P[comm == 3, comm == 1] <- t(P[comm == 1, comm == 3])
  P[comm == 2, comm == 3] <- L[comm == 2, 3, drop = F] %*% t(L[comm == 3, 2, drop = F])
  P[comm == 3, comm == 2] <- t(P[comm == 2, comm == 3])
  diag(P) <- 0
  
  P <- (density/mean(P))*P
  A <- rg_sample(P)
  
  maxiter <- 100
  nstart <- 25
  
  dcbm_vs_pabm(A, k, maxiter, nstart, 200)$pvalue
}

#set.seed(12345)
ff <- pbreplicate(100, expabm2(0.05))
mean(ff)
mean(ff < 0.05)
