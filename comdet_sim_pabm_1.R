library(Rcpp)
library(igraph)
library(irlba)
library(RSpectra)
library(mclust)
library(pbapply)

source("bkm_functions.R")
source("koo_functions.R")

sim_pabm_comp_1 <- function(n, maxiter, nstart) {
  k <- 2
  #assortative, balanced
  L <- cbind(c(rbeta(n/k, 2, 1), rbeta(n/k, 1, 2)),
             c(rbeta(n/k, 1, 2), rbeta(n/k, 2, 1)))
  
  P <- matrix(0, n, n)
  comm <- rep(1:2, c(n/k, n/k))
  
  randomize <- sample(1:n, n) 
  L <- L[randomize, ]
  comm <- comm[randomize]
  
  P[comm == 1, comm == 1] <- tcrossprod(L[comm == 1, 1, drop = F])
  P[comm == 2, comm == 2] <- tcrossprod(L[comm == 2, 2, drop = F])
  P[comm == 1, comm == 2] <- L[comm == 1, 2, drop = F] %*% t(L[comm == 2, 1, drop = F])
  P[comm == 2, comm == 1] <- t(P[comm == 1, comm == 2])
  
  A <- rg_sample(P)
  
  t1 <- system.time(koo_est <- cluster.pabm(A, k, normalize = F))[3]
  conf <- miscorrect(table(factor(comm, levels = 1:k), factor(koo_est, levels = 1:k)))
  
  t2 <- system.time(outQ3 <- projSC(A, k, maxiter, nstart, "PABM"))[3]
  confQ3 <- miscorrect(table(factor(comm, levels = 1:k), factor(outQ3$est_comm, levels = 1:k)))
  
  c(1 - sum(diag(conf))/n, 1 - sum(diag(confQ3))/n, t1, t2)
  }


maxiter <- 100
nstart <- 25

set.seed(12345)
v <- pbreplicate(100, sim_pabm_comp_1(900, maxiter, nstart))
# average mislabelling error: OSC, Q3
mean(v[1,]); mean(v[2,])
# s.d. mislabelling error: OSC, Q3
sd(v[1,]); sd(v[2,])
# average runtime: OSC, Q3
mean(v[3,]); mean(v[4,])
