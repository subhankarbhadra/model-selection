library(pbapply)

source("bkm_functions.R")
source("ECV_modified.R")
source("lei_test.R")

g <- function(n, k, beta, avgdeg) {
  wmat <- (1-beta)*diag(k) + beta*outer(rep(1, k), rep(1, k))
  comm <- rep(1:k, each = n/k)
  M <- wmat[comm, comm]
  diag(M) <- 0
  wmat <- wmat*(avgdeg/n)/mean(M)
  
  maxiter <- 100
  nstart <- 25
  
  A <- as_adj(sample_sbm(n, wmat, rep(n/k, k)))
  
  ecvmod <- ECV.block.modified(A, k)
  lei_out <- lei_eigmax(A, k, maxiter, nstart, 0.05)
  
  c((sbm_vs_dcbm(A, k, maxiter, nstart, 200)$pvalue < 0.05),
    lei_out$wob, lei_out$wb,
    which.min(c(ecvmod$l2, ecvmod$dc.l2)) - 1,
    which.min(c(ecvmod$dev, ecvmod$dc.dev)) - 1)
}

#set.seed(12345)
n <- 600; k <- 5; beta <- 0.2; avgdeg <- 40
gg <- pbreplicate(100, g(n, k, beta, avgdeg))
rowMeans(gg)
