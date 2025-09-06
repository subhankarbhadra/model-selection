library(pbapply)

source("bkm_functions.R")
source("ECV_modified.R")
source("lei_test.R")

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
  
  ecvmod <- ECV.block.modified(A, k)
  lei_out <- lei_eigmax(A, k, maxiter, nstart, 0.05)
  
  c((sbm_vs_dcbm(A, k, maxiter, nstart, 200)$pvalue < 0.05),
    lei_out$wob, lei_out$wb,
    which.min(c(ecvmod$l2, ecvmod$dc.l2)) - 1,
    which.min(c(ecvmod$dev, ecvmod$dc.dev)) - 1)
}

#set.seed(12345)
n <- 600; k <- 5; beta <- 0.2; avgdeg <- 40
ff <- pbreplicate(100, f(n, k, beta, avgdeg))
rowMeans(ff)
