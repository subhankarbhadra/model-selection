library(pbapply)

source("bkm_functions.R")

exsbm <- function(density) {
  n <- 900; k <- 3
  wmat <- (-1/3)*diag(k) + (2/3)*outer(rep(1, k), rep(1, k))
  comm <- rep(1:k, each = n/k)
  M <- wmat[comm, comm]
  diag(M) <- 0
  wmat <- wmat*(density)/mean(M)
  
  maxiter <- 100
  nstart <- 25
  
  A <- as_adjacency_matrix(sample_sbm(n, wmat, rep(n/k, k)))
  
  #ecvmod <- ECV.block.modified(A, k)
  #lei_out <- lei_eigmax(A, k, maxiter, nstart, 0.05)
  
  p1 <- sbm_vs_dcbm(A, k, maxiter, nstart, 200)$pvalue
  if(p1 >= 0.05) {
    return("SBM")
  } else {
    p2 <- dcbm_vs_pabm(A, k, maxiter, nstart, 200)$pvalue
    if(p2 >= 0.05) return("DCBM") else return("PABM")
  }
}

#set.seed(12345)
gg <- pbreplicate(100, exsbm(.05))
table(gg)

exdcbm <- function(density) {
  n <- 900; k <- 3
  wmat <- (-1/3)*diag(k) + (2/3)*outer(rep(1, k), rep(1, k))
  t <- rplcon(n, xmin = 1, alpha = 5)
  comm <- rep(1:k, each = n/k)
  M <- wmat[comm, comm]
  diag(M) <- 0
  wmat <- wmat*(density)/mean(sapply(1:length(t), function(i) mean(t[i]*M[i,]*t)))
  
  maxiter <- 100
  nstart <- 25
  
  A <- DCBM.fast2(n, k, wmat, t, comm)
  
  p1 <- sbm_vs_dcbm(A, k, maxiter, nstart, 200)$pvalue
  if(p1 >= 0.05) {
    return("SBM")
  } else {
    p2 <- dcbm_vs_pabm(A, k, maxiter, nstart, 200)$pvalue
    if(p2 >= 0.05) return("DCBM") else return("PABM")
  }
}

#set.seed(12345)
gg <- pbreplicate(100, exdcbm(.05))
table(gg)

expabm3 <- function(density) {
  n <- 900
  k <- 3
  L <- cbind(c(rbeta(n/k, 1, 2), rbeta(n/k, 2, 1), rbeta(n/k, 2, 1)), 
             c(rbeta(n/k, 2, 1), rbeta(n/k, 1, 2), rbeta(n/k, 2, 1)),
             c(rbeta(n/k, 2, 1), rbeta(n/k, 2, 1), rbeta(n/k, 1, 2)))
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
  
  p1 <- sbm_vs_dcbm(A, k, maxiter, nstart, 200)$pvalue
  if(p1 >= 0.05) {
    return("SBM")
  } else {
    p2 <- dcbm_vs_pabm(A, k, maxiter, nstart, 200)$pvalue
    if(p2 >= 0.05) return("DCBM") else return("PABM")
  }
}

#set.seed(12345)
ff <- pbreplicate(100, expabm3(0.05))
table(ff)
