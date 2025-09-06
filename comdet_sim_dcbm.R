library(Rcpp)
library(igraph)
library(irlba)
library(RSpectra)
library(mclust)
library(pbapply)

source("bkm_functions.R")

#Regularized spectral clustering using laplacian
fast.clustering_DCBM_laplacian <- function(A, k, maxiter, nstart) {
  n <- ncol(A)
  degree <- colSums(A)
  
  D <- sparseMatrix(i = 1:n, j = 1:n, x = 1/sqrt(degree + mean(degree)))
  M <- tcrossprod(crossprod(D, A), D)
  
  e <- irlba(M, nu = k, nv = k)
  S <- e$u
  S_reg <- t(sapply(1:nrow(S), function(i) S[i,]/sqrt(sum(S[i,]^2))))
  c <- kmeans(S_reg, k, maxiter, nstart)$cluster
  
  What <- matrix(0, ncol = k, nrow = k)
  for(i in 1:k) {
    for (j in 1:k) {
      What[i, j] <- mean(A[c == i, c == j])
    }
  }
  list(cluster = c, What = What)
}

sim_dcbm_comp <- function(wmat, t, comm, dt, maxiter, nstart) {
  k <- nrow(wmat)
  n <- length(comm)
  
  M <- wmat[comm, comm]
  diag(M) <- 0
  wmat <- wmat*dt/mean(sapply(1:length(t), function(i) mean(t[i]*M[i,]*t)))
  A <- DCBM.fast2(n, k, wmat, t, comm)
  degree <- colSums(A)
  
  #neglecting zero degree nodes
  sub <- which(degree != 0)
  Asub <- A[sub, sub]
  
  t1 <- system.time(outsub <- fast.clustering_DCBM_laplacian(Asub, k, maxiter, nstart))[3]
  out <- rep(0, n)
  out[sub] <- outsub$cluster
  out[-sub] <- sample.int(k, n-length(sub), replace = T)
  
  conf <- miscorrect(table(factor(comm, levels = 1:k), factor(out, levels = 1:k)))
  
  t2 <- system.time(outQ2 <- projSC(A, k, maxiter, nstart, "DCBM"))[3]
  confQ2 <- miscorrect(table(factor(comm, levels = 1:k), factor(outQ2$est_comm, levels = 1:k)))
  
  c(1 - sum(diag(conf))/n, 1 - sum(diag(confQ2))/n, t1, t2) 
}


#DCBM simulation
wmat <- matrix(c(4, 2, 1, 2, 4, 1, 1, 1, 4), 3, 3)
k <- 3
delta <- 0.05
maxiter <- 100
nstart <- 25

n <- 2000
comm <- rep(1:k, c(n/4, n/4, n/2))

#set.seed(12345)
v <- pbreplicate(100, sim_dcbm_comp(wmat, rbeta(n, 1, 5), comm, delta, maxiter, nstart))
# average mislabelling error: RSC-L, Q2
mean(v[1,]); mean(v[2,])
# s.d. mislabelling error: RSC-L, Q2
sd(v[1,]); sd(v[2,])
# average runtime: RSC, Q2
mean(v[3,]); mean(v[4,])

