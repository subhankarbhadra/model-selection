library(Rcpp)
library(igraph)
library(irlba)
library(RSpectra)
library(mclust)
library(pbapply)

source("bkm_functions.R")

# Performs spectral clustering using laplacian (SC-L)
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  maxiter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  What: k x k numeric matrix, estimated block probability matrix
#
fast.clustering_laplacian <- function(A, k, maxiter, nstart) {
  n <- ncol(A)
  degree <- colSums(A)
  
  D <- sparseMatrix(i = 1:n, j = 1:n, x = 1/sqrt(degree))
  M <- tcrossprod(crossprod(D, A), D)
  
  e <- irlba(M, nu = k, nv = k)
  S <- e$u
  c <- kmeans(S, k, maxiter, nstart)$cluster
  What <- matrix(0, ncol = k, nrow = k)
  for(i in 1:k) {
    for (j in 1:k) {
      What[i, j] <- mean(A[c == i, c == j])
    }
  }
  list(cluster = c, What = What)
}

# Simulates a network from SBM and performs community detection using SC-L and Q_1
#
# Arguments:
#  wmat: symmetric block probability matrix
#  comm: Integer vector of true communities
#  maxiter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#
# Returns: community detection errors from SC-L and Q_1
#
sim_sbm_comp <- function(wmat, comm, maxiter, nstart) {
  n <- length(comm)
  k <- nrow(wmat)
  comm_size <- unname(table(comm))
  A <- as_adjacency_matrix(sample_sbm(n, wmat, comm_size))
  degree <- colSums(A)
  
  #neglecting zero degree nodes
  sub <- which(degree != 0)
  Asub <- A[sub, sub]
  
  outsub <- fast.clustering_laplacian(Asub, k, maxiter, nstart)
  out <- rep(0, n)
  out[sub] <- outsub$cluster
  out[-sub] <- sample.int(k, n-length(sub), replace = T)
  
  conf <- miscorrect(table(factor(comm, levels = 1:k), factor(out, levels = 1:k)))
  
  outQ1 <- projSC(A, k, maxiter, nstart, "SBM")
  confQ1 <- miscorrect(table(factor(comm, levels = 1:k), factor(outQ1$est_comm, levels = 1:k)))
  
  c(1 - sum(diag(conf))/n, 1 - sum(diag(confQ1))/n) 
}

#SBM simulation
maxiter <- 100
nstart <- 25
wmat <- matrix(c(4, 2, 1, 2, 4, 1, 1, 1, 4), 3, 3)
k <- nrow(wmat)
delta <- 0.05

n <- 2000
comm <- rep(1:k, c(n/4, n/4, n/2))
P <- wmat[comm, comm]
diag(P) <- 0
wmat <- (delta/mean(P))*wmat

#set.seed(12345)
v <- pbreplicate(100, sim_sbm_comp(wmat, comm, maxiter, nstart))
# average community detection error: SC-L, Q_1
mean(v[1,]); mean(v[2,])
# s.d. community detection error: SC-L, Q_1
sd(v[1,]); sd(v[2,])

