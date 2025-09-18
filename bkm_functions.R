library(Rcpp) 
library(igraph) 
library(irlba) 
library(RSpectra) 
library(poweRlaw) 

# Rearranges the rows and columns of a matrix so that its largest element is moved to the (1,1) position
#
# Arguments:
#   M : Numeric matrix
#
# Returns:
#   Numeric matrix with the largest element moved to the (1,1) position
#
maxswap <- function(M) {
  n <- nrow(M)
  k <- which.max(M)
  j <- ceiling(k/n)
  i <- k - n*(j-1)
  c(i, j)
  M[c(1, i),] <- M[c(i, 1),]
  M[,c(1, j)] <- M[,c(j, 1)]
  M
}

# Greedy algorithm to correct label permutation
# Iteratively permute the rows and columns of a confusion matrix to align predicted labels with reference labels
#
# Arguments:
#   M : Confusion matrix (numeric)
#
# Returns:
#   Permuted confusion matrix with improved alignment of labels
#
miscorrect <- function(M) {
  n <- nrow(M)
  for (i in 1:(n-1)) M[i:n, i:n] <- maxswap(M[i:n, i:n])
  M
}

# Calls projection-based clustering function
#
sourceCpp("kmeansalt.cpp")

# Generates a network from a DCBM
#
# Arguments:
#  n: Integer, number of nodes
#  k: Integer, number of communities
#  W: k x k numeric matrix, block probability matrix
#  t: Numeric vector of degree parameters
#  m: Integer vector of true communities of nodes
#
# Returns: 
#  n x n symmetric binary adjacency matrix
#
DCBM.fast2 <- function(n, k, W, t, m) {
  on.exit(gc())
  
  stor <- lapply(1:(n-1), function(i) {
    tmp <- rbinom(n-i, 1, pmin(W[m[i],m[(i+1):n]]*t[i]*t[(i+1):n], 1))
    i + which(tmp == 1)
  })
  
  size <- sapply(1:(n-1), function(i) length(stor[[i]]))
  vec <- list('i' = rep(1:(n-1), size), 'j' = unlist(stor))
  
  A <- sparseMatrix(i = c(vec$i, vec$j), j = c(vec$j, vec$i), x = 1, dims = c(n,n))
}

#Generates random graph A ~ Bernoulli(P) from a given probabilty matrix
#
# Arguments:
#   P : n x n numeric matrix of edge probabilities
#
# Returns: n x n symmetric binary adjacency matrix
#
rg_sample <- function(P){
  n <- nrow(P)
  A <- matrix(0, nrow = n, ncol = n)
  A[upper.tri(A, diag = FALSE)] <- rbinom(n*(n-1)/2, 1,
                                          pmin(P[upper.tri(P, diag = FALSE)], 1))
  A <- Matrix(A + t(A))
}

# Computes latent positions under SBM and DCBM
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  d: Integer, number of communities
#
# Returns: n x d numeric matrix of latent positions
#
Latent_positions <- function(A, d) {
  n <- ncol(A)
  e <- irlba(A, nu = d, nv = d)
  Uhat <- e$u
  Uhat
}

# Computes latent positions under PABM
#
# Arguments:
#  A: n x n symmetric binary matrix, adjacency matrix
#  k: Integer, number of communities
#
# Returns: n x k^2 numeric matrix of latent positions
#
Latent_positions_pabm <- function(A, k) {
  n <- ncol(A)
  
  e <- irlba(A, nu = k^2, nv = k^2)
  Uhat <- e$u
  sqrt(n)*Uhat
}

# Minimizes objective functions (Q_1, Q_2 or Q_3) given an adjacency matrix
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  maxiter: Integer, maximum number of iterations for clustering
#  nstart: Integer, number of initializations for clustering
#  method: Character, objective function type - "SBM"(for Q_1), "DCBM"(for Q_2), or "PABM"(for Q_3)
#
# Returns: A list
#  est_comm: Integer vector, estimated community assignments
#  obj: Numeric, minimized value of the objective function
#
projSC <- function(A, k, maxiter, nstart, method) {
  if(method == "SBM") {
    aa <- Latent_positions(A, k)
    out <- kmeans(aa, k, maxiter, nstart)
    out$error <- out$tot.withinss
  } else if(method == "DCBM") {
    aa <- Latent_positions(A, k)
    out <- kmeansAlt(aa, k, 1, maxiter, nstart)
  } else if(method == "PABM") {
    aa <- Latent_positions_pabm(A, k)
    out <- kmeansAltV2(aa, k, k, maxiter, nstart, batch_size = 100)
  }
  
  list(est_comm = c(out$cluster), obj = out$error)
}

# Simulates a network from a given model (SBM, DCBM or PABM) and minimizes an objective function
#
# Arguments:
#  params: List, parameters of the model
#  maxiter: Integer, maximum number of iterations for clustering
#  nstart: Integer, number of initializations for clustering
#  data_str: Character, model type - "SBM", "DCBM" or "PABM"
#  method: Character, objective function type - "SBM"(for Q_1), "DCBM"(for Q_2), or "PABM"(for Q_3)
#
# Returns: A list
#  error: Numeric, community detection error
#  obj: Numeric, minimized value of the objective function
#  est_comm: Integer vector, estimated community assignments
#  adj_mat: n x n adjacency matrix of the generated network
#
sim <- function(params, maxiter, nstart, data_str, method) {
  if(data_str == "SBM") {
    wmat <- params[[1]]
    comm <- params[[2]]
    k <- length(unique(comm))
    n <- length(comm)
    comm_size <- unname(table(comm))
    A <- as_adjacency_matrix(sample_sbm(n, wmat, comm_size))
  } else if(data_str == "DCBM") {
    wmat <- params[[1]]
    t <- params[[2]]
    comm <- params[[3]]
    
    k <- length(unique(comm))
    n <- length(comm)
    
    A <- DCBM.fast2(n, k, wmat, t, comm)
  } else if(data_str == "PABM") {
    P <- params[[1]]
    comm <- params[[2]]
    k <- length(unique(comm))
    n <- length(comm)
    A <- rg_sample(P)
  }
  
  out <- projSC(A, k, maxiter, nstart, method)
  conf <- miscorrect(table(factor(comm, levels = 1:k), factor(out$est_comm, levels = 1:k)))
  
  list(error = 1 - sum(diag(conf))/n, obj = out$obj, est_comm = out$est_comm, adj_mat = A)
}

# Performs test H_0: A ~ SBM vs. H_1: A ~ DCBM using Q_1
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  maxiter: Integer, maximum number of iterations for kmeans 
#  nstart: Integer, number of initializations for kmeans
#  R: Integer, number of bootstrap replicates
#
# Returns: A list
#  pvalue: Numeric, p-value of the test
#  true_obj: Numeric, minimized value of Q_1 for the given network (test statistic)
#  gen_obj: Numeric, minimized values of Q_1 for the bootstrapped networks
#
sbm_vs_dcbm <- function(A, k, maxiter, nstart, R) {
  sbmfit <- projSC(A, k, maxiter, nstart, "SBM")
  
  wmathat <- matrix(0, ncol = k, nrow = k)
  for(i in 1:k) {
    for (j in 1:k) {
      if(i == j) {
        ni <- sum(sbmfit$est_comm == i)
        wmathat[i, i] <- sum(A[sbmfit$est_comm == i, sbmfit$est_comm == i])/(ni*(ni - 1))
      } else wmathat[i, j] <- mean(A[sbmfit$est_comm == i, sbmfit$est_comm == j])
    }
  }
  wmathat <- (wmathat + t(wmathat))/2
  trep <- replicate(R, sim(list(wmathat, sbmfit$est_comm), maxiter, nstart, "SBM", "SBM")$obj)
  list(pvalue = mean(trep > sbmfit$obj), true_obj = sbmfit$obj, gen_obj = trep)
}

# Performs test H_0: A ~ DCBM vs. H_1: A ~ PABM using Q_2
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  maxiter: Integer, maximum number of iterations for clustering 
#  nstart: Integer, number of initializations for clustering
#  R: Integer, number of bootstrap replicates
#
# Returns: A list
#  pvalue: Numeric, p-value of the test
#  true_obj: Numeric, minimized value of Q_2 for the given network (test statistic)
#  gen_obj: Numeric, minimized value of Q_2 for the bootstrapped networks
#
dcbm_vs_pabm <- function(A, k, maxiter, nstart, R) {
  dcbmfit <- projSC(A, k, maxiter, nstart, "DCBM")
  
  wmathat <- matrix(0, ncol = k, nrow = k)
  for(i in 1:k) {
    for (j in 1:k) {
      wmathat[i, j] <- sum(A[dcbmfit$est_comm == i, dcbmfit$est_comm == j])
    }
  }
  wmathat <- (wmathat + t(wmathat))/2
  that <- colSums(A)/colSums(wmathat[,dcbmfit$est_comm])
  
  trep <- replicate(R, sim(list(wmathat, that, dcbmfit$est_comm), maxiter, nstart, "DCBM", "DCBM")$obj)
  list(pvalue = mean(trep > dcbmfit$obj), true_obj = dcbmfit$obj, gen_obj = trep)
}

# dcbm_vs_pabm_2 <- function(A, k, maxiter, nstart, R) {
#   dcbmfit <- projSC(A, k, maxiter, nstart, "DCBM")
#   pabmfit <- projSC(A, k, maxiter, nstart, "PABM")
#   
#   wmathat <- matrix(0, ncol = k, nrow = k)
#   for(i in 1:k) {
#     for (j in 1:k) {
#       wmathat[i, j] <- sum(A[dcbmfit$est_comm == i, dcbmfit$est_comm == j])
#     }
#   }
#   that <- colSums(A)/colSums(wmathat[,dcbmfit$est_comm])
#   
#   trep1 <- replicate(R, sim(list(wmathat, that, dcbmfit$est_comm), maxiter, nstart, "DCBM", "DCBM")$obj)
#   trep2 <- replicate(R, sim(list(wmathat, that, dcbmfit$est_comm), maxiter, nstart, "DCBM", "PABM")$obj)
#   
#   list(pvalue = mean(trep1-trep2 > dcbmfit$obj - pabmfit$obj), true_obj = dcbmfit$obj - pabmfit$obj, gen_obj = trep1-trep2)
# }
