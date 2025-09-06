library(Rcpp) 
library(igraph) 
library(irlba) 
library(RSpectra) 
library(poweRlaw) 

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

miscorrect <- function(M) {
  n <- nrow(M)
  for (i in 1:(n-1)) M[i:n, i:n] <- maxswap(M[i:n, i:n])
  M
}

#Call projection-based kmeans function
sourceCpp("kmeansalt.cpp")

#Generate networks from DCBM
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

#Generate random graph from a probabilty matrix
rg_sample <- function(P){
  n <- nrow(P)
  A <- matrix(0, nrow = n, ncol = n)
  A[upper.tri(A, diag = FALSE)] <- rbinom(n*(n-1)/2, 1,
                                          pmin(P[upper.tri(P, diag = FALSE)], 1))
  A <- Matrix(A + t(A))
}

#Adjacency spectral embedding
Latent_positions <- function(A, d) {
  n <- ncol(A)
  e <- irlba(A, nu = d, nv = d)
  Uhat <- e$u
  Uhat
}

Latent_positions_pabm <- function(A, k) {
  n <- ncol(A)
  
  e <- irlba(A, nu = k^2, nv = k^2)
  Uhat <- e$u
  sqrt(n)*Uhat
}

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
    #dt <- params[[4]]
    
    #M <- wmat[comm, comm]
    #diag(M) <- 0
    #wmat <- wmat*dt/mean(sapply(1:length(t), function(i) mean(t[i]*M[i,]*t)))
    
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

dcbm_vs_pabm_2 <- function(A, k, maxiter, nstart, R) {
  dcbmfit <- projSC(A, k, maxiter, nstart, "DCBM")
  pabmfit <- projSC(A, k, maxiter, nstart, "PABM")
  
  wmathat <- matrix(0, ncol = k, nrow = k)
  for(i in 1:k) {
    for (j in 1:k) {
      wmathat[i, j] <- sum(A[dcbmfit$est_comm == i, dcbmfit$est_comm == j])
    }
  }
  that <- colSums(A)/colSums(wmathat[,dcbmfit$est_comm])
  
  trep1 <- replicate(R, sim(list(wmathat, that, dcbmfit$est_comm), maxiter, nstart, "DCBM", "DCBM")$obj)
  trep2 <- replicate(R, sim(list(wmathat, that, dcbmfit$est_comm), maxiter, nstart, "DCBM", "PABM")$obj)
  
  list(pvalue = mean(trep1-trep2 > dcbmfit$obj - pabmfit$obj), true_obj = dcbmfit$obj - pabmfit$obj, gen_obj = trep1-trep2)
}
