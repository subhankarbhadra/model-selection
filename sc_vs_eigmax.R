library(pbapply)
library(RSpectra)

source("bkm_functions.R")

# Computes runtimes for loss function minimization (Q_1) and test statistic computation for EigMax for SBM networks
# 
time_comp <- function(n) {
  k <- 5; beta <- 0.2; avgdeg <- n*0.1
  wmat <- (1-beta)*diag(k) + beta*outer(rep(1, k), rep(1, k))
  comm <- rep(1:k, each = n/k)
  M <- wmat[comm, comm]
  diag(M) <- 0
  wmat <- wmat*(avgdeg/n)/mean(M)
  
  maxiter <- 100
  nstart <- 25
  
  A <- as_adj(sample_sbm(n, wmat, rep(n/k, k)))
  
  t1 <- system.time(sbmfit <- projSC(A, k, maxiter, nstart, "SBM"))[3]
  
  t2 <- system.time({
    wmathat <- matrix(0, ncol = k, nrow = k);
    for(i in 1:k) {
      for (j in 1:k) {
        if(i == j) {
          ni <- sum(sbmfit$est_comm == i)
          wmathat[i, i] <- sum(A[sbmfit$est_comm == i, sbmfit$est_comm == i])/(ni*(ni - 1))
        } else wmathat[i, j] <- mean(A[sbmfit$est_comm == i, sbmfit$est_comm == j])
      }
    };
    Phat <- wmathat[sbmfit$est_comm, sbmfit$est_comm];
    Atilde <- matrix((A - Phat) / sqrt((n-1)*(Phat*(1 - Phat) + 1e-16)), nrow = n);
    diag(Atilde) <- 0;
    l1l2 <- eigs_sym(Atilde, k = 2, which = "BE")$values;
    #l1l2 <- irlba(Atilde, nu = 1, nv = 1)
    wob <- (n^(2/3))*(max(abs(l1l2)) - 2) > qtw(1 - 0.05/2)})[3]
  c(t1, t2)
}

for (m in seq(1000, 10000, 1000)) {
  print(rowMeans(pbreplicate(100, time_comp(m))))
}

