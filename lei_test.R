library(RSpectra)
library(RMTstat)

lei_eigmax <- function(A, k, maxiter, nstart, alpha = 0.05) {
  sbmfit <- projSC(A, k, maxiter, nstart, "SBM")
  n <- nrow(A)
  
  wmathat <- matrix(0, ncol = k, nrow = k)
  for(i in 1:k) {
    for (j in 1:k) {
      if(i == j) {
        ni <- sum(sbmfit$est_comm == i)
        wmathat[i, i] <- sum(A[sbmfit$est_comm == i, sbmfit$est_comm == i])/(ni*(ni - 1))
      } else wmathat[i, j] <- mean(A[sbmfit$est_comm == i, sbmfit$est_comm == j])
    }
  }
  Phat <- wmathat[sbmfit$est_comm, sbmfit$est_comm]
  Atilde <- matrix((A - Phat) / sqrt((n-1)*(Phat*(1 - Phat) + 1e-16)), nrow = n)
  diag(Atilde) <- 0
  
  l1l2 <- eigs_sym(Atilde, k = 2, which = "BE")$values
  
  wob <- (n^(2/3))*(max(abs(l1l2)) - 2) > qtw(1 - alpha/2)
  
  twmean <- -1.206534
  twsd <- sqrt(1.607781)
  l1l2boots <- replicate(50, sim_lei(wmathat, sbmfit$est_comm, Phat))
  l1mean <- mean(l1l2boots[1,]); l1sd <- sd(l1l2boots[1,])
  l2mean <- mean(l1l2boots[2,]); l2sd <- sd(l1l2boots[2,])
  wb <- twmean + twsd * max(c((l1l2[1] - l1mean)/l1sd, -(l1l2[2] - l2mean)/l2sd)) > qtw(1 - alpha/2)
  
  list(wob = wob, wb = wb)
}

sim_lei <- function(wmathat, est_comm, Phat) {
  k <- length(unique(est_comm))
  n <- length(est_comm)
  comm_size <- unname(table(est_comm))
  Aboot <- as_adjacency_matrix(sample_sbm(n, wmathat, comm_size))
  
  Aboottilde <- matrix((Aboot - Phat) / sqrt((n-1)*(Phat*(1 - Phat) + 1e-16)), nrow = n)
  diag(Aboottilde) <- 0
  
  l1l2 <- eigs_sym(Aboottilde, k = 2, which = "BE")$values
}