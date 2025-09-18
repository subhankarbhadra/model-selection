# Codes for generating networks from NBM (obtained from Noroozi et. al. (2021))

Data_HBM <- function(n, K, rr, q) {
  n_K <- n / K
  K_r <- K / rr
  n_r <- n / rr
  Theta <- matrix(runif(n * K), nrow = n)
  
  theta_1 <- runif(n)
  theta_2 <- runif(n)
  theta_3 <- runif(n)
  theta_4 <- runif(n)
  theta_5 <- runif(n)
  theta_6 <- runif(n)
  theta_7 <- runif(n)
  theta_8 <- runif(n)
  
  if (rr == 1) {
    # Creating r=1 distinct columns for Theta (begins)
    Theta[, 1:6] <- theta_1
    # Creating r=1 distinct columns for Theta (ends)
  } else if (rr == 2) {
    # Creating r=2 distinct columns for Theta (begins)
    for (r in 1:K) {
      theta_1_temp <- theta_1[((r - 1) * n_K + 1):(r * n_K)]
      theta_1[((r - 1) * n_K + 1):(r * n_K)] <- sort(theta_1[((r - 1) * n_K + 1):(r * n_K)])
      theta_2[((r - 1) * n_K + 1):(r * n_K)] <- sort(theta_1_temp, decreasing = TRUE)
    }
    Theta[, 1:6] <- cbind(theta_1, theta_1, theta_1, theta_2, theta_2, theta_2)
    # Creating r=2 distinct columns for Theta (ends)
  } else if (rr == 3) {
    # Test for r=3 Begins
    for (r in 1:K) {
      theta_1_temp <- theta_1[((r - 1) * n_K + 1):(r * n_K)]
      aa <- sort(theta_1[((r - 1) * n_K + 1):(r * n_K)])
      bb <- aa[1:(n_K / 3)]
      cc <- aa[(n_K / 3 + 1):(2 * n_K / 3)]
      dd <- aa[(2 * n_K / 3 + 1):n_K]
      bb <- sort(bb, decreasing = TRUE)
      cc <- sort(cc, decreasing = TRUE)
      dd <- sort(dd, decreasing = TRUE)
      ee <- c(dd, bb, cc)
      ff <- c(cc, dd, bb)
      theta_1[((r - 1) * n_K + 1):(r * n_K)] <- aa
      theta_2[((r - 1) * n_K + 1):(r * n_K)] <- ee
      theta_3[((r - 1) * n_K + 1):(r * n_K)] <- ff
    }
    Theta[, 1:6] <- cbind(theta_1, theta_1, theta_2, theta_2, theta_3, theta_3)
    # Test for r=3 Ends
  } else if (rr == 4) {
    # Test for r=4 Begins
    for (r in 1:K) {
      theta_1_temp <- theta_1[((r - 1) * n_K + 1):(r * n_K)]
      aa <- sort(theta_1[((r - 1) * n_K + 1):(r * n_K)])
      bb <- aa[1:(n_K / 4)]
      cc <- aa[(n_K / 4 + 1):(2 * n_K / 4)]
      dd <- aa[(2 * n_K / 4 + 1):(3 * n_K / 4)]
      ee <- aa[(3 * n_K / 4 + 1):n_K]
      
      bb <- sort(bb, decreasing = TRUE)
      cc <- sort(cc, decreasing = TRUE)
      dd <- sort(dd, decreasing = TRUE)
      ee <- sort(ee, decreasing = TRUE)
      
      ff <- c(ee, bb, cc, dd)
      gg <- c(dd, ee, bb, cc)
      hh <- c(cc, dd, ee, bb)
      
      theta_1[((r - 1) * n_K + 1):(r * n_K)] <- aa
      theta_2[((r - 1) * n_K + 1):(r * n_K)] <- ff
      theta_3[((r - 1) * n_K + 1):(r * n_K)] <- gg
      theta_4[((r - 1) * n_K + 1):(r * n_K)] <- hh
    }
    Theta[, 1:8] <- cbind(theta_1, theta_1, theta_2, theta_2, theta_3, theta_3, theta_4, theta_4)
    # Test for r=4 Ends
  } else if (rr == 6) {
    # Creating r=6 distinct columns for Theta (begins)
    Theta[, 1:6] <- cbind(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6)
    # Creating r=6 distinct columns for Theta (ends)
  }
  
  Epsilon <- 0.01
  B <- matrix(runif(K * K), nrow = K, ncol = K)
  BB <- matrix(runif(K * K), nrow = K, ncol = K)
  one_vector <- rep(1, n_K)
  omega_1 <- 0.001
  
  
  for (r in 1:K) {
    for (s in 1:K) {
      Theta[((r - 1) * n_K + 1):(r * n_K), s] <- n_K * Theta[((r - 1) * n_K + 1):(r * n_K), s] / sum(Theta[((r - 1) * n_K + 1):(r * n_K), s])
    }
  }
  
  Theta_max <- matrix(0, nrow = K, ncol = K)
  
  for (r in 1:K) {
    for (s in 1:K) {
      Theta_max[r, s] <- max(max(Theta[((r - 1) * n_K + 1):(r * n_K), s], Theta[((s - 1) * n_K + 1):(s * n_K), r]))
    }
  }
  
  BB <- 0.65 * B + 0.35
  
  for (i in 1:K) {
    m <- max(BB[, i])
    if (BB[i, i] < m) {
      BB[i, i] <- m + (1 - m) * runif(1)
    }
  }
  
  B <- 1/Theta_max^2
  
  for (i in 1:K) {
    for (j in 1:K) {
      if (i != j || B[i, j] == 1) {
        B[i, j] <- BB[i, j] * B[i, j]
      }
    }
  }
  
  # Original B begins
  if (min(Theta_max) == 1) {
    B <- B * (q^2)
    for (r in 1:K) {
      B[r, r] <- (1 / (q^1)) * B[r, r]
    }
  } else {
    B <- B * (q^1)
    for (r in 1:K) {
      B[r, r] <- (1 / (q^1)) * B[r, r]
    }
  }
  # Original B ends
  
  P <- matrix(0, n, n)
  for (r in 1:K) {
    for (s in 1:K) {
      P[((r - 1) * n_K + 1):(r * n_K), ((s - 1) * n_K + 1):(s * n_K)] <- B[r, s] * Theta[((r - 1) * n_K + 1):(r * n_K), s] %*% t(Theta[((s - 1) * n_K + 1):(s * n_K), r])
    }
  }
  return(P)
}  

Latent_positions <- function(A, d) {
  n <- ncol(A)
  e <- irlba(A, nu = d, nv = d)
  Uhat <- e$u
  Uhat %*% diag(sqrt(e$d))
}

#I <- diag(n)
  #Perm_index_1 <- sample(n)
  #Perm_1 <- matrix(0, nrow = n, ncol = n)
  #for (i in 1:n) {
  #  Perm_1[i, ] <- I[Perm_index_1[i], ]
  #}
  #Z_1 <- matrix(0, nrow = n, ncol = K)
  #for (j in 1:K) {
  #  for (i in 1:n_K) {
  #    Z_1[Perm_index_1[i + (j - 1) * n_K], j] <- 1
  #  }
  #}
  
  #true_assign_1 <- apply(Z_1, 1, function(x) which(x == 1))
  
  #II <- diag(K)
  #Perm_index_2 <- sample(K)
  #Perm_mega <- matrix(0, nrow = K, ncol = K)
  #for (i in 1:K) {
  #  Perm_mega[i, ] <- II[Perm_index_2[i], ]
  #}
  
  #Perm_2 <- matrix(0, nrow = n, ncol = n)
  #for (i in 1:K) {
  #  ind <- which(Perm_mega[i, ] == 1)
  #  Perm_2[(i - 1) * n_K + 1:(i * n_K), ] <- I[(ind - 1) * n_K + 1:(ind * n_K), ]
  #}
  
  #Z_2 <- matrix(0, nrow = K, ncol = rr)
  #for (j in 1:rr) {
  #  for (i in 1:K_r) {
  #    Z_2[Perm_index_2[i + (j - 1) * K_r], j] <- 1
  #  }
  #}
  
  #true_assign_2 <- apply(Z_2, 1, function(x) which(x == 1))
  
  #Perm <- Perm_2 %*% Perm_1
  
  #Perm_index <- apply(Perm, 1, function(x) which(x == 1))
  #Z <- matrix(0, nrow = n, ncol = K)
  #for (j in 1:K) {
  #  for (i in 1:n_K) {
  #    Z[Perm_index[i + (j - 1) * n_K], j] <- 1
  #  }
  #}
  
  #true_assign <- apply(Z, 1, function(x) which(x == 1))
  
  # true_assignment_3 begins
  #Perm_index_3 <- apply(Perm, 1, function(x) which(x == 1))
  #Z_3 <- matrix(0, nrow = n, ncol = rr)
  #for (j in 1:rr) {
  #  for (i in 1:n_r) {
  #    Z_3[Perm_index_3[i + (j - 1) * n_r], j] <- 1
  #  }
  #}
  
  #true_assign_3 <- apply(Z_3, 1, function(x) which(x == 1))
  # true_assignment_3 ends
  
  # true_assign is the true assignment vector
  #P_scrambled <- t(Perm) %*% P %*% Perm
  # P_scrambled is the scrambled matrix of P using the permutation matrix
  
  #A <- matrix(rbinom(n * n, 1, P_scrambled), nrow = n, ncol = n)
  #for (i in 1:n) {
  #  for (j in (i + 1):n) {
  #    A[j, i] <- A[i, j]
  #  }
  #  A[i, i] <- 0
  #}
#}  

# Continue with the rest of the code...
  
  
  
                  