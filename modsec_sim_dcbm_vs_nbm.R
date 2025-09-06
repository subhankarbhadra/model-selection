library(pbapply)

source("bkm_functions.R")
source("nbm.R")

exnbm <- function(n, k, l, w) {
  P <- Data_HBM(n, k, l, w)
  A <- rg_sample(P)
  
  maxiter <- 100
  nstart <- 25
  
  dcbm_vs_pabm(A, k, maxiter, nstart, 200)$pvalue
}

n <- 900
k <- 6

# set.seed(12345)
ff <- pbreplicate(100, exnbm(n, k, 2, 0.6)) 
mean(ff)
mean(ff < 0.05)
