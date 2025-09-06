library(pbapply)

source("bkm_functions.R")

Latent_positions <- function(A, d) {
  n <- ncol(A)
  e <- irlba(A, nu = d, nv = d)
  Uhat <- e$u
  Uhat %*% diag(sqrt(e$d))
}

#Dolphin data
G = read_graph("data/dolphins.gml", format = "gml")
A = as_adjacency_matrix(G, type="both", sparse = FALSE)
dim(A)
sum(A)/2

k <- 2
ff <- pbreplicate(100, sbm_vs_dcbm(A, k, 100, 25, 50)$pvalue)
ff
mean(ff)

gg <- pbreplicate(100, dcbm_vs_pabm(A, k, 100, 25, 50)$pvalue)
gg
mean(gg)
