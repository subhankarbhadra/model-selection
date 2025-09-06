library(pbapply)

source("bkm_functions.R")

Latent_positions <- function(A, d) {
  n <- ncol(A)
  e <- irlba(A, nu = d, nv = d)
  Uhat <- e$u
  Uhat %*% diag(sqrt(e$d))
}

#Dolphin data
G = read_graph("data/polblogs.gml", format = "gml")  # read the graph
G <- as_undirected(G,mode="each") # make undirected
foo <- which(components(G)$membership==1) # node labels for the largest connected component
G = induced_subgraph(G, foo)
vertex_attr(G, "id", index = V(G)) = vertex_attr(G, "label", index = V(G))
A = as_adjacency_matrix(G,sparse=FALSE) 
table(A) # has 0,1,2,3.
sum(diag(A)) # are there non-zero diagonal elements?
diag(A) <- 0
A[A > 0] <- 1
sum(A)/2

k <- 2
ff <- pbreplicate(100, sbm_vs_dcbm(A, k, 100, 25, 50)$pvalue)
ff
mean(ff)

gg <- pbreplicate(100, dcbm_vs_pabm(A, k, 100, 25, 50)$pvalue)
gg
mean(gg)
