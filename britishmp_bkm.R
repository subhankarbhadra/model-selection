library(pbapply)

source("bkm_functions.R")

# British MP dataset
foo = read.table("data/BritishMPcomm.txt", header = T)
mat = read.table("data/politicsuk-retweets.mtx", header = F)
mat[1,]  # no of nodes, and no of edges (all 5 communities)
mat = mat[-1,]
A = matrix(0, nrow(foo), nrow(foo))
for (k in 1:nrow(mat)){ #adjacency matrix with 5 parties
  edgei = mat[k, 1]
  edgej = mat[k, 2]
  i = which(foo[, 1] == edgei)
  j = which(foo[, 1] == edgej)
  A[i, j] <- 1; A[j, i] <- 1
}
b.true <- as.factor(foo[, 2])
levels(b.true) <- 1:5
b.true <- as.numeric(b.true)
table(b.true) # membership sizes of the 5 parties
nodes1 = which(b.true == 1); nodes2 = which(b.true == 2)
nodes = c(nodes1, nodes2)
A <- A[nodes, nodes] # Select only the two big communities
G <- graph_from_adjacency_matrix(A, mode = "undirected")
is_connected(G, mode = "strong") # check if the 2-comm network is connected: nope
foo1 <- components(G, mode = "strong") 
nodes2 <- which(foo1$membership == 1) # largest connected component (basically 31 nodes are disconnected)
A <- A[nodes2, nodes2] # form network with 329 nodes
G <- graph_from_adjacency_matrix(A, mode = "undirected") 
is_connected(G, mode = "strong") # now the network is connected
b.true <- b.true[nodes2]

Latent_positions <- function(A, d) {
  n <- ncol(A)
  e <- irlba(A, nu = d, nv = d)
  Uhat <- e$u
  Uhat %*% diag(sqrt(e$d))
}

k <- 2

# Model selection (SBM vs. DCBM)
ff <- pbreplicate(100, sbm_vs_dcbm(A, k, 100, 25, 50)$pvalue)
ff
mean(ff)

# Model selection (DCBM vs. PABM)
gg <- pbreplicate(100, dcbm_vs_pabm(A, k, 100, 25, 50)$pvalue)
gg
mean(gg)
