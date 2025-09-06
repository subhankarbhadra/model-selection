library(pbapply)

source("bkm_functions.R")

Latent_positions <- function(A, d) {
  n <- ncol(A)
  e <- irlba(A, nu = d, nv = d)
  Uhat <- e$u
  Uhat %*% diag(sqrt(e$d))
}
setwd("data//DBLP_four_area")
a = read.csv("csv/author_label.csv", header = F, col.names = c("a.id","a.lab","a.name"))
a <- a[order(a$a.id),]
a$a.id = as.factor(a$a.id); a$a.lab = as.factor(a$a.lab);
a.list = levels(a$a.id)

ap = read.csv("csv/paper_author.csv", header = F, col.names = c("p.id","a.id"))
pc = read.csv("csv/paper_conf.csv", header = F, col.names = c("p.id","c.id"))
conf = rep(0, dim(ap)[1])
for (i in 1:dim(ap)[1]){
  id = ap[i,]$p.id
  conf[i] = subset(pc, p.id == id)$c.id
}

pac = cbind(ap, conf)
pac1 = subset(pac, a.id %in% a.list, select = c(a.id, conf))
pac1$a.id = as.factor(pac1$a.id); pac1$conf = as.factor(pac1$conf)
n.a = length(levels(pac1$a.id)); n.c = length(levels(pac1$conf))

pac2 = subset(pac, a.id %in% a.list, select = c(a.id, p.id))
pac2$a.id = as.factor(pac2$a.id); pac2$p.id = as.factor(pac2$p.id)
n.p = length(levels(pac2$p.id))

#Adjacency matrices
write.csv(levels(pac1$a.id), "csv/a-levels.csv")
write.csv(levels(pac1$conf), "csv/c-levels.csv")
write.csv(levels(pac2$p.id), "csv/p-levels.csv")
levels(pac1$a.id) <- 1:n.a
levels(pac1$conf) <- 1:n.c
levels(pac2$a.id) <- 1:n.a
levels(pac2$p.id) <- 1:n.p

AC = matrix(0, nrow = n.a, ncol = n.c)
for (i in 1:n.a){
  id = unique(subset(pac1, a.id == i)$conf)
  AC[i, id] = 1
}

ACA = matrix(0, nrow = n.a, ncol = n.a)
for (i in 1:(n.a-1)){
  for (j in (i+1):n.a){
    if (sum(AC[i,]*AC[j,]) > 0) {ACA[i,j] = 1; ACA[j,i] = 1}
  }}

ACA <- as(ACA, "sparseMatrix")
classes <- as.numeric(a$a.lab)
k <- length(unique(classes))

ff <- pbreplicate(100, sbm_vs_dcbm(ACA, k, 100, 25, 50)$pvalue)
ff
mean(ff)

gg <- pbreplicate(100, dcbm_vs_pabm(ACA, k, 100, 25, 50)$pvalue)
gg
mean(gg)
