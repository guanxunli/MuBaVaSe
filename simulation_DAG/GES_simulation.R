library(pcalg)
source("simulation_DAG/graph_generation.R")
p <- 100
n_tol <- 600
K <- 2
n <- n_tol / K
n_graph <- 1
cons <- 2
## define metric function
## True positive rate
TPrate_fun <- function(adj_pre, adj_act) {
  P <- which(adj_act == 1)
  PP <- which(adj_pre == 1)
  return(length(intersect(P, PP)) / length(P))
}
## False positive rate
FPrate_fun <- function(adj_pre, adj_act) {
  N <- which(adj_act == 0)
  PP <- which(adj_pre == 1)
  return(length(intersect(N, PP)) / length(N))
}
#### generate graph
set.seed(2021)
graph_sim <- graph_generation(K = K, n_graph = n_graph, p = p, n_tol = n_tol)
adj_true1 <- t(graph_sim$G[[1]][[1]])
g_true1 <- as(getGraph(adj_true1), "graphNEL")
weight_true1 <- t(graph_sim$A[[1]][[1]])
adj_true2 <- t(graph_sim$G[[1]][[2]])
g_true2 <- as(getGraph(adj_true2), "graphNEL")
weight_true2 <- t(graph_sim$A[[1]][[2]])
#### GES method
## data set 1
# GES graph
# score1 <- new("GaussL0penObsScore", data = t(graph_sim$X[[1]][[1]]), intercept = FALSE,
#               lambda = (cons * log(p) / n)) 
score1 <- new("GaussL0penObsScore", data = graph_sim$X[[1]][[1]], intercept = FALSE) 
ges_fit1 <- ges(score1)
ges_adj1 <- as(ges_fit1$repr, "matrix")
ges_adj1 <- ifelse(ges_adj1 == TRUE, 1, 0)
ges_graph1 <- as(ges_fit1$repr, "graphNEL")
ges_weight1 <- ges_fit1$repr$weight.mat()
# structural Hamming distance (SHD)
shd(g_true1, ges_graph1)
# Mean square error for weight
sum((weight_true1 - ges_weight1)^2)
# TPR & FPR
TPrate_fun(adj_pre = ges_adj1, adj_act = adj_true1)
FPrate_fun(adj_pre = ges_adj1, adj_act = adj_true1)

## data set 2
# GES graph
# score2 <- new("GaussL0penObsScore", data = t(graph_sim$X[[1]][[2]]), intercept = FALSE,
#               lambda = (cons * log(p) / n)) 
score2 <- new("GaussL0penObsScore", data = graph_sim$X[[1]][[2]], intercept = FALSE)
ges_fit2 <- ges(score2)
ges_adj2 <- as(ges_fit2$repr, "matrix")
ges_adj2 <- ifelse(ges_adj2 == TRUE, 1, 0)
ges_graph2 <- as(ges_fit2$repr, "graphNEL")
ges_weight2 <- ges_fit2$repr$weight.mat()
# structural Hamming distance (SHD)
shd(g_true2, ges_graph2)
# Mean square error for weight
sum((weight_true2 - ges_weight2)^2)
# TPR & FPR
TPrate_fun(adj_pre = ges_adj2, adj_act = adj_true2)
FPrate_fun(adj_pre = ges_adj2, adj_act = adj_true2)