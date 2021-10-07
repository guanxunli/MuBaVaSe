library(pcalg)
source("simulation_DAG/graph_generation.R")
p <- 50
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
## remove order edge
check_edge <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  return(sum(abs(adj_pre - adj_act)) / 2)
}

check_weight_l2 <- function(weight_pre, weight_act) {
  weight_pre <- weight_pre + t(weight_pre)
  weight_act <- weight_act + t(weight_act)
  return(sum((weight_pre - weight_act)^2) / 2)
}

check_weight_l1 <- function(weight_pre, weight_act) {
  weight_pre <- weight_pre + t(weight_pre)
  weight_act <- weight_act + t(weight_act)
  return(sum(abs(weight_pre - weight_act)) / 2)
}
#### generate graph
set.seed(202110)
graph_sim <- graph_generation(K = K, n_graph = n_graph, p = p, n_tol = n_tol, e_com = 100, e_pri = 30)
adj_true1 <- t(graph_sim$G[[1]][[1]])
g_true1 <- as(getGraph(adj_true1), "graphNEL")
weight_true1 <- t(graph_sim$A[[1]][[1]])
adj_true2 <- t(graph_sim$G[[1]][[2]])
g_true2 <- as(getGraph(adj_true2), "graphNEL")
weight_true2 <- t(graph_sim$A[[1]][[2]])
#### Given order method
source("Two_dataset/Graph_given_order_two.R")
dta_1 <- graph_sim$X[[1]][[1]]
dta_2 <- graph_sim$X[[1]][[2]]
out_res <- joint_graph_fun_two(dta_1 = dta_1, dta_2 = dta_2)
#### GES method
## data set 1
# score1 <- new("GaussL0penObsScore", data = t(graph_sim$X[[1]][[1]]), intercept = FALSE,
#               lambda = (cons * log(p) / n))
score1 <- new("GaussL0penObsScore", data = graph_sim$X[[1]][[1]], intercept = FALSE)
ges_fit1 <- ges(score1)
ges_adj1 <- as(ges_fit1$repr, "matrix")
ges_adj1 <- ifelse(ges_adj1 == TRUE, 1, 0)
ges_graph1 <- as(ges_fit1$repr, "graphNEL")
ges_weight1 <- ges_fit1$repr$weight.mat()
## our method
adj_1 <- out_res$alpha_res_1
adj_1 <- t(ifelse(adj_1 > 0.5, 1, 0))
g_1 <- as(adj_1, "graphNEL")
weight_1 <- t(out_res$A_res_1)
weight_1[which(adj_1 == 0)] <- 0
# structural Hamming distance (SHD) and undirected edge
print(c(shd(g_true1, ges_graph1), check_edge(adj_true1, ges_adj1)))
print(c(shd(g_true1, g_1), check_edge(adj_true1, adj_1)))
# Mean square error for weight
print(c(round(sum((weight_true1 - ges_weight1)^2), 2), round(check_weight_l2(ges_weight1, weight_true1), 2)))
print(c(round(sum((weight_true1 - weight_1)^2), 2), round(check_weight_l2(weight_1, weight_true1), 2)))
# l1 error
print(c(round(sum(abs(weight_true1 - ges_weight1)), 2), round(check_weight_l1(ges_weight1, weight_true1), 2)))
print(c(round(sum(abs(weight_true1 - weight_1)), 2), round(check_weight_l1(weight_1, weight_true1), 2)))
# TPR & FPR
print(c(round(TPrate_fun(adj_pre = ges_adj1, adj_act = adj_true1), 4), round(FPrate_fun(adj_pre = ges_adj1, adj_act = adj_true1), 4)))
print(c(round(TPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4), round(FPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4)))

## data set 2
## our method
adj_2 <- out_res$alpha_res_2
adj_2 <- t(ifelse(adj_2 > 0.5, 1, 0))
g_2 <- as(adj_2, "graphNEL")
weight_2 <- t(out_res$A_res_2)
weight_2[which(adj_2 == 0)] <- 0
# score2 <- new("GaussL0penObsScore", data = t(graph_sim$X[[1]][[2]]), intercept = FALSE,
#               lambda = (cons * log(p) / n))
score2 <- new("GaussL0penObsScore", data = graph_sim$X[[1]][[2]], intercept = FALSE)
ges_fit2 <- ges(score2)
ges_adj2 <- as(ges_fit2$repr, "matrix")
ges_adj2 <- ifelse(ges_adj2 == TRUE, 1, 0)
ges_graph2 <- as(ges_fit2$repr, "graphNEL")
ges_weight2 <- ges_fit2$repr$weight.mat()
# structural Hamming distance (SHD) and undirected edge
print(c(shd(g_true2, ges_graph2), check_edge(adj_true2, ges_adj2)))
print(c(shd(g_true2, g_2), check_edge(adj_true2, adj_2)))
# Mean square error for weight
print(c(round(sum((weight_true2 - ges_weight2)^2), 2), round(check_weight_l2(ges_weight2, weight_true2), 2)))
print(c(round(sum((weight_true2 - weight_2)^2), 2), round(check_weight_l2(weight_2, weight_true2), 2)))
# l1 error
print(c(round(sum(abs(weight_true2 - ges_weight2)), 2), round(check_weight_l1(ges_weight2, weight_true2), 2)))
print(c(round(sum(abs(weight_true2 - weight_2)), 2), round(check_weight_l1(weight_2, weight_true2), 2)))
# TPR & FPR
print(c(round(TPrate_fun(adj_pre = ges_adj2, adj_act = adj_true2), 4), round(FPrate_fun(adj_pre = ges_adj2, adj_act = adj_true2), 4)))
print(c(round(TPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4), round(FPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4)))

#### output results
cat(
  "GES", "&", shd(g_true1, ges_graph1), "&", check_edge(adj_true1, ges_adj1), "&",
  shd(g_true2, ges_graph2), "&", check_edge(adj_true2, ges_adj2), "&",
  round(sum((weight_true1 - ges_weight1)^2), 2), "&", round(check_weight_l2(ges_weight1, weight_true1), 2), "&",
  round(sum(abs(weight_true1 - ges_weight1)), 2), "&", round(check_weight_l1(ges_weight1, weight_true1), 2), "&",
  round(sum((weight_true2 - ges_weight2)^2), 2), "&", round(check_weight_l2(ges_weight2, weight_true2), 2), "&",
  round(sum(abs(weight_true2 - ges_weight2)), 2), "&", round(check_weight_l1(ges_weight2, weight_true2), 2), "\\\\\n"
)

cat(
  "MCMC", "&", shd(g_true1, g_1), "&", check_edge(adj_true1, adj_1), "&",
  shd(g_true2, g_2), "&", check_edge(adj_true2, adj_2), "&",
  round(sum((weight_true1 - weight_1)^2), 2), "&", round(check_weight_l2(weight_1, weight_true1), 2), "&",
  round(sum(abs(weight_true1 - weight_1)), 2), "&", round(check_weight_l1(weight_1, weight_true1), 2), "&",
  round(sum((weight_true2 - weight_2)^2), 2), "&", round(check_weight_l2(weight_2, weight_true2), 2), "&",
  round(sum(abs(weight_true2 - weight_2)), 2), "&", round(check_weight_l1(weight_2, weight_true2), 2), "\\\\\n"
)