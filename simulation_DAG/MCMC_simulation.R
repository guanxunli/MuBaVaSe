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

#### our method
source("Two_dataset/Graph_MCMC_two.R")
# source("Two_dataset/Graph_MCMC_two_init.R")
dta_1 <- graph_sim$X[[1]][[1]]
dta_2 <- graph_sim$X[[1]][[2]]
out_res <- Graph_MCMC_two(dta_1, dta_2, order_int = NULL, iter_max = 10000, sigma02_int = NULL, sigma2_int = NULL, r = 0.2,
                          q = 0.05, tau = 1.5, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 5000)
# out_res <- Graph_MCMC_two_init(dta_1, dta_2, order_int = NULL, iter_max = 10000, sigma02_int = NULL, sigma2_int = NULL, r = 0.2,
#                           q = 0.05, tau = 1.5, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 5000)

#### analysis results
out_res <- readRDS("simulation_DAG/out_res.rds")
alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
A_mat_1 <- matrix(0, nrow = p, ncol = p)
A_mat_2 <- matrix(0, nrow = p, ncol = p)
for (iter in seq_len(length(out_res$order_list))) {
  order_tmp <- order(out_res$order_list[[iter]])
  alpha_mat_1 <- alpha_mat_1 + out_res$alpha_list_1[[iter]][order_tmp, order_tmp]
  alpha_mat_2 <- alpha_mat_2 + out_res$alpha_list_2[[iter]][order_tmp, order_tmp]
  A_mat_1 <- A_mat_1 + out_res$A_list_1[[iter]][order_tmp, order_tmp]
  A_mat_2 <- A_mat_2 + out_res$A_list_2[[iter]][order_tmp, order_tmp]
}

alpha_mat_1 <- alpha_mat_1 / length(out_res$order_list)
alpha_mat_2 <- alpha_mat_2 / length(out_res$order_list)
A_mat_1 <- A_mat_1 / length(out_res$order_list)
A_mat_2 <- A_mat_2 / length(out_res$order_list)

#### Calculate the error
## data set 1
adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
adj_1 <- t(adj_1)
g_1 <- as(getGraph(adj_1), "graphNEL")
weight_1 <- t(A_mat_1)
# structural Hamming distance (SHD)
shd(g_true1, g_1)
# Mean square error for weight
sum((weight_true1 - weight_1)^2)
# TPR & FPR
TPrate_fun(adj_pre = adj_1, adj_act = adj_true1)
FPrate_fun(adj_pre = adj_1, adj_act = adj_true1)

## data set 2
adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
adj_2 <- t(adj_2)
g_2 <- as(getGraph(adj_2), "graphNEL")
weight_2 <- t(A_mat_2)
# structural Hamming distance (SHD)
shd(g_true2, g_2)
# Mean square error for weight
sum((weight_true2 - weight_2)^2)
# TPR & FPR
TPrate_fun(adj_pre = adj_2, adj_act = adj_true2)
FPrate_fun(adj_pre = adj_2, adj_act = adj_true2)
