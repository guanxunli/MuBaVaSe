library(pcalg)
source("simulation_DAG/graph_generation.R")
# args <- commandArgs()
# init <- args[6]
# p <- as.numeric(args[7])
# n_tol <- as.numeric(args[8])
init <- "noinit"
p <- 100
n_tol <- 2000
K <- 2
n <- n_tol / K
n_graph <- 1

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

# check_weight_l2 <- function(weight_pre, weight_act) {
#   weight_pre <- weight_pre + t(weight_pre)
#   weight_act <- weight_act + t(weight_act)
#   return(sum((weight_pre - weight_act)^2) / 2)
# }
# 
# check_weight_l1 <- function(weight_pre, weight_act) {
#   weight_pre <- weight_pre + t(weight_pre)
#   weight_act <- weight_act + t(weight_act)
#   return(sum(abs(weight_pre - weight_act)) / 2)
# }

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
source("Two_dataset_new/Graph_MCMC_two.R")
dta_1 <- graph_sim$X[[1]][[1]]
dta_2 <- graph_sim$X[[1]][[2]]
## get initialized order
iter_max <- 100000
if (init == "init") {
  out_res <- Graph_MCMC_two(dta_1, dta_2,
                            order_int = seq_len(p), iter_max = 20000, sigma02_int = NULL, sigma2_int = NULL,
                            prior_vec = NULL, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 1
  )
} else if (init == "noinit") {
  out_res <- Graph_MCMC_two(dta_1, dta_2,
                            order_int = NULL, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
                            prior_vec = NULL, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 1
  )
} else if (init == "ges_init") {
  ## get order
  dta <- rbind(dta_1, dta_2)
  score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
  ges_fit <- ges(score_ges)
  ges_adj <- as(ges_fit$repr, "matrix")
  ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
  graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
  order_int <- as.numeric(igraph::topo_sort(graph_i))
  ## Do MCMC
  out_res <- Graph_MCMC_two(dta_1, dta_2,
                            order_int = order_int, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
                            prior_vec = NULL, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 1
  )
}

#### plot log-likelihood
library(ggplot2)
library(gridExtra)
g1 <- ggplot() +
  geom_line(aes(x = seq_len(length(out_res$llike_vec)), y = out_res$llike_vec)) +
  xlab("steps") +
  ylab("log-likelihood") +
  labs(title = "ALL log-likelihood") +
  ylim(c(range(out_res$llike_vec) + c(-5, 5)))
g2 <- ggplot() +
  geom_line(aes(x = seq_len(length(out_res$llike_vec[-seq_len(iter_max - 5000)])), y = out_res$llike_vec[-seq_len(iter_max - 5000)])) +
  xlab("steps") +
  ylab("log-likelihood") +
  labs(title = "Zoom log-likelihood") +
  ylim(c(range(out_res$llike_vec[-seq_len(iter_max - 5000)]) + c(-5, 5)))
layout_matrix <- matrix(c(1, 2), nrow = 2)
pdf(paste0("llikehood_", init, "_", p, "_", n_tol, ".pdf"), width = 10, height = 6.18)
grid.arrange(g1, g2, layout_matrix = layout_matrix)
dev.off()

out_res$alpha_list_1 <- out_res$alpha_list_1[-seq_len(iter_max - 5000)]
out_res$alpha_list_2 <- out_res$alpha_list_2[-seq_len(iter_max - 5000)]
out_res$A_list_1 <- out_res$A_list_1[-seq_len(iter_max - 5000)]
out_res$A_list_2 <- out_res$A_list_2[-seq_len(iter_max - 5000)]
out_res$order_list <- out_res$order_list[-seq_len(iter_max - 5000)]
out_res$llike_vec <- out_res$llike_vec[-seq_len(iter_max - 5000)]
saveRDS(out_res, paste0("simulation_DAG/results/out_res_", init, "_", p, "_", n_tol, ".rds"))

#### analysis results
out_res <- readRDS(paste0("simulation_DAG/results/out_res_", init, "_", p, "_", n_tol, ".rds"))
alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
A_mat_1 <- matrix(0, nrow = p, ncol = p)
A_mat_2 <- matrix(0, nrow = p, ncol = p)
iter_max <- length(out_res[[1]])
for (iter in seq_len(5000)) {
  order_tmp <- order(out_res$order_list[[iter]])
  alpha_mat_1 <- alpha_mat_1 + out_res$alpha_list_1[[iter]][order_tmp, order_tmp]
  alpha_mat_2 <- alpha_mat_2 + out_res$alpha_list_2[[iter]][order_tmp, order_tmp]
  A_mat_1 <- A_mat_1 + out_res$A_list_1[[iter]][order_tmp, order_tmp]
  A_mat_2 <- A_mat_2 + out_res$A_list_2[[iter]][order_tmp, order_tmp]
}

alpha_mat_1 <- alpha_mat_1 / 5000
alpha_mat_2 <- alpha_mat_2 / 5000
A_mat_1 <- A_mat_1 / 5000
A_mat_2 <- A_mat_2 / 5000

#### Calculate the error
## data set 1
adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
adj_1 <- t(adj_1)
g_1 <- as(getGraph(adj_1), "graphNEL")
weight_1 <- t(A_mat_1)
weight_1[which(adj_1 == 0)] <- 0
# structural Hamming distance (SHD) and undirected edge
print(c(shd(g_true1, g_1), check_edge(adj_true1, adj_1)))
# # Mean square error for weight
# print(c(round(sum((weight_true1 - weight_1)^2), 2), round(check_weight_l2(weight_1, weight_true1), 2)))
# # l1 error
# print(c(round(sum(abs(weight_true1 - weight_1)), 2), round(check_weight_l1(weight_1, weight_true1), 2)))
# TPR & FPR
print(c(round(TPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4), 
        round(FPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4)))

## data set 2
adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
adj_2 <- t(adj_2)
g_2 <- as(getGraph(adj_2), "graphNEL")
weight_2 <- t(A_mat_2)
weight_2[which(adj_2 == 0)] <- 0
# structural Hamming distance (SHD) and undirected edge
print(c(shd(g_true2, g_2), check_edge(adj_true2, adj_2)))
# # Mean square error for weight
# print(c(round(sum((weight_true2 - weight_2)^2), 2), round(check_weight_l2(weight_2, weight_true2), 2)))
# # l1 error
# print(c(round(sum(abs(weight_true2 - weight_2)), 2), round(check_weight_l1(weight_2, weight_true2), 2)))
# TPR & FPR
print(c(round(TPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4), 
        round(FPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4)))

## output results
cat(
  "MCMC", "&", shd(g_true1, g_1), "&", check_edge(adj_true1, adj_1), "&",
  shd(g_true2, g_2), "&", check_edge(adj_true2, adj_2), "&", "\\\\\n"
)