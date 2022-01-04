library(pcalg)
source("simulation_DAG/graph_generation.R")
# args <- commandArgs()
# p <- as.numeric(args[6])
# n_tol <- as.numeric(args[7])
# Define parameters
p <- 100
n_tol <- 600
K <- 2
n <- n_tol / K
e_com <- 50
e_pri <- 50
# Define prior
prior_vec_list <- list()
prior_vec_list[[1]] <- c(1 / p^1.5, 1 / p^2)
prior_vec_list[[2]] <- c(1 / (2 * p^1.5), 1 / p^2)
prior_vec_list[[3]] <- c(1 / p^1.5, 1 / p^1.5)
prior_vec_list[[4]] <- c(1 / p^2, 1 / p^2)
# Define MCMC parameters
scale_x <- FALSE
intercept <- TRUE

## define metric function
## remove order edge
check_edge <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  return(sum(abs(adj_pre - adj_act)) / 2)
}
## True positive rate
TPrate_fun <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  P <- which(adj_act == 1)
  PP <- which(adj_pre == 1)
  return(length(intersect(P, PP)) / length(P))
}
## False positive rate
FPrate_fun <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  N <- which(adj_act == 0)
  PP <- which(adj_pre == 1)
  return(length(intersect(N, PP)) / length(N))
}
## False negative rate
FNrate_fun <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  P <- which(adj_act == 1)
  PN <- which(adj_pre == 0)
  return(length(intersect(PN, P)) / length(P))
}
## check adjacency matrix
check_adj_l2 <- function(adj_pre, adj_act) {
  adj_pre <- adj_pre + t(adj_pre)
  adj_act <- adj_act + t(adj_act)
  return(sum((adj_pre - adj_act)^2) / 2)
}
check_adj_l1 <- function(adj_pre, adj_act) {
  adj_pre <- adj_pre + t(adj_pre)
  adj_act <- adj_act + t(adj_act)
  return(sum(abs(adj_pre - adj_act)) / 2)
}

#### generate graph
set.seed(2021)
source("Two_dataset_new/Graph_MCMC_two_sim.R")
n_graph <- 20
graph_sim <- graph_generation(
  K = K, n_graph = n_graph, p = p, n_tol = n_tol,
  e_com = e_com, e_pri = e_pri
)

########################### Check MCMC results ##################
iter_prior <- 2
prior_vec <- prior_vec_list[[iter_prior]]
out_res1 <- readRDS(paste0("simulation_DAG/results/MCMC_results/com", e_com, "_pri", e_pri, "/",
                   prior_vec[1], "_", prior_vec[2], "_", e_com, "_", e_pri, "_", "alpha_mat1.rds"))
out_res2 <- readRDS(paste0("simulation_DAG/results/MCMC_results/com", e_com, "_pri", e_pri, "/",
                           prior_vec[1], "_", prior_vec[2], "_", e_com, "_", e_pri, "_", "alpha_mat2.rds"))
res_1 <- matrix(NA, nrow = n_graph, ncol = 7)
res_2 <- matrix(NA, nrow = n_graph, ncol = 7)
for (iter_graph in seq_len(n_graph)) {
  alpha_mat_1 <- out_res1[[iter_graph]]
  alpha_mat_2 <- out_res2[[iter_graph]]
  ## data set 1
  adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
  adj_1 <- t(adj_1)
  g_1 <- as(getGraph(adj_1), "graphNEL")
  ## data set 2
  adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
  adj_2 <- t(adj_2)
  g_2 <- as(getGraph(adj_2), "graphNEL")
  ## load true value
  adj_true1 <- t(graph_sim$G[[iter_graph]][[1]])
  g_true1 <- as(getGraph(adj_true1), "graphNEL")
  adj_true2 <- t(graph_sim$G[[iter_graph]][[2]])
  g_true2 <- as(getGraph(adj_true2), "graphNEL")
  ## save results
  res_1[iter_graph, ] <- c(
    shd(g_true1, g_1),
    check_edge(adj_true1, adj_1),
    TPrate_fun(adj_pre = adj_1, adj_act = adj_true1),
    FPrate_fun(adj_pre = adj_1, adj_act = adj_true1),
    FNrate_fun(adj_pre = adj_1, adj_act = adj_true1),
    check_adj_l2(adj_pre = alpha_mat_1, adj_act = adj_true1),
    check_adj_l1(adj_pre = alpha_mat_1, adj_act = adj_true1)
  )
  res_2[iter_graph, ] <- c(
    shd(g_true2, g_2),
    check_edge(adj_true2, adj_2),
    TPrate_fun(adj_pre = adj_2, adj_act = adj_true2),
    FPrate_fun(adj_pre = adj_2, adj_act = adj_true2),
    FNrate_fun(adj_pre = adj_2, adj_act = adj_true2),
    check_adj_l2(adj_pre = alpha_mat_2, adj_act = adj_true2),
    check_adj_l1(adj_pre = alpha_mat_2, adj_act = adj_true2)
  )
}

print(round(res_1, 4))
print(round(colMeans(res_1), 4))
print(round(res_2, 4))
print(round(colMeans(res_2), 4))

########################### Check GES results  ##################
out_res <- readRDS("simulation_DAG/results/GES_results/out_res_ges.rds")
iter_lambda <- 2
res_1 <- matrix(NA, nrow = n_graph, ncol = 7)
res_2 <- matrix(NA, nrow = n_graph, ncol = 7)
for (iter_graph in seq_len(n_graph)) {
  ## load true value
  adj_true1 <- t(graph_sim$G[[iter_graph]][[1]])
  g_true1 <- as(getGraph(adj_true1), "graphNEL")
  adj_true2 <- t(graph_sim$G[[iter_graph]][[2]])
  g_true2 <- as(getGraph(adj_true2), "graphNEL")
  ## load results
  adj1 <- out_res[[iter_graph]][[1]][[iter_lambda]]
  g1 <- as(adj1, "graphNEL")
  adj2 <- out_res[[iter_graph]][[2]][[iter_lambda]]
  g2 <- as(adj2, "graphNEL")
  ## save results
  res_1[iter_graph, ] <- c(
    shd(g_true1, g1),
    check_edge(adj_true1, adj1),
    TPrate_fun(adj_pre = adj1, adj_act = adj_true1),
    FPrate_fun(adj_pre = adj1, adj_act = adj_true1),
    FNrate_fun(adj_pre = adj1, adj_act = adj_true1),
    check_adj_l2(adj_pre = adj1, adj_act = adj_true1),
    check_adj_l1(adj_pre = adj1, adj_act = adj_true1)
  )
  res_2[iter_graph, ] <- c(
    shd(g_true2, g2),
    check_edge(adj_true2, adj2),
    TPrate_fun(adj_pre = adj2, adj_act = adj_true2),
    FPrate_fun(adj_pre = adj2, adj_act = adj_true2),
    FNrate_fun(adj_pre = adj2, adj_act = adj_true2),
    check_adj_l2(adj_pre = adj2, adj_act = adj_true2),
    check_adj_l1(adj_pre = adj2, adj_act = adj_true2)
  )
}

print(round(res_1, 4))
print(round(colMeans(res_1), 4))
print(round(res_2, 4))
print(round(colMeans(res_2), 4))

############## Check specific graph ##################
graph_index <- 6
adj_true1 <- t(graph_sim$G[[graph_index]][[1]])
weight_true1 <- t(graph_sim$A[[graph_index]][[1]])
adj_true2 <- t(graph_sim$G[[graph_index]][[2]])
weight_true2 <- t(graph_sim$A[[graph_index]][[2]])
dta_1 <- graph_sim$X[[graph_index]][[1]]
dta_2 <- graph_sim$X[[graph_index]][[2]]

#### check the initial order
# get order
set.seed(2021)
dta <- rbind(dta_1, dta_2)
score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
ges_fit <- ges(score_ges)
ges_adj <- as(ges_fit$repr, "matrix")
ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
order_int <- as.numeric(igraph::topo_sort(graph_i))
# Do MCMC
source("Two_dataset_new/Graph_given_order_two.R")
set.seed(2021)
iter_max <- 100
out_res <- Graph_MCMC_two_sim(dta_1, dta_2, scale_x = scale_x, intercept = intercept,
                              order_int = order_int, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
                              prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
                              burn_in = 1, adj_true1 = adj_true1, adj_true2 = adj_true2
)
# Show likelihood
library(ggplot2)
library(gridExtra)
gl <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$llike_vec)) +
  xlab("Iteration") + ylab("Log likelihood")
gs_1 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat1[1, ])) +
  xlab("Iteration") + ylab("SHD")
gu_1 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat1[2, ])) +
  xlab("Iteration") + ylab("No order")
gs_2 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat2[1, ])) +
  xlab("Iteration") + ylab("SHD")
gu_2 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat2[2, ])) +
  xlab("Iteration") + ylab("No order")
layout_matrix <- matrix(c(1, 2, 3), nrow = 3)
grid.arrange(gl, gs_1, gu_1, layout_matrix = layout_matrix)
grid.arrange(gl, gs_2, gu_2, layout_matrix = layout_matrix)

## save the true index
index1_true <- which(adj_true1 > 0)
index2_true <- which(adj_true2 > 0)
index_intersect_true <- intersect(index1_true, index2_true)

## mcmc results
adj_mcmc1 <- ifelse(out_res1[[graph_index]] > 0.5, 1, 0)
adj_mcmc1 <- ceiling((adj_mcmc1 + t(adj_mcmc1)) / 2)
adj_mcmc1[lower.tri(adj_mcmc1)] <- 0 
index1_mcmc <- which(adj_mcmc1 > 0)
adj_mcmc2 <- ifelse(out_res2[[graph_index]] > 0.5, 1, 0)
adj_mcmc2 <- ceiling((adj_mcmc2 + t(adj_mcmc2)) / 2)
adj_mcmc2[lower.tri(adj_mcmc2)] <- 0 
index2_mcmc <- which(adj_mcmc2 > 0)
index_intersect_mcmc <- intersect(index1_mcmc, index2_mcmc)

## ges results
adj_ges1 <- out_res[[graph_index]][[1]][[2]]
adj_ges1 <- ceiling((adj_ges1 + t(adj_ges1)) / 2)
adj_ges1[lower.tri(adj_ges1)] <- 0 
index1_ges <- which(adj_ges1 > 0)
adj_ges2 <- out_res[[graph_index]][[2]][[2]]
adj_ges2 <- ceiling((adj_ges2 + t(adj_ges2)) / 2)
adj_ges2[lower.tri(adj_ges2)] <- 0 
index2_ges <- which(adj_ges2 > 0)
index_intersect_ges <- intersect(index1_ges, index2_ges)

################## check mcmc and ges ###################
TP1_mcmc <- intersect(index1_true, index1_mcmc)
TP2_mcmc <- intersect(index2_true, index2_mcmc)
TP1_ges <- intersect(index1_true, index1_ges)
TP2_ges <- intersect(index2_true, index2_ges)

FP1_mcmc <- setdiff(index1_mcmc, index1_true)
FP2_mcmc <- setdiff(index2_mcmc, index2_true)
FP1_ges <- setdiff(index1_ges, index1_true)
FP2_ges <- setdiff(index2_ges, index2_true)

FN1_mcmc <- setdiff(index1_true, index1_mcmc)
FN2_mcmc <- setdiff(index2_true, index2_mcmc)
FN1_ges <- setdiff(index1_true, index1_ges)
FN2_ges <- setdiff(index2_true, index2_ges)

################## check mcmc and ges ###################
## check FN 
# data set 1
length(FN1_mcmc)
weight_true1[FN1_mcmc]

length(FN1_ges)
weight_true1[FN1_ges]

# data set 2
length(FN2_mcmc)
weight_true2[FN2_mcmc]

length(FN2_ges)
weight_true2[FN2_ges]

## check FP
# data set 1
length(FP1_mcmc)
weight_true2[FP1_mcmc]

length(FP1_ges)
weight_true2[FP1_ges]

# data set 2
length(FP2_mcmc)
weight_true1[FP2_mcmc]

length(FP2_ges)
weight_true1[FP2_ges]

################## compare mcmc and ges ###################
summary(abs(weight_true1[index1_true]))
summary(abs(weight_true2[index2_true]))

## check TP
# data set 1
length(setdiff(TP1_mcmc, TP1_ges))
weight_true1[setdiff(TP1_mcmc, TP1_ges)]
weight_true2[setdiff(TP1_mcmc, TP1_ges)]

length(setdiff(TP1_ges, TP1_mcmc))
weight_true1[setdiff(TP1_ges, TP1_mcmc)]
weight_true2[setdiff(TP1_ges, TP1_mcmc)]

# data set 2
length(setdiff(TP2_mcmc, TP2_ges))
weight_true2[setdiff(TP2_mcmc, TP2_ges)]
weight_true1[setdiff(TP2_mcmc, TP2_ges)]

length(setdiff(TP2_ges, TP2_mcmc))
weight_true2[setdiff(TP2_ges, TP2_mcmc)]
weight_true1[setdiff(TP2_ges, TP2_mcmc)]

## check FP
# data set 1
length(setdiff(FP1_mcmc, FP1_ges))
weight_true2[setdiff(FP1_mcmc, FP1_ges)]

length(setdiff(FP1_ges, FP1_mcmc))
weight_true2[setdiff(FP1_ges, FP1_mcmc)]

# data set 2
length(setdiff(FP2_mcmc, FP2_ges))
weight_true1[setdiff(FP2_mcmc, FP2_ges)]

length(setdiff(FP2_ges, FP2_mcmc))
weight_true1[setdiff(FP2_ges, FP2_mcmc)]

## check FN
# data set 1
length(setdiff(FN1_mcmc, FN1_ges))
weight_true1[setdiff(FN1_mcmc, FN1_ges)]
weight_true2[setdiff(FN1_mcmc, FN1_ges)]

length(setdiff(FN1_ges, FN1_mcmc))
weight_true1[setdiff(FN1_ges, FN1_mcmc)]
weight_true2[setdiff(FN1_ges, FN1_mcmc)]

# data set 2
length(setdiff(FN2_mcmc, FN2_ges))
weight_true2[setdiff(FN2_mcmc, FN2_ges)]
weight_true1[setdiff(FN2_mcmc, FN2_ges)]

length(setdiff(FN2_ges, FN2_mcmc))
weight_true2[setdiff(FN2_ges, FN2_mcmc)]
weight_true1[setdiff(FN2_ges, FN2_mcmc)]