library(pcalg)
source("simulation_DAG/graph_generation.R")
# args <- commandArgs()
# p <- as.numeric(args[6])
# n_tol <- as.numeric(args[7])
p <- 100
n_tol <- 600
K <- 2
n <- n_tol / K

## define metric function
## True positive rate
TPrate_fun <- function(adj_pre, adj_act) {
  P <- which(adj_act == 1)
  PP <- which(adj_pre == 1)
  return(length(intersect(P, PP)) / length(P))
}
# ## False positive rate
# FPrate_fun <- function(adj_pre, adj_act) {
#   N <- which(adj_act == 0)
#   PP <- which(adj_pre == 1)
#   return(length(intersect(N, PP)) / length(N))
# }
## remove order edge
check_edge <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  return(sum(abs(adj_pre - adj_act)) / 2)
}

eval_fun <- function(adj, adj_true, g_true) {
  g <- as(getGraph(adj), "graphNEL")
  cat(c(
    shd(g_true, g), check_edge(adj_true, adj)
  ), "\n")
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

########################### Do one figure ##################################
# #### generate graph
# set.seed(2021)
# n_graph <- 1
# graph_sim <- graph_generation(K = K, n_graph = n_graph, p = p, n_tol = n_tol)
# adj_true1 <- t(graph_sim$G[[1]][[1]])
# g_true1 <- as(getGraph(adj_true1), "graphNEL")
# weight_true1 <- t(graph_sim$A[[1]][[1]])
# adj_true2 <- t(graph_sim$G[[1]][[2]])
# g_true2 <- as(getGraph(adj_true2), "graphNEL")
# weight_true2 <- t(graph_sim$A[[1]][[2]])

# #### our method
# source("Two_dataset_new/Graph_MCMC_two.R")
# dta_1 <- graph_sim$X[[1]][[1]]
# dta_2 <- graph_sim$X[[1]][[2]]
# prior_vec <- c(1 / (2 * p^1.5), 1 / p^2)
# # prior_vec <- c(1 / p^1.5, 1 / p^1.5)
# # prior_vec <- c(1 / p^1.5, 1 / p^2)

# #### If we know the order
# out_res <- joint_graph_fun_two(dta_1 = dta_1, dta_2 = dta_2, prior_vec = prior_vec, scale_x = TRUE,
#                                intercept = TRUE)
# ## Calculate the error
# ## data set 1
# adj_1 <- ifelse(out_res$alpha_res_1 > 0.5, 1, 0)
# adj_1 <- t(adj_1)
# g_1 <- as(getGraph(adj_1), "graphNEL")
# cat("prior:", prior_vec, "data1:", c(shd(g_true1, g_1), check_edge(adj_true1, adj_1)), "\n")
# ## data set 2
# adj_2 <- ifelse(out_res$alpha_res_2 > 0.5, 1, 0)
# adj_2 <- t(adj_2)
# g_2 <- as(getGraph(adj_2), "graphNEL")
# cat("prior:", prior_vec, "data2:", c(shd(g_true2, g_2), check_edge(adj_true2, adj_2)), "\n")


# prior: 5e-04 1e-04 data1: 12 12
# prior: 5e-04 1e-04 data2: 9 9
# prior: 0.001 0.001 data1: 6 6
# prior: 0.001 0.001 data2: 7 7

# #### Do MCMC
# iter_max <- 50000

# #### with GES Initialization
# # get order
# set.seed(2021)
# dta <- rbind(dta_1, dta_2)
# score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
# ges_fit <- ges(score_ges)
# ges_adj <- as(ges_fit$repr, "matrix")
# ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
# graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
# order_int <- as.numeric(igraph::topo_sort(graph_i))

# # Do MCMC
# set.seed(2021)
# out_res <- Graph_MCMC_two(dta_1, dta_2,
#                           order_int = order_int, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
#                           prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
#                           burn_in = iter_max - 5000
# )
# # analysis
# alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
# alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
# A_mat_1 <- matrix(0, nrow = p, ncol = p)
# A_mat_2 <- matrix(0, nrow = p, ncol = p)
# for (iter in seq_len(5000)) {
#   order_tmp <- order(out_res$order_list[[iter]])
#   alpha_mat_1 <- alpha_mat_1 + out_res$alpha_list_1[[iter]][order_tmp, order_tmp]
#   alpha_mat_2 <- alpha_mat_2 + out_res$alpha_list_2[[iter]][order_tmp, order_tmp]
#   A_mat_1 <- A_mat_1 + out_res$A_list_1[[iter]][order_tmp, order_tmp]
#   A_mat_2 <- A_mat_2 + out_res$A_list_2[[iter]][order_tmp, order_tmp]
# }
# alpha_mat_1 <- alpha_mat_1 / 5000
# alpha_mat_2 <- alpha_mat_2 / 5000
# A_mat_1 <- A_mat_1 / 5000
# A_mat_2 <- A_mat_2 / 5000
# ## data set 1
# adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
# adj_1 <- t(adj_1)
# g_1 <- as(getGraph(adj_1), "graphNEL")
# cat("prior:", prior_vec, "data1:", c(shd(g_true1, g_1), check_edge(adj_true1, adj_1)), "\n")
# ## data set 2
# adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
# adj_2 <- t(adj_2)
# g_2 <- as(getGraph(adj_2), "graphNEL")
# cat("prior:", prior_vec, "data2:", c(shd(g_true2, g_2), check_edge(adj_true2, adj_2)), "\n")

# prior: 5e-04 1e-04 data1: 32 19
# prior: 5e-04 1e-04 data2: 37 25
# prior: 0.001 0.001 data1: 30 17
# prior: 0.001 0.001 data2: 32 22

# ##### without GES initialization
# out_res <- Graph_MCMC_two(dta_1, dta_2,
#                           order_int = NULL, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
#                           prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
#                           burn_in = iter_max - 5000
# )
#
# # analysis
# alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
# alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
# A_mat_1 <- matrix(0, nrow = p, ncol = p)
# A_mat_2 <- matrix(0, nrow = p, ncol = p)
# for (iter in seq_len(5000)) {
#   order_tmp <- order(out_res$order_list[[iter]])
#   alpha_mat_1 <- alpha_mat_1 + out_res$alpha_list_1[[iter]][order_tmp, order_tmp]
#   alpha_mat_2 <- alpha_mat_2 + out_res$alpha_list_2[[iter]][order_tmp, order_tmp]
#   A_mat_1 <- A_mat_1 + out_res$A_list_1[[iter]][order_tmp, order_tmp]
#   A_mat_2 <- A_mat_2 + out_res$A_list_2[[iter]][order_tmp, order_tmp]
# }
#
# alpha_mat_1 <- alpha_mat_1 / 5000
# alpha_mat_2 <- alpha_mat_2 / 5000
# A_mat_1 <- A_mat_1 / 5000
# A_mat_2 <- A_mat_2 / 5000
#
# # Calculate the error
# ## data set 1
# adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
# adj_1 <- t(adj_1)
# eval_fun(adj_1, adj_true1, g_true1)
#
# ## data set 2
# adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
# adj_2 <- t(adj_2)
# eval_fun(adj_2, adj_true2, g_true2)

########################### Do parallel ##################################
#### generate graph
set.seed(2021)
source("Two_dataset_new/Graph_MCMC_two_para.R")
n_graph <- 20
graph_sim <- graph_generation(K = K, n_graph = n_graph, p = p, n_tol = n_tol)
prior_vec_list <- list()
prior_vec_list[[1]] <- c(1 / (2 * p^1.5), 1 / p^2)
prior_vec_list[[2]] <- c(1 / p^1.5, 1 / p^1.5)
prior_vec_list[[3]] <- c(1 / p^1.5, 1 / p^2)
scale_x <- FALSE
intercept <- TRUE
iter_max <- 50000

library(foreach)
library(doParallel)
library(doRNG)

for (iter_prior in seq_len(prior_vec_list)) {
  out_res <- list()
  prior_vec <- prior_vec_list[[iter_prior]]
  ## do parallel 
  cl <- makeCluster(20)
  registerDoParallel(cl)
  out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
    library(pcalg)
    dta_1 <- graph_sim$X[[iter]][[1]]
    dta_2 <- graph_sim$X[[iter]][[2]]
    # get order
    dta <- rbind(dta_1, dta_2)
    score_ges <- new("GaussL0penObsScore", data = dta, intercept = intercept)
    ges_fit <- ges(score_ges)
    ges_adj <- as(ges_fit$repr, "matrix")
    ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
    graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
    order_int <- as.numeric(igraph::topo_sort(graph_i))
    # Do MCMC
    Graph_MCMC_two(dta_1, dta_2, scale_x = scale_x, intercept = intercept,
                   order_int = order_int, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
                   prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
                   left_size = 5000
    )
  }
  stopCluster(cl)
  
  ## check results
  res_1 <- matrix(NA, nrow = n_graph, ncol = 2)
  res_2 <- matrix(NA, nrow = n_graph, ncol = 2)
  for (iter_graph in seq_len(n_graph)) {
    res_tmp <- out_res[[iter_graph]]
    # analysis
    alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
    alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
    A_mat_1 <- matrix(0, nrow = p, ncol = p)
    A_mat_2 <- matrix(0, nrow = p, ncol = p)
    for (iter in seq_len(5000)) {
      order_tmp <- order(res_tmp$order_list[[iter]])
      alpha_mat_1 <- alpha_mat_1 + res_tmp$alpha_list_1[[iter]][order_tmp, order_tmp]
      alpha_mat_2 <- alpha_mat_2 + res_tmp$alpha_list_2[[iter]][order_tmp, order_tmp]
      A_mat_1 <- A_mat_1 + res_tmp$A_list_1[[iter]][order_tmp, order_tmp]
      A_mat_2 <- A_mat_2 + res_tmp$A_list_2[[iter]][order_tmp, order_tmp]
    }
    alpha_mat_1 <- alpha_mat_1 / 5000
    alpha_mat_2 <- alpha_mat_2 / 5000
    A_mat_1 <- A_mat_1 / 5000
    A_mat_2 <- A_mat_2 / 5000
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
      check_edge(adj_true1, adj_1)
    )
    res_2[iter_graph, ] <- c(
      shd(g_true2, g_2),
      check_edge(adj_true2, adj_2)
    )
  }
  cat(
    "prior:", prior_vec,
    "data1:", round(colMeans(res_1), 4),
    "data2:", round(colMeans(res_2), 4), "\n"
  )
  
}

# intercept FALSE prior: 5e-04 1e-04 data1: 31.2 20.75 data2: 29.4 18.45 
# intercept FALSE prior: 0.001 0.001 data1: 37.8 22.65 data2: 37.1 22.05 
# intercept FALSE prior: 0.001 1e-04 data1: 34.8 22.2 data2: 33.55 20.9 