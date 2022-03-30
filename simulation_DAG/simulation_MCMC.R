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
e_com <- 100
e_pri <- 30
# Define prior
prior_vec_list <- list()
prior_vec_list[[1]] <- c(1 / (2 * p^1.25), 1 / p^1.5)
prior_vec_list[[2]] <- c(1 / (2 * p^1.5), 1 / p^2)
prior_vec_list[[3]] <- c(1 / p^2, 1 / p^2.25)
prior_vec_list[[4]] <- c(1 / (2 * p^2), 1 / p^2.25)

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

# ########################### Do one figure ##################################
# #### generate graph
# set.seed(2021)
# n_graph <- 1
# graph_sim <- graph_generation(
#   K = K, n_graph = n_graph, p = p, n_tol = n_tol,
#   e_com = e_com, e_pri = e_pri
# )
# adj_true1 <- t(graph_sim$G[[1]][[1]])
# g_true1 <- as(getGraph(adj_true1), "graphNEL")
# weight_true1 <- t(graph_sim$A[[1]][[1]])
# adj_true2 <- t(graph_sim$G[[1]][[2]])
# g_true2 <- as(getGraph(adj_true2), "graphNEL")
# weight_true2 <- t(graph_sim$A[[1]][[2]])
# 
# #### our method
# source("Two_dataset_new/Graph_given_order_two.R")
# dta_1 <- graph_sim$X[[1]][[1]]
# dta_2 <- graph_sim$X[[1]][[2]]
# 
# #### If we know the order
# for (iter_prior in seq_len(length(prior_vec_list))) {
#   prior_vec <- prior_vec_list[[iter_prior]]
#   out_res <- joint_graph_fun_two(dta_1 = dta_1, dta_2 = dta_2, prior_vec = prior_vec,
#                                  scale_x = scale_x, intercept = intercept)
#   print(round(sum(out_res$llike_1_vec + out_res$llike_2_vec + out_res$llike_penalty_vec), 4))
#   ## Calculate the error
#   ## data set 1
#   adj_1 <- ifelse(out_res$alpha_res_1 > 0.5, 1, 0)
#   adj_1 <- t(adj_1)
#   g_1 <- as(getGraph(adj_1), "graphNEL")
#   cat(
#     "prior = ", prior_vec[1], prior_vec[2],
#     c(shd(g_true1, g_1), check_edge(adj_true1, adj_1)),
#     "TP", round(TPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4),
#     "FP", round(FPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4),
#     "FN", round(FNrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4),
#     "L2", round(check_adj_l2(adj_pre = out_res$alpha_res_1, adj_act = adj_true1), 4),
#     "L1", round(check_adj_l1(adj_pre = out_res$alpha_res_1, adj_act = adj_true1), 4),
#     "\n"
#   )
#   ## data set 2
#   adj_2 <- ifelse(out_res$alpha_res_2 > 0.5, 1, 0)
#   adj_2 <- t(adj_2)
#   g_2 <- as(getGraph(adj_2), "graphNEL")
#   cat(
#     "prior = ", prior_vec[1], prior_vec[2],
#     c(shd(g_true2, g_2), check_edge(adj_true2, adj_2)),
#     "TP", round(TPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4),
#     "FP", round(FPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4),
#     "FN", round(FNrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4),
#     "L2", round(check_adj_l2(adj_pre = out_res$alpha_res_2, adj_act = adj_true2), 4),
#     "L1", round(check_adj_l1(adj_pre = out_res$alpha_res_2, adj_act = adj_true2), 4),
#     "\n"
#   )
# }
# 
# ########################## Do MCMC quick test ############################
# iter_max <- 200
# prior_vec <- prior_vec_list[[2]]
# source("Two_dataset_new/Graph_MCMC_two.R")
# source("Two_dataset_new/Graph_MCMC_two_sim.R")
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
# out_res <- Graph_MCMC_two_sim(dta_1, dta_2, scale_x = scale_x, intercept = intercept,
#                               order_int = order_int, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
#                               prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
#                               burn_in = 1, adj_true1 = adj_true1, adj_true2 = adj_true2
# )
# 
# # Show likelihood
# library(ggplot2)
# library(gridExtra)
# gl <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$llike_vec)) +
#   xlab("Iteration") + ylab("Log likelihood")
# gs_1 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat1[1, ])) +
#   xlab("Iteration") + ylab("SHD")
# gu_1 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat1[2, ])) +
#   xlab("Iteration") + ylab("No order")
# gs_2 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat2[1, ])) +
#   xlab("Iteration") + ylab("SHD")
# gu_2 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat2[2, ])) +
#   xlab("Iteration") + ylab("No order")
# layout_matrix <- matrix(c(1, 2, 3), nrow = 3)
# grid.arrange(gl, gs_1, gu_1, layout_matrix = layout_matrix)
# grid.arrange(gl, gs_2, gu_2, layout_matrix = layout_matrix)
#
# ## analysis
# alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
# alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
# A_mat_1 <- matrix(0, nrow = p, ncol = p)
# A_mat_2 <- matrix(0, nrow = p, ncol = p)
# for (iter in seq_len(length(out_res[[1]]))) {
#   order_tmp <- order(out_res$order_list[[iter]])
#   alpha_mat_1 <- alpha_mat_1 + out_res$alpha_list_1[[iter]][order_tmp, order_tmp]
#   alpha_mat_2 <- alpha_mat_2 + out_res$alpha_list_2[[iter]][order_tmp, order_tmp]
#   A_mat_1 <- A_mat_1 + out_res$A_list_1[[iter]][order_tmp, order_tmp]
#   A_mat_2 <- A_mat_2 + out_res$A_list_2[[iter]][order_tmp, order_tmp]
# }
# alpha_mat_1 <- alpha_mat_1 / length(out_res[[1]])
# alpha_mat_2 <- alpha_mat_2 / length(out_res[[1]])
# A_mat_1 <- A_mat_1 / length(out_res[[1]])
# A_mat_2 <- A_mat_2 / length(out_res[[1]])
# ## data set 1
# adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
# adj_1 <- t(adj_1)
# g_1 <- as(getGraph(adj_1), "graphNEL")
# cat(
#   "prior = ", prior_vec[1], prior_vec[2],
#   c(shd(g_true1, g_1), check_edge(adj_true1, adj_1)),
#   "TP", round(TPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4),
#   "FP", round(FPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4),
#   "FN", round(FNrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4),
#   "L2", round(check_adj_l2(adj_pre = alpha_mat_1, adj_act = adj_true1), 4),
#   "L1", round(check_adj_l1(adj_pre = alpha_mat_1, adj_act = adj_true1), 4),
#   "\n"
# )
# ## data set 2
# adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
# adj_2 <- t(adj_2)
# g_2 <- as(getGraph(adj_2), "graphNEL")
# cat(
#   "prior = ", prior_vec[1], prior_vec[2],
#   c(shd(g_true2, g_2), check_edge(adj_true2, adj_2)),
#   "TP", round(TPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4),
#   "FP", round(FPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4),
#   "FN", round(FNrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4),
#   "L2", round(check_adj_l2(adj_pre = alpha_mat_2, adj_act = adj_true2), 4),
#   "L1", round(check_adj_l1(adj_pre = alpha_mat_2, adj_act = adj_true2), 4),
#   "\n"
# )
#
# ##### without GES initialization
# out_res <- Graph_MCMC_two_sim(dta_1, dta_2, scale_x = scale_x, intercept = intercept,
#                           order_int = NULL, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
#                           prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
#                           burn_in = 1, adj_true1 = adj_true1, adj_true2 = adj_true2
# )
# gl <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$llike_vec)) +
#   xlab("Iteration") + ylab("Log likelihood")
# gs_1 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat1[1, ])) +
#   xlab("Iteration") + ylab("SHD")
# gu_1 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat1[2, ])) +
#   xlab("Iteration") + ylab("No order")
# gs_2 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat2[1, ])) +
#   xlab("Iteration") + ylab("SHD")
# gu_2 <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$error_mat2[2, ])) +
#   xlab("Iteration") + ylab("No order")
# layout_matrix <- matrix(c(1, 2, 3), nrow = 3)
# grid.arrange(gl, gs_1, gu_1, layout_matrix = layout_matrix)
# grid.arrange(gl, gs_2, gu_2, layout_matrix = layout_matrix)
#
# # analysis
# alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
# alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
# A_mat_1 <- matrix(0, nrow = p, ncol = p)
# A_mat_2 <- matrix(0, nrow = p, ncol = p)
# for (iter in seq_len(length(out_res[[1]]))) {
#   order_tmp <- order(out_res$order_list[[iter]])
#   alpha_mat_1 <- alpha_mat_1 + out_res$alpha_list_1[[iter]][order_tmp, order_tmp]
#   alpha_mat_2 <- alpha_mat_2 + out_res$alpha_list_2[[iter]][order_tmp, order_tmp]
#   A_mat_1 <- A_mat_1 + out_res$A_list_1[[iter]][order_tmp, order_tmp]
#   A_mat_2 <- A_mat_2 + out_res$A_list_2[[iter]][order_tmp, order_tmp]
# }
# alpha_mat_1 <- alpha_mat_1 / length(out_res[[1]])
# alpha_mat_2 <- alpha_mat_2 / length(out_res[[1]])
# A_mat_1 <- A_mat_1 / length(out_res[[1]])
# A_mat_2 <- A_mat_2 / length(out_res[[1]])
# ## data set 1
# adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
# adj_1 <- t(adj_1)
# g_1 <- as(getGraph(adj_1), "graphNEL")
# cat(
#   "prior = ", prior_vec[1], prior_vec[2],
#   c(shd(g_true1, g_1), check_edge(adj_true1, adj_1)),
#   "TP", round(TPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4),
#   "FP", round(FPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4),
#   "FN", round(FNrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4),
#   "L2", round(check_adj_l2(adj_pre = alpha_mat_1, adj_act = adj_true1), 4),
#   "L1", round(check_adj_l1(adj_pre = alpha_mat_1, adj_act = adj_true1), 4),
#   "\n"
# )
# ## data set 2
# adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
# adj_2 <- t(adj_2)
# g_2 <- as(getGraph(adj_2), "graphNEL")
# cat(
#   "prior = ", prior_vec[1], prior_vec[2],
#   c(shd(g_true2, g_2), check_edge(adj_true2, adj_2)),
#   "TP", round(TPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4),
#   "FP", round(FPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4),
#   "FN", round(FNrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4),
#   "L2", round(check_adj_l2(adj_pre = alpha_mat_2, adj_act = adj_true2), 4),
#   "L1", round(check_adj_l1(adj_pre = alpha_mat_2, adj_act = adj_true2), 4),
#   "\n"
# )

########################### Do parallel ##################################
#### generate graph
source("Two_dataset_new/Graph_MCMC_two_sim.R")
set.seed(2021)
n_graph <- 20
graph_sim <- graph_generation(
  K = K, n_graph = n_graph, p = p, n_tol = n_tol,
  e_com = e_com, e_pri = e_pri
)
prior_penalty <- TRUE
iter_max <- 1e5

library(foreach)
library(doParallel)
library(doRNG)

for (iter_prior in seq_len(length(prior_vec_list))) {
  out_res <- list()
  prior_vec <- prior_vec_list[[iter_prior]]
  ## do parallel
  cl <- makeCluster(20)
  registerDoParallel(cl)
  set.seed(2021)
  out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
    library(pcalg)
    dta_1 <- graph_sim$X[[iter]][[1]]
    dta_2 <- graph_sim$X[[iter]][[2]]
    adj_true1 <- t(graph_sim$G[[iter]][[1]])
    adj_true2 <- t(graph_sim$G[[iter]][[2]])
    # get order
    dta <- rbind(dta_1, dta_2)
    score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
    ges_fit <- ges(score_ges)
    ges_adj <- as(ges_fit$repr, "matrix")
    ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
    graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
    order_int <- as.numeric(igraph::topo_sort(graph_i))
    # Do MCMC
    Graph_MCMC_two_sim(dta_1, dta_2,
                       prior_penalty = prior_penalty,
                       scale_x = scale_x, intercept = intercept,
                       order_int = order_int, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
                       prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
                       burn_in = iter_max - 5000,
                       adj_true1 = adj_true1, adj_true2 = adj_true2
    )
  }
  stopCluster(cl)
  
  library(ggplot2)
  library(gridExtra)
  ## check results
  res_1 <- matrix(NA, nrow = n_graph, ncol = 7)
  res_2 <- matrix(NA, nrow = n_graph, ncol = 7)
  layout_matrix <- matrix(c(1, 2, 3), nrow = 3)
  alpha_mat_list1 <- list()
  alpha_mat_list2 <- list()
  for (iter_graph in seq_len(n_graph)) {
    res_tmp <- out_res[[iter_graph]]
    # ## plot likelihood and error
    # if (iter_graph %% 10 == 0) {
    #   gl <- ggplot() +
    #     geom_line(aes(x = seq_len(iter_max), y = res_tmp$llike_vec)) +
    #     xlab("Iteration") +
    #     ylab("Log likelihood")
    #   gs_1 <- ggplot() +
    #     geom_line(aes(x = seq_len(iter_max), y = res_tmp$error_mat1[1, ])) +
    #     xlab("Iteration") +
    #     ylab("SHD")
    #   gu_1 <- ggplot() +
    #     geom_line(aes(x = seq_len(iter_max), y = res_tmp$error_mat1[2, ])) +
    #     xlab("Iteration") +
    #     ylab("No order")
    #   gs_2 <- ggplot() +
    #     geom_line(aes(x = seq_len(iter_max), y = res_tmp$error_mat2[1, ])) +
    #     xlab("Iteration") +
    #     ylab("SHD")
    #   gu_2 <- ggplot() +
    #     geom_line(aes(x = seq_len(iter_max), y = res_tmp$error_mat2[2, ])) +
    #     xlab("Iteration") +
    #     ylab("No order")
    #   png(paste0(
    #     "pri", prior_vec[1], "com", prior_vec[2], "graph", iter_graph,
    #     "e_com", e_com, "e_pri", e_pri, "penalty", prior_penalty, "data1_MCMC.png"
    #   ))
    #   grid.arrange(gl, gs_1, gu_1, layout_matrix = layout_matrix)
    #   dev.off()
    #   png(paste0(
    #     "pri", prior_vec[1], "com", prior_vec[2], "graph", iter_graph,
    #     "e_com", e_com, "e_pri", e_pri, "penalty", prior_penalty, "data2_MCMC.png"
    #   ))
    #   grid.arrange(gl, gs_2, gu_2, layout_matrix = layout_matrix)
    #   dev.off()
    # }
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
    alpha_mat_list1[[iter_graph]] <- alpha_mat_1
    alpha_mat_list2[[iter_graph]] <- alpha_mat_2
  }
  # show results
  cat(
    "p:", p, "e_com:", e_com, "e_pri", e_pri, "prior:", round(prior_vec, 4),
    "data1:", round(colMeans(res_1), 4),
    "data2:", round(colMeans(res_2), 4), "\n"
  )
  # ## show tex
  # cat(
  #   "$", prior_vec, "$", "&",
  #   "data1:", "&", round(colMeans(res_1), 4), "&",
  #   "data2:", "&", round(colMeans(res_2), 4), "\\\\\n"
  # )
  # saveRDS(
  #   alpha_mat_list1,
  #   paste0(prior_vec[1], "_", prior_vec[2], "_", e_com, "_", e_pri, "_", "alpha_mat1_MCMC.rds")
  # )
  # saveRDS(
  #   alpha_mat_list2,
  #   paste0(prior_vec[1], "_", prior_vec[2], "_", e_com, "_", e_pri, "_", "alpha_mat2_MCMC.rds")
  # )
}