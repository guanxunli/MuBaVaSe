# ## define parameters
# p <- 100
# n1 <- 300
# n2 <- 400
# p_c <- 100
# p_1 <- 30
# p_2 <- 25
# sigma <- 1
# sigma0 <- 0.6
# A1 <- matrix(0, nrow = p, ncol = p)
# A2 <- matrix(0, nrow = p, ncol = p)
# set.seed(2021)
# # Define the true graph given order
# index_c <- sample(seq_len(p * (p - 1) / 2), size = p_c, replace = FALSE)
# index_1 <- sample(setdiff(seq_len(p * (p - 1) / 2), index_c), size = p_1, replace = FALSE)
# index_2 <- sample(setdiff(seq_len(p * (p - 1) / 2), index_c), size = p_2, replace = FALSE)
#
# A1[lower.tri(A1)][c(index_c, index_1)] <-  rnorm(p_c + p_1, mean = 0, sd = sigma0)
# A2[lower.tri(A2)][c(index_c, index_2)] <-  rnorm(p_c + p_2, mean = 0, sd = sigma0)
#
# alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
# alpha_mat_1[lower.tri(alpha_mat_1)][c(index_c, index_1)] <- 1
# alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
# alpha_mat_2[lower.tri(alpha_mat_2)][c(index_c, index_2)] <- 1
#
# eps_1 <- matrix(rnorm(p * n1), nrow = p, ncol = n1)
# dta_1 <- solve(diag(1, nrow = p) - A1, eps_1)
# dta_1 <- t(dta_1)
# eps_2 <- matrix(rnorm(p * n2), nrow = p, ncol = n2)
# dta_2 <- solve(diag(1, nrow = p) - A2, eps_2)
# dta_2 <- t(dta_2)

## joint inference
# dta_1 and dta_2 are n x p data set
# scale_x : scale the data
# intercept: calculate the mean of Y
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# prior_vec is prior for common part and for single part
# itermax is the maximum iteration
# L_max is the largest number of parents
# tol is the threshold for ELBO
# residual_variance_lowerbound is the lower bound for sigma2

## load variable selection function
source("Two_dataset_v3/sampling/sum_single_effect_two_sampling.R")
joint_graph_fun_two_sampling <- function(dta_1, dta_2, scale_x = FALSE, intercept = TRUE,
                                         sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL,
                                         itermax = 100, L_max = 10, tol = 1e-4,
                                         residual_variance_lowerbound = NULL) {
  ## Initialization
  p <- ncol(dta_1)
  if (p != ncol(dta_2)) stop("The number of features should be same!")
  n1 <- nrow(dta_1)
  n2 <- nrow(dta_2)
  ## define prior vector
  if (is.null(prior_vec)) {
    prior_vec <- c(1 / (2 * p^1.5), 1 / (p^2))
  }
  lprior_vec <- log(prior_vec)
  # probability of the edge exists
  alpha_list <- list()
  sigma02_vec_list <- list()
  sigma2_vec <- rep(NA, p)
  # log likelihood
  llike_vec_1 <- rep(NA, p)
  llike_vec_2 <- rep(NA, p)
  # graph return
  graph_res_1 <- matrix(0, nrow = p, ncol = p)
  graph_res_2 <- matrix(0, nrow = p, ncol = p)
  # graph penalty
  lprior_graph <- rep(0, p)
  lpropose_graph <- rep(0, p)
  # initialization
  if (intercept) {
    mean_1 <- mean(dta_1[, 1])
    mean_2 <- mean(dta_2[, 1])
  } else {
    mean_1 <- 0
    mean_2 <- 0
  }
  alpha_list[[1]] <- NULL
  sigma02_vec_list[[1]] <- NULL
  sigma2_vec[1] <- var(c(dta_1[, 1], dta_2[, 1]))
  llike_vec_1[1] <- sum(dnorm(dta_1[, 1], mean = mean_1, sd = sqrt(sigma2_vec[1]), log = TRUE))
  llike_vec_2[1] <- sum(dnorm(dta_2[, 1], mean = mean_2, sd = sqrt(sigma2_vec[1]), log = TRUE))
  # begin iteration
  for (iter_p in seq_len(p - 1)) {
    X_1 <- dta_1[, seq_len(iter_p), drop = FALSE]
    Y_1 <- dta_1[, iter_p + 1]
    X_2 <- dta_2[, seq_len(iter_p), drop = FALSE]
    Y_2 <- dta_2[, iter_p + 1]
    ## variable selection
    res <- sum_single_effect_two_sampling(
      X_1 = X_1, Y_1 = Y_1, X_2 = X_2, Y_2 = Y_2,
      scale_x = scale_x, intercept = intercept, sigma02_int = sigma02_int,
      sigma2_int = sigma2_int, prior_vec = prior_vec, L = min(iter_p, L_max),
      itermax = itermax, tol = tol,
      residual_variance_lowerbound = residual_variance_lowerbound
    )
    alpha_list[[iter_p + 1]] <- res$alpha_mat
    sigma2_vec[iter_p + 1] <- res$sigma2
    sigma02_vec_list[[iter_p + 1]] <- res$sigma02_vec
    llike_vec_1[iter_p + 1] <- res$loglikelihood_1
    llike_vec_2[iter_p + 1] <- res$loglikelihood_2
    lprior_graph[iter_p + 1] <- res$lprior
    lpropose_graph[iter_p + 1] <- res$lpropose
    if (is.null(res$index1_select) == FALSE) graph_res_1[iter_p + 1, res$index1_select] <- 1
    if (is.null(res$index2_select) == FALSE) graph_res_2[iter_p + 1, res$index2_select] <- 1
  }
  ## return results
  return(list(
    alpha_list = alpha_list, sigma2_vec = sigma2_vec,
    graph_res_1 = graph_res_1, graph_res_2 = graph_res_2,
    llike_vec_1 = llike_vec_1, llike_vec_2 = llike_vec_2,
    lprior_graph = lprior_graph, lpropose_graph = lpropose_graph,
    sigma02_vec_list = sigma02_vec_list
  ))
}

# ################## check results with GES ##################
# set.seed(2022)
# time1 <- Sys.time()
# res_joint_new <- joint_graph_fun_two_sampling(dta_1 = dta_1, dta_2 = dta_2, scale_x = FALSE, intercept = TRUE)
# Sys.time() - time1
# source("Two_dataset_new/Graph_given_order_two.R")
# time1 <- Sys.time()
# res_joint <- joint_graph_fun_two(dta_1 = dta_1, dta_2 = dta_2, scale_x = FALSE, intercept = TRUE)
# Sys.time() - time1
# ## remove order edge
# library(pcalg)
# check_edge <- function(adj_pre, adj_act) {
#   adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
#   adj_act <- ceiling((adj_act + t(adj_act)) / 2)
#   return(sum(abs(adj_pre - adj_act)) / 2)
# }
# ######## data set 1
# #### Define true
# adj_true_1 <- t(alpha_mat_1)
# g_true_1 <- as(adj_true_1, "graphNEL")
# #### our method 1
# adj_1 <- res_joint$alpha_res_1
# adj_1 <- t(ifelse(adj_1 > 0.5, 1, 0))
# g_1 <- as(adj_1, "graphNEL")
# #### our method 2
# adj_1_new <- t(res_joint_new$graph_res_1)
# g_1_new <- as(adj_1_new, "graphNEL")
# #### check results
# # structural Hamming distance (SHD)
# shd(g_true_1, g_1_new)
# shd(g_true_1, g_1)
# # structural Hamming distance (SHD)
# check_edge(adj_1_new, adj_true_1)
# check_edge(adj_1, adj_true_1)
#
# ######## data set 2
# #### Define true
# adj_true_2 <- t(alpha_mat_2)
# g_true_2 <- as(adj_true_2, "graphNEL")
# #### our method 1
# adj_2 <- res_joint$alpha_res_2
# adj_2 <- t(ifelse(adj_2 > 0.5, 1, 0))
# g_2 <- as(adj_2, "graphNEL")
# #### our method 2
# adj_2_new <- t(res_joint_new$graph_res_2)
# g_2_new <- as(adj_2_new, "graphNEL")
# #### check results
# # structural Hamming distance (SHD)
# shd(g_true_2, g_2_new)
# shd(g_true_2, g_2)
# # structural Hamming distance (SHD)
# check_edge(adj_2_new, adj_true_2)
# check_edge(adj_2, adj_true_2)