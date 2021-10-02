# ## define parameters
# p <- 100
# n <- 300
# p_c <- 100
# p_1 <- 20
# p_2 <- 20
# sigma <- 1
# sigma0 <- 0.6
# A1 <- matrix(0, nrow = p, ncol = p)
# A2 <- matrix(0, nrow = p, ncol = p)
# set.seed(202108)
# # Define the true graph given order
# index_c <- sample(seq_len(p * (p - 1) / 2), size = p_c, replace = FALSE)
# index_1 <- sample(setdiff(seq_len(p * (p - 1) / 2), index_c), size = p_1, replace = FALSE)
# index_2 <- sample(setdiff(seq_len(p * (p - 1) / 2), c(index_1, index_c)), size = p_2, replace = FALSE)
# 
# A1[lower.tri(A1)][c(index_c, index_1)] <-  rnorm(p_c + p_1, mean = 0, sd = sigma0)
# A2[lower.tri(A2)][c(index_c, index_2)] <-  rnorm(p_c + p_2, mean = 0, sd = sigma0)
# 
# alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
# alpha_mat_1[lower.tri(alpha_mat_1)][c(index_c, index_1)] <- 1
# alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
# alpha_mat_2[lower.tri(alpha_mat_2)][c(index_c, index_2)] <- 1
# 
# eps_1 <- matrix(rnorm(p * n), nrow = p, ncol = n)
# dta_1 <- solve(diag(1, nrow = p) - A1, eps_1)
# dta_1 <- t(dta_1)
# eps_2 <- matrix(rnorm(p * n), nrow = p, ncol = n)
# dta_2 <- solve(diag(1, nrow = p) - A2, eps_2)
# dta_2 <- t(dta_2)

## joint inference
# dta_1 and dta_2 are n x p data set
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# prior_vecr is prior for common part and for single part
# itermax is the maximum iteration
# L_max is the largest number of parents
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l
# residual_variance_lowerbound is the lower bound for sigma2

## load variable selection function
source("Two_dataset_new/sum_single_effect_two.R")
sigma02_int = NULL
sigma2_int = NULL
prior_vec <- NULL
itermax = 100
L_max = 10
tol = 1e-4
sigma0_low_bd = 1e-8
residual_variance_lowerbound = NULL
joint_graph_fun_two <- function(dta_1, dta_2, sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL, 
                                itermax = 100, L_max = 10, tol = 1e-4, sigma0_low_bd = 1e-8,
                                residual_variance_lowerbound = NULL) {
  ## Initialization
  p <- ncol(dta_1)
  n <- nrow(dta_1)
  ## define prior vector
  if (is.null(prior_vec)) {
    prior_vec <- c(1 / (6 * p^1.5), 2 / (3 * p ^ 1.5))
  }
  ## save matrix
  # probability of the edge exists
  alpha_res_1 <- matrix(0, nrow = p, ncol = p)
  alpha_res_2 <- matrix(0, nrow = p, ncol = p)
  # posterior mean of the edge
  A_res_1 <- matrix(0, nrow = p, ncol = p)
  A_res_2 <- matrix(0, nrow = p, ncol = p)
  # predicted value
  Xb_mat_1 <- matrix(NA, nrow = n, ncol = p)
  Xb_mat_2 <- matrix(NA, nrow = n, ncol = p)
  # log likelihood  
  llike_1_vec <- rep(NA, p)
  llike_2_vec <- rep(NA, p)
  llike_penalty_vec <- rep(0, p)
  sigma2_vec <- rep(NA, p)
  mean_1 <- mean(dta_1[, 1])
  mean_2 <- mean(dta_2[, 1])
  Xb_mat_1[, 1] <- rep(mean_1, n)
  Xb_mat_2[, 1] <- rep(mean_2, n)
  sigma2_vec[1] <- var(c(dta_1[, 1], dta_2[, 1]))
  llike_1_vec[1] <- sum(dnorm(dta_1[, 1], mean = mean_1, sd = sqrt(sigma2_vec[1]), log = TRUE))
  llike_2_vec[1] <- sum(dnorm(dta_2[, 1], mean = mean_2, sd = sqrt(sigma2_vec[1]), log = TRUE))
  # begin iteration
  for (iter_p in seq_len(p - 1)) {
    X_1 <- dta_1[, seq_len(iter_p), drop = FALSE]
    Y_1 <- dta_1[, iter_p + 1]
    X_2 <- dta_2[, seq_len(iter_p), drop = FALSE]
    Y_2 <- dta_2[, iter_p + 1]
    ## variable selection
    res <- sum_single_effect_two(X_1 = X_1, Y_1 = Y_1, X_2 = X_2, Y_2 = Y_2, sigma02_int = sigma02_int,
                                 sigma2_int = sigma2_int, prior_vec = prior_vec, L = min(iter_p, L_max), 
                                 itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd, 
                                 residual_variance_lowerbound = residual_variance_lowerbound)
    # save the matrix we want
    alpha_res_1[iter_p + 1, seq_len(iter_p)] <- res$alpha_1
    alpha_res_2[iter_p + 1, seq_len(iter_p)] <- res$alpha_2
    A_res_1[iter_p + 1, seq_len(iter_p)] <- res$post_mean1
    A_res_2[iter_p + 1, seq_len(iter_p)] <- res$post_mean2
    # calculate the likelihood
    sigma2_vec[iter_p + 1] <- res$sigma2
    Xb_mat_1[, iter_p + 1] <- res$Xb_1
    Xb_mat_2[, iter_p + 1] <- res$Xb_2
    llike_1_vec[iter_p + 1] <- sum(dnorm(x = Y_1, mean = res$Xb_1, sd = sqrt(res$sigma2), log = TRUE))
    llike_2_vec[iter_p + 1] <- sum(dnorm(x = Y_2, mean = res$Xb_2, sd = sqrt(res$sigma2), log = TRUE))
    llike_penalty_vec[iter_p + 1] <- sum(-1.5 * res$alpha * c(rep(prior_vec[1], 2 * iter_p), rep(prior_vec[2], iter_p)))
  }
  ## return results
  return(list(alpha_res_1 = alpha_res_1, alpha_res_2 = alpha_res_2, A_res_1 = A_res_1, A_res_2 = A_res_2, 
              llike_1_vec = llike_1_vec, llike_2_vec = llike_2_vec, llike_penalty_vec = llike_penalty_vec,
              Xb_mat_1 = Xb_mat_1, Xb_mat_2 = Xb_mat_2, sigma2_vec = sigma2_vec))
}

# ################## check results with GES ##################
# res_joint <- joint_graph_fun_two(dta_1 = dta_1, dta_2 = dta_2, r = 0.2, q = 0.05, tau = 1.5)
# ## remove order edge
# check_edge <- function(adj_pre, adj_act) {
#   adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
#   adj_act <- ceiling((adj_act + t(adj_act)) / 2)
#   return(sum(abs(adj_pre - adj_act)) / 2)
# }
# library(pcalg)
# ######## data set 1
# #### Define true
# adj_true_1 <- t(alpha_mat_1)
# g_true_1 <- as(adj_true_1, "graphNEL")
# weight_true_1 <- t(A1)
# #### our method
# adj_1 <- res_joint$alpha_res_1
# adj_1 <- t(ifelse(adj_1 > 0.5, 1, 0))
# g_1 <- as(adj_1, "graphNEL")
# weight_1 <- t(res_joint$A_res_1)
# #### GES method
# score1 <- new("GaussL0penObsScore", data = dta_1, intercept = FALSE) # lambda = sqrt(2 * log(p) / n)
# ges_fit1 <- ges(score1) #  2.2 mins
# ges_adj1 <- as(ges_fit1$repr, "matrix")
# ges_adj1 <- ifelse(ges_adj1 == TRUE, 1, 0)
# ges_graph1 <- as(ges_fit1$repr, "graphNEL")
# ges_weight1 <- ges_fit1$repr$weight.mat()
# #### check results
# # structural Hamming distance (SHD)
# shd(g_true_1, ges_graph1)
# shd(g_true_1, g_1)
# # structural Hamming distance (SHD)
# check_edge(ges_adj1, adj_true_1)
# check_edge(adj_1, adj_true_1)
# # Mean square error for weight
# sum((weight_true_1 - ges_weight1)^2)
# sum((weight_true_1 - weight_1)^2)
# 
# ######## data set 2
# #### Define true
# adj_true_2 <- t(alpha_mat_2)
# g_true_2 <- as(adj_true_2, "graphNEL")
# weight_true_2 <- t(A2)
# #### our method
# adj_2 <- res_joint$alpha_res_2
# adj_2 <- t(ifelse(adj_2 > 0.5, 1, 0))
# g_2 <- as(adj_2, "graphNEL")
# weight_2 <- t(res_joint$A_res_2)
# #### GES method
# score2 <- new("GaussL0penObsScore", data = dta_2, intercept = FALSE)
# ges_fit2 <- ges(score2)
# ges_adj2 <- as(ges_fit2$repr, "matrix")
# ges_adj2 <- ifelse(ges_adj2 == TRUE, 1, 0)
# ges_graph2 <- as(ges_fit2$repr, "graphNEL")
# ges_weight2 <- ges_fit2$repr$weight.mat()
# #### check results
# # structural Hamming distance (SHD)
# shd(g_true_2, ges_graph2)
# shd(g_true_2, g_2)
# # structural Hamming distance (SHD)
# check_edge(ges_adj2, adj_true_2)
# check_edge(adj_2, adj_true_2)
# # Mean square error for weight
# sum((weight_true_2 - ges_weight2)^2)
# sum((weight_true_2 - weight_2)^2)
