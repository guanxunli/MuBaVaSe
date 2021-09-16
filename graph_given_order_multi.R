# ## define parameters
# p <- 100
# n <- 75
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
# 
# dta_list <- list()
# dta_list[[1]] <- dta_1
# dta_list[[2]] <- dta_2

## joint inference
# dta_list are n x p data set
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# prior_vec : prior for different models
# tau is the prior power for null model 1 / (p^tau)
# itermax is the maximum iteration
# L_max is the maximum number of parents
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l
# residual_variance_lowerbound is the lower bound for sigma2

joint_graph_multi <- function(dta_list, sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL, 
                              tau = 1.5, itermax = 100, L_max = 10, tol = 1e-4, sigma0_low_bd = 1e-8,
                              residual_variance_lowerbound = NULL) {
  ## Initialization
  K <- length(dta_list)
  p <- ncol(dta_list[[1]])
  n <- nrow(dta_list[[1]])
  ## save list
  alpha_list <- list()
  A_list <- list()
  Xb_list <- list()
  Y_list <- list()
  mean_list <- list()
  for (iter_K in seq_len(K)) {
    alpha_list[[iter_K]] <- matrix(0, nrow = p, ncol = p)
    A_list[[iter_K]] <- matrix(0, nrow = p, ncol = p)
    Xb_list[[iter_K]] <- matrix(0, nrow = n, ncol = p)
    Y_list[[iter_K]] <- dta_list[[iter_K]][, 1]
    mean_list[[iter_K]] <- mean(Y_list[[iter_K]])
    Xb_list[[iter_K]][, 1] <- mean_list[[iter_K]] 
  }
  # log likelihood
  llike_mat <- matrix(NA, nrow = p, ncol = K)
  sigma2_vec <- rep(NA, p)
  sigma2_vec[1] <- var(unlist(Y_list))
  for (iter_K in seq_len(K)) {
    llike_mat[1, iter_K] <- sum(dnorm(dta_list[[iter_K]][, 1], mean = mean_list[[iter_K]], 
                                 sd = sqrt(sigma2_vec[1]), log = TRUE))
  }
  ## load variable selection function
  source("sum_single_effect_mult.R")
  # begin iteration
  for (iter_p in seq_len(p - 1)) {
    dta_vs_list <- list()
    for (iter_K in seq_len(K)) {
      dta_vs_list[[iter_K]] <- list()
      dta_vs_list[[iter_K]]$X <- dta_list[[iter_K]][, seq_len(iter_p), drop = FALSE]
      dta_vs_list[[iter_K]]$Y <- dta_list[[iter_K]][, iter_p + 1]
    }
    ## variable selection
    res_vs <- sum_single_effect_mult(dta_vs_list, sigma02_int = sigma02_int, sigma2_int = sigma2_int, 
                                     tau = tau, prior_vec = prior_vec, L = min(iter_p, L_max), itermax = itermax, 
                                     tol = tol, sigma0_low_bd = sigma0_low_bd,
                                     residual_variance_lowerbound = residual_variance_lowerbound)
    ## save needed list
    sigma2_vec[iter_p + 1] <- res_vs$sigma2
    for (iter_K in seq_len(K)) {
      alpha_list[[iter_K]][iter_p + 1, seq_len(iter_p)] <- res_vs$res[[iter_K]]$alpha
      A_list[[iter_K]][iter_p + 1, seq_len(iter_p)] <- res_vs$res[[iter_K]]$post_mean
      # calculate the likelihood
      Xb_list[[iter_K]][, iter_p + 1] <- res_vs$res[[iter_K]]$Xb
      llike_mat[iter_p + 1, iter_K] <- sum(dnorm(x = dta_vs_list[[iter_K]]$Y, mean = res_vs$res[[iter_K]]$Xb, 
                                                 sd = res_vs$sigma2, log = TRUE))
    }
  }
  ## return results
  return(list(alpha_list = alpha_list, A_list = A_list, Xb_list = Xb_list, 
              llike_mat = llike_mat, sigma2_vec = sigma2_vec))
}

# ################## check results with two data sets ##################
# time1 <- Sys.time()
# res_multi <- joint_graph_multi(dta_list = dta_list, prior_vec = c(0.05, 0.05, 0.2))
# Sys.time() - time1
# 
# source("Two_dataset/Graph_given_order_two.R")
# time1 <- Sys.time()
# res_two <- joint_graph_fun_two(dta_1, dta_2, r = 0.2, q = 0.05)
# Sys.time() - time1
# 
# sum((res_two$alpha_res_1 - res_multi$alpha_list[[1]])^2)
# sum((res_two$llike_2_vec - res_multi$llike_mat[, 2])^2)
# sum((res_two$Xb_mat_1 - res_multi$Xb_list[[1]])^2)
