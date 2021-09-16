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
# set.seed(2021)
# # Define the true graph given order
# index_c <- sample(seq_len(p * (p - 1) / 2), size = p_c, replace = FALSE)
# index_1 <- sample(setdiff(seq_len(p * (p - 1) / 2), index_c), size = p_1, replace = FALSE)
# index_2 <- sample(setdiff(seq_len(p * (p - 1) / 2), index_1), size = p_2, replace = FALSE)
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

## MCMC method for Graph
# dta_list are n x p data set
# order_int is the initialized order for nodes
# iter_max is the maximun mcmc step
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# prior_vec : prior for different models
# tau is the prior power for null model 1 / (p^tau)
# itermax is the maximum iteration
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l
Graph_MCMC_two <- function(dta_list, order_int = NULL, iter_max = 10000, sigma02_int = NULL, sigma2_int = NULL,
                           prior_vec = NULL, tau = 1.5, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 5000,
                           residual_variance_lowerbound = NULL) {
  ## Initialization
  K <- length(dta_list)
  p <- ncol(dta_list[[1]])
  n <- nrow(dta_list[[1]])
  # Initialize order
  if (is.null(order_int)) {
    order_old <- sample(seq_len(p), p)
  } else {
    order_old <- order_int
  }
  # change data
  dta_old_list <- list()
  for (iter_K in seq_len(K)) {
    dta_old_list[[iter_K]] <- dta_list[[iter_K]][, order_old]
  }
  ## load the main function
  source("graph_given_order_multi.R")
  res_old <- joint_graph_multi(dta_old_list, sigma02_int = sigma02_int, sigma2_int = sigma2_int, prior_vec = prior_vec, 
                               tau = tau, itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
                               residual_variance_lowerbound = residual_variance_lowerbound)
  ## save list
  alpha_res_old <- list()
  A_res_old <- list()
  sigma2_vec_old <- res_old$sigma2_vec
  llike_mat_old <- res_old$llike_mat
  llike_old <- sum(llike_mat_old)
  alpha_list <- list()
  A_list <- list()
  order_list <- list()
  for (iter_K in seq_len(K)) {
    alpha_res_old[[iter_K]] <- res_old$alpha_list[[iter_K]]
    A_res_old[[iter_K]] <- res_old$A_list[[iter_K]]
    alpha_list[[iter_K]] <- list()
    A_list[[iter_K]] <- list()
  }
  ## load the function
  source("sum_single_effect_mult.R")
  ## begin iteration
  dta_pro_list <- list()
  for (iter_MCMC in seq_len(iter_max)) {
    if (iter_MCMC %% 100 == 0) print(iter_MCMC)
    ## Initialize proposal
    for (iter_K in seq_len(K)) {
      dta_pro_list[[iter_K]] <- dta_old_list[[iter_K]]
    }
    order_pro <- order_old
    # log likelihood
    sigma2_vec_pro <- sigma2_vec_old
    llike_mat_pro <- llike_mat_old
    ## propose the new order
    pos_change <- sample(seq_len(p - 1), 1)
    llike_pro <- llike_old - sum(llike_mat_old[c(pos_change, pos_change + 1), ])
    for (iter_K in seq_len(K)) {
      dta_pro_list[[iter_K]][, c(pos_change, pos_change + 1)] <- dta_old_list[[iter_K]][, c(pos_change + 1, pos_change)]
    }
    order_pro[c(pos_change, pos_change + 1)] <- order_old[c(pos_change + 1, pos_change)]
    ## doing variable selection
    if (pos_change == 1) {
      res_pos <- list()
      res_pos$res <- list()
      Y_list <- list()
      for (iter_K in seq_len(K)) {
        Y_list[[iter_K]] <- dta_pro_list[[iter_K]][, 1]
        res_pos$res[[iter_K]] <- list()
        res_pos$res[[iter_K]]$alpha <- rep(0, p)
        res_pos$res[[iter_K]]$post_mean <- rep(0, p)
        res_pos$res[[iter_K]]$Xb <- rep(mean(dta_pro_list[[iter_K]][, 1]), n)
      }
      res_pos$sigma2 <- var(unlist(Y_list))
    } else{
      dta_tmp_list <- list()
      for (iter_K in seq_len(K)) {
        dta_tmp_list[[iter_K]] <- list()
        dta_tmp_list[[iter_K]]$X <- dta_pro_list[[iter_K]][, seq_len(pos_change - 1), drop = FALSE]
        dta_tmp_list[[iter_K]]$Y <- dta_pro_list[[iter_K]][, pos_change]
      }
      res_pos <- sum_single_effect_mult(dta_list = dta_tmp_list, sigma02_int = sigma02_int, 
                                        sigma2_int = sigma2_vec_old[pos_change + 1], prior_vec = prior_vec, 
                                        tau = tau, L = min(pos_change - 1, 10), itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
                                        residual_variance_lowerbound = residual_variance_lowerbound)
    }
    dta_tmp_list <- list()
    for (iter_K in seq_len(K)) {
      dta_tmp_list[[iter_K]] <- list()
      dta_tmp_list[[iter_K]]$X <-  dta_pro_list[[iter_K]][, seq_len(pos_change), drop = FALSE]
      dta_tmp_list[[iter_K]]$Y <- dta_pro_list[[iter_K]][, pos_change + 1]
    }
    res_pos1 <- sum_single_effect_mult(dta_list = dta_tmp_list, sigma02_int = sigma02_int, 
                                      sigma2_int = sigma2_vec_old[pos_change], prior_vec = prior_vec, 
                                      tau = tau, L = min(pos_change, 10), itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
                                      residual_variance_lowerbound = residual_variance_lowerbound)
    # likelihood
    sigma2_vec_pro[c(pos_change, pos_change + 1)] <- c(res_pos$sigma2, res_pos1$sigma2)
    for (iter_K in seq_len(K)) {
      llike_mat_pro[pos_change, iter_K] <- sum(dnorm(x = dta_pro_list[[iter_K]][, pos_change], 
                                                     mean = res_pos$res[[iter_K]]$Xb, sd = sqrt(res_pos$sigma2), log = TRUE))
      llike_mat_pro[pos_change + 1, iter_K] <- sum(dnorm(x = dta_pro_list[[iter_K]][, pos_change + 1], 
                                                     mean = res_pos1$res[[iter_K]]$Xb, sd = sqrt(res_pos1$sigma2), log = TRUE))
    }
    llike_pro <- llike_pro + sum(llike_mat_pro[c(pos_change, pos_change + 1), ])

    # accept or not
    if (llike_pro > llike_old) {
      accept <- TRUE
    } else {
      U <- sample(1)
      thres <- exp(llike_pro - llike_old)
      if (U < thres) {
        accept <- TRUE
      } else {
        accept <- FALSE
      }
    }

    if (accept) {
      for (iter_K in seq_len(K)) {
        alpha_res_old[[iter_K]][c(pos_change, pos_change + 1), ] <- 0
        A_res_old[[iter_K]][c(pos_change, pos_change + 1), ] <- 0
        alpha_res_old[[iter_K]][pos_change, seq_len(pos_change - 1)] <- res_pos$res[[iter_K]]$alpha
        alpha_res_old[[iter_K]][pos_change + 1, seq_len(pos_change)] <- res_pos1$res[[iter_K]]$alpha
        A_res_old[[iter_K]][pos_change,  seq_len(pos_change - 1)] <- res_pos$res[[iter_K]]$post_mean
        A_res_old[[iter_K]][pos_change + 1, seq_len(pos_change)] <- res_pos1$res[[iter_K]]$post_mean
      }
      # likelihood
      sigma2_vec_old <- sigma2_vec_pro
      llike_mat_old <- llike_mat_pro
      llike_old <- llike_pro
      # data and order
      dta_old_list <- dta_pro_list
      order_old <- order_pro
    }
    # save lists
    alpha_list[[iter_MCMC]] <- alpha_res_old
    A_list[[iter_MCMC]] <- A_res_old
    order_list[[iter_MCMC]] <- order_old
  }
  # return results
  return(list(alpha_list = alpha_list[-seq_len(burn_in)], A_list = A_list[-seq_len(burn_in)],
              order_list = order_list[-seq_len(burn_in)]))
}
