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

## MCMC method for Graph
# dta_1 and dta_2 are p x n data set
# scale_x : scale the data
# intercept: calculate the mean of Y
# order_int is the initialized order for nodes
# iter_max is the maximun mcmc step
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# r is for common part and q is for single part
# tau is the prior power for null model 1 / (p^tau)
# itermax is the maximum iteration
# L_max is the largest number of parents
# tol is the threshold for ELBO
# residual_variance_lowerbound is the lower bound for sigma2

source("Two_dataset_v3/Graph_given_order_two_single.R")
source("Two_dataset_v3/sum_single_effect_two_single.R")
Graph_MCMC_two_single <- function(dta_1, dta_2, scale_x = FALSE, intercept = TRUE,
                                  order_int = NULL, iter_max = 50000,
                                  sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL,
                                  itermax = 100, L_max = 10, tol = 1e-4,
                                  burn_in = iter_max - 5000, residual_variance_lowerbound = NULL) {
  ## Initialization
  p <- ncol(dta_1)
  if (p != ncol(dta_2)) stop("The number of features should be same!")
  n1 <- nrow(dta_1)
  n2 <- nrow(dta_2)
  ## define prior vector
  if (is.null(prior_vec)) {
    prior_vec <- c(1 / (2 * p^1.5), 1 / p^2)
  }
  # Initialize order
  if (is.null(order_int)) {
    order_old <- sample(seq_len(p), p)
  } else {
    order_old <- order_int
  }
  # change data
  dta_1_old <- dta_1[, order_old]
  dta_2_old <- dta_2[, order_old]
  ## load the main function
  res_old <- joint_graph_fun_two_single(
    dta_1 = dta_1_old, dta_2 = dta_2_old, scale_x = scale_x, intercept = intercept,
    sigma02_int = sigma02_int, sigma2_int = sigma2_int, prior_vec = prior_vec,
    itermax = itermax, L_max = L_max, tol = tol,
    residual_variance_lowerbound = residual_variance_lowerbound
  )
  # variable selection
  alpha_res_1_old <- res_old$alpha_res_1
  alpha_res_2_old <- res_old$alpha_res_2
  # likelihood
  sigma2_vec_old <- res_old$sigma2_vec
  llike_1_vec_old <- res_old$llike_1_vec
  llike_2_vec_old <- res_old$llike_2_vec
  llike_old <- sum(llike_1_vec_old) + sum(llike_2_vec_old)
  llike_vec <- rep(NA, iter_max)
  ## save lists
  alpha_list_1 <- list()
  alpha_list_2 <- list()
  order_list <- list()
  ## load the function
  for (iter_MCMC in seq_len(iter_max)) {
    if (iter_MCMC %% 1000 == 0) print(iter_MCMC)
    ## Initialize proposal
    dta_1_pro <- dta_1_old
    dta_2_pro <- dta_2_old
    order_pro <- order_old
    # log likelihood
    sigma2_vec_pro <- sigma2_vec_old
    llike_1_vec_pro <- llike_1_vec_old
    llike_2_vec_pro <- llike_2_vec_old
    ## propose the new order
    pos_change <- sample(seq_len(p - 1), 1)
    llike_pro <- llike_old - sum(llike_1_vec_old[c(pos_change, pos_change + 1)]) - sum(llike_2_vec_old[c(pos_change, pos_change + 1)])
    dta_1_pro[, c(pos_change, pos_change + 1)] <- dta_1_old[, c(pos_change + 1, pos_change)]
    dta_2_pro[, c(pos_change, pos_change + 1)] <- dta_2_old[, c(pos_change + 1, pos_change)]
    order_pro[c(pos_change, pos_change + 1)] <- order_old[c(pos_change + 1, pos_change)]
    ## doing variable selection
    if (pos_change == 1) {
      res_pos <- list()
      res_pos$index1_select <- NULL
      res_pos$index2_select <- NULL
      res_pos$sigma2 <- var(c(dta_1_pro[, 1], dta_2_pro[, 1]))
      res_pos$loglikelihood_1 <- sum(dnorm(x = dta_1_pro[, 1], mean = 0, sd = sqrt(res_pos$sigma2), log = TRUE))
      res_pos$loglikelihood_2 <- sum(dnorm(x = dta_2_pro[, 1], mean = 0, sd = sqrt(res_pos$sigma2), log = TRUE))
    } else {
      res_pos <- sum_single_effect_two_single(
        X_1 = dta_1_pro[, seq_len(pos_change - 1), drop = FALSE], Y_1 = dta_1_pro[, pos_change],
        X_2 = dta_2_pro[, seq_len(pos_change - 1), drop = FALSE], Y_2 = dta_2_pro[, pos_change],
        scale_x = scale_x, intercept = intercept,
        sigma02_int = sigma02_int, sigma2_int = sigma2_vec_old[pos_change + 1],
        prior_vec = prior_vec, L = min(pos_change - 1, L_max),
        itermax = itermax, tol = tol,
        residual_variance_lowerbound = residual_variance_lowerbound
      )
    }
    res_pos1 <- sum_single_effect_two_single(
      X_1 = dta_1_pro[, seq_len(pos_change), drop = FALSE], Y_1 = dta_1_pro[, pos_change + 1],
      X_2 = dta_2_pro[, seq_len(pos_change), drop = FALSE], Y_2 = dta_2_pro[, pos_change + 1],
      scale_x = scale_x, intercept = intercept,
      sigma02_int = sigma02_int, sigma2_int = sigma2_vec_old[pos_change],
      prior_vec = prior_vec, L = min(pos_change, L_max),
      itermax = itermax, tol = tol,
      residual_variance_lowerbound = residual_variance_lowerbound
    )
    # likelihood
    sigma2_vec_pro[c(pos_change, pos_change + 1)] <- c(res_pos$sigma2, res_pos1$sigma2)
    llike_1_vec_pro[pos_change] <- res_pos$loglikelihood_1
    llike_2_vec_pro[pos_change] <- res_pos$loglikelihood_2
    llike_1_vec_pro[pos_change + 1] <- res_pos1$loglikelihood_1
    llike_2_vec_pro[pos_change + 1] <- res_pos1$loglikelihood_2
    llike_pro <- llike_pro + sum(llike_1_vec_pro[c(pos_change, pos_change + 1)]) + sum(llike_2_vec_pro[c(pos_change, pos_change + 1)])
    # accept or not
    if (llike_pro > llike_old) {
      accept <- TRUE
    } else {
      U <- runif(1)
      thres <- exp(llike_pro - llike_old)
      if (U < thres) {
        accept <- TRUE
      } else {
        accept <- FALSE
      }
    }
    # update
    if (accept) {
      # change matrix order
      alpha_res_1_old[, c(pos_change, pos_change + 1)] <- alpha_res_1_old[, c(pos_change + 1, pos_change)]
      alpha_res_2_old[, c(pos_change, pos_change + 1)] <- alpha_res_2_old[, c(pos_change + 1, pos_change)]
      # proposed matrix
      alpha_res_1_old[c(pos_change, pos_change + 1), ] <- 0
      alpha_res_2_old[c(pos_change, pos_change + 1), ] <- 0
      # save matrix
      alpha_res_1_old[pos_change, res_pos$index1_select] <- 1
      alpha_res_2_old[pos_change, res_pos$index2_select] <- 1
      alpha_res_1_old[pos_change + 1, res_pos1$index1_select] <- 1
      alpha_res_2_old[pos_change + 1, res_pos1$index2_select] <- 1
      # likelihood
      sigma2_vec_old <- sigma2_vec_pro
      llike_1_vec_old <- llike_1_vec_pro
      llike_2_vec_old <- llike_2_vec_pro
      llike_old <- llike_pro
      # data and order
      dta_1_old <- dta_1_pro
      dta_2_old <- dta_2_pro
      order_old <- order_pro
    }
    ## save lists
    llike_vec[iter_MCMC] <- llike_old
    if (iter_MCMC > burn_in) {
      alpha_list_1[[iter_MCMC - burn_in]] <- alpha_res_1_old
      alpha_list_2[[iter_MCMC - burn_in]] <- alpha_res_2_old
      order_list[[iter_MCMC - burn_in]] <- order_old
    }
  }
  # return results
  return(list(
    alpha_list_1 = alpha_list_1, alpha_list_2 = alpha_list_2,
    order_list = order_list, llike_vec = llike_vec
  ))
}