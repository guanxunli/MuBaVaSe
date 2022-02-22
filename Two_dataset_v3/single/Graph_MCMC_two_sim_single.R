source("Two_dataset_v3/single/Graph_given_order_two_single.R")
source("Two_dataset_v3/single/sum_single_effect_two_single.R")
Graph_MCMC_two_sim_single <- function(dta_1, dta_2, scale_x = FALSE, intercept = FALSE,
                                      order_int = NULL, iter_max = 50000,
                                      sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL,
                                      itermax = 100, L_max = 10, tol = 1e-4,
                                      burn_in = iter_max - 5000, prior_penalty = FALSE,
                                      residual_variance_lowerbound = NULL,
                                      adj_true1 = NULL, adj_true2 = NULL) {
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
  # posterior of parameters
  A_res_1_old <- res_old$A_res_1
  A_res_2_old <- res_old$A_res_2
  # likelihood
  sigma2_vec_old <- res_old$sigma2_vec
  llike_1_vec_old <- res_old$llike_1_vec
  llike_2_vec_old <- res_old$llike_2_vec
  lprior_graph_old <- res_old$lprior_graph
  llike_old <- sum(llike_1_vec_old) + sum(llike_2_vec_old) 
  if (prior_penalty) {
    llike_old <- llike_old + sum(lprior_graph_old)
  }
  llike_vec <- rep(NA, iter_max)
  error_mat1 <- matrix(NA, nrow = 2, ncol = iter_max)
  error_mat2 <- matrix(NA, nrow = 2, ncol = iter_max)
  ## save lists
  alpha_list_1 <- list()
  alpha_list_2 <- list()
  order_list <- list()
  ## load the true results
  # data set 1
  adj_1 <- alpha_res_1_old[order(order_old), order(order_old)]
  adj_1 <- ifelse(adj_1 > 0.5, 1, 0)
  adj_1 <- t(adj_1)
  g_1 <- as(getGraph(adj_1), "graphNEL")
  # data set 2
  adj_2 <- alpha_res_2_old[order(order_old), order(order_old)]
  adj_2 <- ifelse(adj_2 > 0.5, 1, 0)
  adj_2 <- t(adj_2)
  g_2 <- as(getGraph(adj_2), "graphNEL")
  # load true value
  g_true1 <- as(getGraph(adj_true1), "graphNEL")
  g_true2 <- as(getGraph(adj_true2), "graphNEL")
  ## begin mcmc
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
    lprior_graph_pro <- lprior_graph_old
    ## propose the new order
    pos_change <- sample(seq_len(p - 1), 1)
    llike_pro <- llike_old - sum(llike_1_vec_old[c(pos_change, pos_change + 1)]) - 
      sum(llike_2_vec_old[c(pos_change, pos_change + 1)])
    if (prior_penalty) {
      llike_pro <- llike_pro - sum(lprior_graph_old[c(pos_change, pos_change + 1)])
    }
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
      res_pos$lprior_graph <- 0
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
    lprior_graph_pro[pos_change] <- res_pos$lprior_graph
    llike_1_vec_pro[pos_change + 1] <- res_pos1$loglikelihood_1
    llike_2_vec_pro[pos_change + 1] <- res_pos1$loglikelihood_2
    lprior_graph_pro[pos_change + 1] <- res_pos1$lprior_graph
    llike_pro <- llike_pro + sum(llike_1_vec_pro[c(pos_change, pos_change + 1)]) + 
      sum(llike_2_vec_pro[c(pos_change, pos_change + 1)])
    if (prior_penalty) {
      llike_pro <- llike_pro + sum(lprior_graph_pro[c(pos_change, pos_change + 1)])
    }
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
      lprior_graph_old <- lprior_graph_pro
      llike_old <- llike_pro
      # data and order
      dta_1_old <- dta_1_pro
      dta_2_old <- dta_2_pro
      order_old <- order_pro
    }
    ## save lists
    llike_vec[iter_MCMC] <- llike_old
    # check error
    adj_1 <- alpha_res_1_old[order(order_old), order(order_old)]
    adj_1 <- t(adj_1)
    g_1 <- as(getGraph(adj_1), "graphNEL")
    adj_2 <- alpha_res_2_old[order(order_old), order(order_old)]
    adj_2 <- t(adj_2)
    g_2 <- as(getGraph(adj_2), "graphNEL")
    # save results
    error_mat1[, iter_MCMC] <- c(
      pcalg::shd(g_true1, g_1),
      check_edge(adj_true1, adj_1)
    )
    error_mat2[, iter_MCMC] <- c(
      pcalg::shd(g_true2, g_2),
      check_edge(adj_true2, adj_2)
    )
    if (iter_MCMC > burn_in) {
      alpha_list_1[[iter_MCMC - burn_in]] <- alpha_res_1_old
      alpha_list_2[[iter_MCMC - burn_in]] <- alpha_res_2_old
      order_list[[iter_MCMC - burn_in]] <- order_old
    }
  }
  # return results
  return(list(
    alpha_list_1 = alpha_list_1, alpha_list_2 = alpha_list_2,
    order_list = order_list, llike_vec = llike_vec,
    error_mat1 = error_mat1, error_mat2 = error_mat2
  ))
}