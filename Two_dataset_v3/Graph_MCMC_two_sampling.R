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
# set.seed(2022)
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

source("Two_dataset_v3/Graph_given_order_two_sampling.R")
source("Two_dataset_v3/sum_single_effect_two_sampling.R")
Graph_MCMC_two_sampling <- function(dta_1, dta_2, scale_x = FALSE, intercept = TRUE,
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
  lprior_vec <- log(prior_vec)
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
  res_old <- joint_graph_fun_two_sampling(
    dta_1 = dta_1_old, dta_2 = dta_2_old, scale_x = scale_x, intercept = intercept,
    sigma02_int = sigma02_int, sigma2_int = sigma2_int, prior_vec = prior_vec,
    itermax = itermax, L_max = L_max, tol = tol,
    residual_variance_lowerbound = residual_variance_lowerbound
  )
  # graph
  graph_res_1_old <- res_old$graph_res_1
  graph_res_2_old <- res_old$graph_res_2
  # parameters
  alpha_list_old <- res_old$alpha_list
  sigma2_vec_old <- res_old$sigma2_vec
  sigma02_vec_list_old <- res_old$sigma02_vec_list
  lprior_graph_old <- res_old$lprior_graph
  lpropose_graph_old <- res_old$lpropose_graph
  # log likelihood
  llike_vec_1_old <- res_old$llike_vec_1
  llike_vec_2_old <- res_old$llike_vec_2
  ## save lists
  llike_vec <- rep(NA, iter_max)
  graph_list_1 <- list()
  graph_list_2 <- list()
  order_list <- list()
  ## begin MCMC
  for (iter_MCMC in seq_len(iter_max)) {
    if (iter_MCMC %% 1000 == 0) print(iter_MCMC)
    ## Two update methods
    if (sample(c(0, 1), size = 1)) {
      # may change
      lprior_graph_pro <- lprior_graph_old
      lpropose_graph_pro <- lpropose_graph_old
      llike_vec_1_pro <- llike_vec_1_old
      llike_vec_2_pro <- llike_vec_2_old
      graph_res_1_pro <- matrix(0, nrow = p, ncol = p)
      graph_res_2_pro <- matrix(0, nrow = p, ncol = p)
      ## sample one new graph from alpha list
      for (iter_p in seq_len(p - 1)) {
        X_1 <- dta_1_old[, seq_len(iter_p), drop = FALSE]
        Y_1 <- dta_1_old[, iter_p + 1]
        X_2 <- dta_2_old[, seq_len(iter_p), drop = FALSE]
        Y_2 <- dta_2_old[, iter_p + 1]
        # calculate the likelihood
        out_res <- sampling_fun(
          X_1 = X_1, Y_1 = Y_1, X_2 = X_2, Y_2 = Y_2,
          scale_x = scale_x, intercept = intercept,
          lprior_vec = lprior_vec, sigma2 = sigma2_vec_old[iter_p + 1],
          alpha_mat = alpha_list_old[[iter_p + 1]],
          sigma02_vec = sigma02_vec_list_old[[iter_p + 1]]
        )
        llike_vec_1_pro[iter_p + 1] <- out_res$loglikelihood_1
        llike_vec_2_pro[iter_p + 1] <- out_res$loglikelihood_2
        lprior_graph_pro[iter_p + 1] <- out_res$lprior
        lpropose_graph_pro[iter_p + 1] <- out_res$lpropose
        if (is.null(out_res$index1_select) == FALSE) graph_res_1_pro[iter_p + 1, out_res$index1_select] <- 1
        if (is.null(out_res$index2_select) == FALSE) graph_res_2_pro[iter_p + 1, out_res$index2_select] <- 1
      }
      llike_old <- sum(llike_vec_1_old) + sum(llike_vec_2_old) +
        sum(lprior_graph_old) + sum(lpropose_graph_pro)
      llike_pro <- sum(llike_vec_1_pro) + sum(llike_vec_2_pro) +
        sum(lprior_graph_pro) + sum(lpropose_graph_old)
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
        lprior_graph_old <- lprior_graph_pro
        lpropose_graph_old <- lpropose_graph_pro
        llike_vec_1_old <- llike_vec_1_pro
        llike_vec_2_old <- llike_vec_2_pro
        graph_res_1_old <- graph_res_1_pro
        graph_res_2_pro <- graph_res_2_pro
      }
    } else {
      # initialize proposal
      dta_1_pro <- dta_1_old
      dta_2_pro <- dta_2_old
      order_pro <- order_old
      lprior_graph_pro <- lprior_graph_old
      lpropose_graph_pro <- lpropose_graph_old
      llike_vec_1_pro <- llike_vec_1_old
      llike_vec_2_pro <- llike_vec_2_old
      ## propose the new order
      pos_change <- sample(seq_len(p - 1), 1)
      dta_1_pro[, c(pos_change, pos_change + 1)] <- dta_1_old[, c(pos_change + 1, pos_change)]
      dta_2_pro[, c(pos_change, pos_change + 1)] <- dta_2_old[, c(pos_change + 1, pos_change)]
      order_pro[c(pos_change, pos_change + 1)] <- order_old[c(pos_change + 1, pos_change)]
      ## doing variable selection
      if (pos_change == 1) {
        res_pos <- list()
        res_pos$index1_select <- NULL
        res_pos$index2_select <- NULL
        res_pos$sigma2 <- var(c(dta_1_pro[, 1], dta_2_pro[, 1]))
        res_pos$lprior <- 0
        res_pos$lpropose <- 0
        res_pos$alpha_mat <- NULL
        if (intercept) {
          mean_1 <- mean(dta_1_pro[, 1])
          mean_2 <- mean(dta_2_pro[, 1])
        } else {
          mean_1 <- 0
          mean_2 <- 0
        }
        res_pos$loglikelihood_1 <- sum(dnorm(x = dta_1_pro[, 1], mean = mean_1, sd = sqrt(res_pos$sigma2), log = TRUE))
        res_pos$loglikelihood_2 <- sum(dnorm(x = dta_2_pro[, 1], mean = mean_2, sd = sqrt(res_pos$sigma2), log = TRUE))
      } else {
        res_pos <- sum_single_effect_two_sampling(
          X_1 = dta_1_pro[, seq_len(pos_change - 1), drop = FALSE], Y_1 = dta_1_pro[, pos_change],
          X_2 = dta_2_pro[, seq_len(pos_change - 1), drop = FALSE], Y_2 = dta_2_pro[, pos_change],
          scale_x = scale_x, intercept = intercept,
          sigma02_int = sigma02_int, sigma2_int = sigma2_vec_old[pos_change + 1],
          prior_vec = prior_vec, L = min(pos_change - 1, L_max),
          itermax = itermax, tol = tol,
          residual_variance_lowerbound = residual_variance_lowerbound
        )
      }
      res_pos1 <- sum_single_effect_two_sampling(
        X_1 = dta_1_pro[, seq_len(pos_change), drop = FALSE], Y_1 = dta_1_pro[, pos_change + 1],
        X_2 = dta_2_pro[, seq_len(pos_change), drop = FALSE], Y_2 = dta_2_pro[, pos_change + 1],
        scale_x = scale_x, intercept = intercept,
        sigma02_int = sigma02_int, sigma2_int = sigma2_vec_old[pos_change],
        prior_vec = prior_vec, L = min(pos_change, L_max),
        itermax = itermax, tol = tol,
        residual_variance_lowerbound = residual_variance_lowerbound
      )
      lprior_graph_pro[c(pos_change, pos_change + 1)] <- c(res_pos$lprior, res_pos1$lprior)
      lpropose_graph_pro[c(pos_change, pos_change + 1)] <- c(res_pos$lpropose, res_pos1$lpropose)
      llike_vec_1_pro[c(pos_change, pos_change + 1)] <- c(res_pos$loglikelihood_1, res_pos1$loglikelihood_1)
      llike_vec_2_pro[c(pos_change, pos_change + 1)] <- c(res_pos$loglikelihood_2, res_pos1$loglikelihood_2)
      llike_old <- sum(llike_vec_1_old) + sum(llike_vec_2_old) +
        sum(lprior_graph_old) + sum(lpropose_graph_pro)
      llike_pro <- sum(llike_vec_1_pro) + sum(llike_vec_2_pro) +
        sum(lprior_graph_pro) + sum(lpropose_graph_old)
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
        # graphs
        graph_res_1_old[, c(pos_change, pos_change + 1)] <- graph_res_1_old[, c(pos_change + 1, pos_change)]
        graph_res_2_old[, c(pos_change, pos_change + 1)] <- graph_res_2_old[, c(pos_change + 1, pos_change)]
        # proposed matrix
        graph_res_1_old[c(pos_change, pos_change + 1), ] <- 0
        graph_res_2_old[c(pos_change, pos_change + 1), ] <- 0
        # save matrix
        graph_res_1_old[pos_change, res_pos$index1_select] <- 1
        graph_res_2_old[pos_change, res_pos$index2_select] <- 1
        graph_res_1_old[pos_change + 1, res_pos1$index1_select] <- 1
        graph_res_2_old[pos_change + 1, res_pos1$index2_select] <- 1
        # alpha mat
        alpha_list_old[[pos_change]] <- res_pos$alpha_mat
        alpha_list_old[[pos_change + 1]] <- res_pos1$alpha_mat
        if (pos_change + 1 < p) {
          for (iter_p in (pos_change + 2):p) {
            tmp_p <- iter_p - 1
            alpha_list_old[[iter_p]][c(pos_change, pos_change + 1), ] <-
              alpha_list_old[[iter_p]][c(pos_change + 1, pos_change), ]
            alpha_list_old[[iter_p]][c(pos_change + tmp_p, pos_change + tmp_p + 1), ] <-
              alpha_list_old[[iter_p]][c(pos_change + tmp_p + 1, pos_change + tmp_p), ]
            alpha_list_old[[iter_p]][c(pos_change + 2 * tmp_p, pos_change + 2 * tmp_p + 1), ] <-
              alpha_list_old[[iter_p]][c(pos_change + 2 * tmp_p + 1, pos_change + 2 * tmp_p), ]
          }
        }
        # hyper parameters
        sigma2_vec_old[c(pos_change, pos_change + 1)] <- c(res_pos$sigma2, res_pos1$sigma2)
        sigma02_vec_list_old[[pos_change]] <- res_pos$sigma02_vec
        sigma02_vec_list_old[[pos_change + 1]] <- res_pos1$sigma02_vec
        # likelihood
        llike_vec_1_old <- llike_vec_1_pro
        llike_vec_2_old <- llike_vec_2_pro
        lprior_graph_old <- lprior_graph_pro
        lpropose_graph_old <- lpropose_graph_pro
        # data and order
        dta_1_old <- dta_1_pro
        dta_2_old <- dta_2_pro
        order_old <- order_pro
      }
    }
    # save lists
    llike_vec[iter_MCMC] <- sum(llike_vec_1_old) + sum(llike_vec_2_old) +
      sum(lprior_graph_old)
    if (iter_MCMC > burn_in) {
      graph_list_1[[iter_MCMC - burn_in]] <- graph_res_1_old
      graph_list_2[[iter_MCMC - burn_in]] <- graph_res_2_old
      order_list[[iter_MCMC - burn_in]] <- order_old
    }
  }
  # return results
  return(list(
    graph_list_1 = graph_list_1, graph_list_2 = graph_list_2,
    order_list = order_list, llike_vec = llike_vec
  ))
}