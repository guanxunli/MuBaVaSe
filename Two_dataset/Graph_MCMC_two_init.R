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

## MCMC method for Graph
# dta_1 and dta_2 are p x n data set
# order_int is the initialized order for nodes
# iter_max is the maximun mcmc step
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# r is for common part and q is for single part
# tau is the prior power for null model 1 / (p^tau)
# itermax is the maximum iteration
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l
# residual_variance_lowerbound is the lower bound for sigma2

Graph_MCMC_two_init <- function(dta_1, dta_2, order_int = NULL, iter_max = 10000, sigma02_int = NULL, sigma2_int = NULL, r = 0.2, 
                                q = 0.05, tau = 1.5, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 5000, 
                                residual_variance_lowerbound = NULL) {
  ## Initialization
  p <- ncol(dta_1)
  n <- nrow(dta_1)
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
  source("Two_dataset/Graph_given_order_two_init.R")
  res_old <- joint_graph_fun_two_init(dta_1 = dta_1_old, dta_2 = dta_2_old, sigma02_int = sigma02_int, sigma2_int = sigma2_int, r = r, 
                                      q = q, tau = tau, itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
                                      residual_variance_lowerbound = residual_variance_lowerbound)
  # variable selection
  alpha_res_1_old <- res_old$alpha_res_1
  alpha_res_2_old <- res_old$alpha_res_2
  # posterior of parameters
  A_res_1_old <- res_old$A_res_1
  A_res_2_old <- res_old$A_res_2
  # coefficient matrix
  b_list_1_old <- res_old$b_list_1
  b_list_2_old <- res_old$b_list_2
  # likelihood
  sigma2_vec_old <- res_old$sigma2_vec
  llike_1_vec_old <- res_old$llike_1_vec
  llike_2_vec_old <- res_old$llike_2_vec
  llike_old <- sum(llike_1_vec_old) + sum(llike_2_vec_old)
  ## save lists
  alpha_list_1 <- list()
  alpha_list_2 <- list()
  A_list_1 <- list()
  A_list_2 <- list()
  order_list <- list()
  ## load the function
  source("Two_dataset/sum_single_effect_two_init.R")
  for (iter_MCMC in seq_len(iter_max)) {
    if (iter_MCMC %% 100 == 0) print(iter_MCMC)
    ## Initialize proposal
    dta_1_pro <- dta_1_old
    dta_2_pro <- dta_2_old
    order_pro <- order_old
    b_list_1_pro <- b_list_1_old
    b_list_2_pro <- b_list_2_old
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
      res_pos$sigma2 <- var(c(dta_1_pro[, 1], dta_2_pro[, 1]))
      res_pos$alpha_1 <- rep(0, p)
      res_pos$post_mean1 <- rep(0, p)
      res_pos$Xb_1 <- rep(mean(dta_1_pro[, 1]), n)
      res_pos$alpha_2 <- rep(0, p)
      res_pos$post_mean2 <- rep(0, p)
      res_pos$Xb_2 <- rep(mean(dta_2_pro[, 1]), n)
    } else {
      b_list_1_pro[[pos_change]] <- b_list_1_old[[pos_change + 1]][-pos_change, seq_len(min(pos_change - 1, 10))]
      b_list_2_pro[[pos_change]] <- b_list_2_old[[pos_change + 1]][-pos_change, seq_len(min(pos_change - 1, 10))]
      res_pos <- sum_single_effect_two_init(X_1 = dta_1_pro[, seq_len(pos_change - 1), drop = FALSE], Y_1 = dta_1_pro[, pos_change],
                                            X_2 = dta_2_pro[, seq_len(pos_change - 1), drop = FALSE], Y_2 = dta_2_pro[, pos_change],
                                            sigma02_int = sigma02_int, sigma2_int = sigma2_vec_old[pos_change + 1], 
                                            r = r, q = q, tau = tau, L = min(pos_change - 1, 10), 
                                            itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
                                            residual_variance_lowerbound = residual_variance_lowerbound,
                                            b_mat_1 = b_list_1_pro[[pos_change]], b_mat_2 = b_list_2_pro[[pos_change]])
    }
    if (pos_change < 11) {
      b_list_1_pro[[pos_change + 1]] <- cbind(rbind(b_list_1_old[[pos_change]], 0), 0)
      b_list_2_pro[[pos_change + 1]] <- cbind(rbind(b_list_2_old[[pos_change]], 0), 0)
    } else{
      b_list_1_pro[[pos_change + 1]] <- rbind(b_list_1_old[[pos_change]], 0)
      b_list_2_pro[[pos_change + 1]] <- rbind(b_list_2_old[[pos_change]], 0)
    }
    res_pos1 <- sum_single_effect_two_init(X_1 = dta_1_pro[, seq_len(pos_change), drop = FALSE], Y_1 = dta_1_pro[, pos_change + 1],
                                           X_2 = dta_2_pro[, seq_len(pos_change), drop = FALSE], Y_2 = dta_2_pro[, pos_change + 1],
                                           sigma02_int = sigma02_int, sigma2_int = sigma2_vec_old[pos_change], 
                                           r = r, q = q, tau = tau, L = min(pos_change, 10), 
                                           itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
                                           residual_variance_lowerbound = residual_variance_lowerbound,
                                           b_mat_1 = b_list_1_pro[[pos_change + 1]], b_mat_2 = b_list_2_pro[[pos_change + 1]] )
    # likelihood
    sigma2_vec_pro[c(pos_change, pos_change + 1)] <- c(res_pos$sigma2, res_pos1$sigma2)
    llike_1_vec_pro[pos_change] <- sum(dnorm(x = dta_1_pro[, pos_change], mean = res_pos$Xb_1, sd = sqrt(res_pos$sigma2), log = TRUE))
    llike_2_vec_pro[pos_change] <- sum(dnorm(x = dta_2_pro[, pos_change], mean = res_pos$Xb_2, sd = sqrt(res_pos$sigma2), log = TRUE))
    llike_1_vec_pro[pos_change + 1] <- sum(dnorm(x = dta_1_pro[, pos_change + 1], mean = res_pos1$Xb_1, sd = sqrt(res_pos1$sigma2), log = TRUE))
    llike_2_vec_pro[pos_change + 1] <- sum(dnorm(x = dta_2_pro[, pos_change + 1], mean = res_pos1$Xb_2, sd = sqrt(res_pos1$sigma2), log = TRUE))
    llike_pro <- llike_pro + sum(llike_1_vec_pro[c(pos_change, pos_change + 1)]) + sum(llike_2_vec_pro[c(pos_change, pos_change + 1)])
    b_list_1_pro[[pos_change]] <- res_pos$b_mat_1
    b_list_2_pro[[pos_change]] <- res_pos$b_mat_2
    b_list_1_pro[[pos_change + 1]] <- res_pos1$b_mat_1
    b_list_2_pro[[pos_change + 1]] <- res_pos1$b_mat_2
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
      # proposed matrix
      alpha_res_1_old[c(pos_change, pos_change + 1), ] <- 0
      alpha_res_2_old[c(pos_change, pos_change + 1), ] <- 0
      A_res_1_old[c(pos_change, pos_change + 1), ] <- 0
      A_res_2_old[c(pos_change, pos_change + 1), ] <- 0
      # save matrix
      alpha_res_1_old[pos_change, seq_len(pos_change - 1)] <- res_pos$alpha_1
      alpha_res_2_old[pos_change, seq_len(pos_change - 1)] <- res_pos$alpha_2
      alpha_res_1_old[pos_change + 1, seq_len(pos_change)] <- res_pos1$alpha_1
      alpha_res_2_old[pos_change + 1, seq_len(pos_change)] <- res_pos1$alpha_2
      A_res_1_old[pos_change, seq_len(pos_change - 1)] <- res_pos$post_mean1
      A_res_2_old[pos_change, seq_len(pos_change - 1)] <- res_pos$post_mean2
      A_res_1_old[pos_change + 1, seq_len(pos_change)] <- res_pos1$post_mean1
      A_res_2_old[pos_change + 1, seq_len(pos_change)] <- res_pos1$post_mean2
      b_list_1_old <- b_list_1_pro
      b_list_2_old <- b_list_2_pro
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
    alpha_list_1[[iter_MCMC]] <- alpha_res_1_old
    alpha_list_2[[iter_MCMC]] <- alpha_res_2_old
    A_list_1[[iter_MCMC]] <- A_res_1_old
    A_list_2[[iter_MCMC]] <- A_res_2_old
    order_list[[iter_MCMC]] <- order_old
  }
  # return results
  return(list(alpha_list_1 = alpha_list_1[-seq_len(burn_in)], alpha_list_2 = alpha_list_2[-seq_len(burn_in)], 
              A_list_1 = A_list_1[-seq_len(burn_in)], A_list_2 = A_list_2[-seq_len(burn_in)],
              order_list = order_list[-seq_len(burn_in)]))
}

# # ## MCMC
# time1 <- Sys.time()
# res <- Graph_MCMC_two_init(dta_1 = dta_1, dta_2 = dta_2, iter_max = 200, burn_in = 1)
# Sys.time() - time1 # 
