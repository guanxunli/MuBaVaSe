## load variable selection function
joint_graph_fun_two_parallel <- function(dta_1, dta_2, scale_x = FALSE, intercept = TRUE,
                                sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL,
                                itermax = 100, L_max = 10, tol = 1e-4, sigma0_low_bd = 1e-8,
                                residual_variance_lowerbound = NULL, n_cores = 4) {
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
  ## save matrix
  # probability of the edge exists
  alpha_res_1 <- matrix(0, nrow = p, ncol = p)
  alpha_res_2 <- matrix(0, nrow = p, ncol = p)
  # posterior mean of the edge
  A_res_1 <- matrix(0, nrow = p, ncol = p)
  A_res_2 <- matrix(0, nrow = p, ncol = p)
  # predicted value
  Xb_mat_1 <- matrix(NA, nrow = n1, ncol = p)
  Xb_mat_2 <- matrix(NA, nrow = n2, ncol = p)
  # log likelihood
  llike_1_vec <- rep(NA, p)
  llike_2_vec <- rep(NA, p)
  llike_penalty_vec <- rep(0, p)
  sigma2_vec <- rep(NA, p)
  if (intercept) {
    mean_1 <- mean(dta_1[, 1])
    mean_2 <- mean(dta_2[, 1])
  } else {
    mean_1 <- 0
    mean_2 <- 0
  }
  Xb_mat_1[, 1] <- rep(mean_1, n1)
  Xb_mat_2[, 1] <- rep(mean_2, n2)
  sigma2_vec[1] <- var(c(dta_1[, 1], dta_2[, 1]))
  llike_1_vec[1] <- sum(dnorm(dta_1[, 1], mean = mean_1, sd = sqrt(sigma2_vec[1]), log = TRUE))
  llike_2_vec[1] <- sum(dnorm(dta_2[, 1], mean = mean_2, sd = sqrt(sigma2_vec[1]), log = TRUE))
  # begin iteration
  library(foreach)
  library(doParallel)
  library(doRNG)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  set.seed(2021)
  out_res <- foreach(iter_p = seq_len(p - 1)) %dorng% {
    source("Two_dataset_new/sum_single_effect_two.R")
    X_1 <- dta_1[, seq_len(iter_p), drop = FALSE]
    Y_1 <- dta_1[, iter_p + 1]
    X_2 <- dta_2[, seq_len(iter_p), drop = FALSE]
    Y_2 <- dta_2[, iter_p + 1]
    ## variable selection
    res <- sum_single_effect_two(
      X_1 = X_1, Y_1 = Y_1, X_2 = X_2, Y_2 = Y_2,
      scale_x = scale_x, intercept = intercept, sigma02_int = sigma02_int,
      sigma2_int = sigma2_int, prior_vec = prior_vec, L = min(iter_p, L_max),
      itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
      residual_variance_lowerbound = residual_variance_lowerbound
    )
  }
  stopCluster(cl)
  for (iter_p in seq_len(p - 1)) {
    Y_1 <- dta_1[, iter_p + 1]
    Y_2 <- dta_2[, iter_p + 1]
    # save the matrix we want
    alpha_res_1[iter_p + 1, seq_len(iter_p)] <- out_res[[iter_p]]$alpha_1
    alpha_res_2[iter_p + 1, seq_len(iter_p)] <- out_res[[iter_p]]$alpha_2
    A_res_1[iter_p + 1, seq_len(iter_p)] <- out_res[[iter_p]]$post_mean1
    A_res_2[iter_p + 1, seq_len(iter_p)] <- out_res[[iter_p]]$post_mean2
    # calculate the likelihood
    sigma2_vec[iter_p + 1] <- out_res[[iter_p]]$sigma2
    Xb_mat_1[, iter_p + 1] <- out_res[[iter_p]]$Xb_1
    Xb_mat_2[, iter_p + 1] <- out_res[[iter_p]]$Xb_2
    llike_1_vec[iter_p + 1] <- sum(dnorm(x = Y_1, mean = out_res[[iter_p]]$Xb_1, sd = sqrt(out_res[[iter_p]]$sigma2), log = TRUE))
    llike_2_vec[iter_p + 1] <- sum(dnorm(x = Y_2, mean = out_res[[iter_p]]$Xb_2, sd = sqrt(out_res[[iter_p]]$sigma2), log = TRUE))
    llike_penalty_vec[iter_p + 1] <- sum(out_res[[iter_p]]$alpha * c(rep(lprior_vec[1], 2 * iter_p), rep(lprior_vec[2], iter_p)))
  }
  ## return results
  return(list(
    alpha_res_1 = alpha_res_1, alpha_res_2 = alpha_res_2, A_res_1 = A_res_1, A_res_2 = A_res_2,
    llike_1_vec = llike_1_vec, llike_2_vec = llike_2_vec, llike_penalty_vec = llike_penalty_vec,
    Xb_mat_1 = Xb_mat_1, Xb_mat_2 = Xb_mat_2, sigma2_vec = sigma2_vec
  ))
}