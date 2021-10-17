# ## Define parameters
# n <- 500
# p <- 1000
# p_c <- 25
# p_1 <- 5
# p_2 <- 5
# sigma <- 1
# sigma0 <- 0.6
# set.seed(2021)
# ## Generate data
# index_c <- sample(seq_len(p), size = p_c, replace = FALSE)
# index_1 <- sample(setdiff(seq_len(p), index_c), size = p_1, replace = FALSE)
# index_2 <- sample(setdiff(seq_len(p), index_c), size = p_2, replace = FALSE)
# 
# b_1 <- rep(0, p)
# b_1[c(index_c, index_1)] <- rnorm(p_c + p_1, mean = 0, sd = sigma0)
# # b_1[c(index_c, index_1)] <- c(rep(1,15), rep(0.05, 10), rep(0.1, 5))
# b_2 <- rep(0, p)
# b_2[c(index_c, index_2)] <- rnorm(p_c + p_2, mean = 0, sd = sigma0)
# # b_2[c(index_c, index_2)] <- c(rep(0.05,15), rep(1, 10), rep(0.1, 5))
# 
# alpha_1 <- rep(0, p)
# alpha_1[c(index_c, index_1)] <- 1
# alpha_2 <- rep(0, p)
# alpha_2[c(index_c, index_2)] <- 1
# 
# X_1 <- matrix(rnorm(p * n), nrow = n, ncol = p)
# X_2 <- matrix(rnorm(p * n), nrow = n, ncol = p)
# Y_1 <- X_1 %*% b_1 + rnorm(n, sd = sigma)
# Y_2 <- X_2 %*% b_2 + rnorm(n, sd = sigma)
# 
# dta_list <- list()
# dta_list[[1]] <- list(X = X_1, Y = Y_1)
# dta_list[[2]] <- list(X = X_2, Y = Y_2)

## main function with null model
# dta_list: dta list has X and Y
# scale_x : scale the data
# intercept: calculate the mean of Y
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# L is the effict size
# prior_vec : prior for different models
# itermax is the maximum iteration
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l
source("multi_utility.R")
sum_single_effect_mult <- function(dta_list, scale_x = TRUE, intercept = TRUE,
                                   sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL,
                                   L = NULL, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
                                   com_mat = NULL, com_list = NULL,
                                   residual_variance_lowerbound = NULL) {
  ## Initialization
  K <- length(dta_list)
  n_group <- 2^K - 1
  p <- ncol(dta_list[[1]]$X)
  n <- nrow(dta_list[[1]]$X)
  # combinatorics matrix
  if (is.null(com_mat)) {
    com_list <- list()
    com_mat <- matrix(c(0, 1), ncol = 1)
    for (iter in 2:K) {
      com_mat_copy <- com_mat
      com_mat <- cbind(1, com_mat)
      com_mat_copy <- cbind(0, com_mat_copy)
      com_mat <- rbind(com_mat_copy, com_mat)
    }
    com_mat <- com_mat[-1, ]
    
    for (iter_com in seq_len(n_group)) {
      com_list[[iter_com]] <- which(com_mat[iter_com, ] == 1)
    }
  }
  
  # initialize list
  if (is.null(L)) L <- min(10, p)
  Y_list <- list()
  b_list <- list()
  b2_list <- list()
  alpha_list <- list()
  alpha_mat <- matrix(0, nrow = n_group * p, ncol = L)
  for (iter_K in seq_len(K)) {
    Y_list[[iter_K]] <- dta_list[[iter_K]]$Y
    # Initialize b and alpha
    b_list[[iter_K]] <- matrix(0, nrow = p, ncol = L)
    b2_list[[iter_K]] <- matrix(0, nrow = p, ncol = L)
    alpha_list[[iter_K]] <- matrix(0, nrow = p, ncol = L)
    # Initialize data set
    if (intercept) {
      dta_list[[iter_K]]$mean_Y <- mean(dta_list[[iter_K]]$Y)
    } else {
      dta_list[[iter_K]]$mean_Y <- 0
    }
    dta_list[[iter_K]]$Y <- dta_list[[iter_K]]$Y - dta_list[[iter_K]]$mean_Y
    # Initialize X
    if (scale_x) {
      dta_list[[iter_K]]$X_scale <- scale(dta_list[[iter_K]]$X)
    } else {
      dta_list[[iter_K]]$X_scale <- dta_list[[iter_K]]$X
    }
    dta_list[[iter_K]]$X_scale2 <- dta_list[[iter_K]]$X_scale * dta_list[[iter_K]]$X_scale
    dta_list[[iter_K]]$X2 <- colSums(dta_list[[iter_K]]$X_scale2)
  }
  
  # Initialize sigma
  if (is.null(sigma2_int)) sigma2_int <- as.numeric(var(unlist(Y_list)))
  if (is.null(sigma02_int)) sigma02_int <- 0.2 * sigma2_int
  if (is.null(residual_variance_lowerbound)) residual_variance_lowerbound <- sigma2_int / 1e4
  
  # Initialize ELBO
  ELBO <- rep(NA, itermax + 1)
  ELBO[1] <- -Inf
  sigma2 <- sigma2_int
  sigma02_vec <- rep(sigma02_int, L)
  
  # Initialize prior
  if (is.null(prior_vec)) {
    prior_vec <- rep(1 / (n_group * p^1.5))
  }
  prior_pi <- rep(prior_vec, each = p)
  prior_pi <- c(prior_pi, 1 - sum(prior_pi))
  ## storage list
  res_list <- matrix(NA, nrow = n, ncol = K)
  res_tmp_list <- matrix(NA, nrow = n, ncol = K)
  XtY_list <- matrix(NA, nrow = p, ncol = K)
  b_hat_list <- matrix(NA, nrow = p, ncol = K)
  s2_list <- matrix(NA, nrow = p, ncol = K)
  z2_list <- matrix(NA, nrow = p, ncol = K)
  for (iter in seq_len(itermax)) {
    ## calculate residuals
    for (iter_K in seq_len(K)) {
      res_list[, iter_K] <- dta_list[[iter_K]]$Y - dta_list[[iter_K]]$X_scale %*% rowSums(b_list[[iter_K]])
    }
    KL_div <- 0
    ## begin iteration
    for (l in seq_len(L)) {
      # update parameters
      for (iter_K in seq_len(K)) {
        res_tmp_list[, iter_K] <- res_list[, iter_K] + dta_list[[iter_K]]$X_scale %*% b_list[[iter_K]][, l]
        XtY_list[, iter_K] <- crossprod(dta_list[[iter_K]]$X_scale, res_tmp_list[, iter_K])
        b_hat_list[, iter_K] <- XtY_list[, iter_K] / dta_list[[iter_K]]$X2
        s2_list[, iter_K] <- sigma2 / dta_list[[iter_K]]$X2
        z2_list[, iter_K] <- b_hat_list[, iter_K]^2 / s2_list[, iter_K]
      }
      # calculate sigma0
      lsigma02_int <- max(log(sigma02_vec[l]), -30)
      sigma02 <- sigma0_opt_multi(lsigma02_int, prior_pi, z2_list, s2_list, b_hat_list, n_group, p, com_list)
      sigma02_vec[l] <- sigma02
      ## Get Bayesian Factor
      tmp1 <- log(sqrt(s2_list / (sigma02 + s2_list)))
      tmp2 <- z2_list / 2 * sigma02 / (sigma02 + s2_list)
      lBF_list <- tmp1 + tmp2
      # get bayesian factor
      lBF <- rep(0, n_group * p + 1)
      for (iter_com in seq_len(n_group)) {
        if (length(com_list[[iter_com]]) == 1) {
          lBF[(1 + p * (iter_com - 1)):(p * iter_com)] <- lBF_list[, com_list[[iter_com]]]
        } else {
          lBF[(1 + p * (iter_com - 1)):(p * iter_com)] <- rowSums(lBF_list[, com_list[[iter_com]], drop = FALSE])
        }
      }
      maxlBF <- max(lBF)
      wBF <- exp(lBF - maxlBF)
      wBF_sum <- sum(prior_pi * wBF)
      lBF_model <- maxlBF + log(wBF_sum)
      ## Get posterior
      post_alpha <- prior_pi * wBF / wBF_sum
      alpha_mat[, l] <- post_alpha[-length(post_alpha)]
      post_sigma2_list <- 1 / (1 / s2_list + 1 / sigma02)
      post_mu_list <- post_sigma2_list / s2_list * b_hat_list
      ## Calculate posterior mean
      post_alpha_mat <- matrix(post_alpha[-length(post_alpha)], nrow = p, ncol = n_group)
      for (iter_K in seq_len(K)) {
        alpha_list[[iter_K]][, l] <- rowSums(post_alpha_mat[, which(com_mat[, iter_K] == 1), drop = FALSE])
        b_list[[iter_K]][, l] <- alpha_list[[iter_K]][, l] * post_mu_list[, iter_K]
        b2_list[[iter_K]][, l] <- alpha_list[[iter_K]][, l] * (post_mu_list[, iter_K]^2 + post_sigma2_list[, iter_K])
      }
      ## Calculate the KL-divergence
      KL_div <- KL_div + KL_fun_multi(dta_list, res_tmp_list, sigma2, b_list, b2_list, lBF_model, K, n, l)
      for (iter_K in seq_len(K)) {
        res_list[, iter_K] <- res_tmp_list[, iter_K] - dta_list[[iter_K]]$X_scale %*% b_list[[iter_K]][, l]
      }
    }
    ERSS <- 0
    for (iter_K in seq_len(K)) {
      ERSS <- ERSS + ERSS_fun_single(dta_list[[iter_K]]$X_scale, dta_list[[iter_K]]$X_scale2,
                                     Y = dta_list[[iter_K]]$Y, b_mat = b_list[[iter_K]], b2_mat = b2_list[[iter_K]]
      )
    }
    ELBO[iter + 1] <- -K * n / 2 * log(2 * pi * sigma2) - 1 / (2 * sigma2) * ERSS + KL_div
    # estimate sigma2
    sigma2 <- max(ERSS / (K * n), residual_variance_lowerbound)
    if (ELBO[iter + 1] - ELBO[iter] < 1e-4) break
  }
  
  ELBO <- as.numeric(na.omit(ELBO[-1]))
  # select effect index
  index_L <- which(sigma02_vec > sigma0_low_bd)
  ## return results
  out_res <- list()
  out_res$ELBO <- ELBO
  out_res$sigma2 <- sigma2
  out_res$sigma02_vec <- sigma02_vec
  out_res$res <- list()
  
  if (length(index_L) > 0) {
    out_res$alpha <- 1 - matrixStats::rowProds(1 - alpha_mat[, index_L, drop = FALSE])
    for (iter_K in seq_len(K)) {
      out_res$res[[iter_K]] <- list()
      out_res$res[[iter_K]]$alpha <- 1 - matrixStats::rowProds(1 - alpha_list[[iter_K]][, index_L, drop = FALSE])
      out_res$res[[iter_K]]$post_mean <- rowSums(b_list[[iter_K]][, index_L, drop = FALSE])
      out_res$res[[iter_K]]$Xb <- dta_list[[iter_K]]$mean_Y + dta_list[[iter_K]]$X_scale %*% out_res$res[[iter_K]]$post_mean
    }
  } else {
    out_res$alpha <- rep(0, n_group * p)
    for (iter_K in seq_len(K)) {
      out_res$res[[iter_K]] <- list()
      out_res$res[[iter_K]]$alpha <- rep(0, p)
      out_res$res[[iter_K]]$post_mean <- rep(0, p)
      out_res$res[[iter_K]]$Xb <- rep(dta_list[[iter_K]]$mean_Y, n)
    }
  }
  # return results
  return(out_res)
}

# ## general function
# time1 <- Sys.time()
# out_res <- sum_single_effect_mult(dta_list,
#                                   L = p_c + p_1 + p_2,
#                                   prior_vec = c(1 / (6 * p^1.5), 1 / (6 * p^1.5), 2 / (3 * p^1.5))
# )
# Sys.time() - time1
# 
# ## original function for twp data sets
# source("Two_dataset_new/sum_single_effect_two.R")
# time1 <- Sys.time()
# res <- sum_single_effect_two(X_1, Y_1, X_2, Y_2, L = p_c + p_1 + p_2, prior_vec = c(1 / (6 * p^1.5), 2 / (3 * p^1.5)))
# Sys.time() - time1
# 
# ## compare two results
# sum((out_res$ELBO - res$ELBO)^2)