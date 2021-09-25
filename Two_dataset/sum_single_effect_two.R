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

## main function with null model
# X_1 and X_2 are regressors, n x p matrix, each column is one feature
# Y_1 and Y_2 are response, n x 1 vector
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# r is for common part and q is for single part
# tau is the prior power for null model 1 / (p^tau)
# L is the effict size
# itermax is the maximum iteration
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l
# residual_variance_lowerbound is the lower bound for sigma2

source("Two_dataset/utility_two.R")
sum_single_effect_two <- function(X_1, Y_1, X_2, Y_2, sigma02_int = NULL, sigma2_int = NULL, 
                                 r = 0.2, q = 0.05, tau = 1.5, L = NULL, itermax = 100, 
                                 tol = 1e-4, sigma0_low_bd = 1e-8, residual_variance_lowerbound = NULL) {
  ## Initialization
  p <- ncol(X_1)
  n <- nrow(X_1)
  
  # Initialize sigma
  if (is.null(sigma2_int)) sigma2_int <- as.numeric(var(c(Y_1, Y_2)))
  if (is.null(sigma02_int)) sigma02_int <- 0.2 * sigma2_int
  if (is.null(L)) L <- min(10, p)
  if(is.null(residual_variance_lowerbound)) residual_variance_lowerbound <- sigma2_int / 1e4
  
  # data set 1
  X_scale_1 <- scale(X_1)
  X_scale2_1 <- X_scale_1 * X_scale_1
  X2_1 <- colSums(X_scale2_1)
  mean_Y_1 <- mean(Y_1)
  Y_1 <- Y_1 - mean_Y_1
  # data set 2
  X_scale_2 <- scale(X_2)
  X_scale2_2 <- X_scale_2 * X_scale_2
  X2_2 <- colSums(X_scale2_2)
  mean_Y_2 <- mean(Y_2)
  Y_2 <- Y_2 - mean_Y_2
  
  # Initialize prior
  prior_pi <- c(rep(q, 2 * p), rep(r, p))
  prior_pi <- prior_pi / sum(prior_pi)
  pnull <- 1 - p ^ (1 - tau)
  prior_pi <- c((1 - pnull) * prior_pi, pnull)
  
  ELBO <- rep(NA, itermax + 1)
  ELBO[1] <- -Inf
  sigma2 <- sigma2_int
  sigma02_vec <- rep(sigma02_int, L)
  
  # Save matrix
  b_mat_1 <- matrix(0, nrow = p, ncol = L)
  b2_mat_1 <- matrix(0, nrow = p, ncol = L)
  b_mat_2 <- matrix(0, nrow = p, ncol = L)
  b2_mat_2 <- matrix(0, nrow = p, ncol = L)
  alpha_mat_1 <- matrix(0, nrow = p, ncol = L)
  alpha_mat_2 <- matrix(0, nrow = p, ncol = L)
  
  # Begin iteration
  for (iter in seq_len(itermax)) {
    res_1 <- Y_1 - X_scale_1 %*% rowSums(b_mat_1)
    res_2 <- Y_2 - X_scale_2 %*% rowSums(b_mat_2)
    KL_div <- 0
    for (l in seq_len(L)) {
      ## data set 1
      res_tmp_1 <- res_1 + X_scale_1 %*% b_mat_1[, l]
      # update parameters
      XtY_1 <- crossprod(X_scale_1, res_tmp_1)
      b_hat_1 <- XtY_1 / X2_1
      s2_1 <- sigma2 / X2_1
      z2_1 <- b_hat_1^2 / s2_1
      ## data set 2
      res_tmp_2 <- res_2 + X_scale_2 %*% b_mat_2[, l]
      # update parameters
      XtY_2 <- crossprod(X_scale_2, res_tmp_2)
      b_hat_2 <- XtY_2 / X2_2
      s2_2 <- sigma2 / X2_2
      z2_2 <- b_hat_2^2 / s2_2
      # calculate sigma0
      lsigma02_int <- max(log(sigma02_vec[l]), -30)
      sigma02 <- sigma0_opt_two(lsigma02_int, prior_pi, z2_1, s2_1, z2_2, s2_2, b_hat_1, b_hat_2)
      sigma02_vec[l] <- sigma02
      ## Get Bayesian Factor
      # data set 1
      tmp1_1 <- log(sqrt(s2_1 / (sigma02 + s2_1))) 
      tmp2_1 <- z2_1 / 2 * sigma02 / (sigma02 + s2_1)
      lBF_1 <- tmp1_1 + tmp2_1
      # data set 2
      tmp1_2 <- log(sqrt(s2_2 / (sigma02 + s2_2))) 
      tmp2_2 <- z2_2 / 2 * sigma02 / (sigma02 + s2_2)
      lBF_2 <- tmp1_2 + tmp2_2
      # get bayesian factor
      lBF <- c(lBF_1, lBF_2, lBF_1 + lBF_2, 0)
      maxlBF <- max(lBF)
      wBF <- exp(lBF - maxlBF)
      wBF_sum <- sum(prior_pi * wBF)
      lBF_model <- maxlBF + log(wBF_sum)
      ## Get posterior
      post_alpha <- prior_pi * wBF / wBF_sum
      # data set 1
      post_sigma2_1 <- 1 / (1/s2_1 + 1/sigma02)
      post_mu_1 <- post_sigma2_1 / s2_1 * b_hat_1
      # data set 2
      post_sigma2_2 <- 1 / (1/s2_2 + 1/sigma02)
      post_mu_2 <- post_sigma2_2 / s2_2 * b_hat_2
      ## Calculate posterior mean
      # data set 1
      alpha_mat_1[, l] <- post_alpha[1 : p] + post_alpha[(2 * p + 1) : (3 * p)]
      alpha_mat_2[, l] <- post_alpha[(p + 1) : (2 * p)] + post_alpha[(2 * p + 1) : (3 * p)]
      b_mat_1[, l] <-  alpha_mat_1[, l] * post_mu_1
      b2_mat_1[, l] <- alpha_mat_1[, l] * (post_mu_1^2 + post_sigma2_1)
      # data set 2
      b_mat_2[, l] <- alpha_mat_2[, l] * post_mu_2
      b2_mat_2[, l] <- alpha_mat_2[, l] * (post_mu_2^2 + post_sigma2_2)
      ## calculate the KL divergence
      KL_div <- KL_div + KL_fun_two(X_scale_1 = X_scale_1, Y_1 = res_tmp_1, X_scale_2 = X_scale_2, Y_2 = res_tmp_2, 
                                    X_scale2_1 = X_scale2_1, X_scale2_2 = X_scale2_2,sigma2 = sigma2, 
                                    b_1 = b_mat_1[, l], b2_1 = b2_mat_1[, l], b_2 = b_mat_2[, l], 
                                    b2_2 = b2_mat_2[, l], lBF = lBF_model)
      res_1 <- res_tmp_1 - X_scale_1 %*% b_mat_1[, l]
      res_2 <- res_tmp_2 - X_scale_2 %*% b_mat_2[, l]
    }
    # calculate ELBO
    ERSS_1 <- ERSS_fun_single(X_scale = X_scale_1, X_scale2 = X_scale2_1, Y = Y_1, b_mat = b_mat_1, b2_mat = b2_mat_1)
    ERSS_2 <- ERSS_fun_single(X_scale = X_scale_2, X_scale2 = X_scale2_2, Y = Y_2, b_mat = b_mat_2, b2_mat = b2_mat_2)
    ERSS <- ERSS_1 + ERSS_2
    ELBO[iter + 1] <- - n * log(2 * pi * sigma2) - 1 / (2 * sigma2) * ERSS + KL_div
    # estimate sigma2
    sigma2 <- max(ERSS / (2 * n), residual_variance_lowerbound)
    if (ELBO[iter + 1] -   ELBO[iter] < 1e-4) break
  }
  ELBO <- as.numeric(na.omit(ELBO[-1]))
  # select effect index
  index_L <- which(sigma02_vec > sigma0_low_bd)
  ## return results
  res <- list()
  res$ELBO <- ELBO
  res$sigma2 <- sigma2
  res$sigma02_vec <- sigma02_vec
  
  if (length(index_L) > 0) {
    # data set 1
    res$alpha_1 <- 1 - matrixStats::rowProds(1 - alpha_mat_1[, index_L, drop = FALSE])
    res$post_mean1 <- rowSums(b_mat_1[, index_L, drop = FALSE])
    res$Xb_1 <- mean_Y_1 + X_scale_1 %*% res$post_mean1
    # data set 2
    res$alpha_2 <- 1 - matrixStats::rowProds(1 - alpha_mat_2[, index_L, drop = FALSE])
    res$post_mean2 <- rowSums(b_mat_2[, index_L, drop = FALSE])
    res$Xb_2 <- mean_Y_2 + X_scale_2 %*% res$post_mean2
  } else{
    # data set 1
    res$alpha_1 <- rep(0, p)
    res$post_mean1 <- rep(0, p)
    res$Xb_1 <- rep(mean_Y_1, n)
    # data set 2
    res$alpha_2 <- rep(0, p)
    res$post_mean2 <- rep(0, p)
    res$Xb_2 <- rep(mean_Y_2, n)
  }
  # return results
  return(res)
}

# res <- sum_single_effect_two(X_1, Y_1, X_2, Y_2, sigma02_int = NULL, sigma2_int = NULL,
#                              r = 0.2, q = 0.05, tau = 1.5, L = p_c + p_1 + p_2, itermax = 100,
#                              tol = 1e-4, sigma0_low_bd = 1e-8)