# ## Define parameters
# n <- 100
# p <- 1000
# p_c <- 17
# p_1 <- 3
# p_2 <- 3
# sigma <- 0.1
# sigma0 <- 0.6
# r <- 0.2
# q <- 0.05
# set.seed(2021)
# ## Generate data
# index_c <- sample(seq_len(p), size = p_c, replace = FALSE)
# index_1 <- sample(setdiff(seq_len(p), index_c), size = p_1, replace = FALSE)
# index_2 <- sample(setdiff(seq_len(p), c(index_1, index_c)), size = p_2, replace = FALSE)
# 
# b_1 <- rep(0, p)
# b_1[c(index_c, index_1)] <- rnorm(p_c + p_1, mean = 0, sd = sigma0)
# b_2 <- rep(0, p)
# b_2[c(index_c, index_1)] <- rnorm(p_c + p_2, mean = 0, sd = sigma0)
# 
# X_1 <- matrix(rnorm(p * n), nrow = n, ncol = p)
# X_2 <- matrix(rnorm(p * n), nrow = n, ncol = p)
# Y_1 <- X_1 %*% b_1 + rnorm(n, sd = sigma)
# Y_2 <- X_2 %*% b_2 + rnorm(n, sd = sigma)

#### Define functions
## get sigma0
lBF_model_multi <- function(sigma0, prior_pi, z2_1, s2_1, z2_2, s2_2) {
  # Data set 1
  tmp1_1 <- log(sqrt(s2_1 / (sigma0^2 + s2_1))) 
  tmp2_1 <- z2_1 / 2 * sigma0^2 / (sigma0^2 + s2_1)
  lBF_1 <- tmp1_1 + tmp2_1
  # Data set 2
  tmp1_2 <- log(sqrt(s2_2 / (sigma0^2 + s2_2))) 
  tmp2_2 <- z2_2 / 2 * sigma0^2 / (sigma0^2 + s2_2)
  lBF_2 <- tmp1_2 + tmp2_2
  # combine them all
  lBF <- c(lBF_1, lBF_2, lBF_1 + lBF_2)
  maxlBF <- max(lBF)
  wBF <- exp(lBF - maxlBF)
  wBF_sum <- sum(prior_pi * wBF)
  return(- maxlBF - log(wBF_sum))
}

sigma0_opt_multi <- function(sigma0_opt_fun, sigma0_int, prior_pi, z2_1, s2_1, z2_2, s2_2) {
  tmp1 <- sigma0_opt_fun(sigma0 = sigma0_int, prior_pi = prior_pi, 
                         z2_1 = z2_1, s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2)
  sigma0 <- optim(sigma0_int, sigma0_opt_fun, method = "L-BFGS-B", lower = 0, 
                  prior_pi = prior_pi,  z2_1 = z2_1, s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2)$par
  tmp2 <- sigma0_opt_fun(sigma0 = sigma0, prior_pi = prior_pi, 
                         z2_1 = z2_1, s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2)
  if (tmp2 < tmp1) {
    return(sigma0)
  } else{
    return(sigma0_int)
  }
}

## Calculate the KL divergence
KL_fun_multi <- function(X_1, Y_1, X_2, Y_2, sigma, b1_1, b2_1, b1_2, b2_2, 
                   sigma0, prior_pi, z2_1, s2_1, z2_2, s2_2) {
  n <- length(Y_1)
  tmp1_1 <- sum(dnorm(Y_1, mean = 0, sd = sigma, log = TRUE))
  tmp1_2 <- sum(dnorm(Y_2, mean = 0, sd = sigma, log = TRUE))
  tmp2 <- -lBF_model_multi(sigma0 = sigma0, prior_pi = prior_pi, z2_1 = z2_1, s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2)
  tmp3 <- n * log(2 * pi * sigma^2)
  tmp4_1 <- 1 / (2 * sigma^2) * (crossprod(Y_1) - 2 * crossprod(Y_1, X_1 %*% b1_1) + sum(X_1^2 %*% b1_2))
  tmp4_2 <- 1 / (2 * sigma^2) * (crossprod(Y_2) - 2 * crossprod(Y_2, X_2 %*% b2_1) + sum(X_2^2 %*% b2_2))
  return(tmp1_1 + tmp1_2 + tmp2 + tmp3 + tmp4_1 + tmp4_2)
}

## Calculate ERSS
ERSS_fun_multi <- function(X, Y, b_mat, b2_mat) {
  mu_lmat <- X %*% b_mat
  mu2_lmat <- X^2 %*% b2_mat
  mu_pred <- rowSums(mu_lmat)
  res_tmp <- sum((Y - mu_pred)^2)
  var_sum <- sum(mu2_lmat - mu_lmat^2)
  return(res_tmp + var_sum)
}

sum_single_effect_multi <- function(X_1, Y_1, X_2, Y_2, sigma0_int = 1, sigma_int = 1, 
                              r = 0.2, q = 0.05, L = 10, itermax = 100, tol = 1e-4) {
  # Initialization
  p <- ncol(X_1)
  n <- nrow(X_1)
  prior_pi <- c(rep(q, 2 * p), rep(r, p))
  prior_pi <- prior_pi / sum(prior_pi)
  XtX_1 <- crossprod(X_1)
  XtX_2 <- crossprod(X_2)
  sigma <- sigma_int
  sigma0 <- sigma0_int
  ELBO <- rep(NA, itermax + 1)
  ELBO[1] <- -Inf
  # Save matrix
  b_mat_1 <- matrix(0, nrow = p, ncol = L)
  b2_mat_1 <- matrix(0, nrow = p, ncol = L)
  b_mat_2 <- matrix(0, nrow = p, ncol = L)
  b2_mat_2 <- matrix(0, nrow = p, ncol = L)
  b_mat_1_list <- list()
  b_mat_2_list <- list()
  # Begin iteration
  for (iter in seq_len(itermax)) {
    res_1 <- Y_1 - X_1 %*% rowSums(b_mat_1)
    res_2 <- Y_2 - X_2 %*% rowSums(b_mat_2)
    KL_div <- 0
    for (l in seq_len(L)) {
      res_tmp_1 <- res_1 + X_1 %*% b_mat_1[, l]
      res_tmp_2 <- res_2 + X_2 %*% b_mat_2[, l]
      # update parameters
      XtY_1 <- crossprod(X_1, res_tmp_1)
      XtY_2 <- crossprod(X_2, res_tmp_2)
      b_hat_1 <- XtY_1 / diag(XtX_1)
      b_hat_2 <- XtY_2 / diag(XtX_2)
      s2_1 <- sigma^2 / diag(XtX_1)
      s2_2 <- sigma^2 / diag(XtX_2)
      z2_1 <- b_hat_1^2 / s2_1
      z2_2 <- b_hat_2^2 / s2_2
      # calculate sigma0
      sigma0_int <- max(c(b_hat_1^2 - s2_1, b_hat_2^2 - s2_2, 1))
      sigma0 <- sigma0_opt_multi(lBF_model_multi, sigma0_int, prior_pi = prior_pi, 
                           z2_1 = z2_1, s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2)
      ## Get Bayesian Factor
      # Data set 1
      tmp1_1 <- log(sqrt(s2_1 / (sigma0^2 + s2_1))) 
      tmp2_1 <- z2_1 / 2 * sigma0^2 / (sigma0^2 + s2_1)
      lBF_1 <- tmp1_1 + tmp2_1
      # Data set 2
      tmp1_2 <- log(sqrt(s2_2 / (sigma0^2 + s2_2))) 
      tmp2_2 <- z2_2 / 2 * sigma0^2 / (sigma0^2 + s2_2)
      lBF_2 <- tmp1_2 + tmp2_2
      # combine them all
      lBF <- c(lBF_1, lBF_2, lBF_1 + lBF_2)
      maxlBF <- max(lBF)
      wBF <- exp(lBF - maxlBF)
      wBF_sum <- sum(prior_pi * wBF)
      # Get posterior
      post_alpha <- prior_pi * wBF / wBF_sum
      post_sigma2_1 <- 1 / (1/s2_1 + 1/sigma0^2)
      post_mu_1 <- post_sigma2_1 / s2_1 * b_hat_1
      post_sigma2_2 <- 1 / (1/s2_2 + 1/sigma0^2)
      post_mu_2 <- post_sigma2_2 / s2_2 * b_hat_2
      # Calculate posterior mean
      b_mat_1[, l] <- (post_alpha[1 : p] + post_alpha[(2 * p + 1) : (3 * p)]) * post_mu_1
      b2_mat_1[, l] <- (post_alpha[1 : p] + post_alpha[(2 * p + 1) : (3 * p)]) * (post_mu_1^2 + post_sigma2_1)
      b_mat_2[, l] <- (post_alpha[(p + 1) : (2 * p)] + post_alpha[(2 * p + 1) : (3 * p)]) * post_mu_2
      b2_mat_2[, l] <- (post_alpha[(p + 1) : (2 * p)] + post_alpha[(2 * p + 1) : (3 * p)]) * (post_mu_2^2 + post_sigma2_2)
      KL_div <- KL_div + KL_fun_multi(X_1 = X_1, Y_1 = res_tmp_1, X_2 = X_2, Y_2 = res_tmp_2, sigma = sigma, 
                                b1_1 = b_mat_1[, l], b2_1 = b_mat_2[, l], b1_2 = b2_mat_1[, l], b2_2 = b2_mat_2[, l], 
                                sigma0 = sigma0, prior_pi = prior_pi, z2_1 = z2_1, s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2)
      res_1 <- res_tmp_1 - X_1 %*% b_mat_1[, l]
      res_2 <- res_tmp_2 - X_2 %*% b_mat_2[, l]
    }
    ERSS_1 <- ERSS_fun_multi(X = X_1, Y = Y_1, b_mat = b_mat_1, b2_mat = b2_mat_1)
    ERSS_2 <- ERSS_fun_multi(X = X_2, Y = Y_2, b_mat = b_mat_2, b2_mat = b2_mat_2)
    sigma <- (ERSS_1 + ERSS_2) / (2 * n)
    ELBO[iter + 1] <- -n * log(2 * pi * sigma) - 1 / (2 * sigma) * (ERSS_1 + ERSS_2) + KL_div
    b_mat_1_list[[iter]] <- b_mat_1
    b_mat_2_list[[iter]] <- b_mat_2
    if (abs(ELBO[iter + 1] -   ELBO[iter]) < 1e-4)  break
  }
  ELBO <- as.numeric(na.omit(ELBO[-1]))
  b_mat_1 <- b_mat_1_list[[which.max(ELBO)]]
  b_mat_2 <- b_mat_2_list[[which.max(ELBO)]]
  res <- list()
  res$ELBO <- ELBO
  res$b_mat_1 <- b_mat_1
  res$b_mat_2 <- b_mat_2
  return(res)
}

# res <- sum_single_effect_multi(X_1, Y_1, X_2, Y_2, sigma0_int = 1, sigma_int = 1, 
#                          r = 0.2, q = 0.05, L = 22, itermax = 100, tol = 1e-4)
# print(sort(c(index_1, index_c)))
# which(abs(round(rowSums(res$b_mat_1), 4)) > 0)
# print(sort(c(index_2, index_c)))
# which(abs(round(rowSums(res$b_mat_2), 4)) > 0)
