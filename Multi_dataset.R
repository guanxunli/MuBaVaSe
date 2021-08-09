# ## Define parameters
# n <- 500
# p <- 1000
# p_c <- 30
# p_1 <- 5
# p_2 <- 5
# sigma <- 1
# sigma0 <- 0.6
# r <- 0.2
# q <- 0.05
# set.seed(1234)
# ## Generate data
# index_c <- sample(seq_len(p), size = p_c, replace = FALSE)
# index_1 <- sample(setdiff(seq_len(p), index_c), size = p_1, replace = FALSE)
# index_2 <- sample(setdiff(seq_len(p), c(index_1, index_c)), size = p_2, replace = FALSE)
# 
# b_1 <- rep(0, p)
# b_1[c(index_c, index_1)] <- rnorm(p_c + p_1, mean = 0, sd = sigma0)
# b_2 <- rep(0, p)
# b_2[c(index_c, index_2)] <- rnorm(p_c + p_2, mean = 0, sd = sigma0)
# 
# X_1 <- matrix(rnorm(p * n), nrow = n, ncol = p)
# X_2 <- matrix(rnorm(p * n), nrow = n, ncol = p)
# Y_1 <- X_1 %*% b_1 + rnorm(n, sd = sigma)
# Y_2 <- X_2 %*% b_2 + rnorm(n, sd = sigma)

#### Define functions
## get sigma0
lBF_model_multi <- function(lsigma02, prior_pi, z2_1, s2_1, z2_2, s2_2) {
  sigma02 <- exp(lsigma02)
  # data set 1
  tmp1_1 <- log(sqrt(s2_1 / (sigma02 + s2_1)))
  tmp2_1 <- z2_1 / 2 * sigma02 / (sigma02 + s2_1)
  lBF_1 <- tmp1_1 + tmp2_1
  # data set 2
  tmp1_2 <- log(sqrt(s2_2 / (sigma02 + s2_2)))
  tmp2_2 <- z2_2 / 2 * sigma02 / (sigma02 + s2_2)
  lBF_2 <- tmp1_2 + tmp2_2
  # combine
  lBF <- c(lBF_1, lBF_2, lBF_1 + lBF_2)
  maxlBF <- max(lBF)
  wBF <- exp(lBF - maxlBF)
  wBF_sum <- sum(prior_pi * wBF)
  return(- maxlBF - log(wBF_sum))
}

sigma0_opt_multi <- function(lsigma02_int, prior_pi, z2_1, s2_1, z2_2, s2_2, b_hat_1, b_hat_2) {
  tmp1 <- lBF_model_multi(lsigma02 = lsigma02_int, prior_pi = prior_pi, z2_1 = z2_1, 
                          s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2)
  lsigma02 <- optim(par = log(max(c(b_hat_1^2 - s2_1, 1, b_hat_2^2 - s2_2))), fn = lBF_model_multi, 
                    method = "Brent", lower = -30, upper = 15, prior_pi = prior_pi, z2_1 = z2_1, 
                    s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2)$par
  tmp2 <- lBF_model_multi(lsigma02 = lsigma02, prior_pi = prior_pi, z2_1 = z2_1, 
                          s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2)
  if (tmp2 < tmp1) {
    return(exp(lsigma02))
  } else{
    return(exp(lsigma02_int))
  }
}

## Calculate the KL divergence
KL_fun_multi <- function(X_scale_1, X_scale2_1, Y_1, X_scale_2, Y_2, X_scale2_2,
                         sigma2, b_1, b2_1, b_2, b2_2, lBF) {
  n <- length(Y_1)
  tmp1_1 <- sum(dnorm(Y_1, mean = 0, sd = sqrt(sigma2), log = TRUE))
  tmp1_2 <- sum(dnorm(Y_2, mean = 0, sd = sqrt(sigma2), log = TRUE))
  tmp3 <- n * log(2 * pi * sigma2)
  tmp4_1 <- 1 / (2 * sigma2) * (crossprod(Y_1) - 2 * crossprod(Y_1, X_scale_1%*% b_1) + sum(X_scale2_1 %*% b2_1))
  tmp4_2 <- 1 / (2 * sigma2) * (crossprod(Y_2) - 2 * crossprod(Y_2, X_scale_2%*% b_2) + sum(X_scale2_2 %*% b2_2))
  return(tmp1_1 + tmp1_2 + lBF + tmp3 + tmp4_1 + tmp4_2)
}

## Calculate ERSS
ERSS_fun_single <- function(X_scale, X_scale2, Y, b_mat, b2_mat) {
  mu_lmat <- X_scale %*% b_mat
  mu2_lmat <- X_scale2 %*% b2_mat
  mu_pred <- rowSums(mu_lmat)
  res_tmp <- sum((Y - mu_pred)^2)
  var_sum <- sum(mu2_lmat - mu_lmat^2)
  return(res_tmp + var_sum)
}

## Select effect set
effset_fun <- function(index_L, alpha_mat, Xcor, cor_low_bd = 0.5, sigma0_low_bd = 1e-8) {
  if (length(index_L) == 0) {
    return(NULL)
  } else{
    index_eff <- NULL
    for (iter_l in seq_len(length(index_L))) {
      index_sel <- which(alpha_mat[, iter_l] > 1/p)
      if (length(index_sel) == 1) {
        index_eff <- c(index_eff, index_sel)
      } else if (length(index_sel) > 1) {
        cor_tmp <- max(abs(Xcor[index_sel, index_sel]))
        if (cor_tmp > cor_low_bd) index_eff <- c(index_eff, index_sel)
      }
    }
  }
  return(index_eff)
}

# ## main function
# sigma02_int = NULL
# sigma2_int = NULL
# r = 0.2
# q = 0.05
# L = NULL
# itermax = 100
# tol = 1e-4
# cor_low_bd = 0.5
# sigma0_low_bd = 1e-8
sum_single_effect_multi <- function(X_1, Y_1, X_2, Y_2, sigma02_int = NULL, sigma2_int = NULL, 
                              r = 0.2, q = 0.05, L = NULL, itermax = 100, tol = 1e-4, 
                              cor_low_bd = 0.5, sigma0_low_bd = 1e-8) {
  ## Initialization
  p <- ncol(X_1)
  n <- nrow(X_1)
  # data set 1
  mean_Y_1 <- mean(Y_1)
  Y_1 <- Y_1 - mean_Y_1
  X_scale_1 <- scale(X_1)
  X2_1 <- colSums(X_scale_1 * X_scale_1)
  X_scale2_1 <- X_scale_1 * X_scale_1
  Xcor_1 <- cor(X_1)
  diag(Xcor_1) <- 0
  # data set 2
  mean_Y_2 <- mean(Y_2)
  Y_2 <- Y_2 - mean_Y_2
  X_scale_2 <- scale(X_2)
  X2_2 <- colSums(X_scale_2 * X_scale_2)
  X_scale2_2 <- X_scale_2 * X_scale_2
  Xcor_2 <- cor(X_2)
  diag(Xcor_2) <- 0
  
  if (is.null(sigma2_int)) sigma2_int <- as.numeric(var(c(Y_1, Y_2)))
  if (is.null(sigma02_int)) sigma02_int <- 0.2 * sigma2_int
  if (is.null(L)) L <- min(10, p)
  prior_pi <- c(rep(q, 2 * p), rep(r, p))
  prior_pi <- prior_pi / sum(prior_pi)
  
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
      sigma02 <- sigma0_opt_multi(lsigma02_int, prior_pi, z2_1, s2_1, z2_2, s2_2, b_hat_1, b_hat_2)
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
      lBF <- c(lBF_1, lBF_2, lBF_1 + lBF_2)
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
      KL_div <- KL_div + KL_fun_multi(X_scale_1 = X_scale_1, Y_1 = res_tmp_1, X_scale_2 = X_scale_2, Y_2 = res_tmp_2, 
                                      X_scale2_1 = X_scale2_1, X_scale2_2 = X_scale2_2,sigma2 = sigma2, 
                                      b_1 = b_mat_1[, l], b2_1 = b2_mat_1[, l], b_2 = b_mat_2[, l], 
                                      b2_2 = b2_mat_2[, l], lBF = lBF_model)
      res_1 <- res_tmp_1 - X_scale_1 %*% b_mat_1[, l]
      res_2 <- res_tmp_2 - X_scale_2 %*% b_mat_2[, l]
    }
    ERSS_1 <- ERSS_fun_single(X_scale = X_scale_1, X_scale2 = X_scale2_1, Y = Y_1, b_mat = b_mat_1, b2_mat = b2_mat_1)
    ERSS_2 <- ERSS_fun_single(X_scale = X_scale_2, X_scale2 = X_scale2_2, Y = Y_2, b_mat = b_mat_2, b2_mat = b2_mat_2)
    ERSS <- ERSS_1 + ERSS_2
    ELBO[iter + 1] <- - n * log(2 * pi * sigma2) - 1 / (2 * sigma2) * ERSS + KL_div
    sigma2 <- ERSS / (2 * n)
    if (ELBO[iter + 1] -   ELBO[iter] < 1e-4) break
  }
  ELBO <- as.numeric(na.omit(ELBO[-1]))
  # select effect index
  index_L <- which(sigma02_vec > sigma0_low_bd)
  index_eff_1 <- effset_fun(index_L, alpha_mat_1, Xcor_1, cor_low_bd, sigma0_low_bd)
  index_eff_2 <- effset_fun(index_L, alpha_mat_2, Xcor_2, cor_low_bd, sigma0_low_bd)
  # return results
  res <- list()
  res$ELBO <- ELBO
  res$sigma2 <- sigma2
  res$sigma02_vec <- sigma02_vec
  
  res$alpha_mat_1 <- alpha_mat_1
  res$post_mean1 <- rowSums(b_mat_1[, index_L, drop = FALSE])
  res$index_eff_1 <- index_eff_1
  
  res$alpha_mat_2 <- alpha_mat_2
  res$post_mean2 <- rowSums(b_mat_2[, index_L, drop = FALSE])
  res$index_eff_2 <- index_eff_2

  return(res)
}

# #### check results
# ## package
# # data set 1
# res <- susieR::susie(X = X_1, y = Y_1, L = 35)
# res1 <- as.numeric(res$sets$cs)
# length(intersect(res1, c(index_1, index_c))) / (p_1 + p_c)
# length(intersect(res1, c(index_1, index_c))) / length(res1)
# sum((colSums(res$alpha * res$mu) - b_1)^2)
# # data set 2
# res <- susieR::susie(X = X_2, y = Y_2, L = 35)
# res2 <- as.numeric(res$sets$cs)
# length(intersect(res2, c(index_2, index_c))) / (p_2 + p_c)
# length(intersect(res2, c(index_2, index_c))) / length(res2)
# sum((colSums(res$alpha * res$mu) - b_2)^2)
# # new method
# res <- sum_single_effect_multi(X_1, Y_1, X_2, Y_2, L = 40)
# res1 <- res$index_eff_1
# length(intersect(res1, c(index_1, index_c))) / (p_1 + p_c)
# length(intersect(res1, c(index_1, index_c))) / length(res1)
# sum((res$post_mean1- b_1)^2)
# res2 <- res$index_eff_2
# length(intersect(res2, c(index_2, index_c))) / (p_2 + p_c)
# length(intersect(res2, c(index_2, index_c))) / length(res2)
# sum((res$post_mean2 - b_2)^2)
