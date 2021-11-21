n <- 500
p <- 1000
sigma <- 1
sigma0 <- 0.6
L <- 20
set.seed(2021)
## Generate data
index_t <- sample(seq_len(p), size = L, replace = FALSE)
b <- rep(0, p)
b[index_t] <- rnorm(L, mean = 0, sd = sigma0)
# b[index_t] <- 100
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- X %*% b + rnorm(n, sd = sigma)

## main function with null model
# X are regressors, n x p matrix, each column is one feature
# Y is response, n x 1 vector
# scale_x : scale the data
# intercept: calculate the mean of Y
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# prior_null is prior for the null model
# L is the effect size
# itermax is the maximum iteration
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l
# residual_variance_lowerbound is the lower bound for sigma2

source("single_dataset/utility_single.R")
sum_single_effect_single <- function(X, Y, scale_x = TRUE, intercept = TRUE,
                                  sigma02_int = NULL, sigma2_int = NULL, prior_null = NULL,
                                  L = NULL, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
                                  residual_variance_lowerbound = NULL) {
  ## Initialization
  p <- ncol(X)
  n <- nrow(X)
  
  # Initialize sigma
  if (is.null(sigma2_int)) sigma2_int <- as.numeric(var(Y))
  if (is.null(sigma02_int)) sigma02_int <- 0.2 * sigma2_int
  if (is.null(L)) L <- min(10, p)
  if (is.null(residual_variance_lowerbound)) residual_variance_lowerbound <- sigma2_int / 1e4
  
  ## data preprocess
  # scale
  if (scale_x) {
    X_scale <- scale(X)
  } else {
    X_scale <- X_1
  }
  # intercept
  if (intercept) {
    mean_Y <- mean(Y)
  } else {
    mean_Y <- 0
  }
  # pre-calculate
  X_scale2 <- X_scale * X_scale
  X2 <- colSums(X_scale2)
  Y <- Y- mean_Y

  # Initialize prior
  if (is.null(prior_null)) {
    prior_null <- 1 - 1 / (p ^ 1.5)
  } 
  prior_pi <- c(rep((1 - prior_null) / p, p), prior_null)
  
  # initialize ELBO
  ELBO <- rep(NA, itermax + 1)
  ELBO[1] <- -Inf
  sigma2 <- sigma2_int
  sigma02_vec <- rep(sigma02_int, L)
  
  # Save matrix
  b_mat <- matrix(0, nrow = p, ncol = L)
  b2_mat <- matrix(0, nrow = p, ncol = L)
  alpha_mat <- matrix(0, nrow = p, ncol = L)
  alpha_null <- rep(0, L)
  
  # Begin iteration
  for (iter in seq_len(itermax)) {
    res <- Y - X_scale %*% rowSums(b_mat)
    KL_div <- 0
    # # keep old results
    # sigma02_vec_old <- sigma02_vec
    # sigma2_old <- sigma2
    # alpha_null_old <- alpha_null
    # alpha_mat_old <- alpha_mat
    # b_mat_old <- b_mat
    for (l in seq_len(L)) {
      # residuals
      res_tmp <- res + X_scale %*% b_mat[, l]
      # update parameters
      XtY <- crossprod(X_scale, res_tmp)
      b_hat <- XtY / X2
      s2 <- sigma2 / X2
      z2 <- b_hat^2 / s2
      # calculate sigma0
      lsigma02_int <- max(log(sigma02_vec[l]), -30)
      sigma02 <- sigma0_opt_single(lsigma02_int, prior_pi, z2, s2, b_hat)
      sigma02_vec[l] <- sigma02
      ## Get Bayesian Factor
      tmp1 <- log(sqrt(s2 / (sigma02 + s2)))
      tmp2 <- z2 / 2 * sigma02 / (sigma02 + s2)
      lBF <- c(tmp1 + tmp2, 0)
      maxlBF <- max(lBF)
      wBF <- exp(lBF - maxlBF)
      wBF_sum <- sum(prior_pi * wBF)
      lBF_model <- maxlBF + log(wBF_sum)
      ## Get posterior
      post_alpha <- prior_pi * wBF / wBF_sum
      alpha_null[l] <- post_alpha[length(post_alpha)]
      alpha_mat[, l] <- post_alpha[-length(post_alpha)]
      post_sigma2 <- 1 / (1 / s2 + 1 / sigma02)
      post_mu <- post_sigma2 / s2 * b_hat
      ## Calculate posterior mean
      b_mat[, l] <- alpha_mat[, l] * post_mu
      b2_mat[, l] <- alpha_mat[, l] * (post_mu^2 + post_sigma2)
      ## calculate the KL divergence
      KL_div <- KL_div + KL_fun_single(
        X_scale = X_scale, Y = res_tmp, X_scale2 = X_scale2, sigma2 = sigma2,
        b = b_mat[, l], b2 = b2_mat[, l], lBF = lBF_model
      )
      res <- res_tmp - X_scale %*% b_mat[, l]
    }
    # calculate ELBO
    ERSS <- ERSS_fun_single(X_scale = X_scale, X_scale2 = X_scale2, Y = Y, b_mat = b_mat, b2_mat = b2_mat)
    ELBO[iter + 1] <- -n / 2 * log(2 * pi * sigma2) - 1 / (2 * sigma2) * ERSS + KL_div
    # estimate sigma2
    sigma2 <- max(ERSS / (2 * n), residual_variance_lowerbound)
    if (ELBO[iter + 1] - ELBO[iter] < 1e-4) break
  }
  ELBO <- as.numeric(na.omit(ELBO[-1]))
  # if (ELBO[length(ELBO)] < ELBO[length(ELBO) - 1]) {
  #   sigma02_vec <- sigma02_vec_old
  #   sigma2 <- sigma2_old
  #   alpha_null <- alpha_null_old
  #   alpha_mat <- alpha_mat_old
  #   b_mat <- b_mat_old
  # }
  # select effect index
  index_L <- which(sigma02_vec > sigma0_low_bd)
  ## return results
  res <- list()
  res$ELBO <- ELBO
  res$sigma2 <- sigma2
  res$sigma02_vec <- sigma02_vec
  res$alpha_null <- alpha_null
  res$index_L <- index_L
  
  if (length(index_L) > 0) {
    res$alpha <- 1 - matrixStats::rowProds(1 - alpha_mat[, index_L, drop = FALSE])
    res$post_mean <- rowSums(b_mat[, index_L, drop = FALSE])
    res$Xb <- mean_Y + X_scale %*% res$post_mean
  } else {
    # data set 1
    res$alpha <- rep(0, p)
    res$post_mean <- rep(0, p)
    res$Xb <- rep(mean_Y, n)
  }
  # return results
  return(res)
}

## check results
res <- sum_single_effect_single(X = X, Y = Y, L = L)
res$ELBO
res1 <- which(res$alpha > 0.5)
length(intersect(res1, index_t)) / L
length(intersect(res1, index_t)) / length(res1)
sum((res$post_mean- b)^2)
