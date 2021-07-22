# n <- 100
# p <- 1000
# sigma <- 0.1
# sigma0 <- 0.6
# L <- 10
# set.seed(2021)
# ## Generate data
# index_t <- sample(seq_len(p), size = L, replace = FALSE)
# b <- rep(0, p)
# b[index_t] <- rnorm(L, mean = 0, sd = sigma0)
# X <- matrix(rnorm(n * p), nrow = n, ncol = p)
# Y<- X %*% b + rnorm(n, sd = sigma)

#### Define functions
## get sigma0
lBF_model_single <- function(sigma0, prior_pi, z2, s2) {
  tmp1 <- log(sqrt(s2 / (sigma0^2 + s2))) 
  tmp2 <- z2 / 2 * sigma0^2 / (sigma0^2 + s2)
  lBF <- tmp1 + tmp2
  maxlBF <- max(lBF)
  wBF <- exp(lBF - maxlBF)
  wBF_sum <- sum(prior_pi * wBF)
  return(- maxlBF - log(wBF_sum))
}

sigma0_opt_single <- function(sigma0_opt_fun, sigma0_int, prior_pi, z2, s2) {
  tmp1 <- sigma0_opt_fun(sigma0 = sigma0_int, prior_pi = prior_pi, z2 = z2, s2 = s2)
  sigma0 <- optim(sigma0_int, sigma0_opt_fun, method = "L-BFGS-B", lower = 0, 
                  prior_pi = prior_pi, z2 = z2, s2 = s2)$par
  tmp2 <- sigma0_opt_fun(sigma0 = sigma0, prior_pi = prior_pi, z2 = z2, s2 = s2)
  if (tmp2 < tmp1) {
    return(sigma0)
  } else{
    return(sigma0_int)
  }
}

## Calculate the KL divergence
KL_fun_single <- function(X, Y, sigma, b, b2, sigma0, prior_pi, z2, s2) {
  n <- length(Y)
  tmp1 <- sum(dnorm(Y, mean = 0, sd = sigma, log = TRUE))
  tmp2 <- -lBF_model_single(sigma0 = sigma0, prior_pi = prior_pi, z2 = z2, s2 = s2)
  tmp3 <- n / 2 * log(2 * pi * sigma^2)
  tmp4 <- 1 / (2 * sigma^2) * (crossprod(Y) - 2 * crossprod(Y, X %*% b) + sum(X^2 %*% b2))
  return(tmp1 + tmp2 + tmp3 + tmp4)
}

## Calculate ERSS
ERSS_fun_single <- function(X, Y, b_mat, b2_mat) {
  mu_lmat <- X %*% b_mat
  mu2_lmat <- X^2 %*% b2_mat
  mu_pred <- rowSums(mu_lmat)
  res_tmp <- sum((Y - mu_pred)^2)
  var_sum <- sum(mu2_lmat - mu_lmat^2)
  return(res_tmp + var_sum)
}

## Sum of single effect model
sum_single_effect_single <- function(X, Y, sigma_int = 1, sigma0_int = 1, prior_pi, L, itermax = 100, tol = 1e-4) {
  # Initialization
  p <- ncol(X)
  n <- nrow(X)
  XtX <- crossprod(X)
  sigma <- sigma_int
  sigma0 <- sigma0_int
  ELBO <- rep(NA, itermax + 1)
  ELBO[1] <- -Inf
  # Save matrix
  b_mat <- matrix(0, nrow = p, ncol = L)
  b2_mat <- matrix(0, nrow = p, ncol = L)
  b_mat_list <- list()
  # Begin iteration
  for (iter in seq_len(itermax)) {
    res <- Y - X %*% rowSums(b_mat)
    KL_div <- 0
    for (l in seq_len(L)) {
      res_tmp <- res + X %*% b_mat[, l]
      # update parameters
      XtY <- crossprod(X, res_tmp)
      b_hat <- XtY / diag(XtX)
      s2 <- sigma^2 / diag(XtX)
      z2 <- b_hat^2 / s2
      # calculate sigma0
      sigma0_int <- max(c(b_hat^2 - s2, 1))
      sigma0 <- sigma0_opt_single(lBF_model_single, sigma0_int, prior_pi = prior_pi, z2 = z2, s2 = s2)
      # Get Bayesian Factor
      tmp1 <- log(sqrt(s2 / (sigma0^2 + s2))) 
      tmp2 <- z2 / 2 * sigma0^2 / (sigma0^2 + s2)
      lBF <- tmp1 + tmp2
      maxlBF <- max(lBF)
      wBF <- exp(lBF - maxlBF)
      wBF_sum <- sum(prior_pi * wBF)
      # Get posterior
      post_alpha <- prior_pi * wBF / wBF_sum
      post_sigma2 <- 1 / (1/s2 + 1/sigma0^2)
      post_mu <- post_sigma2 / s2 * b_hat
      # Calculate posterior mean
      b_mat[, l] <- post_alpha * post_mu
      b2_mat[, l] <- post_alpha * (post_mu^2 + post_sigma2)
      KL_div <- KL_div + KL_fun_single(X = X, Y = res_tmp, sigma = sigma, b = b_mat[, l], b2 = b2_mat[, l],
                                sigma0 = sigma0, prior_pi = prior_pi, z2 = z2, s2 = s2)
      res <- res_tmp - X %*% b_mat[, l]
    }
    ERSS <- ERSS_fun_single(X = X, Y = Y, b_mat = b_mat, b2_mat = b2_mat)
    sigma <- ERSS / n
    ELBO[iter + 1] <- -n / 2 * log(2 * pi * sigma) - 1 / (2 * sigma) * ERSS + KL_div
    b_mat_list[[iter]] <- b_mat
    if (abs(ELBO[iter + 1] -   ELBO[iter]) < 1e-4) break
  }
  ELBO <- as.numeric(na.omit(ELBO[-1]))
  b_mat <- b_mat_list[[which.max(ELBO)]]
  res <- list()
  res$ELBO <- ELBO
  res$b_mat <- b_mat
  return(res)
}

# #### check results
# prior_pi <- rep(1/p, p)
# res <- sum_single_effect_single(X = X, Y = Y, prior_pi = prior_pi, L = L)
# print(sort(index_t))
# which(abs(round(rowSums(res$b_mat), 4)) > 0)

# res <- susieR::susie(X = X, y = Y)
# print(sort(index_t))
# which(abs(round(colSums(res$alpha * res$mu), 4)) > 0)
