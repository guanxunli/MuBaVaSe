#### Define functions
## get sigma0
lBF_model_single <- function(lsigma02, prior_pi, z2, s2) {
  sigma02 <- exp(lsigma02)
  tmp1 <- log(sqrt(s2 / (sigma02 + s2)))
  tmp2 <- z2 / 2 * sigma02 / (sigma02 + s2)
  lBF <- c(tmp1 + tmp2, 0)
  maxlBF <- max(lBF)
  wBF <- exp(lBF - maxlBF)
  wBF_sum <- sum(prior_pi * wBF)
  return(-maxlBF - log(wBF_sum))
}

sigma0_opt_single <- function(lsigma02_int, prior_pi, z2, s2, b_hat) {
  tmp1 <- lBF_model_single(lsigma02 = lsigma02_int, prior_pi = prior_pi, z2 = z2, s2 = s2)
  lsigma02 <- optim(
    par = log(max(c(b_hat^2 - s2, 1))), fn = lBF_model_single, method = "Brent", lower = -30, upper = 15,
    prior_pi = prior_pi, z2 = z2, s2 = s2
  )$par
  tmp2 <- lBF_model_single(lsigma02 = lsigma02, prior_pi = prior_pi, z2 = z2, s2 = s2)
  if (tmp2 < tmp1) {
    return(exp(lsigma02))
  } else {
    return(exp(lsigma02_int))
  }
}

## Calculate the minus KL divergence
KL_fun_single <- function(X_scale, X_scale2, Y, sigma2, b, b2, lBF) {
  n <- length(Y)
  tmp1 <- sum(dnorm(Y, mean = 0, sd = sqrt(sigma2), log = TRUE))
  tmp3 <- n / 2 * log(2 * pi * sigma2)
  tmp4 <- 1 / (2 * sigma2) * (crossprod(Y) - 2 * crossprod(Y, X_scale %*% b) + sum(X_scale2 %*% b2))
  return(tmp1 + lBF + tmp3 + tmp4)
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