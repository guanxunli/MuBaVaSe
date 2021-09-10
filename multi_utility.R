#### Define functions
## get sigma0
lBF_model_multi <- function(lsigma02, prior_pi, z2_list, s2_list, K, p, com_mat) {
  sigma02 <- exp(lsigma02)
  # calculate the lBF
  tmp1 <- log(sqrt(s2_list / (sigma02 + s2_list)))
  tmp2 <- z2_list / 2 * sigma02 / (sigma02 + s2_list)
  lBF_list <- tmp1 + tmp2
  # combine
  lBF <- rep(0, (2 ^ K - 1) * p + 1)
  for (iter_com in seq_len(nrow(com_mat))) {
    lBF[(1 + p * (iter_com - 1)) : (p * iter_com)] <- 
      rowSums(lBF_list[, which(com_mat[iter_com, ] == 1), drop = FALSE])
  }
  maxlBF <- max(lBF)
  wBF <- exp(lBF - maxlBF)
  wBF_sum <- sum(prior_pi * wBF)
  return(- maxlBF - log(wBF_sum))
}

sigma0_opt_multi <- function(lsigma02_int, prior_pi, z2_list, s2_list, b_hat_list, K, p, com_mat) {
  tmp1 <- lBF_model_multi(lsigma02 = lsigma02_int, prior_pi, z2_list, s2_list, K, p, com_mat)
  lsigma02 <- optim(par = log(max(b_hat_list ^ 2 - s2_list, 1)), fn = lBF_model_multi, 
                    method = "Brent", lower = -30, upper = 15, prior_pi = prior_pi, z2_list = z2_list, 
                    s2_list = s2_list, K = K, p = p, com_mat = com_mat)$par
  tmp2 <- lBF_model_multi(lsigma02 = lsigma02, prior_pi, z2_list, s2_list, K, p, com_mat)
  if (tmp2 < tmp1) {
    return(exp(lsigma02))
  } else{
    return(exp(lsigma02_int))
  }
}

## Calculate the KL divergence
KL_fun_multi <- function(dta_list, Y_list, sigma2, b_list, b2_list, lBF, K, n, l) {
  tmp1 <- 0
  tmp4 <- 0
  tmp3 <- K * n / 2 * log(2 * pi * sigma2)
  for (iter_K in seq_len(K)) {
    tmp1 <- tmp1 + sum(dnorm(Y_list[, iter_K], mean = 0, sd = sqrt(sigma2), log = TRUE))
    tmp4 <- tmp4 + 1 / (2 * sigma2) * 
      (crossprod(Y_list[, iter_K]) - 2 * crossprod(Y_list[, iter_K], dta_list[[iter_K]]$X_scale %*% b_list[[iter_K]][, l]) + 
         sum(dta_list[[iter_K]]$X_scale2 %*% b2_list[[iter_K]][, l]))
  }
  return(tmp1 + tmp3 + tmp4 + lBF)
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


