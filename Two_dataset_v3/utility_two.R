library(Matrix)
#### Define functions
## get sigma0
lBF_model_two <- function(lsigma02, prior_pi, z2_1, s2_1, z2_2, s2_2) {
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
  lBF <- c(lBF_1, lBF_2, lBF_1 + lBF_2, 0)
  maxlBF <- max(lBF)
  wBF <- exp(lBF - maxlBF)
  wBF_sum <- sum(prior_pi * wBF)
  return(-maxlBF - log(wBF_sum))
}

sigma0_opt_two <- function(lsigma02_int, prior_pi, z2_1, s2_1, z2_2, s2_2, b_hat_1, b_hat_2) {
  tmp1 <- lBF_model_two(
    lsigma02 = lsigma02_int, prior_pi = prior_pi, z2_1 = z2_1,
    s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2
  )
  lsigma02 <- optim(
    par = log(max(c(b_hat_1^2 - s2_1, 1, b_hat_2^2 - s2_2))), fn = lBF_model_two,
    method = "Brent", lower = -30, upper = 15, prior_pi = prior_pi, z2_1 = z2_1,
    s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2
  )$par
  tmp2 <- lBF_model_two(
    lsigma02 = lsigma02, prior_pi = prior_pi, z2_1 = z2_1,
    s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2
  )
  if (tmp2 < tmp1) {
    return(exp(lsigma02))
  } else {
    return(exp(lsigma02_int))
  }
}

sigma0_opt_two_test <- function(lsigma02_int, prior_pi, z2_1, s2_1, z2_2, s2_2, b_hat_1, b_hat_2) {
  tmp1 <- lBF_model_two(
    lsigma02 = lsigma02_int, prior_pi = prior_pi, z2_1 = z2_1,
    s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2
  )
  lsigma02 <- optimize(
    lBF_model_two,
    lower = -30, upper = 15, prior_pi = prior_pi, z2_1 = z2_1,
    s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2
  )$minimum
  tmp2 <- lBF_model_two(
    lsigma02 = lsigma02, prior_pi = prior_pi, z2_1 = z2_1,
    s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2
  )
  if (tmp2 < tmp1) {
    return(exp(lsigma02))
  } else {
    return(exp(lsigma02_int))
  }
}

## Calculate the -KL divergence
KL_fun_two <- function(X_scale_1, X_scale2_1, Y_1, X_scale_2, Y_2, X_scale2_2,
                       sigma2, b_1, b2_1, b_2, b2_2, lBF) {
  n1 <- length(Y_1)
  n2 <- length(Y_2)
  tmp1_1 <- sum(dnorm(Y_1, mean = 0, sd = sqrt(sigma2), log = TRUE))
  tmp1_2 <- sum(dnorm(Y_2, mean = 0, sd = sqrt(sigma2), log = TRUE))
  tmp3_1 <- n1 / 2 * log(2 * pi * sigma2)
  tmp3_2 <- n2 / 2 * log(2 * pi * sigma2)
  tmp4_1 <- 1 / (2 * sigma2) * (crossprod(Y_1) - 2 * crossprod(Y_1, X_scale_1 %*% b_1) + sum(X_scale2_1 %*% b2_1))
  tmp4_2 <- 1 / (2 * sigma2) * (crossprod(Y_2) - 2 * crossprod(Y_2, X_scale_2 %*% b_2) + sum(X_scale2_2 %*% b2_2))
  return(tmp1_1 + tmp1_2 + lBF + tmp3_1 + tmp3_2 + tmp4_1 + tmp4_2)
}

KL_fun_two_graph <- function(XtY_1, XtXbeta_use_1, X_scale2_1,
                             XtY_2, XtXbeta_use_2, X_scale2_2,
                             sigma2, b_1, b2_1, b_2, b2_2, lBF) {
  tmp1_1 <- -2 * (crossprod(b_1, XtY_1) - crossprod(b_1, XtXbeta_use_1))
  tmp1_2 <- -2 * (crossprod(b_2, XtY_2) - crossprod(b_2, XtXbeta_use_2))
  tmp2_1 <- sum(X_scale2_1 %*% b2_1)
  tmp2_2 <- sum(X_scale2_2 %*% b2_2)
  return(lBF + 1 / (2 * sigma2) * (tmp1_1 + tmp1_2 + tmp2_1 + tmp2_2))
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

## sampling function
sampling_fun <- function(X_1, Y_1, X_2, Y_2, scale_x, intercept,
                         lprior_vec, sigma2, alpha_mat, sigma02_vec) {
  ## Initialization
  p <- ncol(X_1)
  if (p != ncol(X_2)) stop("The number of features should be same!")
  n1 <- nrow(X_1)
  n2 <- nrow(X_2)
  L <- ncol(alpha_mat)
  ## data preprocess
  X_scale_1 <- scale(X_1, center = intercept, scale = scale_x)
  X_scale_2 <- scale(X_2, center = intercept, scale = scale_x)
  if (intercept) {
    mean_Y_1 <- mean(Y_1)
    mean_Y_2 <- mean(Y_2)
  } else {
    mean_Y_1 <- 0
    mean_Y_2 <- 0
  }
  Y_1 <- Y_1 - mean_Y_1
  Y_2 <- Y_2 - mean_Y_2
  ## return list
  res <- list()
  res$alpha_mat <- alpha_mat
  res$sigma2 <- sigma2
  res$sigma02_vec <- sigma02_vec
  res$index1_select <- NULL
  res$loglikelihood_1 <- sum(dnorm(Y_1, sd = sqrt(res$sigma2), log = TRUE))
  res$index2_select <- NULL
  res$loglikelihood_2 <- sum(dnorm(Y_2, sd = sqrt(res$sigma2), log = TRUE))
  res$lprior <- 0
  ## variable selection
  index_all <- apply(res$alpha_mat, 2, function(x) {
    index_use <- sample(seq_len(length(x)), 1, prob = x)
    return(c(index_use, x[index_use]))
  })
  res$lpropose <- sum(log(index_all[2, ]))
  index_select <- which(index_all[1, ] < (3 * p + 1))
  # calculate the likelihood and index
  if (length(index_select) > 0) {
    sigma02_select <- res$sigma02_vec[index_select]
    index_select <- index_all[1, ][index_select]
    # index
    index_1 <- which(index_select < p + 1)
    index_c <- which(index_select > 2 * p)
    index_2 <- intersect(which(index_select > p), which(index_select < 2 * p + 1))
    res$lprior <- log(choose(L, length(index_1))) + log(choose(L - length(index_1), length(index_2))) +
      log(choose(L - length(index_1) - length(index_2), length(index_c))) +
      (length(index_1) + length(index_2)) * lprior_vec[1] + length(index_c) * lprior_vec[2]
    # select index
    sigma02_select1 <- sigma02_select[c(index_1, index_c)]
    index1_select <- c(index_select[index_1], index_select[index_c] - 2 * p)
    sigma02_select2 <- sigma02_select[c(index_2, index_c)]
    index2_select <- c(index_select[index_2] - p, index_select[index_c] - 2 * p)
    # return likelihood
    if (length(index1_select) > 0) {
      Sigma_1_inverse <- diag(1, n1) - X_scale_1[, index1_select, drop = FALSE] %*%
        solve(
          crossprod(X_scale_1[, index1_select, drop = FALSE]) +
            res$sigma2 * diag(1 / sigma02_select1, nrow = length(sigma02_select1)),
          t(X_scale_1[, index1_select, drop = FALSE])
        )
      loglikelihood_1 <- -n1 / 2 * log(2 * pi) + 1 / 2 * log(det(Sigma_1_inverse)) -
        1 / 2 * crossprod(Y_1, Sigma_1_inverse %*% Y_1)
      res$index1_select <- index1_select
      res$loglikelihood_1 <- loglikelihood_1
    }
    if (length(index2_select) > 0) {
      Sigma_2_inverse <- diag(1, n2) - X_scale_2[, index2_select, drop = FALSE] %*%
        solve(
          crossprod(X_scale_2[, index2_select, drop = FALSE]) +
            res$sigma2 * diag(1 / sigma02_select2, nrow = length(sigma02_select2)),
          t(X_scale_2[, index2_select, drop = FALSE])
        )
      loglikelihood_2 <- -n2 / 2 * log(2 * pi) + 1 / 2 * log(det(Sigma_2_inverse)) -
        1 / 2 * crossprod(Y_2, Sigma_2_inverse %*% Y_2)
      res$index2_select <- index2_select
      res$loglikelihood_2 <- loglikelihood_2
    }
  }
  return(res)
}