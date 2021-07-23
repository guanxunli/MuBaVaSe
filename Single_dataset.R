n <- 100
p <- 1000
sigma <- 0.1
sigma0 <- 0.6
L <- 10
set.seed(202107)
## Generate data
index_t <- sample(seq_len(p), size = L, replace = FALSE)
b <- rep(0, p)
b[index_t] <- rnorm(L, mean = 0, sd = sigma0)
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- X %*% b + rnorm(n, sd = sigma)

#### Define functions
## get sigma0
lBF_model_single <- function(lsigma02, prior_pi, z2, s2) {
  sigma02 <- exp(lsigma02)
  tmp1 <- log(sqrt(s2 / (sigma02 + s2)))
  tmp2 <- z2 / 2 * sigma02 / (sigma02 + s2)
  lBF <- tmp1 + tmp2
  maxlBF <- max(lBF)
  wBF <- exp(lBF - maxlBF)
  wBF_sum <- sum(prior_pi * wBF)
  return(- maxlBF - log(wBF_sum))
}

sigma0_opt_single <- function(lsigma02_int, prior_pi, z2, s2, b_hat) {
  tmp1 <- lBF_model_single(lsigma02 = lsigma02_int, prior_pi = prior_pi, z2 = z2, s2 = s2)
  lsigma02 <- optim(par = log(max(c(b_hat^2 - s2, 1))), fn = lBF_model_single, method = "Brent", lower = -30, upper = 15,
                    prior_pi = prior_pi, z2 = z2, s2 = s2)$par
  tmp2 <- lBF_model_single(lsigma02 = lsigma02, prior_pi = prior_pi, z2 = z2, s2 = s2)
  if (tmp2 < tmp1) {
    return(exp(lsigma02))
  } else{
    return(exp(lsigma02_int))
  }
}

## Calculate the KL divergence
KL_fun_single <- function(X, Y, sigma2, b, b2, lsigma02, prior_pi, z2, s2) {
  n <- length(Y)
  tmp1 <- sum(dnorm(Y, mean = 0, sd = sqrt(sigma2), log = TRUE))
  tmp2 <- -lBF_model_single(lsigma02 = lsigma02, prior_pi = prior_pi, z2 = z2, s2 = s2)
  tmp3 <- n / 2 * log(2 * pi * sigma2)
  tmp4 <- 1 / (2 * sigma2) * (crossprod(Y) - 2 * crossprod(Y, scale(X) %*% b) + sum(scale(X)^2 %*% b2))
  return(tmp1 + tmp2 + tmp3 + tmp4)
}

## Calculate ERSS
ERSS_fun_single <- function(X, Y, b_mat, b2_mat) {
  mu_lmat <- scale(X) %*% b_mat
  mu2_lmat <- scale(X)^2 %*% b2_mat
  mu_pred <- rowSums(mu_lmat)
  res_tmp <- sum((Y - mu_pred)^2)
  var_sum <- sum(mu2_lmat - mu_lmat^2)
  return(res_tmp + var_sum)
}

## Main functions
sum_single_effect_single <- function(X, Y, sigma2_int = NULL, sigma02_int = NULL, prior_pi = NULL, 
                                     L = NULL, itermax = 100, tol = 1e-4) {
  # Initialization
  p <- ncol(X)
  n <- nrow(X)
  mean_Y <- mean(Y)
  Y <- Y - mean_Y
  X_scale <- scale(X)
  XtX <- crossprod(X_scale)
  
  if (is.null(sigma2_int)) sigma2_int <- as.numeric(var(Y))
  if (is.null(sigma02_int)) sigma02_int <- 0.2 * sigma2_int
  if (is.null(prior_pi)) prior_pi <- rep(1/p ,p)
  if (is.null(L)) L <- min(10, ncol(X))
  
  ELBO <- rep(NA, itermax + 1)
  ELBO[1] <- -Inf
  sigma2 <- sigma2_int
  sigma02_vec <- rep(sigma02_int, L)
  # Save matrix
  b_mat <- matrix(0, nrow = p, ncol = L)
  b2_mat <- matrix(0, nrow = p, ncol = L)
  # Begin iteration
  for (iter in seq_len(itermax)) {
    res <- Y - X_scale %*% rowSums(b_mat)
    KL_div <- 0
    for (l in seq_len(L)) {
      res_tmp <- res + X_scale %*% b_mat[, l]
      # update parameters
      XtY <- crossprod(X_scale, res_tmp)
      b_hat <- XtY / diag(XtX)
      s2 <- sigma2 / diag(XtX)
      z2 <- b_hat^2 / s2
      # calculate sigma0
      lsigma02_int <- sigma02_vec[l]
      sigma02 <- sigma0_opt_single(lsigma02_int, prior_pi = prior_pi, z2 = z2, s2 = s2, b_hat = b_hat)
      sigma02_vec[l] <- sigma02
      # Get Bayesian Factor
      tmp1 <- log(sqrt(s2 / (sigma02 + s2))) 
      tmp2 <- z2 / 2 * sigma02 / (sigma02 + s2)
      lBF <- tmp1 + tmp2
      maxlBF <- max(lBF)
      wBF <- exp(lBF - maxlBF)
      wBF_sum <- sum(prior_pi * wBF)
      # Get posterior
      post_alpha <- prior_pi * wBF / wBF_sum
      post_sigma2 <- 1 / (1/s2 + 1/sigma02)
      post_mu <- post_sigma2 / s2 * b_hat
      # Calculate posterior mean
      b_mat[, l] <- post_alpha * post_mu
      b2_mat[, l] <- post_alpha * (post_mu^2 + post_sigma2)
      KL_div <- KL_div + KL_fun_single(X = X, Y = res_tmp, sigma2 = sigma2, b = b_mat[, l], b2 = b2_mat[, l],
                                       lsigma02 = log(sigma02), prior_pi = prior_pi, z2 = z2, s2 = s2)
      res <- res_tmp - scale(X) %*% b_mat[, l]
    }
    ERSS <- ERSS_fun_single(X = X, Y = Y, b_mat = b_mat, b2_mat = b2_mat)
    ELBO[iter + 1] <- -n / 2 * log(2 * pi * sigma2) - 1 / (2 * sigma2) * ERSS + KL_div
    sigma2 <- ERSS / n
    if (abs(ELBO[iter + 1] -   ELBO[iter]) < 1e-4) break
  }
  ELBO <- as.numeric(na.omit(ELBO[-1]))
  # return results
  res <- list()
  res$ELBO <- ELBO
  res$b_mat <- b_mat
  return(res)
}

#### check results
## package
res <- susieR::susie(X = X, y = Y, L = L)
res1 <- which(abs(round(colSums(res$alpha * res$mu), 4)) > 0)
length(intersect(res1, index_t)) / L
length(intersect(res1, index_t)) / length(res1)

## My code
res <- sum_single_effect_single(X = X, Y = Y, L = L)
res1 <- which(abs(round(rowSums(res$b_mat), 4)) > 0)
length(intersect(res1, index_t)) / L
length(intersect(res1, index_t)) / length(res1)