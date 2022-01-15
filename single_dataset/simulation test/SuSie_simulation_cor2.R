########################### Do simulations ###########################
library(parallel)
library(foreach)
library(doParallel)
library(doRNG)
library(ggplot2)
## load functions
source("single_dataset/sum_single_effect_single_null.R")
# scale_x = FALSE
# intercept = FALSE
# sigma02_int = NULL
# sigma2_int = NULL
# prior_null = NULL
# L = 3
# itermax = 100
# tol = 1e-4
# sigma0_low_bd = 1e-8
# residual_variance_lowerbound = NULL
sum_single_effect_single_null_simu <- function(X, Y, scale_x = TRUE, intercept = TRUE,
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
  if (is.null(residual_variance_lowerbound)) residual_variance_lowerbound <- 1e-4
  
  ## data preprocess
  # scale
  if (scale_x) {
    X_scale <- scale(X)
  } else {
    X_scale <- X
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
  Y <- Y - mean_Y
  
  # Initialize prior
  if (is.null(prior_null)) {
    prior_null <- 1 - 1 / (p^0.5)
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
    KL_div_vec <- rep(0, L)
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
      ## calculate the minus KL divergence
      KL_div_vec[l] <- KL_fun_single(
        X_scale = X_scale, Y = res_tmp, X_scale2 = X_scale2, sigma2 = sigma2,
        b = b_mat[, l], b2 = b2_mat[, l], lBF = lBF_model
      )
      res <- res_tmp - X_scale %*% b_mat[, l]
      cat(
        "Iteration:", iter, "l:", l, "Posterior:", round(post_alpha[1:3], 4),
        "NULL model:", round(post_alpha[length(post_alpha)], 4), "\n"
      )
    }
    # calculate ELBO
    ERSS <- ERSS_fun_single(X_scale = X_scale, X_scale2 = X_scale2, Y = Y, b_mat = b_mat, b2_mat = b2_mat)
    ELBO[iter + 1] <- -n / 2 * log(2 * pi * sigma2) - 1 / (2 * sigma2) * ERSS + sum(KL_div_vec)
    # estimate sigma2
    sigma2 <- max(ERSS / n, residual_variance_lowerbound)
    if (ELBO[iter + 1] - ELBO[iter] < 1e-4) break
  }
  ELBO <- as.numeric(na.omit(ELBO[-1]))
  # select effect index
  index_L <- which(sigma02_vec > sigma0_low_bd)
  ## return results
  res <- list()
  res$ELBO <- ELBO
  res$sigma2 <- sigma2
  res$sigma02_vec <- sigma02_vec
  res$alpha_null <- alpha_null
  res$index_L <- index_L
  res$KL_div_vec <- -KL_div_vec
  
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

#### Define parameters
sigma <- 1
n <- 500
p <- 1000
L <- 2
index_t <- seq_len(L)
Sigma <- diag(1, nrow = p)
Sigma[1, 3] <- Sigma[3, 1] <- 0.99
Sigma_chol <- chol(Sigma)

## Do simulation
n_iter <- 100
cl <- makeCluster(5)
registerDoParallel(cl)
set.seed(2021)
out_res <- foreach(iter = seq_len(n_iter)) %dorng% {
  ## Generate data
  b <- rep(0, p)
  b_tmp <- sample(c(-1, 1), size = L, replace = TRUE) *
    runif(L, min = 0.1, max = 1)
  b[1] <- b_tmp[which.max(abs(b_tmp))]
  b[2] <- b_tmp[-which.max(abs(b_tmp))]
  # b[1] <- b_tmp[-which.max(abs(b_tmp))]
  # b[2] <- b_tmp[which.max(abs(b_tmp))]
  X <- t(Sigma_chol) %*% matrix(rnorm(n * p), nrow = p, ncol = n)
  X <- t(X)
  X_scale <- scale(X, scale = FALSE)
  Y <- X %*% b + rnorm(n, sd = sigma)
  ## Do variable selection
  res <- sum_single_effect_single_null(
    X = X_scale, Y = Y, L = L + 1,
    scale_x = FALSE, intercept = FALSE
  )
  res1 <- which(res$alpha > 0.5)
  list(b = b[1:3], X = X_scale, Y = Y,
       acc = length(intersect(res1, index_t)) / L,
       index_find = res1, alpha = res$alpha[1:3])
}
stopCluster(cl)

res_df <- data.frame("acc" = rep(NA, n_iter), "min_coef" = rep(NA, n_iter),
                     "corXY1" = rep(NA, n_iter),
                     "corXY2" = rep(NA, n_iter),
                     "corXY3" = rep(NA, n_iter),
                     "corXY_small" = rep(NA, n_iter))
for (iter in seq_len(n_iter)) {
  res_tmp <- out_res[[iter]]
  res_df$acc[iter] <- ifelse(res_tmp$acc == 1, 1, 0)
  res_df$min_coef[iter] <- min(abs(res_tmp$b[1:2]))
  res_df$corXY1[iter] <- cor(res_tmp$Y, res_tmp$X[, 1])
  res_df$corXY2[iter] <- cor(res_tmp$Y, res_tmp$X[, 2])
  res_df$corXY3[iter] <- cor(res_tmp$Y, res_tmp$X[, 3])
  res_df$corXY_small[iter] <- min(abs(c(res_df$corXY1[iter],
                                        res_df$corXY2[iter],
                                        res_df$corXY3[iter])))
}
sum(res_df$acc)

## Show figures
res_df$acc <- as.factor(res_df$acc)
# the minium coefficients
ggplot(data = res_df, aes(x = min_coef, color = acc, fill = acc)) +
  geom_histogram(position = "dodge", bins = 30) +
  scale_colour_manual(values = c("red", "green4")) +
  theme_bw()

index_wrong <- which(res_df$acc == 0)
check_df <- apply(res_df[index_wrong, -1], 2, round, digits = 4)
check_df <- as.data.frame(check_df)
check_df$find <- NA
check_df$small_coef_index <- NA
check_df$small_cor_index <- NA
for (iter in seq_len(length(index_wrong))) {
  res_tmp <- out_res[[index_wrong[iter]]]
  if (length(res_tmp$index_find) == 0) {
    check_df$find[iter] <- 0
  } else if (length(res_tmp$index_find) == 2) {
    check_df$find[iter] <- 3
  } else {
    check_df$find[iter] <- res_tmp$index_find
  }
  b_tmp <- res_tmp$b
  check_df$small_coef_index[iter] <- which.min(abs(b_tmp[1:2]))
  cor_tmp <- c(
    cor(res_tmp$Y, res_tmp$X[, 1]),
    cor(res_tmp$Y, res_tmp$X[, 2]),
    cor(res_tmp$Y, res_tmp$X[, 3])
  )
  check_df$small_cor_index[iter] <- which.min(abs(cor_tmp))
}
print(check_df)

# the correlation figure
ggplot(data = res_df, aes(x = abs(corXY1), y = abs(corXY2), color = acc, fill = acc)) +
  geom_point() +
  scale_colour_manual(values = c("red", "green4")) +
  theme_bw()
# correlation 1 and 3
ggplot(data = res_df, aes(x = abs(corXY_small), y = abs(corXY3), color = acc, fill = acc)) +
  geom_point() +
  scale_colour_manual(values = c("red", "green4")) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()
ggplot(data = res_df, aes(x = abs(corXY3), y = abs(corXY1), color = acc, fill = acc)) +
  geom_point() +
  scale_colour_manual(values = c("red", "green4")) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()

########################### check specific examples###########################
b <- out_res[[83]]$b
print(b)
X <- out_res[[83]]$X
Y <- out_res[[83]]$Y
round(c(cor(Y, X[, 1]), cor(Y, X[, 2]), cor(Y, X[, 3])), 4)
res <- sum_single_effect_single_null_simu(X = X, Y = Y, L = L + 1, 
                                          scale_x = FALSE, intercept = FALSE)
round(res$alpha[1:3], 4)
