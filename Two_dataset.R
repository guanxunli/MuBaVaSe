#### Define functions
## For one data set
source("One_dataset.R")
## posterior for b
single_effect2 <- function(X1, Y1, X2, Y2, sigma1, sigma2, sigma0, pi) {
  p <- ncol(X1)
  # update b for data set 1
  b_hat1 <- crossprod(X1, Y1) / diag(crossprod(X1))
  s_square1 <- diag(sigma1^2 / crossprod(X1))
  sigma_square_p1 <- 1 / (1/s_square1 + 1/sigma0^2)
  mu_p1 <- sigma_square_p1 / s_square1 * b_hat1
  # update b for data set 2
  b_hat2 <- crossprod(X2, Y2) / diag(crossprod(X2))
  s_square2 <- diag(sigma2^2 / crossprod(X2))
  sigma_square_p2 <- 1 / (1/s_square2 + 1/sigma0^2)
  mu_p2 <- sigma_square_p2 / s_square2 * b_hat2
  # update BF
  z_square1 <- mu_p1^2 / s_square1
  z_square2 <- mu_p2^2 / s_square2
  BF1 <- exp(log(sqrt(s_square1 / (sigma0^2 + s_square1))) + 
              z_square1 / 2 * sigma0^2 / (sigma0^2 + s_square1))
  BF2 <- exp(log(sqrt(s_square2 / (sigma0^2 + s_square2))) + 
               z_square2 / 2 * sigma0^2 / (sigma0^2 + s_square2))
  index_inf1 <- which(BF1 == Inf)
  index_inf2 <- which(BF2 == Inf)
  if (length(intersect(index_inf1, index_inf2)) > 0) {
    BF <- rep(0, 3 * p)
    BF[intersect(index_inf1, index_inf2) + 2 * p] <- 1
  } else if (length(index_inf1) > 0) {
    BF <- rep(0, 3 * p)
    BF[index_inf1] <- 1
    if (length(index_inf2) > 0) {
      BF[index_inf2 + p] <- 1
    }
  } else if (length(index_inf2) > 0) {
    BF <- rep(0, 3 * p)
    BF[index_inf2 + p] <- 1
  } else{
    BF <- c(BF1, BF2, BF1 * BF2)
    index_inf <- which(BF == Inf)
    if (length(index_inf) > 0) {
      BF <- rep(0, 3 * p)
      BF[index_inf] <- 1
    }
  }
  alpha_p <- pi * BF / sum(pi * BF)
  # return results
  res <- list()
  res$mu1 <- mu_p1
  res$mu2 <- mu_p2
  res$sigma1 <- sigma_square_p1
  res$sigma2 <- sigma_square_p2
  res$alpha <- alpha_p
  return(res)
}

sum_single_effect2 <- function(X1, Y1, X2, Y2, sigma1, sigma2, sigma0, pi, L, itermax = 10) {
  p <- ncol(X1)
  # Initialization
  b_mat_new1 <- matrix(0, nrow = p, ncol = L)
  b_mat_old1 <- matrix(0, nrow = p, ncol = L)
  mu_mat1 <- matrix(0, nrow = p, ncol = L)
  sigma_mat1 <- matrix(0, nrow = p, ncol = L)
  b_mat_new2 <- matrix(0, nrow = p, ncol = L)
  b_mat_old2 <- matrix(0, nrow = p, ncol = L)
  mu_mat2 <- matrix(0, nrow = p, ncol = L)
  sigma_mat2 <- matrix(0, nrow = p, ncol = L)
  # alpha matrix
  alpha_mat <- matrix(0, nrow = 3 * p, ncol = L)
  # Iteration
  for (iter in seq_len(itermax)) {
    for (l in seq_len(L)) {
      res_l1 <- Y1 - X1 %*% rowSums(b_mat_new1[, -l])
      res_l2 <- Y2 - X2 %*% rowSums(b_mat_new2[, -l])
      # Single effect
      res_pos <- single_effect2(X1 = X1, Y1 = res_l1, X2 = X2, Y2 = res_l2, 
                               sigma1 = sigma1, sigma2 = sigma2, sigma0 = sigma0, pi = pi)
      mu_mat1[, l] <- res_pos$mu1
      sigma_mat1[, l] <- res_pos$sigma1
      mu_mat2[, l] <- res_pos$mu2
      sigma_mat2[, l] <- res_pos$sigma2
      alpha_mat[, l] <- res_pos$alpha
      b_mat_new1[, l] <- (res_pos$alpha[1:p] + res_pos$alpha[(2 * p + 1) : (3 * p)]) * 
        res_pos$mu1
      b_mat_new2[, l] <- (res_pos$alpha[(p + 1) : (2 * p)] + res_pos$alpha[(2 * p + 1) : (3 * p)]) * 
        res_pos$mu2
    }
    if (sum((b_mat_old1 - b_mat_new1)^2 + (b_mat_old2 - b_mat_new2)^2) < 1e-4) {
      break
    } else {
      b_mat_old1 <- b_mat_new1
      b_mat_old2 <- b_mat_new2
    }
  }
  res <- list()
  res$mu1 <- mu_mat1
  res$sigma1 <- sigma_mat1
  res$mu2 <- mu_mat2
  res$sigma2 <- sigma_mat2
  res$alpha <- alpha_mat
  return(res)
}

#### Initialization
## Define parameters
iter_sim_max <- 20
out_res <- data.frame("dat1_md_sens" = rep(NA, iter_sim_max), "dat1_md_prec" = rep(NA, iter_sim_max),
                      "dat2_md_sens" = rep(NA, iter_sim_max), "dat2_md_prec" = rep(NA, iter_sim_max),
                      "dat1_sd_sens" = rep(NA, iter_sim_max), "dat1_sd_prec" = rep(NA, iter_sim_max),
                      "dat2_sd_sens" = rep(NA, iter_sim_max), "dat1_sd_prec" = rep(NA, iter_sim_max))
n <- 100
p <- 1000
p_c <- 45
p_1 <- 5
p_2 <- 5
sigma1 <- 0.1
sigma2 <- 0.1
sigma0 <- 0.6
r <- 0.2
q <- 0.05
set.seed(2021)
for (iter_sim in seq_len(iter_sim_max)) {
  if (iter_sim %% 10 == 0) print(iter_sim)
  ## Generate data
  index_c <- sample(seq_len(p), size = p_c, replace = FALSE)
  index_1 <- sample(setdiff(seq_len(p), index_c), size = p_1, replace = FALSE)
  index_2 <- sample(setdiff(seq_len(p), c(index_1, index_c)), size = p_2, replace = FALSE)
  
  b_1 <- rep(0, p)
  b_1[c(index_c, index_1)] <- rnorm(p_c + p_1, mean = 0, sd = sigma0)
  b_2 <- rep(0, p)
  b_2[c(index_c, index_1)] <- rnorm(p_c + p_2, mean = 0, sd = sigma0)
  
  X_1 <- matrix(rnorm(p * n), nrow = n, ncol = p)
  X_2 <- matrix(rnorm(p * n), nrow = n, ncol = p)
  Y_1 <- X_1 %*% b_1 + rnorm(n, sd = sigma1)
  Y_2 <- X_2 %*% b_2 + rnorm(n, sd = sigma2)
  
  #### two data set at the same times
  pi <- c(rep(q, 2 * p), rep(r, p))
  pi <- pi / sum(pi)
  res <- sum_single_effect2(X1 = X_1, X2 = X_2, Y1 = Y_1, Y2 = Y_2, 
                           sigma1 = sigma1, sigma2 = sigma2, sigma0 = sigma0, pi = pi,
                           L = (p_1 + p_2 + p_c), itermax = 1000)
  res1 <- which(abs(round(rowSums((res$alpha[1:p, ] + res$alpha[(2 * p + 1) : (3 * p), ]) * 
                                    res$mu1),4)) > 0)
  res2 <- which(abs(round(rowSums((res$alpha[(p + 1): (2 * p), ] + res$alpha[(2 * p + 1) : (3 * p), ]) * 
                                    res$mu2),4)) > 0)
  out_res[iter_sim, 1] <- length(intersect(res1, c(index_1, index_c))) / (p_c + p_1)
  out_res[iter_sim, 2] <- length(intersect(res1, c(index_1, index_c))) / length(res1)
  out_res[iter_sim, 3] <- length(intersect(res2, c(index_2, index_c))) / (p_c + p_2)
  out_res[iter_sim, 4] <- length(intersect(res2, c(index_2, index_c))) / length(res2)
  
  #### Single data set
  pi <- rep(1 / p , p)
  ## data set 1
  res <- sum_single_effect(X = X_1, Y = Y_1, sigma = sigma1, sigma0 = sigma0, 
                           pi = pi, L = (p_1 + p_c), itermax = 1000)
  res1 <- which(abs(round(rowSums(res$alpha * res$mu), 4)) > 0)
  out_res[iter_sim, 5] <- length(intersect(res1, c(index_1, index_c))) / (p_c + p_1)
  out_res[iter_sim, 6] <- length(intersect(res1, c(index_1, index_c))) / length(res1)
  ## data set 2
  res <- sum_single_effect(X = X_2, Y = Y_2, sigma = sigma2, sigma0 = sigma0, 
                           pi = pi, L = (p_2 + p_c), itermax = 1000)
  res2 <- which(abs(round(rowSums(res$alpha * res$mu), 4)) > 0)
  out_res[iter_sim, 7] <- length(intersect(res2, c(index_1, index_c))) / (p_c + p_1)
  out_res[iter_sim, 8] <- length(intersect(res2, c(index_1, index_c))) / length(res2)
}

