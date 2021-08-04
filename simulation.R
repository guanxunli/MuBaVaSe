source("Multi_dataset.R")

#### Initialization
## Define parameters
iter_sim_max <- 500
out_res <- data.frame("dat1_md_sens" = rep(NA, iter_sim_max), "dat1_md_prec" = rep(NA, iter_sim_max),
                      "dat2_md_sens" = rep(NA, iter_sim_max), "dat2_md_prec" = rep(NA, iter_sim_max),
                      "dat1_sd_sens" = rep(NA, iter_sim_max), "dat1_sd_prec" = rep(NA, iter_sim_max),
                      "dat2_sd_sens" = rep(NA, iter_sim_max), "dat2_sd_prec" = rep(NA, iter_sim_max))
n <- 500
p <- 1000
p_c <- 30
p_1 <- 5
p_2 <- 5
sigma <- 1
sigma0 <- 0.6
set.seed(2021)
for (iter_sim in seq_len(iter_sim_max)) {
  if (iter_sim %% 50 == 0) print(iter_sim)
  ## Generate data
  index_c <- sample(seq_len(p), size = p_c, replace = FALSE)
  index_1 <- sample(setdiff(seq_len(p), index_c), size = p_1, replace = FALSE)
  index_2 <- sample(setdiff(seq_len(p), c(index_1, index_c)), size = p_2, replace = FALSE)
  
  b_1 <- rep(0, p)
  b_1[c(index_c, index_1)] <- rnorm(p_c + p_1, mean = 0, sd = sigma0)
  b_2 <- rep(0, p)
  b_2[c(index_c, index_2)] <- rnorm(p_c + p_2, mean = 0, sd = sigma0)
  
  X_1 <- matrix(rnorm(p * n), nrow = n, ncol = p)
  X_2 <- matrix(rnorm(p * n), nrow = n, ncol = p)
  Y_1 <- X_1 %*% b_1 + rnorm(n, sd = sigma)
  Y_2 <- X_2 %*% b_2 + rnorm(n, sd = sigma)
  
  #### two data set at the same times
  res <- sum_single_effect_multi(X_1 = X_1, X_2 = X_2, Y_1 = Y_1, Y_2 = Y_2)
  res1 <- which(abs(round(rowSums(res$b_mat_1),4)) > 0)
  res2 <- which(abs(round(rowSums(res$b_mat_2),4)) > 0)
  out_res[iter_sim, 1] <- length(intersect(res1, c(index_1, index_c))) / (p_c + p_1)
  out_res[iter_sim, 2] <- length(intersect(res1, c(index_1, index_c))) / length(res1)
  out_res[iter_sim, 3] <- length(intersect(res2, c(index_2, index_c))) / (p_c + p_2)
  out_res[iter_sim, 4] <- length(intersect(res2, c(index_2, index_c))) / length(res2)
  
  #### Single data set
  ## data set 1
  res <- susieR::susie(X = X_1, y = Y_1)
  res1 <- which(abs(round(colSums(res$alpha * res$mu),4)) > 0)
  out_res[iter_sim, 5] <- length(intersect(res1, c(index_1, index_c))) / (p_c + p_1)
  out_res[iter_sim, 6] <- length(intersect(res1, c(index_1, index_c))) / length(res1)
  ## data set 2
  res <- susieR::susie(X = X_2, y = Y_2)
  res2 <- which(abs(round(colSums(res$alpha * res$mu),4)) > 0)
  out_res[iter_sim, 7] <- length(intersect(res2, c(index_1, index_c))) / (p_c + p_1)
  out_res[iter_sim, 8] <- length(intersect(res2, c(index_1, index_c))) / length(res2)
}