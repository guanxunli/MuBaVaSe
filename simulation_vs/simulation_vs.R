#### Initialization
## Define parameters
source("sum_single_effect_mult.R")
source("single_dataset/sum_single_effect_single_null.R")
## simulation function
simu_vs_fun <- function(K, n, p, p_c, p_s, sigma, sigma0) {
  ## Generate data
  index_c <- sample(seq_len(p), size = p_c, replace = FALSE)
  b <- list()
  index_s <- list()
  dta_list <- list()
  for (k in seq_len(K)) {
    dta_list[[k]] <- list()
    index_s[[k]] <- sample(setdiff(seq_len(p), index_c), size = p_s, replace = FALSE)
    b[[k]] <- rep(0, p)
    b[[k]][c(index_c, index_s[[k]])] <- rnorm(p_c + p_s, mean = 0, sd = sigma0)
    # b[[k]][c(index_c, index_s[[k]])] <- sample(c(-1, 1), size = p_c + p_s, replace = TRUE) *
    #   runif(p_c + p_s, min = 0.1, max = 1)
    dta_list[[k]]$X <- scale(matrix(rnorm(p * n), nrow = n, ncol = p))
    y_tmp <- dta_list[[k]]$X %*% b[[k]]
    dta_list[[k]]$Y <- y_tmp + rnorm(n, sd = sigma)
  }
  #### multiple data sets at the same times
  res_c <- sum_single_effect_mult(dta_list,
    scale_x = TRUE, intercept = TRUE,
    sigma02_int = NULL, sigma2_int = NULL,
    L = p_c + K * p_s + K,
    prior_vec = c(1 / (2 * p^1.1), 1 / (2 * p^1.1), 1 / p^1.25)
  )
  index_res_c <- list()
  sens_res_c <- rep(NA, K)
  prec_res_c <- rep(NA, K)
  for (k in seq_len(K)) {
    index_res_c[[k]] <- which(res_c$res[[k]]$alpha > 0.5)
    if (length(index_res_c[[k]]) == 0) {
      sens_res_c[k] <- prec_res_c[k] <- 0
    } else {
      sens_res_c[k] <- length(intersect(
        index_res_c[[k]],
        c(index_s[[k]], index_c)
      )) / (p_c + p_s)
      prec_res_c[k] <- length(intersect(
        index_res_c[[k]],
        c(index_s[[k]], index_c)
      )) / length(index_res_c[[k]])
    }
  }
  #### Single data set
  index_res_s <- list()
  sens_res_s <- rep(NA, K)
  prec_res_s <- rep(NA, K)
  ## data set 1
  for (k in seq_len(K)) {
    res <- susieR::susie(X = dta_list[[k]]$X, y = dta_list[[k]]$Y, L = p_c + p_s + 1)
    index_res_s[[k]] <- as.numeric(res$sets$cs)
    if (length(index_res_s[[k]]) == 0) {
      sens_res_s[k] <- prec_res_s[k] <- 0
    } else {
      sens_res_s[k] <- length(intersect(
        index_res_s[[k]],
        c(index_s[[k]], index_c)
      )) / (p_c + p_s)
      prec_res_s[k] <- length(intersect(
        index_res_s[[k]],
        c(index_s[[k]], index_c)
      )) / length(index_res_s[[k]])
    }
  }
  return(list(
    sens_res_c = mean(sens_res_c), sens_res_s = mean(sens_res_s),
    prec_res_c = mean(prec_res_c), prec_res_s = mean(prec_res_s)
  ))
}

## define parameters
K <- 2
n <- 500
p <- 1000
sigma <- 1
sigma0 <- 0.6
p_c_vec <- c(10, 25)
p_s_vec <- c(2, 5)
iter_sim_max <- 500

library(doParallel)
library(foreach)
library(doRNG)

for (iter_pc in seq_len(length(p_c_vec))) {
  p_c <- p_c_vec[iter_pc]
  for (iter_ps in seq_len(length(p_s_vec))) {
    p_s <- p_s_vec[iter_ps]
    cl <- makeCluster(25)
    registerDoParallel(cl)
    set.seed(2022)
    out_res <- foreach(iter = seq_len(iter_sim_max)) %dorng% {
      simu_vs_fun(
        K = K, n = n, p = p, p_c = p_c, p_s = p_s,
        sigma = sigma, sigma0 = sigma0
      )
    }
    stopCluster(cl)
    print(paste0(
      "Finish K", K, "n", n, "p", p, "sigma", sigma,
      "pc", p_c, "ps", p_s
    ))
    saveRDS(out_res, paste0(
      "Finish K", K, "n", n, "p", p, "sigma", sigma, "pc", p_c,
      "ps", p_s, ".rds"
    ))
  }
}