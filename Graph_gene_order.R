## define parameters
p <- 100
n <- 75
p_c <- 100
p_1 <- 20
p_2 <- 20
sigma <- 1
sigma0 <- 0.6
A1 <- matrix(0, nrow = p, ncol = p)
A2 <- matrix(0, nrow = p, ncol = p)
set.seed(202108)
# Define the true graph given order
index_c <- sample(seq_len(p * (p - 1) / 2), size = p_c, replace = FALSE)
index_1 <- sample(setdiff(seq_len(p * (p - 1) / 2), index_c), size = p_1, replace = FALSE)
index_2 <- sample(setdiff(seq_len(p * (p - 1) / 2), c(index_1, index_c)), size = p_2, replace = FALSE)

A1[lower.tri(A1)][c(index_c, index_1)] <-  rnorm(p_c + p_1, mean = 0, sd = sigma0)
A2[lower.tri(A2)][c(index_c, index_2)] <-  rnorm(p_c + p_2, mean = 0, sd = sigma0)

alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
alpha_mat_1[lower.tri(alpha_mat_1)][c(index_c, index_1)] <- 1
alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
alpha_mat_2[lower.tri(alpha_mat_2)][c(index_c, index_2)] <- 1

eps_1 <- matrix(rnorm(p * n), nrow = p, ncol = n)
dta_1 <- solve(diag(1, nrow = p) - A1, eps_1)
eps_2 <- matrix(rnorm(p * n), nrow = p, ncol = n)
dta_2 <- solve(diag(1, nrow = p) - A2, eps_2)

## joint inference
joint_graph_fun <- function(dta_1, dta_2, sigma02_int = NULL, sigma2_int = NULL, r = 1, 
                            q = 1, tau = 1.5, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8) {
  source("Multi_dataset_null.R")
  p <- nrow(dta_1)
  alpha_res_1 <- matrix(0, nrow = p, ncol = p)
  alpha_res_2 <- matrix(0, nrow = p, ncol = p)
  A_res_1 <- matrix(0, nrow = p, ncol = p)
  A_res_2 <- matrix(0, nrow = p, ncol = p)
  for (iter_p in seq_len(p - 1)) {
    X_1 <- t(dta_1[seq_len(iter_p), , drop = FALSE])
    Y_1 <- dta_1[iter_p + 1, ]
    X_2 <- t(dta_2[seq_len(iter_p), , drop = FALSE])
    Y_2 <- dta_2[iter_p + 1, ]
    
    res <- sum_single_effect_multi_null(X_1 = X_1, Y_1 = Y_1, X_2 = X_2, Y_2 = Y_2, sigma02_int = sigma02_int,
                                        sigma2_int = sigma2_int, r = r, q = q, tau = tau, L = iter_p, 
                                        itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd)
    alpha_res_1[iter_p + 1, seq_len(iter_p)] <- res$alpha_1
    alpha_res_2[iter_p + 1, seq_len(iter_p)] <- res$alpha_2
    A_res_1[iter_p + 1, seq_len(iter_p)] <- res$post_mean1
    A_res_2[iter_p + 1, seq_len(iter_p)] <- res$post_mean2
  }
  ## return results
  return(list(alpha_res_1 = alpha_res_1, alpha_res_2 = alpha_res_2,
              A_res_1 = A_res_1, A_res_2 = A_res_2))
}


## single inference
single_graph_fun <- function(dta) {
  p <- nrow(dta)
  alpha_res <- matrix(0, nrow = p, ncol = p)
  A_res<- matrix(0, nrow = p, ncol = p)
  for (iter_p in seq_len(p - 1)) {
    X <- t(dta[seq_len(iter_p), , drop = FALSE])
    Y <- dta[iter_p + 1, ]
    res <- susieR::susie(X = X_1, y = Y_1, L = iter_p)
    alpha_res[iter_p + 1, seq_len(iter_p)] <- 1 - apply(1 - res$alpha[,,drop = FALSE], 1, prod)
    A_res[iter_p + 1, seq_len(iter_p)] <- rowSums(res$mu * res$alpha)
  }
  ## return results
  return(list(alpha_res = alpha_res, A_res = A_res))
}

res_joint <- joint_graph_fun(dta_1 = dta_1, dta_2 = dta_2, r = 1, q = 1, tau = 1.5)
res_single1 <- single_graph_fun(dta_1)
res_single2 <- single_graph_fun(dta_2)

#### check error
## data set 1
# l1 error
sum(abs(res_single1$alpha_res - alpha_mat_1))
sum(abs(res_joint$alpha_res_1 - alpha_mat_1))
# l2 error
sum((res_single1$alpha_res - alpha_mat_1)^2)
sum((res_joint$alpha_res_1 - alpha_mat_1)^2)
## data set 2
# l1 error
sum(abs(res_single2$alpha_res - alpha_mat_2))
sum(abs(res_joint$alpha_res_2 - alpha_mat_2))
# l2 error
sum((res_single2$alpha_res - alpha_mat_2)^2)
sum((res_joint$alpha_res_2 - alpha_mat_2)^2)

#### check posterior
## data set 1
# l1 error
sum(abs(res_single1$A_res - A1))
sum(abs(res_joint$A_res_1 - A1))
# l2 error
sum((res_single1$A_res - A1)^2)
sum((res_joint$A_res_1 - A1)^2)
## data set 2
# l1 error
sum(abs(res_single2$A_res - A2))
sum(abs(res_joint$A_res_2 - A2))
# l2 error
sum((res_single2$A_res - A2)^2)
sum((res_joint$A_res_2 - A2)^2)

