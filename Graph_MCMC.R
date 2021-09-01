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
set.seed(2021)
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

order_tmp <- sample(1:p, p)
dta_1 <- dta_1[order_tmp, ]
dta_2 <- dta_2[order_tmp, ]

## filp the order
order_change <- function(order_old) {
  order_new <- order_old
  index_change <- sample(1:p, 2)
  order_new[index_change] <- order_old[c(index_change[2], index_change[1])]
  return(order_new)
}

## MCMC method for Graph
# dta_1 and dta_2 are p x n data set
# order_int is the initialized order for nodes
# iter_max is the maximun mcmc step
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# r is for common part and q is for single part
# tau is the prior power for null model 1 / (p^tau)
# itermax is the maximum iteration
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l
Graph_MCMC <- function(dta_1, dta_2, order_int = NULL, iter_max = 10000, sigma02_int = NULL, sigma2_int = NULL, r = 1, 
                       q = 1, tau = 1.5, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8) {
  ## load the main function
  source("Graph_given_order.R")
  ## Initialization
  p <- nrow(dta_1)
  if (is.null(order_int)) {
    order_old <- sample(1:p, size = p)
  } else {
    order_old <- order_int
  }
  res_old <- joint_graph_fun(dta_1 = dta_1, dta_2 = dta_2, sigma02_int = sigma02_int, sigma2_int = sigma2_int, r = r, 
                             q = q, tau = tau, itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
                             order_dta = order_int)
  llike_old <- res_old$llike
  ## save lists
  alpha_list_1 <- list()
  alpha_list_2 <- list()
  A_list_1 <- list()
  A_list_2 <- list()
  order_list <- list()
  ## begin iteration
  for (iter_MCMC in seq_len(iter_max)) {
    # propose the new order
    order_pro <- order_change(order_old)
    # calculate the new graph
    res_pro <- joint_graph_fun(dta_1 = dta_1, dta_2 = dta_2, sigma02_int = sigma02_int, sigma2_int = sigma2_int, r = r, 
                               q = q, tau = tau, itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
                               order_dta = order_pro)
    llike_pro <- res_pro$llike
    # accept or reject
    if (llike_pro > llike_old) {
      # accept
      llike_old <- llike_pro
      res_old <- res_pro
      order_old <- order_pro
    } else {
      p_ratio <- exp(llike_pro - llike_old) # accept prob
      U <- runif(1)
      if (U < p_ratio) {
        # accept
        llike_old <- llike_pro
        res_old <- res_pro
        order_old <- order_pro
      }
    }
    # save results
    alpha_list_1[[iter_MCMC]] <- res_old$alpha_res_1
    alpha_list_2[[iter_MCMC]] <- res_old$alpha_res_2
    A_list_1[[iter_MCMC]] <- res_old$A_res_1
    A_list_2[[iter_MCMC]] <- res_old$A_res_2
    order_list[[iter_MCMC]] <- order_old 
  }
  # return results
  return(list(alpha_list_1 = alpha_list_1, alpha_list_2 = alpha_list_2, A_list_1 = A_list_1, A_list_2 = A_list_2,
              order_list = order_list))
}

# ## MCMC
# time1 <- Sys.time()
# res <- Graph_MCMC(dta_1 = dta_1, dta_2 = dta_2, iter_max = 10) 
# Sys.time() - time1 (13.72 min)