## load data
load("real_data/ovarian.rda")
dta_1 <- data[[1]]
dta_2 <- data[[2]]
p <- ncol(data[[1]])

## generate graph
source("Two_dataset_new/Graph_MCMC_two.R")
prior_vec_list <- list()
prior_vec_list[[1]] <- c(1 / (2 * p^1.5), 1 / p^2)
prior_vec_list[[2]] <- c(1 / p^1.5, 1 / p^1.5)
prior_vec_list[[3]] <- c(1 / p^1.5, 1 / p^2)
scale_x <- TRUE
intercept <- TRUE
iter_max <- 50000

#### Different priors test
for (iter_prior in seq_len(length(prior_vec_list))) {
  prior_vec <- prior_vec_list[[iter_prior]]
  ## get order
  dta <- rbind(dta_1, dta_2)
  set.seed(2021)
  library(pcalg)
  score_ges <- new("GaussL0penObsScore", data = dta, intercept = intercept)
  ges_fit <- ges(score_ges)
  ges_adj <- as(ges_fit$repr, "matrix")
  ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
  graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
  order_int <- as.numeric(igraph::topo_sort(graph_i))
  ## Do MCMC
  out_res <- Graph_MCMC_two(dta_1, dta_2,
                            scale_x = scale_x, intercept = intercept,
                            order_int = order_int, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
                            prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
                            burn_in = 1
  )
  ## save results
  out_res$alpha_list_1 <- out_res$alpha_list_1[-seq_len(iter_max - 5001)]
  out_res$alpha_list_2 <- out_res$alpha_list_2[-seq_len(iter_max - 5001)]
  out_res$A_list_1 <- out_res$A_list_1[-seq_len(iter_max - 5001)]
  out_res$A_list_2 <- out_res$A_list_2[-seq_len(iter_max - 5001)]
  out_res$order_list <- out_res$order_list[-seq_len(iter_max - 5001)]
  out_res$llike_vec <- out_res$llike_vec[-seq_len(iter_max - 5001)]
  ## check results
  alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
  alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
  A_mat_1 <- matrix(0, nrow = p, ncol = p)
  A_mat_2 <- matrix(0, nrow = p, ncol = p)
  for (iter in seq_len(5000)) {
    order_tmp <- order(out_res$order_list[[iter]])
    alpha_mat_1 <- alpha_mat_1 + out_res$alpha_list_1[[iter]][order_tmp, order_tmp]
    alpha_mat_2 <- alpha_mat_2 + out_res$alpha_list_2[[iter]][order_tmp, order_tmp]
    A_mat_1 <- A_mat_1 + out_res$A_list_1[[iter]][order_tmp, order_tmp]
    A_mat_2 <- A_mat_2 + out_res$A_list_2[[iter]][order_tmp, order_tmp]
  }
  alpha_mat_1 <- alpha_mat_1 / 5000
  alpha_mat_2 <- alpha_mat_2 / 5000
  A_mat_1 <- A_mat_1 / 5000
  A_mat_2 <- A_mat_2 / 5000
  
  #### Compare results
  ## data set 1
  adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
  adj_1 <- t(adj_1)
  ## data set 2
  adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
  adj_2 <- t(adj_2)
  ## check results
  cat(
    "scale_x:", scale_x, "intercept", intercept, "prior_vec", prior_vec, "\n",
    sum(adj_1), sum(adj_2), length(intersect(which(adj_1 == 1), which(adj_2 == 1))), "\n"
  )
}

# scale_x: TRUE intercept TRUE prior_vec 0.000754657 0.0001731302 40 97 39 
# scale_x: TRUE intercept TRUE prior_vec 0.001509314 0.001509314 63 111 62 
# scale_x: TRUE intercept TRUE prior_vec 0.001509314 0.0001731302 34 101 33

# scale_x: FALSE intercept TRUE prior_vec 0.000754657 0.0001731302 39 97 39 
# scale_x: FALSE intercept TRUE prior_vec 0.001509314 0.001509314 70 113 69 
# scale_x: FALSE intercept TRUE prior_vec 0.001509314 0.0001731302 35 106 34 

