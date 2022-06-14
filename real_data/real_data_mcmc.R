## load data
library(ggplot2)
load("real_data/ovarian.rda")
dta_1 <- data[[1]]
dta_2 <- data[[2]]
p <- ncol(data[[1]])

## generate graph
source("two_data_sets/Graph_MCMC_two.R")
prior_vec_list <- list()
prior_vec_list[[1]] <- c(1 / (2 * p^1.25), 1 / p^1.5)
prior_vec_list[[2]] <- c(1 / (2 * p^1.5), 1 / p^2)
prior_vec_list[[3]] <- c(1 / p^2, 1 / p^2.25)
prior_vec_list[[4]] <- c(1 / (2 * p^2), 1 / p^2.25)

scale_x <- FALSE
intercept <- TRUE
iter_max <- 1e5

# ## A quick test
# set.seed(2021)
# prior_vec <- prior_vec_list[[1]]
# ## get order
# dta <- rbind(dta_1, dta_2)
# set.seed(2021)
# library(pcalg)
# score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
# ges_fit <- ges(score_ges)
# ges_adj <- as(ges_fit$repr, "matrix")
# ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
# graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
# order_int <- as.numeric(igraph::topo_sort(graph_i))
# out_res <- Graph_MCMC_two(dta_1, dta_2,
#                           scale_x = scale_x, intercept = intercept,
#                           order_int = order_int, iter_max = 200, sigma02_int = NULL, sigma2_int = NULL,
#                           prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
#                           burn_in = 1
# )
#
# ## plot results
# library(ggplot2)
# library(gridExtra)
# ggplot() + geom_line(aes(x = seq_len(length(out_res$llike_vec)), y = out_res$llike_vec))
# alpha_mat_1 <- out_res$alpha_mat_1
# alpha_mat_2 <- out_res$alpha_mat_2
# A_mat_1 <- out_res$A_mat_1
# A_mat_2 <- out_res$A_mat_2
# #### Compare results
# ## data set 1
# adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
# adj_1 <- t(adj_1)
# ## data set 2
# adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
# adj_2 <- t(adj_2)
# # intersection
# inter_order <- length(intersect(which(adj_1 == 1), which(adj_2 == 1)))
# ## check results
# cat(
#   "scale_x:", scale_x, "intercept", intercept,
#   sum(adj_1), sum(adj_2), inter_order, "\n"
# )

#### Do MCMC with order
## get order
dta <- rbind(dta_1, dta_2)
set.seed(2021)
library(pcalg)
# order_int <- NULL
score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
ges_fit <- ges(score_ges)
ges_adj <- as(ges_fit$repr, "matrix")
ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
order_int <- as.numeric(igraph::topo_sort(graph_i))
## Do MCMC

library(foreach)
library(doParallel)
library(doRNG)

cl <- makeCluster(4)
registerDoParallel(cl)
out_res <- foreach(iter_prior = seq_len(length(prior_vec_list))) %dorng% {
  prior_vec <- prior_vec_list[[iter_prior]]
  Graph_MCMC_two(dta_1, dta_2,
    scale_x = scale_x, intercept = intercept,
    order_int = order_int, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
    prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
    burn_in = iter_max - 5000
  )
}
stopCluster(cl)
saveRDS(out_res, "real_data/results/out_mcmc.rds")
## save results
out_res <- readRDS("real_data/results/out_mcmc.rds")
for (iter_prior in seq_len(length(prior_vec_list))) {
  res_tmp <- out_res[[iter_prior]]
  # png(paste0("prior", iter_prior, "realdata_mcmctwo.png"))
  # ggplot() +
  #   geom_line(aes(x = seq_len(iter_max), y = res_tmp$llike_vec)) +
  #   xlab("Iteration") +
  #   ylab("Log likelihood")
  # dev.off()
  alpha_mat_1 <- res_tmp$alpha_mat_1
  alpha_mat_2 <- res_tmp$alpha_mat_2
  A_mat_1 <- res_tmp$A_mat_1
  A_mat_2 <- res_tmp$A_mat_1
  #### Compare results
  ## data set 1
  adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
  adj_1 <- t(adj_1)
  ## data set 2
  adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
  adj_2 <- t(adj_2)
  # intersection
  inter_order <- length(intersect(which(adj_1 == 1), which(adj_2 == 1)))
  ## check results
  cat(
    "scale_x:", scale_x, "intercept", intercept, "prior", iter_prior,
    sum(adj_1), sum(adj_2), inter_order, "\n"
  )
}