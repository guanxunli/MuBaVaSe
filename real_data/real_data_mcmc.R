## load data
load("real_data/ovarian.rda")
dta_1 <- data[[1]]
dta_2 <- data[[2]]
p <- ncol(data[[1]])

## generate graph
source("real_data/Graph_MCMC_two.R")

## Do MCMC without order
out_res <- Graph_MCMC_two(dta_1, dta_2,
                          scale_x = FALSE, intercept = FALSE,
                          order_int = order_int, iter_max = 50000, sigma02_int = NULL, sigma2_int = NULL,
                          prior_vec = c(1 / (2 * p^1.5), 1 / p ^ 2), 
                          itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 1
)

library(ggplot2)
library(gridExtra)
p1 <- ggplot() +
  geom_line(aes(x = seq_len(length(out_res[[1]])), y = out_res$llike_vec))
p2 <- ggplot() +
  geom_line(aes(x = seq_len(5000), y = out_res$llike_vec[(length(out_res[[1]]) - 4999):length(out_res[[1]])]))
layout_matrix <- matrix(c(1, 2), nrow = 2)
pdf(file = "real_data/llikehoold_without.pdf")
grid.arrange(p1, p2)
dev.off()

## save results
iter_max <- length(out_res[[1]])
out_res$alpha_list_1 <- out_res$alpha_list_1[-seq_len(iter_max - 5000)]
out_res$alpha_list_2 <- out_res$alpha_list_2[-seq_len(iter_max - 5000)]
out_res$A_list_1 <- out_res$A_list_1[-seq_len(iter_max - 5000)]
out_res$A_list_2 <- out_res$A_list_2[-seq_len(iter_max - 5000)]
out_res$order_list <- out_res$order_list[-seq_len(iter_max - 5000)]
out_res$llike_vec <- out_res$llike_vec[-seq_len(iter_max - 5000)]
saveRDS(out_res, "real_data/results/out_mcmc_without.rds")
## check results
out_res <- readRDS("real_data/results/out_mcmc_without.rds")
p <- ncol(dta)
alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
A_mat_1 <- matrix(0, nrow = p, ncol = p)
A_mat_2 <- matrix(0, nrow = p, ncol = p)
iter_max <- length(out_res[[1]])
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
sum(adj_1) #
## data set 2
adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
adj_2 <- t(adj_2)
sum(adj_2) #
## check results
length(intersect(which(adj_1 == 1), which(adj_2 == 1)))

#### Do MCMC with order
## get order
dta <- rbind(dta_1, dta_2)
set.seed(2021)
library(pcalg)
score_ges <- new("GaussL0penObsScore", data = dta, intercept = TRUE)
ges_fit <- ges(score_ges)
ges_adj <- as(ges_fit$repr, "matrix")
ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
order_int <- as.numeric(igraph::topo_sort(graph_i))
## Do MCMC
out_res <- Graph_MCMC_two(dta_1, dta_2,
                          scale_x = FALSE, intercept = FALSE,
                          order_int = order_int, iter_max = 50000, sigma02_int = NULL, sigma2_int = NULL,
                          prior_vec = c(1 / (2 * p^1.5), 1 / p ^ 2), 
                          itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 1
)

library(ggplot2)
library(gridExtra)
p1 <- ggplot() +
  geom_line(aes(x = seq_len(length(out_res[[1]])), y = out_res$llike_vec))
p2 <- ggplot() +
  geom_line(aes(x = seq_len(5000), y = out_res$llike_vec[(length(out_res[[1]]) - 4999):length(out_res[[1]])]))
layout_matrix <- matrix(c(1, 2), nrow = 2)
pdf(file = "real_data/llikehoold.pdf")
grid.arrange(p1, p2)
dev.off()

## save results
iter_max <- length(out_res[[1]])
out_res$alpha_list_1 <- out_res$alpha_list_1[-seq_len(iter_max - 5000)]
out_res$alpha_list_2 <- out_res$alpha_list_2[-seq_len(iter_max - 5000)]
out_res$A_list_1 <- out_res$A_list_1[-seq_len(iter_max - 5000)]
out_res$A_list_2 <- out_res$A_list_2[-seq_len(iter_max - 5000)]
out_res$order_list <- out_res$order_list[-seq_len(iter_max - 5000)]
out_res$llike_vec <- out_res$llike_vec[-seq_len(iter_max - 5000)]
saveRDS(out_res, "real_data/results/out_mcmc.rds")
## check results
out_res <- readRDS("real_data/results/out_mcmc.rds")
p <- ncol(dta)
alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
A_mat_1 <- matrix(0, nrow = p, ncol = p)
A_mat_2 <- matrix(0, nrow = p, ncol = p)
iter_max <- length(out_res[[1]])
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
sum(adj_1) #
## data set 2
adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
adj_2 <- t(adj_2)
sum(adj_2) #
## check results
length(intersect(which(adj_1 == 1), which(adj_2 == 1)))

## 73 73 73