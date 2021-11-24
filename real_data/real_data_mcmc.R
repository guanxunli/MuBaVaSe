## load data
# dta <- readRDS("dta_use.rds")
# dta_group <- read.csv("raw_data/dta.csv")
# dta_group <- dta_group[, -1]
# rownames(dta_group) <- dta_group$AOCSID
# sample_group <- dta_group[rownames(dta), 8]
# sample_group <- ifelse(sample_group == "1", 1, 0)
# dta_1 <- dta[which(sample_group == 1), ]
# dta_2 <- dta[which(sample_group == 0), ]
load("real_data/ovarian.rda")
dta_1 <- data[[1]]
dta_2 <- data[[2]]

## generate graph
source("real_data/method_code/Graph_MCMC_two.R")
## get order
dta <- rbind(dta_1, dta_2)
library(pcalg)
score_ges <- new("GaussL0penObsScore", data = dta, intercept = TRUE)
ges_fit <- ges(score_ges)
ges_adj <- as(ges_fit$repr, "matrix")
ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
order_int <- as.numeric(igraph::topo_sort(graph_i))
## Do MCMC
out_res <- Graph_MCMC_two(dta_1, dta_2,
                          scale_x = FALSE, intercept = TRUE,
                          order_int = order_int, iter_max = 50000, sigma02_int = NULL, sigma2_int = NULL,
                          prior_vec = NULL, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 40000
)
saveRDS(out_res, "out_res.rds")

## check results
out_res <- readRDS("out_res.rds")
p <- ncol(dta)
alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
A_mat_1 <- matrix(0, nrow = p, ncol = p)
A_mat_2 <- matrix(0, nrow = p, ncol = p)
iter_max <- length(out_res[[1]])
for (iter in seq_len(iter_max)) {
  order_tmp <- order(out_res$order_list[[iter]])
  alpha_mat_1 <- alpha_mat_1 + out_res$alpha_list_1[[iter]][order_tmp, order_tmp]
  alpha_mat_2 <- alpha_mat_2 + out_res$alpha_list_2[[iter]][order_tmp, order_tmp]
  A_mat_1 <- A_mat_1 + out_res$A_list_1[[iter]][order_tmp, order_tmp]
  A_mat_2 <- A_mat_2 + out_res$A_list_2[[iter]][order_tmp, order_tmp]
}
alpha_mat_1 <- alpha_mat_1 / iter_max
alpha_mat_2 <- alpha_mat_2 / iter_max
A_mat_1 <- A_mat_1 / iter_max
A_mat_2 <- A_mat_2 / iter_max

#### Compare results
## data set 1
adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
adj_1 <- t(adj_1)
sum(adj_1) # 22
## data set 2
adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
adj_2 <- t(adj_2)
sum(adj_2) # 22
## check results
length(intersect(which(adj_1 == 1), which(adj_2 == 1))) # 22

#### check centers
## get center
mcmc_adj <- ceiling((adj_1 + adj_2) / 2)
c_row <- rowSums(mcmc_adj)
c_col <- colSums(mcmc_adj)
c_sum <- c_row + c_col
names(c_sum) <- colnames(dta_1)
head(sort(c_sum, decreasing = TRUE))