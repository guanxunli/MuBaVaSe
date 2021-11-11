## load data
setwd("~/Programming/R/MuBaVaSe/real_data")
dta <- readRDS("dta_use.rds")
dta_group <- read.csv("raw_data/dta.csv")
dta_group <- dta_group[, -1]
rownames(dta_group) <- dta_group$AOCSID
sample_group <- dta_group[rownames(dta), 8]
sample_group <- ifelse(sample_group == "1", 1, 0)
dta_1 <- dta[which(sample_group == 1), ]
dta_2 <- dta[which(sample_group == 0), ]

## generate graph
source("Graph_MCMC_two.R")
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
out_res <- Graph_MCMC_two(dta_1, dta_2, scale_x = FALSE, intercept = TRUE,
                          order_int = order_int, iter_max = 10, sigma02_int = NULL, sigma2_int = NULL,
                          prior_vec = NULL, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 1
)