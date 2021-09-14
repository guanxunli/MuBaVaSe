library(pcalg)
source("simulation_DAG/graph_generation.R")
#### generate graph
set.seed(2021)
graph_sim <- graph_generation(K = 2, n_graph = 1, p = 100, n_tol = 600)
#### GES method
## data set 1
time1 <- Sys.time()
g_true1 <- as(getGraph(t(graph_sim$G[[1]][[1]])), "graphNEL")
score1 <- new("GaussL0penObsScore", data = t(graph_sim$X[[1]][[1]]), intercept = FALSE)
ges_fit1 <- ges(score1)
ges_adj1 <- as(ges_fit1$repr, "matrix")
ges_graph1 <- as(ges_fit1$repr, "graphNEL")
ges_weight1 <- ges_fit1$repr$weight.mat()
## structural Hamming distance (SHD)
shd(g_true1, ges_graph1)
## Mean square error for weight
sum((t(graph_sim$A[[1]][[1]]) - ges_weight1)^2)
## data set 2
g_true2 <- as(getGraph(t(graph_sim$G[[1]][[2]])), "graphNEL")
score2 <- new("GaussL0penObsScore", data = graph_sim$X[[1]][[2]], intercept = FALSE)
ges_fit2 <- ges(score2)
ges_adj2 <- as(ges_fit2$repr, "matrix")
ges_graph2 <- as(ges_fit2$repr, "graphNEL")
ges_weight2 <- ges_fit2$repr$weight.mat()
## structural Hamming distance (SHD)
shd(g_true2, ges_graph2)
## Mean square error for weight
sum((t(graph_sim$A[[1]][[2]]) - ges_weight2)^2)
Sys.time() - time1

#### our method
source("Two_dataset/Graph_MCMC_two.R")
dta_1 <- graph_sim$X[[1]][[1]]
dta_2 <- graph_sim$X[[1]][[2]]
out_res <- Graph_MCMC_two(dta_1, dta_2, order_int = NULL, iter_max = 10000, sigma02_int = NULL, sigma2_int = NULL, r = 0.2, 
                          q = 0.05, tau = 1.5, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 5000)
