library(pcalg)
library(parallel)
source("simulation_DAG/graph_generation.R")
# args <- commandArgs()
# p <- as.numeric(args[6])
# n_tol <- as.numeric(args[7])
p <- 100
n_tol <- 600
K <- 2
n <- n_tol / K

## define metric function
# ## True positive rate
# TPrate_fun <- function(adj_pre, adj_act) {
#   P <- which(adj_act == 1)
#   PP <- which(adj_pre == 1)
#   return(length(intersect(P, PP)) / length(P))
# }
# ## False positive rate
# FPrate_fun <- function(adj_pre, adj_act) {
#   N <- which(adj_act == 0)
#   PP <- which(adj_pre == 1)
#   return(length(intersect(N, PP)) / length(N))
# }
## remove order edge
check_edge <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  return(sum(abs(adj_pre - adj_act)) / 2)
}

########################### Do one figure ##################################
#### generate graph
set.seed(2021)
n_graph <- 1
graph_sim <- graph_generation(K = K, n_graph = n_graph, p = p, n_tol = n_tol)
adj_true1 <- t(graph_sim$G[[1]][[1]])
g_true1 <- as(getGraph(adj_true1), "graphNEL")
weight_true1 <- t(graph_sim$A[[1]][[1]])
adj_true2 <- t(graph_sim$G[[1]][[2]])
g_true2 <- as(getGraph(adj_true2), "graphNEL")
weight_true2 <- t(graph_sim$A[[1]][[2]])
data <- graph_sim$X[[1]]

## PC method
pc_fun <- function(dta, alphas = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)) {
  p <- ncol(dta)
  dta_cor <- cor(dta)
  dag_list <- list()
  for (iter_alpha in seq_len(length(alphas))) {
    alpha <- alphas[iter_alpha]
    pc_fit <- pc(
      suffStat = list(C = dta_cor, n = dim(dta)[1]),
      indepTest = gaussCItest, alpha = alpha,
      labels = sapply(1:p, toString)
    )
    dag <- as(pc_fit@graph, "matrix")
    dag_list[[iter_alpha]] <- ifelse(dag == TRUE, 1, 0)
  }
  return(dag_list)
}

dag_list1 <- pc_fun(data[[1]])
dag_list2 <- pc_fun(data[[2]])

#### check results
eval_fun <- function(dag_list, g_true, adj_true, alphas = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)) {
  for (iter in seq_len(length(dag_list))) {
    adj <- dag_list[[iter]]
    g <- as(adj, "graphNEL")
    cat(
      "alpha = ", alphas[iter], c(shd(g_true, g), check_edge(adj_true, adj)), "\n"
    )
  }
}

## data set 1
eval_fun(dag_list1, g_true = g_true1, adj_true = adj_true1)

# alpha =  1e-04 69 43
# alpha =  5e-04 64 36
# alpha =  0.001 60 33
# alpha =  0.005 59 32
# alpha =  0.01 55 30
# alpha =  0.05 60 40

## data set 2
eval_fun(dag_list2, g_true = g_true2, adj_true = adj_true2)

# alpha =  1e-04 69 40
# alpha =  5e-04 61 35
# alpha =  0.001 58 34
# alpha =  0.005 49 26
# alpha =  0.01 50 26
# alpha =  0.05 53 32

########################### Do parallel ##################################
#### generate graph
set.seed(2021)
n_graph <- 20
graph_sim <- graph_generation(K = K, n_graph = n_graph, p = p, n_tol = n_tol)
alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)

pc_fun <- function(dta, alphas = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)) {
  p <- ncol(dta)
  dta_cor <- cor(dta)
  dag_list <- list()
  for (iter_alpha in seq_len(length(alphas))) {
    alpha <- alphas[iter_alpha]
    pc_fit <- pc(
      suffStat = list(C = dta_cor, n = dim(dta)[1]),
      indepTest = gaussCItest, alpha = alpha,
      labels = sapply(1:p, toString)
    )
    dag <- as(pc_fit@graph, "matrix")
    dag_list[[iter_alpha]] <- ifelse(dag == TRUE, 1, 0)
  }
  return(dag_list)
}

library(foreach)
library(doParallel)
library(doRNG)
cl <- makeCluster(20)
registerDoParallel(cl)
out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
  library(pcalg)
  ## load data
  data <- graph_sim$X[[iter]]
  ## Do GES
  dag_list1 <- pc_fun(data[[1]], alphas = alphas)
  dag_list2 <- pc_fun(data[[2]], alphas = alphas)
  list(dag_list1 = dag_list1, dag_list2 = dag_list2)
}
stopCluster(cl)

## check results
res_1 <- list()
res_2 <- list()
for (iter_alpha in seq_len(length(alphas))) {
  res_1[[iter_alpha]] <- matrix(NA, nrow = n_graph, ncol = 2)
  res_2[[iter_alpha]] <- matrix(NA, nrow = n_graph, ncol = 2)
  for (iter_graph in seq_len(n_graph)) {
    ## load true value
    adj_true1 <- t(graph_sim$G[[iter_graph]][[1]])
    g_true1 <- as(getGraph(adj_true1), "graphNEL")
    adj_true2 <- t(graph_sim$G[[iter_graph]][[2]])
    g_true2 <- as(getGraph(adj_true2), "graphNEL")
    ## load results
    adj1 <- out_res[[iter_graph]][[1]][[iter_alpha]]
    g1 <- as(adj1, "graphNEL")
    adj2 <- out_res[[iter_graph]][[2]][[iter_alpha]]
    g2 <- as(adj2, "graphNEL")
    ## save results
    res_1[[iter_alpha]][iter_graph, ] <- c(
      shd(g_true1, g1),
      check_edge(adj_true1, adj1)
    )
    res_2[[iter_alpha]][iter_graph, ] <- c(
      shd(g_true2, g2),
      check_edge(adj_true2, adj2)
    )
  }
  cat(
    "alpha:", alphas[iter_alpha],
    "data1:", round(colMeans(res_1[[iter_alpha]]), 4),
    "data2:", round(colMeans(res_2[[iter_alpha]]), 4), "\n"
  )
}

# alpha: 1e-04 data1: 68.4 43.85 data2: 65.95 41.3
# alpha: 5e-04 data1: 62.1 38.5 data2: 59.4 36.4
# alpha: 0.001 data1: 59.65 36.35 data2: 56.7 34.15
# alpha: 0.005 data1: 53.9 32.05 data2: 51.45 29.55
# alpha: 0.01 data1: 52.8 32 data2: 50.4 29.45
# alpha: 0.05 data1: 61.4 45.8 data2: 56.75 42