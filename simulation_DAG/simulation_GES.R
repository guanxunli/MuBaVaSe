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

#### GES method
ges_fun <- function(dta, lambdas = c(1, 2, 3, 4, 5)) {
  p <- ncol(dta)
  dag_list <- list()
  for (iter_lambda in seq_len(length(lambdas))) {
    lambda <- lambdas[iter_lambda]
    l0score <- new("GaussL0penObsScore", data = dta, lambda = lambda * log(p), intercept = FALSE)
    ges_fit <- ges(l0score)
    dag <- as(ges_fit$repr, "matrix")
    dag_list[[iter_lambda]] <- ifelse(dag == TRUE, 1, 0)
  }
  return(dag_list)
}

dag_list1 <- ges_fun(data[[1]])
dag_list2 <- ges_fun(data[[2]])

#### check results
eval_fun <- function(dag_list, g_true, adj_true, lambdas = c(1, 2, 3, 4, 5)) {
  for (iter in seq_len(length(dag_list))) {
    adj <- dag_list[[iter]]
    g <- as(adj, "graphNEL")
    cat(
      "lambda = ", lambdas[iter], c(shd(g_true, g), check_edge(adj_true, adj)), "\n"
    )
  }
}

## data set 1
eval_fun(dag_list1, g_true = g_true1, adj_true = adj_true1)

# lambda =  1 48 39
# lambda =  2 21 19
# lambda =  3 26 25
# lambda =  4 31 30
# lambda =  5 40 39

## data set 2
eval_fun(dag_list2, g_true = g_true2, adj_true = adj_true2)

# lambda =  1 37 29
# lambda =  2 20 17
# lambda =  3 21 19
# lambda =  4 24 22
# lambda =  5 28 25

########################### Do parallel ##################################
#### generate graph
set.seed(2021)
n_graph <- 20
graph_sim <- graph_generation(K = K, n_graph = n_graph, p = p, n_tol = n_tol)
lambdas <- c(1, 2, 3, 4, 5)

ges_fun <- function(dta, lambdas = c(1, 2, 3, 4, 5)) {
  p <- ncol(dta)
  dag_list <- list()
  for (iter_lambda in seq_len(length(lambdas))) {
    lambda <- lambdas[iter_lambda]
    l0score <- new("GaussL0penObsScore", data = dta, lambda = lambda * log(p), intercept = FALSE)
    ges_fit <- ges(l0score)
    dag <- as(ges_fit$repr, "matrix")
    dag_list[[iter_lambda]] <- ifelse(dag == TRUE, 1, 0)
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
  dag_list1 <- ges_fun(data[[1]], lambdas = lambdas)
  dag_list2 <- ges_fun(data[[2]], lambdas = lambdas)
  list(dag_list1 = dag_list1, dag_list2 = dag_list2)
}
stopCluster(cl)

## check results
res_1 <- list()
res_2 <- list()
for (iter_lambda in seq_len(length(lambdas))) {
  res_1[[iter_lambda]] <- matrix(NA, nrow = n_graph, ncol = 2)
  res_2[[iter_lambda]] <- matrix(NA, nrow = n_graph, ncol = 2)
  for (iter_graph in seq_len(n_graph)) {
    ## load true value
    adj_true1 <- t(graph_sim$G[[iter_graph]][[1]])
    g_true1 <- as(getGraph(adj_true1), "graphNEL")
    adj_true2 <- t(graph_sim$G[[iter_graph]][[2]])
    g_true2 <- as(getGraph(adj_true2), "graphNEL")
    ## load results
    adj1 <- out_res[[iter_graph]][[1]][[iter_lambda]]
    g1 <- as(adj1, "graphNEL")
    adj2 <- out_res[[iter_graph]][[2]][[iter_lambda]]
    g2 <- as(adj2, "graphNEL")
    ## save results
    res_1[[iter_lambda]][iter_graph, ] <- c(
      shd(g_true1, g1),
      check_edge(adj_true1, adj1)
    )
    res_2[[iter_lambda]][iter_graph, ] <- c(
      shd(g_true2, g2),
      check_edge(adj_true2, adj2)
    )
  }
  cat(
    "lambda:", lambdas[iter_lambda],
    "data1:", round(colMeans(res_1[[iter_lambda]]), 4),
    "data2:", round(colMeans(res_2[[iter_lambda]]), 4), "\n"
  )
}

# lambda: 1 data1: 44.9 34.8 data2: 40.05 31
# lambda: 2 data1: 26.5 22.3 data2: 22.95 18.95
# lambda: 3 data1: 32.35 28.85 data2: 28.8 25.35
# lambda: 4 data1: 38.75 36.05 data2: 34.25 31.35
# lambda: 5 data1: 44.55 42.35 data2: 40.2 37.85