library(pcalg)
library(stabs)
source("simulation_DAG/graph_generation.R")
# args <- commandArgs()
# p <- as.numeric(args[6])
# n_tol <- as.numeric(args[7])
p <- 100
n_tol <- 600
K <- 2
n <- n_tol / K
n_graph <- 1

## define metric function
## True positive rate
TPrate_fun <- function(adj_pre, adj_act) {
  P <- which(adj_act == 1)
  PP <- which(adj_pre == 1)
  return(length(intersect(P, PP)) / length(P))
}
## False positive rate
FPrate_fun <- function(adj_pre, adj_act) {
  N <- which(adj_act == 0)
  PP <- which(adj_pre == 1)
  return(length(intersect(N, PP)) / length(N))
}
## remove order edge
check_edge <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  return(sum(abs(adj_pre - adj_act)) / 2)
}

#### generate graph
set.seed(2021)
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
    pc_fit <- pc(suffStat = list(C = dta_cor, n = dim(dta)[1]), 
                 indepTest = gaussCItest, alpha = alpha, 
                 labels = sapply(1:p, toString))
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
      "alpha = ", alphas[iter], c(shd(g_true, g), check_edge(adj_true, adj),
                                    round(TPrate_fun(adj_pre = adj, adj_act = adj_true), 4),
                                    round(FPrate_fun(adj_pre = adj, adj_act = adj_true), 4)
      ), "\n"
    )
  }
}

## data set 1
eval_fun(dag_list1, g_true = g_true1, adj_true = adj_true1)

# alpha =  1e-04 69 43 0.5538 0.0026 
# alpha =  5e-04 64 36 0.5769 0.0028 
# alpha =  0.001 60 33 0.6 0.0027 
# alpha =  0.005 59 32 0.6538 0.0029 
# alpha =  0.01 55 30 0.7154 0.0029 
# alpha =  0.05 60 40 0.7615 0.004 

## data set 2
eval_fun(dag_list2, g_true = g_true2, adj_true = adj_true2)

# alpha =  1e-04 69 40 0.6231 0.0029 
# alpha =  5e-04 61 35 0.6615 0.0026 
# alpha =  0.001 58 34 0.7 0.0025 
# alpha =  0.005 49 26 0.7615 0.0024 
# alpha =  0.01 50 26 0.7692 0.0026 
# alpha =  0.05 53 32 0.8 0.0034 