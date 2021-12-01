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

#### joint GES method the first step
ges_joint_fun <- function(data, lambdas = c(2, 3, 4, 5)) {
  source("simulation_DAG/newclass.R")
  p <- ncol(data[[1]])
  dag_list <- list()
  for (iter_lambda in seq_len(length(lambdas))) {
    lambda <- lambdas[iter_lambda]
    l0score <- new("MultiGaussL0pen",
                   data = data, lambda = lambda * log(p),
                   intercept = TRUE, use.cpp = FALSE
    )
    ges_fit <- ges(l0score)
    dag <- as(ges_fit$essgraph, "matrix")
    dag_list[[iter_lambda]] <- ifelse(dag == TRUE, 1, 0)
  }
  return(dag_list)
}

dag_list <- ges_joint_fun(data)

## Joint GES the second step
subset <- function(y, x, data) {
  t <- rep(0, ncol(data))
  if (length(x) <= 1) {
    t[x] <- 1
  } else {
    model <- glmnet::cv.glmnet(as.matrix(data[, x]), data[, y], family = "gaussian", intercept = FALSE)
    nonz <- which(as.vector(coef(model)) != 0) - 1
    t[x[nonz]] <- 1
  }
  return(t)
}

# do joint estimation given single data
ges_alg <- function(dag_list, dta) {
  adj_list <- list()
  for (iter in seq_len(length(dag_list))) {
    in_mat <- dag_list[[iter]]
    joint_mat <- sapply(seq_len(ncol(dta)), function(i) subset(i, which(in_mat[, i] != 0), dta))
    adj_list[[iter]] <- joint_mat
  }
  return(adj_list)
}

dag_list1 <- ges_alg(dag_list, data[[1]])
dag_list2 <- ges_alg(dag_list, data[[2]])

#### check results
eval_fun <- function(dag_list, g_true, adj_true, lambdas = c(2, 3, 4, 5)) {
  for (iter in seq_len(length(dag_list))) {
    adj <- dag_list[[iter]]
    g <- as(adj, "graphNEL")
    cat(
      "lambda = ", lambdas[iter], c(shd(g_true, g), check_edge(adj_true, adj),
                                    round(TPrate_fun(adj_pre = adj, adj_act = adj_true), 4),
                                    round(FPrate_fun(adj_pre = adj, adj_act = adj_true), 4)
      ), "\n"
    )
  }
}

## data set 1
eval_fun(dag_list1, g_true = g_true1, adj_true = adj_true1)

# lambda =  2 37 16 0.8769 0.0027 
# lambda =  3 35 15 0.9 0.0025 
# lambda =  4 38 16 0.8615 0.0025 
# lambda =  5 44 21 0.8077 0.0027 

## data set 2
eval_fun(dag_list2, g_true = g_true2, adj_true = adj_true2)

# lambda =  2 49 25 0.8615 0.0041 
# lambda =  3 32 10 0.9077 0.0025 
# lambda =  4 49 23 0.8462 0.0039 
# lambda =  5 52 23 0.7846 0.0039
