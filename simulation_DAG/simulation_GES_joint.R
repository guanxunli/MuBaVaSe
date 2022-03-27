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
e_com <- 50
e_pri <- 50

## define metric function
## remove order edge
check_edge <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  return(sum(abs(adj_pre - adj_act)) / 2)
}
## True positive rate
TPrate_fun <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  P <- which(adj_act == 1)
  PP <- which(adj_pre == 1)
  return(length(intersect(P, PP)) / length(P))
}
## False positive rate
FPrate_fun <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  N <- which(adj_act == 0)
  PP <- which(adj_pre == 1)
  return(length(intersect(N, PP)) / length(N))
}
## False negative rate
FNrate_fun <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  P <- which(adj_act == 1)
  PN <- which(adj_pre == 0)
  return(length(intersect(PN, P)) / length(P))
}
## check adjacency matrix
check_adj_l2 <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  return(sum((adj_pre - adj_act)^2) / 2)
}
check_adj_l1 <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  return(sum(abs(adj_pre - adj_act)) / 2)
}

# ########################### Do one figure ##################################
# #### generate graph
# set.seed(2021)
# n_graph <- 1
# graph_sim <- graph_generation(
#   K = K, n_graph = n_graph, p = p, n_tol = n_tol,
#   e_com = e_com, e_pri = e_pri
# )
# adj_true <- list()
# g_true <- list()
# for (iter_K in seq_len(K)) {
#   adj_true[[iter_K]] <- t(graph_sim$G[[1]][[iter_K]])
#   g_true[[iter_K]] <- as(getGraph(adj_true[[iter_K]]), "graphNEL")
# }
# data <- graph_sim$X[[1]]
# 
# #### joint GES method the first step
# ges_joint_fun <- function(data, lambdas = c(0.5, 1, 2, 3, 4, 5)) {
#   source("simulation_DAG/newclass.R")
#   p <- ncol(data[[1]])
#   dag_list <- list()
#   for (iter_lambda in seq_len(length(lambdas))) {
#     lambda <- lambdas[iter_lambda]
#     l0score <- new("MultiGaussL0pen",
#                    data = data, lambda = lambda * log(p),
#                    intercept = TRUE, use.cpp = FALSE
#     )
#     ges_fit <- ges(l0score)
#     dag <- as(ges_fit$essgraph, "matrix")
#     dag_list[[iter_lambda]] <- ifelse(dag == TRUE, 1, 0)
#   }
#   return(dag_list)
# }
# 
# dag_list_com <- ges_joint_fun(data)
# 
# ## Joint GES the second step
# subset <- function(y, x, data) {
#   t <- rep(0, ncol(data))
#   if (length(x) <= 1) {
#     t[x] <- 1
#   } else {
#     model <- glmnet::cv.glmnet(as.matrix(data[, x]), data[, y], family = "gaussian", intercept = FALSE)
#     nonz <- which(as.vector(coef(model)) != 0) - 1
#     t[x[nonz]] <- 1
#   }
#   return(t)
# }
# 
# # do joint estimation given single data
# ges_alg <- function(dag_list, dta) {
#   adj_list <- list()
#   for (iter in seq_len(length(dag_list))) {
#     in_mat <- dag_list[[iter]]
#     joint_mat <- sapply(seq_len(ncol(dta)), function(i) subset(i, which(in_mat[, i] != 0), dta))
#     adj_list[[iter]] <- joint_mat
#   }
#   return(adj_list)
# }
# 
# dag_list <- list()
# for (iter_K in seq_len(K)) {
#   dag_list[[iter_K]] <- ges_alg(dag_list_com, data[[iter_K]])
# }
# 
# #### check results
# eval_fun <- function(dag_list, g_true, adj_true,
#                      lambdas = c(0.5, 1, 2, 3, 4, 5)) {
#   for (iter in seq_len(length(dag_list))) {
#     adj <- dag_list[[iter]]
#     g <- as(adj, "graphNEL")
#     cat(
#       "lambda = ", lambdas[iter], c(shd(g_true, g), check_edge(adj_true, adj)),
#       "TP", round(TPrate_fun(adj_pre = adj, adj_act = adj_true), 4),
#       "FP", round(FPrate_fun(adj_pre = adj, adj_act = adj_true), 4),
#       "FN", round(FNrate_fun(adj_pre = adj, adj_act = adj_true), 4),
#       "L2", round(check_adj_l2(adj_pre = adj, adj_act = adj_true), 4),
#       "L1", round(check_adj_l1(adj_pre = adj, adj_act = adj_true), 4), "\n"
#     )
#   }
# }
# 
# ## data set 1
# for (iter_K in seq_len(K)) {
#   cat("data set", iter_K, "\n")
#   eval_fun(dag_list[[iter_K]], g_true = g_true[[iter_K]], adj_true = adj_true[[iter_K]])
# }

########################### Do parallel ##################################
#### generate graph
set.seed(2021)
n_graph <- 20
graph_sim <- graph_generation(
  K = K, n_graph = n_graph, p = p, n_tol = n_tol,
  e_com = e_com, e_pri = e_pri
)
lambdas <- c(0.5, 1, 2, 3, 4, 5)

#### joint GES method the first step
ges_joint_fun <- function(data, lambdas = c(0.5, 1, 2, 3, 4, 5)) {
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

library(foreach)
library(doParallel)
library(doRNG)
cl <- makeCluster(20)
registerDoParallel(cl)
out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
  library(pcalg)
  ## load data
  data <- graph_sim$X[[iter]]
  dag_list <- list()
  ## Do joint GES the first step
  dag_list_com <- ges_joint_fun(data)
  ## Do joint GES the second step
  for (iter_K in seq_len(K)) {
    dag_list[[iter_K]] <- ges_alg(dag_list_com,  data[[iter_K]])
  }
  return(dag_list)
}
stopCluster(cl)

## check results
res <- list()
for (iter_K in seq_len(K)) {
  res[[iter_K]] <- list()
}
for (iter_lambda in seq_len(length(lambdas))) {
  for (iter_K in seq_len(K)) {
    res[[iter_K]][[iter_lambda]] <- matrix(NA, nrow = n_graph, ncol = 7)
  }
  for (iter_graph in seq_len(n_graph)) {
    for (iter_K in seq_len(K)) {
      ## load true value
      adj_true <- t(graph_sim$G[[iter_graph]][[iter_K]])
      g_true <- as(getGraph(adj_true), "graphNEL")
      ## load results
      adj <- out_res[[iter_graph]][[iter_K]][[iter_lambda]]
      g <- as(adj, "graphNEL")
      res[[iter_K]][[iter_lambda]][iter_graph, ] <- c(
        shd(g_true, g),
        check_edge(adj_true, adj),
        TPrate_fun(adj_pre = adj, adj_act = adj_true),
        FPrate_fun(adj_pre = adj, adj_act = adj_true),
        FNrate_fun(adj_pre = adj, adj_act = adj_true),
        check_adj_l2(adj_pre = adj, adj_act = adj_true),
        check_adj_l1(adj_pre = adj, adj_act = adj_true)
      )
    }
  }
  cat("lambda:", lambdas[iter_lambda], "p:", p, "e_com:", e_com, "e_pri", e_pri, "\n")
  for (iter_K in seq_len(K)){
    cat("data", iter_K, round(colMeans(res[[iter_K]][[iter_lambda]]), 4), "\n")
  }
}

res_ave <- list()
for (iter_lambda in seq_len(length(lambdas))) {
  res_ave[[iter_lambda]] <- matrix(0, nrow = n_graph, ncol = 7)
  for (iter_K in seq_len(K)) {
    res_ave[[iter_lambda]] <- res_ave[[iter_lambda]] + res[[iter_K]][[iter_lambda]]
  }
  res_ave[[iter_lambda]] <- res_ave[[iter_lambda]] / K
  cat("lambda:", lambdas[iter_lambda], "p:", p, "e_com:", e_com, "e_pri", e_pri, "\n")
  cat(round(colMeans(res_ave[[iter_lambda]]), 4), "\n")
}