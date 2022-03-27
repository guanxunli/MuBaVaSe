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
# ## PC method
# pc_fun <- function(dta,
#                    alphas = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)) {
#   p <- ncol(dta)
#   dta_cor <- cor(dta)
#   dag_list <- list()
#   for (iter_alpha in seq_len(length(alphas))) {
#     alpha <- alphas[iter_alpha]
#     pc_fit <- pc(
#       suffStat = list(C = dta_cor, n = dim(dta)[1]),
#       indepTest = gaussCItest, alpha = alpha,
#       labels = sapply(1:p, toString)
#     )
#     dag <- as(pc_fit@graph, "matrix")
#     dag_list[[iter_alpha]] <- ifelse(dag == TRUE, 1, 0)
#   }
#   return(dag_list)
# }
# 
# dag_list <- list()
# for (iter_K in seq_len(K)) {
#   dag_list[[iter_K]] <- ges_fun(data[[iter_K]])
# }
# 
# #### check results
# eval_fun <- function(dag_list, g_true, adj_true,
#                      alphas = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)) {
#   for (iter in seq_len(length(dag_list))) {
#     adj <- dag_list[[iter]]
#     g <- as(adj, "graphNEL")
#     cat(
#       "alpha = ", alphas[iter], c(shd(g_true, g), check_edge(adj_true, adj)),
#       "TP", round(TPrate_fun(adj_pre = adj, adj_act = adj_true), 4),
#       "FP", round(FPrate_fun(adj_pre = adj, adj_act = adj_true), 4),
#       "FN", round(FNrate_fun(adj_pre = adj, adj_act = adj_true), 4),
#       "L2", round(check_adj_l2(adj_pre = adj, adj_act = adj_true), 4),
#       "L1", round(check_adj_l1(adj_pre = adj, adj_act = adj_true), 4), "\n"
#     )
#   }
# }
# 
# 
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
  dag_list <- list()
  ## Do PC
  for (iter_K in seq_len(K)) {
    dag_list[[iter_K]] <- pc_fun(data[[iter_K]],  alphas = alphas)
  }
  return(dag_list)
}
stopCluster(cl)

## check results
res <- list()
for (iter_K in seq_len(K)) {
  res[[iter_K]] <- list()
}
for (iter_alpha in seq_len(length(alphas))) {
  for (iter_K in seq_len(K)) {
    res[[iter_K]][[iter_alpha]] <- matrix(NA, nrow = n_graph, ncol = 7)
  }
  for (iter_graph in seq_len(n_graph)) {
    for (iter_K in seq_len(K)) {
      ## load true value
      adj_true <- t(graph_sim$G[[iter_graph]][[iter_K]])
      g_true <- as(getGraph(adj_true), "graphNEL")
      ## load results
      adj <- out_res[[iter_graph]][[iter_K]][[iter_alpha]]
      g <- as(adj, "graphNEL")
      res[[iter_K]][[iter_alpha]][iter_graph, ] <- c(
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
  cat("alpha:", alphas[iter_alpha], "p:", p, "e_com:", e_com, "e_pri", e_pri, "\n")
  for (iter_K in seq_len(K)){
    cat("data", iter_K, round(colMeans(res[[iter_K]][[iter_alpha]]), 4), "\n")
  }
}

res_ave <- list()
for (iter_alpha in seq_len(length(alphas))) {
  res_ave[[iter_alpha]] <- matrix(0, nrow = n_graph, ncol = 7)
  for (iter_K in seq_len(K)) {
    res_ave[[iter_alpha]] <- res_ave[[iter_alpha]] + res[[iter_K]][[iter_alpha]]
  }
  res_ave[[iter_alpha]] <- res_ave[[iter_alpha]] / K
  cat("alpha:", alphas[iter_alpha], "p:", p, "e_com:", e_com, "e_pri", e_pri, "\n")
  cat(round(colMeans(res_ave[[iter_alpha]]), 4), "\n")
}