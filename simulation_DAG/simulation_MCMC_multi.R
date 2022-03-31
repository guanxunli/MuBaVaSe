library(pcalg)
source("simulation_DAG/graph_generation.R")
# args <- commandArgs()
# p <- as.numeric(args[6])
# n_tol <- as.numeric(args[7])
# Define parameters
p <- 100
n_tol <- 600
K <- 2
n <- n_tol / K
e_com <- 100
e_pri <- 20
# Define prior
prior_vec_list <- list()
prior_vec_list[[1]] <- c(1 / (2 * p^1.25), 1 / (2 * p^1.25), 1 / p^1.5)
prior_vec_list[[2]] <- c(1 / (2 * p^1.5), 1 / (2 * p^1.5), 1 / p^2)
prior_vec_list[[3]] <- c(1 / p^2, 1 / p^2, 1 / p^2.25)
prior_vec_list[[4]] <- c(1 / (2 * p^2), 1 / (2 * p^2), 1 / p^2.25)

# Define MCMC parameters
scale_x <- FALSE
intercept <- TRUE

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
  adj_pre <- adj_pre + t(adj_pre)
  adj_act <- adj_act + t(adj_act)
  return(sum((adj_pre - adj_act)^2) / 2)
}
check_adj_l1 <- function(adj_pre, adj_act) {
  adj_pre <- adj_pre + t(adj_pre)
  adj_act <- adj_act + t(adj_act)
  return(sum(abs(adj_pre - adj_act)) / 2)
}

########################### Do one figure ##################################
#### generate graph
set.seed(2021)
n_graph <- 1
graph_sim <- graph_generation(
  K = K, n_graph = n_graph, p = p, n_tol = n_tol,
  e_com = e_com, e_pri = e_pri
)
adj_true <- list()
g_true <- list()
for (iter_K in seq_len(K)) {
  adj_true[[iter_K]] <- t(graph_sim$G[[1]][[iter_K]])
  g_true[[iter_K]] <- as(getGraph(adj_true[[iter_K]]), "graphNEL")
}

#### our method
source("graph_given_order_multi.R")
dta_list <- graph_sim$X[[1]]

#### If we know the order
for (iter_prior in seq_len(length(prior_vec_list))) {
  prior_vec <- prior_vec_list[[iter_prior]]
  out_res <- joint_graph_multi(dta_list = dta_list, prior_vec = prior_vec,
                                 scale_x = scale_x, intercept = intercept)
  print(round(sum(out_res$llike_mat) + sum(out_res$llike_penalty), 4))
  ## Calculate the error
  for (iter_K in seq_len(K)) {
    adj <-  ifelse(out_res$alpha_list[[iter_K]] > 0.5, 1, 0)
    adj <- t(adj)
    g <- as(getGraph(adj),  "graphNEL")
    cat(
      "prior = ", prior_vec[1], prior_vec[3],
      c(shd(g_true[[iter_K]], g), check_edge(adj_true[[iter_K]], adj)),
      "TP", round(TPrate_fun(adj_pre = adj, adj_act = adj_true[[iter_K]]), 4),
      "FP", round(FPrate_fun(adj_pre = adj, adj_act = adj_true[[iter_K]]), 4),
      "FN", round(FNrate_fun(adj_pre = adj, adj_act = adj_true[[iter_K]]), 4),
      "\n"
    )
  }
}

########################## Do MCMC quick test ############################
iter_max <- 200
prior_vec <- prior_vec_list[[2]]
source("graph_mcmc_multi_sim.R")
#### with GES Initialization
# get order
set.seed(2022)
dta <- matrix(NA, nrow = K * n, ncol = p)
for (iter_K in seq_len(K)) {
  dta[(1 + (iter_K - 1) * n):(iter_K * n), ] <- graph_sim$X[[1]][[iter_K]]
}
score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
ges_fit <- ges(score_ges)
ges_adj <- as(ges_fit$repr, "matrix")
ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
order_int <- as.numeric(igraph::topo_sort(graph_i))
# Do MCMC
set.seed(2021)
out_res <- Graph_MCMC_multi(dta_list, scale_x = scale_x, intercept = intercept,
                            order_int = order_int, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
                            prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
                            burn_in = 1, adj_true = adj_true
)
# Show likelihood
library(ggplot2)
library(gridExtra)
gl <- ggplot() + geom_line(aes(x = seq_len(iter_max), y = out_res$llike_vec)) +
  xlab("Iteration") + ylab("Log likelihood")
gs_1 <- ggplot() + geom_line(aes(x = seq_len(iter_max),
                                 y = out_res$error_mat_list[[1]][1, ])) +
  xlab("Iteration") + ylab("SHD")
gu_1 <- ggplot() + geom_line(aes(x = seq_len(iter_max),
                                 y = out_res$error_mat_list[[1]][2, ])) +
  xlab("Iteration") + ylab("No order")
gs_2 <- ggplot() + geom_line(aes(x = seq_len(iter_max),
                                 y = out_res$error_mat_list[[2]][1, ])) +
  xlab("Iteration") + ylab("SHD")
gu_2 <- ggplot() + geom_line(aes(x = seq_len(iter_max),
                                 y = out_res$error_mat_list[[2]][2, ])) +
  xlab("Iteration") + ylab("No order")
layout_matrix <- matrix(c(1, 2, 3), nrow = 3)
grid.arrange(gl, gs_1, gu_1, layout_matrix = layout_matrix)
grid.arrange(gl, gs_2, gu_2, layout_matrix = layout_matrix)

########################### Do parallel ##################################
#### generate graph
source("graph_mcmc_multi.R")
set.seed(2022)
n_graph <- 50
graph_sim <- graph_generation(
  K = K, n_graph = n_graph, p = p, n_tol = n_tol,
  e_com = e_com, e_pri = e_pri
)
iter_max <- 1e5

library(foreach)
library(doParallel)
library(doRNG)

for (iter_prior in seq_len(length(prior_vec_list))) {
  out_res <- list()
  prior_vec <- prior_vec_list[[iter_prior]]
  ## do parallel
  cl <- makeCluster(25)
  registerDoParallel(cl)
  out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
    library(pcalg)
    # get order
    dta <- matrix(NA, nrow = K * n, ncol = p)
    adj_true <- list()
    for (iter_K in seq_len(K)) {
      dta[(1 + (iter_K - 1) * n):(iter_K * n), ] <- graph_sim$X[[iter]][[iter_K]]
      adj_true[[iter_K]] <- t(graph_sim$G[[iter]][[iter_K]])
    }
    score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
    ges_fit <- ges(score_ges)
    ges_adj <- as(ges_fit$repr, "matrix")
    ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
    graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
    order_int <- as.numeric(igraph::topo_sort(graph_i))
    # Do MCMC
    Graph_MCMC_multi(
      dta_list = graph_sim$X[[iter]], scale_x = scale_x, intercept = intercept,
      order_int = order_int, iter_max = iter_max,
      prior_vec = prior_vec, itermax = 100, L_max = 10,
      burn_in = iter_max - 5000
    )
  }
  stopCluster(cl)
  saveRDS(out_res, paste0("prior", iter_prior, "K", K, "e_com", e_com, "e_pri", e_pri, "mcmc.rds"))
  ## check results
  res <- list()
  for (iter_K in seq_len(K)) {
    res[[iter_K]] <- matrix(NA, nrow = n_graph, ncol = 7)
  }
  res_ave <- matrix(0, nrow = n_graph, ncol = 7)
  
  for (iter_graph in seq_len(n_graph)) {
    res_tmp <- out_res[[iter_graph]]
    A_mat_list <- res_tmp$A_list
    alpha_mat_list <- res_tmp$alpha_list
    # analysis
    # for (iter_K in seq_len(K)) {
    #   alpha_mat_list[[iter_K]] <- matrix(0, nrow = p, ncol = p)
    #   A_mat_list[[iter_K]] <- matrix(0, nrow = p, ncol = p)
    # }
    # 
    # for (iter in seq_len(5000)) {
    #   order_tmp <- order(res_tmp$order_list[[iter]])
    #   for (iter_K in seq_len(K)) {
    #     alpha_mat_list[[iter_K]] <- alpha_mat_list[[iter_K]] +
    #       res_tmp$alpha_list[[iter]][[iter_K]][order_tmp, order_tmp]
    #     A_mat_list[[iter_K]] <- A_mat_list[[iter_K]] +
    #       res_tmp$A_list[[iter]][[iter_K]][order_tmp, order_tmp]
    #   }
    # }
    for (iter_K in seq_len(K)) {
      alpha_mat_list[[iter_K]] <- alpha_mat_list[[iter_K]] / 5000
      A_mat_list[[iter_K]] <- A_mat_list[[iter_K]] / 5000
      adj <- ifelse(alpha_mat_list[[iter_K]] > 0.5, 1, 0)
      adj <- t(adj)
      g <- as(getGraph(adj), "graphNEL")
      adj_true <- t(graph_sim$G[[iter_graph]][[iter_K]])
      g_true <- as(getGraph(adj_true), "graphNEL")
      ## save results
      res[[iter_K]][iter_graph, ] <- c(
        shd(g_true, g),
        check_edge(adj_true, adj),
        TPrate_fun(adj_pre = adj, adj_act = adj_true),
        FPrate_fun(adj_pre = adj, adj_act = adj_true),
        FNrate_fun(adj_pre = adj, adj_act = adj_true),
        check_adj_l2(adj_pre = alpha_mat_list[[iter_K]], adj_act = adj_true),
        check_adj_l1(adj_pre = alpha_mat_list[[iter_K]], adj_act = adj_true)
      )
    }
  }
  # show results
  cat("p:", p, "e_com:", e_com, "e_pri", e_pri, "prior:", round(prior_vec, 4), "\n")
  for (iter_K in seq_len(K)) {
    cat("data", iter_K, round(colMeans(res[[iter_K]]), 4), "\n")
    res_ave <- res_ave + res[[iter_K]]
  }
  # average results
  res_ave <- res_ave / K
  # cat("prior:", round(prior_vec, 4), "p:", p, "e_com:", e_com, "e_pri", e_pri, "\n")
  # cat(round(colMeans(res_ave), 4), "\n")
  cat(
    K, "&", "prior", iter_prior, "&", e_com, "&", e_pri, "&",
    res_ave[2], "&", res_ave[3], "&", res_ave[4], "&", res_ave[6], "\\\\", "\n"
  )
}