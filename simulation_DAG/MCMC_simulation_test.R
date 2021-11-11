library(pcalg)
source("simulation_DAG/graph_generation.R")

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

check_weight_l2 <- function(weight_pre, weight_act) {
  weight_pre <- weight_pre + t(weight_pre)
  weight_act <- weight_act + t(weight_act)
  return(sum((weight_pre - weight_act)^2) / 2)
}

check_weight_l1 <- function(weight_pre, weight_act) {
  weight_pre <- weight_pre + t(weight_pre)
  weight_act <- weight_act + t(weight_act)
  return(sum(abs(weight_pre - weight_act)) / 2)
}

## define parameters
p_vec <- c(5, 10, 20, 50)
n_vec <- 6 * p
e_com_vec <- c(2, 6, 10, 25)
e_pri_vec <- c(1, 2, 3, 4)

for (iter in seq_len(4)) {
  p <- p_vec[iter]
  n_tol <- n_vec[iter]
  e_com <- e_com_vec[iter]
  e_pri <- e_pri_vec[iter]
  K <- 2
  n_graph <- 1
  
  #### generate graph
  set.seed(2021)
  graph_sim <- graph_generation(
    K = K, n_graph = n_graph, p = p,
    n_tol = n_tol, e_com = e_com, e_pri = e_pri
  )
  adj_true1 <- t(graph_sim$G[[1]][[1]])
  g_true1 <- as(getGraph(adj_true1), "graphNEL")
  weight_true1 <- t(graph_sim$A[[1]][[1]])
  adj_true2 <- t(graph_sim$G[[1]][[2]])
  g_true2 <- as(getGraph(adj_true2), "graphNEL")
  weight_true2 <- t(graph_sim$A[[1]][[2]])
  
  #### our method
  # source("Two_dataset/Graph_MCMC_two.R")
  source("Two_dataset_new/Graph_MCMC_two.R")
  dta_1 <- graph_sim$X[[1]][[1]]
  dta_2 <- graph_sim$X[[1]][[2]]
  ## GES method
  score1 <- new("GaussL0penObsScore", data = dta_1, intercept = FALSE)
  ges_fit1 <- ges(score1)
  ges_adj1 <- as(ges_fit1$repr, "matrix")
  ges_adj1 <- ifelse(ges_adj1 == TRUE, 1, 0)
  ges_graph1 <- as(ges_fit1$repr, "graphNEL")
  ges_weight1 <- ges_fit1$repr$weight.mat()
  
  score2 <- new("GaussL0penObsScore", data = graph_sim$X[[1]][[2]], intercept = FALSE)
  ges_fit2 <- ges(score2)
  ges_adj2 <- as(ges_fit2$repr, "matrix")
  ges_adj2 <- ifelse(ges_adj2 == TRUE, 1, 0)
  ges_graph2 <- as(ges_fit2$repr, "graphNEL")
  ges_weight2 <- ges_fit2$repr$weight.mat()
  
  ## our method
  out_res <- Graph_MCMC_two(dta_1, dta_2,
                            order_int = NULL, iter_max = 100000, sigma02_int = NULL, sigma2_int = NULL,
                            prior_vec = NULL, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8, burn_in = 95000
  )
  #### analysis results
  alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
  alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
  A_mat_1 <- matrix(0, nrow = p, ncol = p)
  A_mat_2 <- matrix(0, nrow = p, ncol = p)
  out_res <- res[[i]]
  for (iter in seq_len(5000)) {
    order_tmp <- order(out_res$order_list[[iter]])
    alpha_mat_1 <- alpha_mat_1 + out_res$alpha_list_1[[iter]][order_tmp, order_tmp]
    alpha_mat_2 <- alpha_mat_2 + out_res$alpha_list_2[[iter]][order_tmp, order_tmp]
    A_mat_1 <- A_mat_1 + out_res$A_list_1[[iter]][order_tmp, order_tmp]
    A_mat_2 <- A_mat_2 + out_res$A_list_2[[iter]][order_tmp, order_tmp]
  }
  
  alpha_mat_1 <- alpha_mat_1 / 5000
  alpha_mat_2 <- alpha_mat_2 / 5000
  A_mat_1 <- A_mat_1 / 5000
  A_mat_2 <- A_mat_2 / 5000
  
  #### Calculate the error
  ## data set 1
  adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
  adj_1 <- t(adj_1)
  g_1 <- as(getGraph(adj_1), "graphNEL")
  weight_1 <- t(A_mat_1)
  weight_1[which(adj_1 == 0)] <- 0
  
  ## data set 2
  adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
  adj_2 <- t(adj_2)
  g_2 <- as(getGraph(adj_2), "graphNEL")
  weight_2 <- t(A_mat_2)
  weight_2[which(adj_2 == 0)] <- 0
  
  ## output results
  cat(
    "MCMC", p, "&", shd(g_true1, g_1), "&", check_edge(adj_true1, adj_1), "&",
    shd(g_true2, g_2), "&", check_edge(adj_true2, adj_2), "&",
    round(sum((weight_true1 - weight_1)^2), 2), "&", round(check_weight_l2(weight_1, weight_true1), 2), "&",
    round(sum(abs(weight_true1 - weight_1)), 2), "&", round(check_weight_l1(weight_1, weight_true1), 2), "&",
    round(sum((weight_true2 - weight_2)^2), 2), "&", round(check_weight_l2(weight_2, weight_true2), 2), "&",
    round(sum(abs(weight_true2 - weight_2)), 2), "&", round(check_weight_l1(weight_2, weight_true2), 2), "\\\\\n"
  )
  
  cat(
    "GES", "&", shd(g_true1, ges_graph1), "&", check_edge(adj_true1, ges_adj1), "&",
    shd(g_true2, ges_graph2), "&", check_edge(adj_true2, ges_adj2), "&",
    round(sum((weight_true1 - ges_weight1)^2), 2), "&", round(check_weight_l2(ges_weight1, weight_true1), 2), "&",
    round(sum(abs(weight_true1 - ges_weight1)), 2), "&", round(check_weight_l1(ges_weight1, weight_true1), 2), "&",
    round(sum((weight_true2 - ges_weight2)^2), 2), "&", round(check_weight_l2(ges_weight2, weight_true2), 2), "&",
    round(sum(abs(weight_true2 - ges_weight2)), 2), "&", round(check_weight_l1(ges_weight2, weight_true2), 2), "\\\\\n"
  )
}