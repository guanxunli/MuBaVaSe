library(pcalg)
library(stabs)
source("simulation_DAG/graph_generation.R")
# args <- commandArgs()
# p <- as.numeric(args[6])
# n_tol <- as.numeric(args[7])
p <- 100
n_tol <- 2000
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

## PC input
stab_input <- function(i) {
  p <- ncol(data[[i]])
  dt <- cbind(as.matrix(data[[i]]), matrix(0, nrow=nrow(data[[i]]), ncol=p * (p-1)))
  return(list(x=dt, y=rep(i, nrow(data[[i]]))))
}
stab_input_list <- lapply(seq_len(length(data)), stab_input)

## learn causal networks
stabs_pc <- function(x, y, q, ...){
  #Y is the label of the classes, X is the input matrix
  idx <- y[1]
  totcol <- nrow(data[[idx]])
  dt <- data[[idx]][sample(1:totcol, as.integer(0.9 * totcol), replace=FALSE), ]
  p <- ncol(dt)
  
  #train the model
  alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
  model_alpha <- function(alpha){
    pc_fit <- pc(suffStat = list(C = cor(x), n = dim(x)[1]), 
                 indepTest = gaussCItest, alpha=alpha, 
                 labels = sapply(1:p, toString))
    dag <- as(pc_fit@graph, "matrix")
    as.vector(dag != 0)
  }
  
  #get the path and selected variables
  path <- sapply(alphas, model_alpha)
  selected <- rowSums(path) != 0
  return(list(selected = selected, path = path))
}

cutoff <- 0.5
gesdag_list <- lapply(stab_input_list, 
                      function(stab_input) stabsel(x = stab_input$x, y = stab_input$y, fitfun = stabs_pc, cutoff = cutoff, PFER = 1))

#### Calculate the error
## data set 1
adj_1 <- matrix(as.vector(gesdag_list[[1]]$max > cutoff), nrow = p, ncol = p)
adj_1 <- ifelse(adj_1 == TRUE, 1, 0)
g_1 <- as(getGraph(ges_adj1), "graphNEL")
# structural Hamming distance (SHD) and undirected edge
print(c(shd(g_true1, g_1), check_edge(adj_true1, adj_1)))
# TPR & FPR
print(c(round(TPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4), 
        round(FPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4)))

## data set 2
adj_2 <- matrix(as.vector(gesdag_list[[2]]$max > cutoff), nrow = p, ncol = p)
adj_3 <- ifelse(adj_3 == TRUE, 1, 0)
g_2 <- as(getGraph(adj_2), "graphNEL")
# structural Hamming distance (SHD) and undirected edge
print(c(shd(g_true2, g_2), check_edge(adj_true2, adj_2)))
# TPR & FPR
print(c(round(TPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4), 
        round(FPrate_fun(adj_pre = adj_2, adj_act = adj_true2), 4)))

## output results
cat(
  "GES", "&", shd(g_true1, g_1), "&", check_edge(adj_true1, adj_1), "&",
  shd(g_true2, g_2), "&", check_edge(adj_true2, adj_2), "&", "\\\\\n"
)