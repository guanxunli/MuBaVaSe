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

## learn causal networks
stabs_ges <- function(x, y, q, ...){
  sample_data <- function(sing_dt) {
    totcol <- nrow(sing_dt)
    sing_dt[sample(1:totcol, as.integer(0.9 * totcol), replace=FALSE), ]
  }
  #Y is the label of the classes, X is the input matrix
  dt <- lapply(data, sample_data)
  lambdas <- c(2,3,4,5)
  model_lambda <- function(lambda){
    l0score <- new("MultiGaussL0pen", data = dt, lambda = lambda * log(ncol(dt[[1]])), intercept = TRUE, use.cpp = FALSE)
    ges_fit <- ges(l0score)
    dag <- as(ges_fit$essgraph, "matrix")
    as.vector(dag != 0)
  }
  path <- sapply(lambdas, model_lambda)
  selected <- rowSums(path) != 0
  return(list(selected = selected, path = path))
}

## joint GES the first step
cutoff <- 0.6
#construct x
x <- do.call(rbind, data)
x <- cbind(x, matrix(0, nrow=nrow(x), ncol=p * (p-1)))
#construct y
y <- c()
for(i in 1:length(data)){
  y <- c(y, rep(i, nrow(data[[i]])))
}
## stable joint GES
stab_result <- stabsel(x = x, y = y, fitfun = stabs_ges, cutoff = cutoff, PFER = 1)
p <- ncol(data[[1]])
dag <- matrix(as.vector(stab_result$max > cutoff), nrow = p, ncol = p)

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
ges_alg <- function(data, dag) {
  in_mat <- as(pdag2dag(dag)$graph, "matrix")
  joint_mat <- lapply(data, function(dt) sapply(seq_len(ncol(dt)), function(i) subset(i, which(in_mat[, i] != 0), dt)))
  return(lapply(joint_mat, function(sing_mat) dag2cpdag(as(sing_mat, "graphNEL"))))
}

gesdag <- ges_alg(data, dag)

#### Calculate the error
## data set 1
adj_1 <- gesdag[[1]]
adj_1 <- ifelse(adj_1 == TRUE, 1, 0)
g_1 <- as(getGraph(ges_adj1), "graphNEL")
# structural Hamming distance (SHD) and undirected edge
print(c(shd(g_true1, g_1), check_edge(adj_true1, adj_1)))
# TPR & FPR
print(c(round(TPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4), 
        round(FPrate_fun(adj_pre = adj_1, adj_act = adj_true1), 4)))

## data set 2
adj_2 <- gesdag[[2]]
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