library(graph)
library(pcalg)
library(stabs)
source("simulation_DAG/graph_generation.R")
p <- 50
n_tol <- 600
K <- 2
n <- n_tol / K
n_graph <- 1
cons <- 2
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
#### generate graph
set.seed(2021)
graph_sim <- graph_generation(K = K, n_graph = n_graph, p = p, n_tol = n_tol, e_com = 100, e_pri = 30)
adj_true1 <- t(graph_sim$G[[1]][[1]])
g_true1 <- as(getGraph(adj_true1), "graphNEL")
weight_true1 <- t(graph_sim$A[[1]][[1]])
adj_true2 <- t(graph_sim$G[[1]][[2]])
g_true2 <- as(getGraph(adj_true2), "graphNEL")
weight_true2 <- t(graph_sim$A[[1]][[2]])
# generate data
data <- graph_sim$X[[1]]
#### GES method
cutoff <- 0.6
stabs_ges <- function(x, y, q, ...) {
  # Y is the label of the classes, X is the input matrix
  idx <- y[1]
  totcol <- nrow(data[[idx]])
  dt <- data[[idx]][sample(1:totcol, as.integer(0.9 * totcol), replace = FALSE), ]
  
  # train the model
  lambdas <- c(1, 2, 3, 4, 5)
  model.lambda <- function(lambda) {
    l0score <- new("GaussL0penObsScore", data = dt, lambda = lambda * log(ncol(dt)), intercept = FALSE, use.cpp = TRUE)
    ges.fit <- ges(l0score)
    dag <- as(ges.fit$essgraph, "matrix")
    as.vector(dag != 0)
  }
  
  # get the path and selected variables
  path <- sapply(lambdas, model.lambda)
  selected <- rowSums(path) != 0
  return(list(selected = selected, path = path))
}

# run jobs on cluster
stab_inputlist <- list()
for (iter in seq_len(length(data))) {
  p <- ncol(data[[iter]])
  dt <- cbind(as.matrix(data[[iter]]), matrix(0, nrow = nrow(data[[iter]]), ncol = p * (p - 1)))
  stab_inputlist[[iter]] <- list(x = dt, y = rep(iter, nrow(data[[iter]])))
}

gesdag_list <- list()
for (iter in seq_len(length(stab_inputlist))) {
  stab.input <- stab_inputlist[[iter]]
  gesdag_list[[iter]] <- stabsel(x = stab.input$x, y = stab.input$y, fitfun = stabs_ges, cutoff = cutoff, PFER = 1)
}

#### check results
cutoff <- 0.75
## data set 1 results
ges_adj1 <- matrix(as.vector(gesdag_list[[1]]$max > cutoff), nrow = p, ncol = p)
ges_adj1 <- ifelse(ges_adj1 == TRUE, 1, 0)
ges_graph1 <- as(getGraph(ges_adj1), "graphNEL")
# structural Hamming distance (SHD) and undirected edge
print(c(shd(g_true1, ges_graph1), check_edge(adj_true1, ges_adj1)))
# TPR & FPR
print(c(round(TPrate_fun(adj_pre = ges_adj1, adj_act = adj_true1), 4), round(FPrate_fun(adj_pre = ges_adj1, adj_act = adj_true1), 4)))

## data set 2
ges_adj2 <- matrix(as.vector(gesdag_list[[2]]$max > cutoff), nrow = p, ncol = p)
ges_adj2 <- ifelse(ges_adj2 == TRUE, 1, 0)
ges_graph2 <- as(getGraph(ges_adj2), "graphNEL")
# structural Hamming distance (SHD) and undirected edge
print(c(shd(g_true2, ges_graph2), check_edge(adj_true2, ges_adj2)))
# TPR & FPR
print(c(round(TPrate_fun(adj_pre = ges_adj2, adj_act = adj_true2), 4), round(FPrate_fun(adj_pre = ges_adj2, adj_act = adj_true2), 4)))

#### joint GES method
cutoff <- 0.6
## The first step
source("simulation_DAG/newclass.R")
p <- ncol(data[[1]])
# learn causal networks
stabs_ges <- function(x, y, q, ...) {
  # Y is the label of the classes, X is the input matrix
  dt <- lapply(data, function(sing.dt) {
    totcol <- nrow(sing.dt)
    sing.dt[sample(1:totcol, as.integer(0.9 * totcol), replace = FALSE), ]
  })

  lambdas <- c(2, 3, 4, 5)
  model.lambda <- function(lambda) {
    l0score <- new("MultiGaussL0pen", data = dt, lambda = lambda * log(ncol(dt[[1]])), intercept = TRUE, use.cpp = FALSE)
    ges.fit <- ges(l0score)
    dag <- as(ges.fit$essgraph, "matrix")
    as.vector(dag != 0)
  }
  path <- sapply(lambdas, model.lambda)
  selected <- rowSums(path) != 0
  return(list(selected = selected, path = path))
}

# construct x
x <- do.call(rbind, data)
x <- cbind(x, matrix(0, nrow = nrow(x), ncol = p * (p - 1)))

# construct y
y <- c()
for (i in seq_len(length(data))) {
  y <- c(y, rep(i, nrow(data[[i]])))
}

stab_result <- stabsel(x = x, y = y, fitfun = stabs_ges, cutoff = cutoff, PFER = 1)

