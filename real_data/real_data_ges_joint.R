# # load data
load("real_data/ovarian.rda")
p <- ncol(data[[1]])

## GES method
library(pcalg)
library(stabs)

#### joint GSE method
source("real_data/method_code/newclass.R")

set.seed(2021)
## learn causal networks
stabs_ges <- function(x, y, q, ...) {
  sample_data <- function(sing_dt) {
    totcol <- nrow(sing_dt)
    sing_dt[sample(1:totcol, as.integer(0.9 * totcol), replace = FALSE), ]
  }
  # Y is the label of the classes, X is the input matrix
  dt <- lapply(data, sample_data)
  # dt <- data
  lambdas <- c(2, 3, 4, 5)
  model_lambda <- function(lambda) {
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
# construct x
x <- do.call(rbind, data)
x <- cbind(x, matrix(0, nrow = nrow(x), ncol = p * (p - 1)))
# construct y
y <- c()
for (i in seq_len(length(data))) {
  y <- c(y, rep(i, nrow(data[[i]])))
}
## stable joint GES
cutoff <- 0.75
stab_result <- stabsel(x = x, y = y, fitfun = stabs_ges, cutoff = cutoff, PFER = 1)
saveRDS(stab_result, "real_data/results/out_ges_joint.rds")

stab_result <- readRDS("real_data/results/out_ges_joint.rds")
cutoff <- 0.6
dag <- matrix(as.vector(stab_result$max > cutoff), nrow = p, ncol = p)
dag <- as(dag, "graphNEL")

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
## data set 1 results
ges_joint_graph1 <- gesdag[[1]]
ges_joint_graph1 <- as(ges_joint_graph1, "matrix")
ges_joint_graph1 <- ifelse(ges_joint_graph1 == 1, TRUE, FALSE)
ges_joint_graph1 <- ges_joint_graph1 | t(ges_joint_graph1)
sum(ges_joint_graph1) / 2

## data set 2
ges_joint_graph2 <- gesdag[[2]]
ges_joint_graph2 <- as(ges_joint_graph2, "matrix")
ges_joint_graph2 <- ifelse(ges_joint_graph2 == 1, TRUE, FALSE)
ges_joint_graph2 <- ges_joint_graph2 | t(ges_joint_graph2)
sum(ges_joint_graph2) / 2

## intersections
ges_joint_graph <- ges_joint_graph1 & ges_joint_graph2
sum(ges_joint_graph) / 2

## 51 52 48
## 55 59 49