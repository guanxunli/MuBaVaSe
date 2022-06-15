# # load data
load("real_data/ovarian.rda")
p <- ncol(data[[1]])
library(pcalg)
library(graph)

################################ with out stable selection ########################
lambdas <- c(1, 2, 3, 4, 5)

## Joint GES the first step
ges_joint_fun <- function(data, lambda) {
  source("real_data/newclass.R")
  p <- ncol(data[[1]])
  dag_list <- list()
  l0score <- new("MultiGaussL0pen",
    data = data, lambda = lambda * log(p),
    intercept = TRUE, use.cpp = FALSE
  )
  ges_fit <- ges(l0score)
  dag <- as(ges_fit$essgraph, "matrix")
  return(ifelse(dag == TRUE, 1, 0))
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
ges_alg <- function(dag_com, dta) {
  in_mat <- dag_com
  joint_mat <- sapply(seq_len(ncol(dta)), function(i) subset(i, which(in_mat[, i] != 0), dta))
  adj_pri <- joint_mat
  return(adj_pri)
}

for (iter in seq_len(length(lambdas))) {
  lambda_use <- lambdas[iter]
  dag_com <- ges_joint_fun(data, lambda_use)
  dag1 <- ges_alg(dag_com, data[[1]])
  dag2 <- ges_alg(dag_com, data[[2]])
  ## data set 1 results
  ges_joint_graph1 <- as(dag1, "matrix")
  ges_joint_graph1 <- ifelse(ges_joint_graph1 == 1, TRUE, FALSE)
  ges_joint_graph1 <- ges_joint_graph1 | t(ges_joint_graph1)
  n1 <- sum(ges_joint_graph1) / 2
  ## data set 2
  ges_joint_graph2 <- as(dag2, "matrix")
  ges_joint_graph2 <- ifelse(ges_joint_graph2 == 1, TRUE, FALSE)
  ges_joint_graph2 <- ges_joint_graph2 | t(ges_joint_graph2)
  n2 <- sum(ges_joint_graph2) / 2
  ## intersections
  ges_joint_graph <- ges_joint_graph1 & ges_joint_graph2
  n_com <- sum(ges_joint_graph) / 2
  n_total <- n1 + n2 - n_com
  n_ratio <- n_com / n_total
  ## check results
  cat("joint GES &$\\lambda = ", lambda_use, "$&", n1, "&", n2, "&", n_com,
      "&", n_total, "&", round(n_ratio, 4), "\\\\\n")
  # cat("lambda: ", lambda_use, c(
  #   sum(ges_joint_graph1), sum(ges_joint_graph2),
  #   sum(ges_joint_graph)
  # ) / 2, "\n")
}

################################ with stable selection ########################
library(stabs)
#### joint GSE method
source("real_data/newclass.R")
cutoff_vec <- seq(0.6, 0.9, by = 0.05)
set.seed(1)

## learn causal networks
stabs_ges <- function(x, y, q, ...) {
  sample_data <- function(sing_dt) {
    totcol <- nrow(sing_dt)
    sing_dt[sample(1:totcol, as.integer(0.9 * totcol), replace = FALSE), ]
  }
  # Y is the label of the classes, X is the input matrix
  dt <- lapply(data, sample_data)
  # dt <- data
  lambdas <- c(1, 2, 3, 4, 5)
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
stab_fun <- function(cutoff) {
  return(stabsel(x = x, y = y, fitfun = stabs_ges, cutoff = cutoff, PFER = 1))
}

stab_result_list <- mclapply(cutoff_vec, stab_fun, mc.cores = length(cutoff_vec))
saveRDS(stab_result_list, "out_ges_joint.rds")

## Joint GES the second step
stab_result_list <- readRDS("real_data/results/out_ges_joint.rds")
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

cutoff_vec2 <- seq(0.5, 0.9, by = 0.05)

for (iter in seq_len(length(cutoff_vec))) {
  stab_result <- stab_result_list[[iter]]
  for (iter2 in seq_len(length(cutoff_vec2))) {
    cutoff <- cutoff_vec2[iter2]
    dag <- matrix(as.vector(stab_result$max > cutoff), nrow = p, ncol = p)
    dag <- as(dag, "graphNEL")

    gesdag <- ges_alg(data, dag)
    ## data set 1 results
    ges_joint_graph1 <- gesdag[[1]]
    ges_joint_graph1 <- as(ges_joint_graph1, "matrix")
    ges_joint_graph1 <- ifelse(ges_joint_graph1 == 1, TRUE, FALSE)
    ges_joint_graph1 <- ges_joint_graph1 | t(ges_joint_graph1)
    n1 <- sum(ges_joint_graph1) / 2

    ## data set 2
    ges_joint_graph2 <- gesdag[[2]]
    ges_joint_graph2 <- as(ges_joint_graph2, "matrix")
    ges_joint_graph2 <- ifelse(ges_joint_graph2 == 1, TRUE, FALSE)
    ges_joint_graph2 <- ges_joint_graph2 | t(ges_joint_graph2)
    n2 <- sum(ges_joint_graph2) / 2

    ## intersections
    ges_joint_graph <- ges_joint_graph1 & ges_joint_graph2
    n_com <- sum(ges_joint_graph) / 2
    n_total <- n1 + n2 - n_com
    n_ratio <- n_com / n_total
    ## check results
    cat("joint GES &", cutoff_vec[iter], "&", cutoff, "&", n1, "&", n2, "&", n_com,
        "&", n_total, "&", round(n_ratio, 4), "\\\\\n")
    # cat(
    #   "cutoff1: ", cutoff_vec[iter], "cutoff2: ", cutoff,
    #   c(sum(ges_joint_graph1), sum(ges_joint_graph2), sum(ges_joint_graph)) / 2, "\n"
    # )
  }
}