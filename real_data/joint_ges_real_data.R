library(pcalg)
library(graph)
library(stabs)

set.seed(1)

cutoff <- 0.6
source("simulation_DAG/newclass.R")

#load the drosophila data
load("real_data/ovarian.rda")
p <- ncol(data[[1]])

#learn causal networks
stabs.ges <- function(x, y, q, ...){
  #Y is the label of the classes, X is the input matrix
  dt <- lapply(data, function(sing.dt) {totcol <- nrow(sing.dt); sing.dt[sample(1:totcol, as.integer(0.9 * totcol), replace=FALSE), ]})
  #dt <- lapply(1:max(y), function(i) x[y==i, 1:p])
  print(ncol(x))
  lambdas <- c(2,3,4,5)
  model.lambda <- function(lambda){
    l0score <- new("MultiGaussL0pen", data = dt, lambda = lambda * log(ncol(dt[[1]])), intercept = TRUE, use.cpp = FALSE)
    ges.fit <- ges(l0score)
    dag <- as(ges.fit$essgraph, "matrix")
    as.vector(dag != 0)
  }
  path <- sapply(lambdas, model.lambda)
  selected <- rowSums(path) != 0
  return(list(selected=selected, path=path))
}

#construct x
x <- do.call(rbind, data)
x <- cbind(x, matrix(0, nrow=nrow(x), ncol=p * (p-1)))

#construct y
y <- c()
for(i in 1:length(data)){
  y <- c(y, rep(i, nrow(data[[i]])))
}

stab.result <- stabsel(x = x, y = y, fitfun = stabs.ges, cutoff = cutoff, PFER = 1)
save(stab.result, "out_joint_ges.rds")

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
cutoff_vec <- seq(0.5, 0.9, by = 0.05)

for (iter in seq_len(length(cutoff_vec))) {
  cutoff <- cutoff_vec[iter]
  dag <- matrix(as.vector(stab_result$max > cutoff), nrow = p, ncol = p)
  dag <- as(dag, "graphNEL")
  
  gesdag <- ges_alg(data, dag)
  ## data set 1 results
  ges_joint_graph1 <- gesdag[[1]]
  ges_joint_graph1 <- as(ges_joint_graph1, "matrix")
  ges_joint_graph1 <- ifelse(ges_joint_graph1 == 1, TRUE, FALSE)
  ges_joint_graph1 <- ges_joint_graph1 | t(ges_joint_graph1)
  
  ## data set 2
  ges_joint_graph2 <- gesdag[[2]]
  ges_joint_graph2 <- as(ges_joint_graph2, "matrix")
  ges_joint_graph2 <- ifelse(ges_joint_graph2 == 1, TRUE, FALSE)
  ges_joint_graph2 <- ges_joint_graph2 | t(ges_joint_graph2)
  
  ## intersections
  ges_joint_graph <- ges_joint_graph1 & ges_joint_graph2
  
  ## check results
  cat("joint GES &", cutoff, "&", sum(ges_joint_graph1) / 2, "&", 
      sum(ges_joint_graph2) / 2, "&",  sum(ges_joint_graph) / 2, "\\\\\n")
  # cat("cutoff: ", cutoff, c(sum(ges_joint_graph1), sum(ges_joint_graph2), 
  #                           sum(ges_joint_graph)) / 2, "\n")
}