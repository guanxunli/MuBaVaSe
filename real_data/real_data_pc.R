# # load data
load("real_data/ovarian.rda")

## PC method
library(pcalg)
library(stabs)

set.seed(2021)
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

#### check results
cutoff <- 0.5
## data set 1 results
ges_adj1 <- matrix(as.vector(gesdag_list[[1]]$max > cutoff), nrow = p, ncol = p)
ges_adj1 <- ges_adj1 | t(ges_adj1)
sum(ges_adj1) / 2

## data set 2
ges_adj2 <- matrix(as.vector(gesdag_list[[2]]$max > cutoff), nrow = p, ncol = p)
ges_adj2 <- ges_adj2 | t(ges_adj2)
sum(ges_adj2) / 2

## intersections
ges_adj <- ges_adj1 & ges_adj2
sum(ges_adj)