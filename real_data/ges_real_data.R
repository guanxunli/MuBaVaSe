library(pcalg)
library(graph)
library(stabs)

set.seed(1)

load("ovarian/data/ovarian.rda")
cutoff <- 0.6

#learn causal networks
stabs.ges <- function(x, y, q, ...){
  #Y is the label of the classes, X is the input matrix
  idx <- y[1]
  totcol <- nrow(data[[idx]])
  dt <- data[[idx]][sample(1:totcol, as.integer(0.9 * totcol), replace=FALSE), ]
  
  #train the model
  lambdas <- c(1,2,3,4,5)
  model.lambda <- function(lambda){
    l0score <- new("GaussL0penObsScore", data=dt, lambda = lambda * log(ncol(dt)), intercept = FALSE, use.cpp = TRUE)
    ges.fit <- ges(l0score)
    dag <- as(ges.fit$essgraph, "matrix")
    as.vector(dag != 0)
  }
  
  #get the path and selected variables
  path <- sapply(lambdas, model.lambda)
  selected <- rowSums(path) != 0
  return(list(selected=selected, path=path))
}

#run jobs on cluster
stab.inputlist <- lapply(1:length(data), function(i) {p <- ncol(data[[i]]); 
dt <- cbind(as.matrix(data[[i]]), matrix(0, nrow=nrow(data[[i]]), ncol=p * (p-1))); 
list(x=dt, y=rep(i, nrow(data[[i]])))
})
gesdag.list <- lapply(stab.inputlist, function(stab.input) stabsel(x = stab.input$x, y = stab.input$y, fitfun = stabs.ges, cutoff = cutoff, PFER = 1))

#save to file
save(gesdag.list, "out_ges.rds")

cutoff_vec <- seq(0.5, 0.9, by = 0.05)
for (iter in seq_len(length(cutoff_vec))) {
  cutoff <- cutoff_vec[iter]
  ## data set 1 results
  ges_adj1 <- matrix(as.vector(gesdag_list[[1]]$max > cutoff), nrow = p, ncol = p)
  ges_adj1 <- ges_adj1 | t(ges_adj1)
  ## data set 2
  ges_adj2 <- matrix(as.vector(gesdag_list[[2]]$max > cutoff), nrow = p, ncol = p)
  ges_adj2 <- ges_adj2 | t(ges_adj2)
  ## intersections
  ges_adj <- ges_adj1 & ges_adj2
  ## check results
  cat("GES &", cutoff, "&", sum(ges_adj1) / 2, "&", sum(ges_adj2) / 2, "&",  sum(ges_adj) / 2, "\\\\\n")
  # cat("cutoff: ", cutoff, c(sum(ges_adj1), sum(ges_adj2), sum(ges_adj)) / 2, "\n")
}
