library(pcalg)
library(graph)
library(stabs)

set.seed(1)

load("real_data/ovarian.rda")
genenames <- colnames(data[[1]])
p <- length(genenames)

cutoff <- 0.6
#learn causal networks
stabs.pc <- function(x, y, q, ...){
  #Y is the label of the classes, X is the input matrix
  idx <- y[1]
  totcol <- nrow(data[[idx]])
  dt <- data[[idx]][sample(1:totcol, as.integer(0.9 * totcol), replace=FALSE), ]
  p <- ncol(dt)
  
  #train the model
  alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
  model.alpha <- function(alpha){
    pc.fit <- pc(suffStat=list(C=cor(x), n=dim(x)[1]), indepTest = gaussCItest, alpha=alpha, labels = sapply(1:p, toString))
    dag <- as(pc.fit@graph, "matrix")
    as.vector(dag != 0)
  }
  
  #get the path and selected variables
  path <- sapply(alphas, model.alpha)
  selected <- rowSums(path) != 0
  return(list(selected=selected, path=path))
}

#run jobs on cluster
stab.inputlist <- lapply(1:length(data), function(i) {p <- ncol(data[[i]]); 
dt <- cbind(as.matrix(data[[i]]), matrix(0, nrow=nrow(data[[i]]), ncol=p * (p-1))); 
list(x=dt, y=rep(i, nrow(data[[i]])))
})
gesdag.list <- lapply(stab.inputlist, function(stab.input) stabsel(x = stab.input$x, y = stab.input$y, fitfun = stabs.pc, cutoff = cutoff, PFER = 1))

## save to file
save(gesdag.list, "out_pc.rds")
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
  cat("PC &", cutoff, "&", sum(ges_adj1) / 2, "&", sum(ges_adj2) / 2, "&",  sum(ges_adj) / 2, "\\\\\n")
  # cat("cutoff: ", cutoff, c(sum(ges_adj1), sum(ges_adj2), sum(ges_adj)) / 2, "\n")
}
