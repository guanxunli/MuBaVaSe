# # load data
load("real_data/ovarian.rda")
p <- ncol(data[[1]])

## GES method
library(pcalg)
library(stabs)
library(parallel)

set.seed(2021)
#### GES method
## GES input
stab_input <- function(i) {
  p <- ncol(data[[i]])
  dt <- cbind(as.matrix(data[[i]]), matrix(0, nrow = nrow(data[[i]]), ncol = p * (p - 1)))
  return(list(x = dt, y = rep(i, nrow(data[[i]]))))
}
stab_input_list <- lapply(seq_len(length(data)), stab_input)

# stable GES function
stab_ges <- function(x, y, q, ...) {
  # Y is the label of the classes, X is the input matrix
  idx <- y[1]
  totcol <- nrow(data[[idx]])
  # dt <- data[[idx]][sample(1:totcol, as.integer(0.9 * totcol), replace = FALSE), ]
  dt <- data[[idx]]
  
  # train the model
  lambdas <- c(1, 2, 3, 4, 5)
  model_lambda <- function(lambda) {
    l0score <- new("GaussL0penObsScore", data = dt, lambda = lambda * log(ncol(dt)), intercept = FALSE)
    ges_fit <- ges(l0score)
    dag <- as(ges_fit$essgraph, "matrix")
    as.vector(dag != 0)
  }
  
  # get the path and selected variables
  path <- sapply(lambdas, model_lambda)
  selected <- rowSums(path) != 0
  return(list(selected = selected, path = path))
}

cutoff <- 0.6
gesdag_list <- lapply(
  stab_input_list,
  function(stab_input) stabsel(x = stab_input$x, y = stab_input$y, fitfun = stab_ges, cutoff = cutoff, PFER = 1)
)
saveRDS(gesdag_list, "real_data/results/out_ges.rds")

#### check results
gesdag_list <- readRDS("real_data/results/out_ges.rds")
cutoff <- 0.75
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
sum(ges_adj) / 2

## 66 102 32