# # load data
load("real_data/ovarian.rda")
p <- ncol(data[[1]])
library(pcalg)
library(graph)
################################ with out stable selection ########################
alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)

pc_fun <- function(dta, alpha) {
  p <- ncol(dta)
  dta_cor <- cor(dta)
  alpha <- alphas[iter_alpha]
  pc_fit <- pc(
    suffStat = list(C = dta_cor, n = dim(dta)[1]),
    indepTest = gaussCItest, alpha = alpha,
    labels = sapply(1:p, toString)
  )
  dag <- as(pc_fit@graph, "matrix")
  return(ifelse(dag == TRUE, 1, 0))
}

for (iter_alpha in seq_len(length(alphas))) {
  alpha_use <- alphas[iter_alpha]
  dag1 <- pc_fun(data[[1]], alpha_use)
  dag2 <- pc_fun(data[[2]], alpha_use)
  pc_adj1 <- dag1 | t(dag1)
  n1 <- sum(pc_adj1) / 2
  pc_adj2 <- dag2 | t(dag2)
  n2 <- sum(pc_adj2) / 2
  pc_adj <- pc_adj1 & pc_adj2
  n_com <- sum(pc_adj) / 2
  n_total <- n1 + n2 - n_com
  n_ratio <- n_com / n_total
  ## check results
  cat("PC & $\\alpha = ", alpha_use, "$&", n1, "&", n2, "&", n_com,
      "&", n_total, "&", round(n_ratio, 4), "\\\\\n")
  # cat("alpha: ", alpha_use, c(sum(pc_adj1), sum(pc_adj2), sum(pc_adj)) / 2, "\n")
}

################################ with stable selection ########################
library(stabs)
set.seed(1)
cutoff_vec <- seq(0.6, 0.9, by = 0.05)
## PC input
stab_input <- function(i) {
  p <- ncol(data[[i]])
  dt <- cbind(as.matrix(data[[i]]), matrix(0, nrow = nrow(data[[i]]), ncol = p * (p - 1)))
  return(list(x = dt, y = rep(i, nrow(data[[i]]))))
}
stab_input_list <- lapply(seq_len(length(data)), stab_input)

## learn causal networks
stabs_pc <- function(x, y, q, ...) {
  # Y is the label of the classes, X is the input matrix
  idx <- y[1]
  totcol <- nrow(data[[idx]])
  dt <- data[[idx]][sample(1:totcol, as.integer(0.9 * totcol), replace = FALSE), ]
  # dt <- data[[idx]]
  p <- ncol(dt)

  # train the model
  alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
  model_alpha <- function(alpha) {
    pc_fit <- pc(
      suffStat = list(C = cor(dt), n = nrow(dt)),
      indepTest = gaussCItest, alpha = alpha,
      labels = sapply(1:p, toString)
    )
    dag <- as(pc_fit@graph, "matrix")
    as.vector(dag != 0)
  }

  # get the path and selected variables
  path <- sapply(alphas, model_alpha)
  selected <- rowSums(path) != 0
  return(list(selected = selected, path = path))
}

library(parallel)

pcdag_fun <- function(cutoff) {
  res <- lapply(
    stab_input_list,
    function(stab_input) stabsel(x = stab_input$x, y = stab_input$y, fitfun = stabs_pc, cutoff = cutoff, PFER = 1)
  )
  return(res)
}

pcdag_list <- mclapply(cutoff_vec, pcdag_fun,
  mc.cores = length(cutoff_vec)
)
saveRDS(pcdag_list, "out_pc.rds")

#### check results
pcdag_list <- readRDS("real_data/results/out_pc.rds")
cutoff_vec2 <- seq(0.5, 0.9, by = 0.05)
for (iter in seq_len(length(cutoff_vec))) {
  pcdag_list_tmp <- pcdag_list[[iter]]
  for (iter2 in seq_len(length(cutoff_vec2))) {
    cutoff <- cutoff_vec2[iter2]
    ## data set 1 results
    pc_adj1 <- matrix(as.vector(pcdag_list_tmp[[1]]$max > cutoff), nrow = p, ncol = p)
    pc_adj1 <- pc_adj1 | t(pc_adj1)
    n1 <- sum(pc_adj1) / 2
    ## data set 2
    pc_adj2 <- matrix(as.vector(pcdag_list_tmp[[2]]$max > cutoff), nrow = p, ncol = p)
    pc_adj2 <- pc_adj2 | t(pc_adj2)
    n2 <-  sum(pc_adj2) / 2
    ## intersections
    pc_adj <- pc_adj1 & pc_adj2
    n_com <- sum(pc_adj) / 2
    n_total <- n1 + n2 - n_com
    n_ratio <- n_com / n_total
    cat("PC &", cutoff_vec[iter], "&", cutoff, "&", n1, "&", n2, "&", n_com,
        "&", n_total, "&", round(n_ratio, 4), "\\\\\n")
    # cat(
    #   "cutoff1: ", cutoff_vec[iter], "cutoff2: ", cutoff,
    #   c(sum(pc_adj1), sum(pc_adj2), sum(pc_adj)) / 2, "\n"
    # )
  }
}