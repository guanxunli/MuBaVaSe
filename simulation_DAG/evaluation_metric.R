# calculate the Hamming distance
hammingDistance <- function(G_true, G_estimate) {
  return(sum(abs(G_true - G_estimate)))
}


# recall= TP / (TP + FN)
# precison = TP / (TP + FP)
# False discovery rate (FDR) = FP / (FP + TP)
calculate_metric <- function(G_true, G_estimate) {
  # actual edge
  P <- which(G_true == 1)
  # predicted edge
  PP <- which(G_estimate == 1)
  # true positive
  TP <- length(intersect(P, PP))
  # calculate metric
  recall <- TP / length(P)
  precision <- TP / length(PP)
  fdr <- 1- precision
  return(list(recall = recall, precision = precision, fdr = fdr))
}


