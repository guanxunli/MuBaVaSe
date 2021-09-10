## generate joint graph
# p : the number of features
# n : the number of data points
# n_graph : number of graphs
# pro_c : probability of edges for the common part
# pro_1, pro_2 : probability of edges for the data set 1 and data set 2, respectively 
# signal_level : scale for the coefficients matrix
# sig : sd for the error term
graph_generate_unif <- function(p, n, n_graph = 1, pro_c = 0.1, pro_1 = 0.02, pro_2 = 0.02, 
                                signal_level = 0.2, sig = 1) {
  # G is the true graph
  # X is the data set from the true graph
  # A is the coefficients matrix
  X <- list()
  G <- list()
  A <- list()
  
  for (iter_graph in seq_len(n_graph)) {
    G[[iter_graph]] <- list()
    # Generate graph first
    G_1 <- matrix(0, p, p)
    G_2 <- matrix(0, p, p)
    edge_c <- rbinom(p * (p - 1) / 2, 1, prob = pro_c)
    edge_1 <- edge_c
    edge_1[which(edge_1 == 0)] <- rbinom(length(which(edge_1 == 0)), 1, prob = pro_1)
    edge_2 <- edge_c
    edge_2[which(edge_2 == 0)] <- rbinom(length(which(edge_2 == 0)), 1, prob = pro_2)
    G_1[lower.tri(G_1)] <- edge_1
    G_2[lower.tri(G_2)] <- edge_2
    G[[iter_graph]]$G_1 <- G_1
    G[[iter_graph]]$G_2 <- G_2
    
    # Generate coefficient matrix
    A[[iter_graph]] <- list()
    A_1 <- matrix(0, nrow = p, ncol = p)
    A_2 <- matrix(0, nrow = p, ncol = p)
    n_1 <- sum(G_1)
    n_2 <- sum(G_2)
    A_1[which(G_1 == 1)] <-  2 * rep(signal_level, n_1) *  (rbinom(n_1, 1, 0.5) - 1/2)
    A_2[which(G_2 == 1)] <-  2 * rep(signal_level, n_2) *  (rbinom(n_2, 1, 0.5) - 1/2)
    A[[iter_graph]]$A_1 <- A_1
    A[[iter_graph]]$A_2 <- A_2
    
    # Generate data
    X[[iter_graph]] <- list()
    eps_1 <- matrix(rnorm(n * p, mean = 0, sd = sig), ncol = n)
    X_1 <- solve(diag(1, nrow = p) - A_1, eps_1)
    eps_2 <- matrix(rnorm(n * p, mean = 0, sd = sig), ncol = n)
    X_2 <- solve(diag(1, nrow = p) - A_2, eps_2)
    X[[iter_graph]]$X_1 <- X_1   
    X[[iter_graph]]$X_2 <- X_2
  }
  # return results
  return(list(X = X, G = G, A = A))
}

## generate joint graph
# p : the number of features
# n : the number of data points
# n_graph : number of graphs
# pro_c : probability of edges for the common part
# pro_1, pro_2 : probability of edges for the data set 1 and data set 2, respectively 
# signal_sig : sd for the coefficient matrix
# sig : sd for the error term
graph_generate_gaussian <- function(p, n, n_graph = 1, pro_c = 0.1, pro_1 = 0.02, pro_2 = 0.02, 
                                    signal_sig = 1, sig = 0.2) {
  # G is the true graph
  # X is the data set from the true graph
  # A is the coefficients matrix
  X <- list()
  G <- list()
  A <- list()
  
  for (iter_graph in seq_len(n_graph)) {
    G[[iter_graph]] <- list()
    # Generate graph first
    G_1 <- matrix(0, p, p)
    G_2 <- matrix(0, p, p)
    edge_c <- rbinom(p * (p - 1) / 2, 1, prob = pro_c)
    edge_1 <- edge_c
    edge_1[which(edge_1 == 0)] <- rbinom(length(which(edge_1 == 0)), 1, prob = pro_1)
    edge_2 <- edge_c
    edge_2[which(edge_2 == 0)] <- rbinom(length(which(edge_2 == 0)), 1, prob = pro_2)
    G_1[lower.tri(G_1)] <- edge_1
    G_2[lower.tri(G_2)] <- edge_2
    G[[iter_graph]]$G_1 <- G_1
    G[[iter_graph]]$G_2 <- G_2
    
    # Generate coefficient matrix
    A[[iter_graph]] <- list()
    A_1 <- matrix(0, nrow = p, ncol = p)
    A_2 <- matrix(0, nrow = p, ncol = p)
    n_1 <- sum(G_1)
    n_2 <- sum(G_2)
    A_1[which(G_1 == 1)] <-  rnorm(n_1, mean = 0, sd = signal_sig)
    A_2[which(G_2 == 1)] <-  rnorm(n_2, mean = 0, sd = signal_sig)
    A[[iter_graph]]$A_1 <- A_1
    A[[iter_graph]]$A_2 <- A_2
    
    # Generate data
    X[[iter_graph]] <- list()
    eps_1 <- matrix(rnorm(n * p, mean = 0, sd = sig), ncol = n)
    X_1 <- solve(diag(1, nrow = p) - A_1, eps_1)
    eps_2 <- matrix(rnorm(n * p, mean = 0, sd = sig), ncol = n)
    X_2 <- solve(diag(1, nrow = p) - A_2, eps_2)
    X[[iter_graph]]$X_1 <- X_1   
    X[[iter_graph]]$X_2 <- X_2
  }
  # return results
  return(list(X = X, G = G, A = A))
}
