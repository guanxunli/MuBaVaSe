## Define parameters
n <- 500
p <- 1000
p_c <- 25
p_1 <- 5
p_2 <- 5
sigma <- 1
sigma0 <- 0.6
set.seed(2021)
## Generate data
index_c <- sample(seq_len(p), size = p_c, replace = FALSE)
index_1 <- sample(setdiff(seq_len(p), index_c), size = p_1, replace = FALSE)
index_2 <- sample(setdiff(seq_len(p), c(index_1, index_c)), size = p_2, replace = FALSE)

b_1 <- rep(0, p)
b_1[c(index_c, index_1)] <- rnorm(p_c + p_1, mean = 0, sd = sigma0)
# b_1[c(index_c, index_1)] <- c(rep(1,15), rep(0.05, 10), rep(0.1, 5))
b_2 <- rep(0, p)
b_2[c(index_c, index_2)] <- rnorm(p_c + p_2, mean = 0, sd = sigma0)
# b_2[c(index_c, index_2)] <- c(rep(0.05,15), rep(1, 10), rep(0.1, 5))

alpha_1 <- rep(0, p)
alpha_1[c(index_c, index_1)] <- 1
alpha_2 <- rep(0, p)
alpha_2[c(index_c, index_2)] <- 1

X_1 <- matrix(rnorm(p * n), nrow = n, ncol = p)
X_2 <- matrix(rnorm(p * n), nrow = n, ncol = p)
Y_1 <- X_1 %*% b_1 + rnorm(n, sd = sigma)
Y_2 <- X_2 %*% b_2 + rnorm(n, sd = sigma)

## Variable selection seperatly
res1 <- susieR::susie(X = X_1, y = Y_1, L = p_1 + p_c)
res2 <- susieR::susie(X = X_2, y = Y_2, L = p_2 + p_c)

## Joint inference
source("original code//Multi_dataset.R")
res <- sum_single_effect_multi(X_1, Y_1, X_2, Y_2, L = p_1 + p_c + p_2, r = 0.2, q = 0.05)
source("Multi_dataset_null.R")
res_null <- sum_single_effect_multi_null(X_1, Y_1, X_2, Y_2, L = p_1 + p_c + p_2, r = 1, q = 1, tau = 1.5)

#### check results
## data set 1
# l1 error 
print(paste0("Single method l1 error for data set 1 is ", round(sum(abs(1 - apply(1 - res1$alpha, 2, prod) - alpha_1)), 4)))
print(paste0("Joint method l1 error for data set 1 without null model is ", round(sum(abs(res$alpha_1 - alpha_1)), 4)))
print(paste0("Joint method l1 error for data set 1 with null model is ", round(sum(abs(res_null$alpha_1 - alpha_1)), 4)))

# l2 error
print(paste0("Single method l2 error for data set 1 is ", round(sum((1 - apply(1 - res1$alpha, 2, prod) - alpha_1)^2), 4)))
print(paste0("Joint method l2 error for data set 1 without null model is ", round(sum((res$alpha_1 - alpha_1)^2), 4)))
print(paste0("Joint method l2 error for data set 1 with null model is ", round(sum((res_null$alpha_1 - alpha_1)^2), 4)))

## data set 2
# l1 error 
print(paste0("Single method l1 error for data set 2 is ", round(sum(abs(1 - apply(1 - res2$alpha, 2, prod) - alpha_2)), 4)))
print(paste0("Joint method l1 error for data set 2 without null model is ", round(sum(abs(res$alpha_2 - alpha_2)), 4)))
print(paste0("Joint method l1 error for data set 2 with null model is ", round(sum(abs(res_null$alpha_2 - alpha_2)), 4)))

# l2 error
print(paste0("Single method l2 error for data set 2 is ", round(sum((1 - apply(1 - res2$alpha, 2, prod) - alpha_2)^2), 4)))
print(paste0("Joint method l2 error for data set 2 without null model is ", round(sum((res$alpha_2 - alpha_2)^2), 4)))
print(paste0("Joint method l2 error for data set 2 with null model is ", round(sum((res_null$alpha_2 - alpha_2)^2), 4)))

#### check posterior mean
## data set 1
# l1 error 
print(paste0("Single method l1 error for data set 1 is ", round(sum(abs(colSums(res1$mu * res1$alpha) - b_1)), 4)))
print(paste0("Joint method l1 error for data set 1 without null model is ", round(sum(abs(res$post_mean1 - b_1)), 4)))
print(paste0("Joint method l1 error for data set 1 with null model is ", round(sum(abs(res_null$post_mean1 - b_1)), 4)))

# l2 error
print(paste0("Single method l2 error for data set 1 is ", round(sum((colSums(res1$mu * res1$alpha) - b_1)^2), 4)))
print(paste0("Joint method l2 error for data set 1 without null model is ", round(sum((res$post_mean1 - b_1)^2), 4)))
print(paste0("Joint method l2 error for data set 1 with null model is ", round(sum((res_null$post_mean1 - b_1)^2), 4)))

## data set 2
# l1 error 
print(paste0("Single method l1 error for data set 2 is ", round(sum(abs(colSums(res2$mu * res2$alpha) - b_2)), 4)))
print(paste0("Joint method l1 error for data set 2 without null model is ", round(sum(abs(res$post_mean2 - b_2)), 4)))
print(paste0("Joint method l1 error for data set 2 with null model is ", round(sum(abs(res_null$post_mean2 - b_2)), 4)))

# l2 error
print(paste0("Single method l2 error for data set 2 is ", round(sum((colSums(res2$mu * res2$alpha) - b_2)^2), 4)))
print(paste0("Joint method l2 error for data set 2 without null model is ", round(sum((res$post_mean2 - b_2)^2), 4)))
print(paste0("Joint method l2 error for data set 2 with null model is ", round(sum((res_null$post_mean2 - b_2)^2), 4)))
