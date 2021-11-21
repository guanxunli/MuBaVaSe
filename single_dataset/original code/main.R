n <- 500
p <- 1000
sigma <- 1
sigma0 <- 0.6
L <- 20
set.seed(2021)
## Generate data
index_t <- sample(seq_len(p), size = L, replace = FALSE)
b <- rep(0, p)
b[index_t] <- rnorm(L, mean = 0, sd = sigma0)
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- X %*% b + rnorm(n, sd = sigma)

source("single_dataset/original code/initialization.R")
source("single_dataset/original code/single_effect.R")
source("single_dataset/original code/update_effect.R")
source("single_dataset/original code/elbo.R")
susie = function (X,y,L = min(10,ncol(X)),
                  scaled_prior_variance = 0.2,
                  residual_variance = NULL,
                  prior_weights = NULL,
                  null_weight = NULL,
                  standardize = TRUE,
                  intercept = TRUE,
                  estimate_residual_variance = TRUE,
                  estimate_prior_variance = TRUE,
                  estimate_prior_method = c("optim", "EM", "simple"),
                  check_null_threshold = 0,
                  prior_tol = 1e-9,
                  residual_variance_upperbound = Inf,
                  s_init = NULL,
                  coverage = 0.95,
                  min_abs_corr = 0.5,
                  compute_univariate_zscore = FALSE,
                  na.rm = FALSE,
                  max_iter = 100,
                  tol = 1e-3,
                  verbose = FALSE,
                  track_fit = FALSE,
                  residual_variance_lowerbound = var(drop(y))/1e4,
                  refine = FALSE) {
  
  # Process input estimate_prior_method.
  estimate_prior_method = "optim"
  
  # Check input X.
  if (!(is.double(X) & is.matrix(X)) &
      !inherits(X,"CsparseMatrix") &
      is.null(attr(X,"matrix.type")))
    stop("Input X must be a double-precision matrix, or a sparse matrix, or ",
         "a trend filtering matrix")
  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight) && is.null(attr(X,"matrix.type"))) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(X) * (1 - null_weight),ncol(X)),null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight),null_weight)
    X = cbind(X,0)
  }
  if (any(is.na(X)))
    stop("Input X must not contain missing values")
  if (any(is.na(y))) {
    if (na.rm) {
      samples_kept = which(!is.na(y))
      y = y[samples_kept]
      X = X[samples_kept,]
    } else
      stop("Input y must not contain missing values")
  }
  
  # Check input y.
  p = ncol(X)
  n = nrow(X)
  mean_y = mean(y)
  
  # Center and scale input.
  if (intercept)
    y = y - mean_y
  X = set_X_attributes(X,center = intercept,scale = standardize)
  
  # Initialize susie fit.
  s = init_setup(n,p,L,scaled_prior_variance,residual_variance,prior_weights,
                 null_weight,as.numeric(var(y)),standardize)
  if (!missing(s_init) && !is.null(s_init)) {
    if (!inherits(s_init,"susie"))
      stop("s_init should be a susie object")
    if (max(s_init$alpha) > 1 || min(s_init$alpha) < 0)
      stop("s_init$alpha has invalid values outside range [0,1]; please ",
           "check your input")
    
    # First, remove effects with s_init$V = 0
    s_init = susie_prune_single_effects(s_init)
    num_effects = nrow(s_init$alpha)
    if(missing(L)){
      L = num_effects
    }else if(L < num_effects){
      warning(paste("Specified number of effects L =",L,
                    "is smaller than the number of effects",num_effects,
                    "in input SuSiE model. The SuSiE model will have",
                    num_effects,"effects."))
      L = num_effects
    }
    # expand s_init if L > num_effects.
    s_init = susie_prune_single_effects(s_init, L, s$V)
    s = modifyList(s,s_init)
    s = init_finalize(s,X = X)
  } else {
    s = init_finalize(s)
  }
  # Initialize elbo to NA.
  elbo = rep(as.numeric(NA),max_iter + 1)
  elbo[1] = -Inf;
  tracking = list()
  
  for (i in 1:max_iter) {
    s = update_each_effect(X,y,s,estimate_prior_variance,estimate_prior_method,
                           check_null_threshold)
    if (verbose)
      print(paste0("objective:",get_objective(X,y,s)))
    
    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance before the update.
    elbo[i+1] = get_objective(X,y,s)
    if ((elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }
    if (estimate_residual_variance) {
      s$sigma2 = pmax(residual_variance_lowerbound,
                      estimate_residual_variance_fun(X,y,s))
      if (s$sigma2 > residual_variance_upperbound)
        s$sigma2 = residual_variance_upperbound
      if (verbose)
        print(paste0("objective:",get_objective(X,y,s)))
    }
  }
  
  # Remove first (infinite) entry, and trailing NAs.
  elbo = elbo[2:(i+1)]
  s$elbo = elbo
  s$niter = i
  
  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }
  
  if (intercept) {
    
    # Estimate unshrunk intercept.
    s$intercept = mean_y - sum(attr(X,"scaled:center") *
                                 (colSums(s$alpha * s$mu)/attr(X,"scaled:scale")))
    s$fitted = s$Xr + mean_y
  } else {
    s$intercept = 0
    s$fitted = s$Xr
  }
  s$fitted = drop(s$fitted)
  names(s$fitted) = `if`(is.null(names(y)),rownames(X),names(y))
  
  # For prediction.
  s$X_column_scale_factors = attr(X,"scaled:scale")
  return(s)
}

#### check results
## package
res <- susieR::susie(X = X, y = Y, L = L)
res$elbo
res1 <- as.numeric(res$sets$cs)
length(intersect(res1, index_t)) / L
length(intersect(res1, index_t)) / length(res1)
sum((colSums(res$alpha * res$mu) - b)^2)

## My code
res <- susie(X = X, y = Y, L = L)
res$elbo
