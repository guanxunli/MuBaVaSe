set_X_attributes = function(X, center = TRUE, scale = TRUE) {
  
  # if X is a trend filtering matrix
  if (!is.null(attr(X,"matrix.type"))) {
    order = attr(X,"order")
    n = ncol(X)
    
    # Set three attributes for X.
    attr(X,"scaled:center") = compute_tf_cm(order,n)
    attr(X,"scaled:scale") = compute_tf_csd(order,n)
    attr(X,"d") = compute_tf_d(order,n,attr(X,"scaled:center"),
                               attr(X,"scaled:scale"),scale,center)
    if (!center)
      attr(X,"scaled:center") = rep(0,n)
    if (!scale)
      attr(X,"scaled:scale") = rep(1,n)
  } else {
    
    # If X is either a dense or sparse ordinary matrix.
    # Get column means.
    cm = colMeans(X,na.rm = TRUE)
    
    # Get column standard deviations.
    csd = compute_colSds(X)
    
    # Set sd = 1 when the column has variance 0.
    csd[csd == 0] = 1
    if (!center)
      cm = rep(0,length = length(cm))
    if (!scale) 
      csd = rep(1,length = length(cm))
    X.std = (t(X) - cm)/csd
    
    # Set three attributes for X.
    attr(X,"d") = rowSums(X.std * X.std)
    attr(X,"scaled:center") = cm
    attr(X,"scaled:scale") = csd
  }
  return(X)
}

compute_colSds = function(X) {
  n = nrow(X)
  return(sqrt((colSums(X^2)/n - (colSums(X)/n)^2)*(n/(n-1))))
}

# Set default susie initialization.
init_setup = function (n, p, L, scaled_prior_variance, residual_variance,
                       prior_weights, null_weight, varY, standardize) {
  if (!is.numeric(scaled_prior_variance) || scaled_prior_variance < 0)
    stop("Scaled prior variance should be positive number")
  if (scaled_prior_variance > 1 && standardize)
    stop("Scaled prior variance should be no greater than 1 when ",
         "standardize = TRUE")
  if(is.null(residual_variance))
    residual_variance = varY
  if(is.null(prior_weights))
    prior_weights = rep(1/p,p)
  else
    prior_weights = prior_weights / sum(prior_weights)
  if(length(prior_weights) != p)
    stop("Prior weights must have length p")
  if (p < L)
    L = p
  s = list(alpha  = matrix(1/p,nrow = L,ncol = p),
           mu     = matrix(0,nrow = L,ncol = p),
           mu2    = matrix(0,nrow = L,ncol = p),
           Xr     = rep(0,n),
           KL     = rep(as.numeric(NA),L),
           lbf    = rep(as.numeric(NA),L),
           lbf_variable = matrix(as.numeric(NA),L,p),
           sigma2 = residual_variance,
           V      = scaled_prior_variance*varY,
           pi     = prior_weights)
  if (is.null(null_weight))
    s$null_index = 0
  else
    s$null_index = p
  class(s) = "susie"
  return(s)
}

init_finalize = function (s, X = NULL, Xr = NULL) {
  if(length(s$V) == 1)
    s$V = rep(s$V,nrow(s$alpha))
  
  # Check sigma2.
  if (!is.numeric(s$sigma2))
    stop("Input residual variance sigma2 must be numeric")
  
  # Avoid problems with dimension if input is a 1 x 1 matrix.
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (s$sigma2 <= 0)
    stop("Residual variance sigma2 must be positive (is your var(Y) zero?)")
  
  # check prior variance
  if (!is.numeric(s$V))
    stop("Input prior variance must be numeric")
  if (!all(s$V >= 0))
    stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (!all(dim(s$mu) == dim(s$alpha)))
    stop("dimension of mu and alpha in input object do not match")
  if (nrow(s$alpha) != length(s$V))
    stop("Input prior variance V must have length of nrow of alpha in ",
         "input object")
  
  # Update Xr.
  if (!missing(Xr))
    s$Xr = Xr
  if (!missing(X))
    s$Xr = compute_Xb(X,colSums(s$mu * s$alpha))
  
  # Reset KL and lbf.
  s$KL = rep(as.numeric(NA),nrow(s$alpha))
  s$lbf = rep(as.numeric(NA),nrow(s$alpha))
  class(s) = "susie"
  return(s)
}

compute_Xb = function (X, b) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")
  
  # Scale Xb.
  if (!is.null(attr(X,"matrix.type")))
    
    # When X is a trend filtering matrix.
    scaled.Xb = compute_tf_Xb(attr(X,"order"),b/csd)
  else
    
    # When X is an ordinary sparse/dense matrix.
    scaled.Xb = tcrossprod(X,t(b/csd))
  
  # Center Xb.
  Xb = scaled.Xb - sum(cm*b/csd)
  return(as.numeric(Xb))
}

compute_Xty = function (X, y) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")
  ytX = crossprod(y,X)
  
  # Scale Xty.
  if (!is.null(attr(X,"matrix.type")))
    
    # When X is a trend filtering matrix.
    scaled.Xty = compute_tf_Xty(attr(X,"order"),y)/csd
  else
    
    # When X is an ordinary sparse/dense matrix.
    scaled.Xty = t(ytX/csd)
  
  # Center Xty.
  centered.scaled.Xty = scaled.Xty - cm/csd * sum(y)
  return(as.numeric(centered.scaled.Xty))
}

compute_MXt = function (M, X) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")
  
  if (!is.null(attr(X,"matrix.type")))
    
    # When X is a trend filtering matrix.
    return(as.matrix(t(apply(M,1,function(b) compute_Xb(X,b)))))
  else
    
    # When X is an ordinary sparse/dense matrix.
    return(as.matrix(tcrossprod(M,sweep(X,2,csd,"/")) -
                       drop(tcrossprod(M,t(cm/csd)))))
  
}