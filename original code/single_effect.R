single_effect_regression =
  function (y, X, V, residual_variance = 1, prior_weights = NULL,
            optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
            check_null_threshold = 0) {
    optimize_V = "optim"
    Xty = compute_Xty(X,y)
    betahat = (1/attr(X,"d")) * Xty
    shat2 = residual_variance/attr(X,"d")
    if (is.null(prior_weights))
      prior_weights = rep(1/ncol(X),ncol(X))
    if (optimize_V != "EM" && optimize_V != "none")
      V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
                                  alpha = NULL,post_mean2 = NULL,V_init = V,
                                  check_null_threshold = check_null_threshold)
    
    # log(bf) for each SNP
    lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
      dnorm(betahat,0,sqrt(shat2),log = TRUE)
    
    # Deal with special case of infinite shat2 (e.g., happens if X does
    # not vary).
    lbf[is.infinite(shat2)] = 0 
    maxlbf = max(lbf)
    
    # w is proportional to BF, but subtract max for numerical stability.
    w = exp(lbf - maxlbf)
    
    # Posterior prob for each SNP.
    w_weighted = w * prior_weights
    weighted_sum_w = sum(w_weighted)
    alpha = w_weighted / weighted_sum_w
    post_var = (1/V + attr(X,"d")/residual_variance)^(-1) # Posterior variance.
    post_mean = (1/residual_variance) * post_var * Xty
    post_mean2 = post_var + post_mean^2 # Second moment.
    
    # BF for single effect model.
    lbf_model = maxlbf + log(weighted_sum_w)
    loglik = lbf_model + sum(dnorm(y,0,sqrt(residual_variance),log = TRUE))
    
    if(optimize_V == "EM")
      V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
                                  alpha,post_mean2,
                                  check_null_threshold = check_null_threshold)
    
    return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
                lbf_model = lbf_model,V = V,loglik = loglik))
  }

optimize_prior_variance = function (optimize_V, betahat, shat2, prior_weights,
                                    alpha = NULL, post_mean2 = NULL,
                                    V_init = NULL, check_null_threshold = 0) {
  V = V_init
  if (optimize_V != "simple") {
    if(optimize_V == "optim") {
      lV = optim(par = log(max(c(betahat^2-shat2,1),na.rm = TRUE)),
                 fn = neg.loglik.logscale,betahat = betahat,shat2 = shat2,
                 prior_weights = prior_weights,method = "Brent",lower = -30,
                 upper = 15)$par
      ## if the estimated one is worse than current one, don't change it.
      if(neg.loglik.logscale(lV, betahat = betahat,shat2 = shat2,prior_weights = prior_weights) > 
         neg.loglik.logscale(log(V), betahat = betahat,
                             shat2 = shat2,prior_weights = prior_weights)){
        lV = log(V)
      }
      V = exp(lV)
    } else if (optimize_V == "uniroot")
      V = est_V_uniroot(betahat,shat2,prior_weights)
    else if (optimize_V == "EM")
      V = sum(alpha * post_mean2)
    else
      stop("Invalid option for optimize_V method")
  }
  
  if (loglik(0,betahat,shat2,prior_weights) +
      check_null_threshold >= loglik(V,betahat,shat2,prior_weights))
    V = 0
  return(V)
}


loglik = function (V, betahat, shat2, prior_weights) {
  
  #log(bf) for each SNP
  lbf = dnorm(betahat,0,sqrt(V+shat2),log = TRUE) -
    dnorm(betahat,0,sqrt(shat2),log = TRUE)
  
  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0 
  
  maxlbf = max(lbf)
  w = exp(lbf - maxlbf) # w = BF/BFmax
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  return(log(weighted_sum_w) + maxlbf)
}

neg.loglik.logscale = function(lV,betahat,shat2,prior_weights)
  -loglik(exp(lV),betahat,shat2,prior_weights)

loglik.grad = function(V, betahat, shat2, prior_weights) {
  
  # log(bf) for each SNP.
  lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
    dnorm(betahat,0,sqrt(shat2),log = TRUE)
  
  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0 
  
  maxlbf = max(lbf)
  w = exp(lbf - maxlbf) # w = BF/BFmax
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  return(sum(alpha * lbf.grad(V,shat2,betahat^2/shat2)))
}

# Define loglikelihood and gradient as function of lV:=log(V)
# to improve numerical optimization
negloglik.grad.logscale = function (lV, betahat, shat2, prior_weights)
  -exp(lV) * loglik.grad(exp(lV),betahat,shat2,prior_weights)

# Vector of gradients of logBF_j for each j, with respect to prior
# variance V.
lbf.grad = function (V, shat2, T2) {
  l = 0.5*(1/(V + shat2)) * ((shat2/(V + shat2))*T2 - 1)
  l[is.nan(l)] = 0
  return(l)
}

lbf = function (V, shat2, T2) {
  l = 0.5*log(shat2/(V + shat2)) + 0.5*T2*(V/(V + shat2))
  l[is.nan(l)] = 0
  return(l)
}
