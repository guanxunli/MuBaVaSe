get_objective = function (X, Y, s) {
  return(Eloglik(X,Y,s) - sum(s$KL))
}

# Expected loglikelihood for a susie fit.
Eloglik = function (X, Y, s) {
  n = nrow(X)
  return(-(n/2) * log(2*pi*s$sigma2) - (1/(2*s$sigma2)) * get_ER2(X,Y,s))
}

# Expected squared residuals.
# s$Xr is column sum of Xr_L
get_ER2 = function (X, Y, s) {
  Xr_L = compute_MXt(s$alpha * s$mu,X) # L by N matrix
  postb2 = s$alpha * s$mu2 # Posterior second moment.
  return(sum((Y - s$Xr)^2) - sum(Xr_L^2) + sum(attr(X,"d") * t(postb2)))
}

SER_posterior_e_loglik = function (X, Y, s2, Eb, Eb2) {
  n = nrow(X)
  return(-0.5*n*log(2*pi*s2) - 0.5/s2*(sum(Y*Y)
                                       - 2*sum(Y*compute_Xb(X,Eb))
                                       + sum(attr(X,"d") * Eb2)))
}

estimate_residual_variance_fun = function (X, y, s) {
  n = nrow(X)
  return((1/n)*get_ER2(X,y,s))
}