require(quantreg)

### Utility functions

muvarToAlpha <- function(mu, var) {
  # Mean/variance to alpha
  ((1 - mu) / var - 1 / mu) * mu ^ 2
}

muvarToBeta <- function(mu, var){
  # Mean/variance to beta
  alpha = muvarToAlpha(mu, var)
  alpha * (1 / mu - 1)
}

muphiToVar <- function(mu, phi){
  # Mean/precision to variance
  mu * (1 - mu) / (1 + phi)
}

muphiToAlpha <- function(mu, phi){
  # Mean/precision to alpha
  mu * phi
}
muphiToBeta <- function(mu, phi){
  # Mean/precision to beta
  phi * (1-mu)
}

abToMu <- function(alpha, beta){
  # Alpha, beta to mean
  alpha / (alpha + beta)
}

abToPhi <- function(alpha, beta){
  # Alpha, beta to mean
  alpha + beta
}


# Probability density auxillary functions
dbeta.loglik <- function(y, mu, phi){
  # Log-likelihood beta for one observation
  sum(lgamma(phi) - lgamma(phi * mu) - lgamma(phi - phi * mu) +
        (phi * mu - 1) * log(y) + (phi * (1 - mu) - 1) * log(1-y))
}


### Data conversion
intToZo <- function(x){
  # Arbitrary interval to zero-one [0, 1] ; idempotent
  (x - min(x)) / (max(x) - min(x))
}

zoToint <- function(x, a, b){
  # Interval [0, 1] to [a, b] ; idempotent
  ((b - a) * x) + a
}

zoSqueeze <- function(x){
  # [0, 1] to (0, 1) ; not idempotent
  N = length(x)
  ((N - 1) * x + 0.5) / N
}


### Parameter fitting

estimateBetaPrecisions <- function(X, method="qreg", tau=0.05){
  # Estimate precisions for PSI-distributed data based on sample variance
  # Model is of form log(phi) ~ 1 + v + v^(-1)
  # X: expression matrix in a bounded interval, (rows: features, cols: conditions)
  # method: lm - linear model, qreg - quantile regression
  # tau: target quantile in the qreg model
  # ...: Parameters for lm/qreg fitting

  # TODO: add min/max polynomial degree (allowing negatives)

  # Return values
  #   means: sample means
  #   vars: sample variances
  #   phi: sample precisions
  #   phi.fit: fitted precisions
  #   phi.fit.min: min(phi.fit, phi)
  #   phi.fit.weights: inverse probability density (related to sample variance)
  LM = "lm"
  QREG = "qreg"
  stopifnot(method %in% c(LM, QREG))
  
  # Convert to (0, 1)
  Xo = zoSqueeze(intToZo(X))
  
  # Compute sample variances and precisions
  means = apply(Xo, 1, mean)
  vars = apply(Xo, 1, var)
  alpha = muvarToAlpha(means, vars)
  beta = muvarToBeta(means, vars)
  phi = abToPhi(alpha, beta)
  phi.fit = rep(0, length(phi))
  phi.fit.min = rep(0, length(phi))
  
  #  Name vectors
  names(means) = row.names(Xo)
  names(vars) = row.names(Xo)
  names(phi) = row.names(Xo)
  names(phi.fit) = row.names(Xo)
  names(phi.fit.min) = row.names(Xo)
  
  # Select non-zero entries
  inxs = which(vars > 0)
  nms = row.names(Xo)[inxs]
  y = log(phi[inxs])
  v = matrix(vars[inxs], nrow = length(inxs), ncol=1)
  z = cbind(v^-1, v^0, v^1)
  
  # Density model to compute weights 
  message(sprintf("Fitting variance proabaility density"))
  df = approxfun(density(v))
  weights = 1.0 / df(v)

  # Fit precisions and take minimum(observed, fit)
  if(method == LM){
    message(sprintf("Fitting model with lm.wfit"))
    model = lm.wfit(x=z, y=y, w=weights)  
    pf = exp(model$fitted.values)
    pfm = exp(pmin(y, model$fitted.values))
  } else if (method == QREG){
    message(sprintf("Fitting model with quantreg::rq.wfit (tau=%f)", tau))
    qmodel = rq.wfit(x=z, y=y, weights=weights, tau=tau)
    pf = exp(qmodel$fitted.values)
    pfm = exp(pmin(y, qmodel$fitted.values))
  }

  # Ensure correct order and fill-in missing values (zeros)
  phi.fit[nms] = pf
  phi.fit.min[nms] = pfm
  phi.fit[phi.fit == 0] = min(pf)
  phi.fit.min[phi.fit.min == 0] = min(pfm)
  
  # Return
  data.frame(mean=means, 
       var=vars, phi=phi, phi.fit=phi.fit, phi.fit.min=phi.fit.min,
       row.names=row.names(X))
}


plotPrecisionFit <- function(var, phi, phi.fit, phi.fit.min){
  # Plot observed an fit precisions
  ord = order(var)
  smoothScatter(var[ord], log(phi)[ord], xlab="Variance", ylab="Log precision")
  lines(var[ord], log(phi.fit)[ord], col="black")
  lines(var[ord], log(phi.fit.min)[ord], col="red")
  legend("topright",  
         legend = c("model", "min(data, model)"), 
         fill = c("black", "red"))
}