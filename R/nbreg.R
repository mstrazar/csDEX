# Iterative reweighted least-squares fitting of Negative-Binomial general linear models


nbreg.fit <- function (X, y, phi, beta.init=NULL, offset=0,
    max.iter=100, tol=1e-6, lambda.0=0.0001, lambda.max=100, verbose=FALSE){

    # A General method to perform Iterative reweighted least-squares for a negative-binomial GLM.

    # Data
    #   X: Design matrix.
    #   y: Target vector.

    # Prior parameters for a previously fitted model. Not to be confused with the null model.
    #   offset: Expected offset from NB distributed data

    # Hyperparameters
    #   phi: Precision.
    #   tol: Convergence criterion tolerance.
    #   max.iter: Allowed number of iterations.
    #   seed: Random initialization seed.
    #   lambda.0: Initial Levenberg-Marquardt fudging coefficient.
    #   lambda.max: Upper-bound on lambda.

    # Return:
    #   beta: parameters fit to current data.
  
    nb.log.likelihood <- function(y, mu, phi) {
      # Log likelihood for the NB
      l = sum(log(dnbinom(y, mu=mu, size=phi)))
      return(l)
    }
    
    # NB model auxilliary variables
    b.der <- function(mu, phi){
      return(phi / (mu * (mu+phi)))
    }
    
    b.der.2 <- function(mu, phi){
      return( - (2 * mu * phi + phi^2)/(mu^2 * (mu + phi)^2))
    }
    
    c.der <- function(mu, phi){
      return( - phi / (mu + phi))
    }
    
    c.der.2 <- function(mu, phi){
      return(phi / (mu + phi)^2)
    }
    
    var.nb <- function(mu, phi){
      # When phi -> inf, var.nb -> mu, which is consistent with the Poisson model (no overdispersion,
      # variance equal to mean)
      v = (b.der.2(mu, phi) * c.der(mu, phi) - c.der.2(mu, phi) * b.der(mu, phi)) / b.der(mu,phi)^3
      v[is.na(v) | is.infinite(v)] = max(v[!(is.na(v) | is.infinite(v))])
      stopifnot(v >= 0)
      return(v)
    }
    
    U.nb <- function(X, y, beta, phi, inv.link=exp, inv.link.der=exp){
      # Derivative of NB likelihood function
      p = ncol(X)
      mu = inv.link(X %*% beta)
      v =  var.nb(mu, phi)
      U = t(X) %*% ((y - mu) / v * inv.link.der(X %*% beta))
      return(U)
    }
    
    J.nb <- function(X, beta, phi, inv.link=exp, inv.link.der=exp){
      # Matrix of second-order derivatives for NB likelihood function
      p = ncol(X)
      mu = inv.link(X %*% beta)
      v =  var.nb(mu, phi)
      D = diag(as.vector(1 / v * inv.link.der(X %*% beta)^2))
      J = (t(X) %*% D) %*% X
      return(J)
    }

    # Enables providing a prior beta parameters, where some parameters are readily fit, and others are
    # not (indicated by NA).
    n = nrow(X)
    p = ncol(X)

    link = log                              #   mu
    inv.link = exp                          #   nu
    log.likelihood = nb.log.likelihood      #   y, mu, phi
    score.function = U.nb                   #   X, y, beta, phi
    inf.function = J.nb                     #   X, beta, phi
  
    # Initial values for beta ; Very good initialization that lands upon values
    # producing sensible order of magnitude.
    # Basically, it fits the betas to mean log y value.
    beta = beta.init
    if(is.null(beta.init)){
        N = inv.link(offset)
        rate = y/N                                          # Targets minus offset
        w = 1 * N/(1 + phi * N)                             # weights are inversely-proportional to precision
        weighted.avg.y = log(sum(w * rate) / sum(w))        # weighted average of target signal
        beta = qr.coef(qr(X),rep(weighted.avg.y,length=n))  # Fit features to weighted average
    }
    mu = drop(inv.link(X %*% beta + offset))                # Fitted means
    
    jnxs = 1:p
    converged = FALSE
    lambda = lambda.0

    for (itr in 1:max.iter){

        # Current parameters become previous ; new ones to be computed
        beta.prev = beta
        beta.trial = rep(0, p)

        ymax = max(y)
        mu.prev = inv.link(X %*% beta.prev)

        l.prev = log.likelihood(y, mu.prev, phi)
        if(verbose) message(sprintf("Likelihood (start): %f", l.prev))

        # Compute parameter updatedes
        U = score.function(X, y, beta, phi)
        J = inf.function(X, beta, phi)

        repeat {

            # Fast solution x = inv(J)U -> Jx = U -> RRTx = U -> x = inv(R) (inv(RT) U )
            # Use pivoting just to obtain rank parameter
            G = chol(J + lambda * diag(p), pivot=TRUE)
            while(is.null(attr(G, "rank")) || attr(G, "rank") < p) {
                # Conditional assignment
                lambda = lambda * 10
                G = chol(J + lambda * diag(p), pivot=TRUE)
            }
            jnxs = attr(G, "pivot")

            b = backsolve(G, backsolve(G, U[jnxs], transpose=TRUE))
            beta.trial[jnxs] = beta.prev[jnxs] + b
            mu.trial = inv.link(X %*% beta.trial)
            l.trial = log.likelihood(y, mu.trial, phi)

            # Print debug info
            if(verbose) message(sprintf("Likelihood (lambda=%f): %f", lambda, l.trial))

            # Levenberg-Marquardt tradeoff between gradient step and direct optimization
            if(l.trial > l.prev) {
                if(verbose) message(sprintf("Likelihood improvement achieved with lambda %f", lambda))
                beta = beta.trial
                mu = mu.trial
                lambda = lambda / 10
                break
            }
            if(lambda > lambda.max){
                if(verbose) message(sprintf("Lambda too large: %f", lambda))
                beta = beta.prev
                break
            } else {
                lambda = lambda * 10
                if(verbose) message(sprintf("Repeating with lambda %f", lambda))
            }
        }

        l.current = log.likelihood(y, mu, phi)
        if(verbose) message(sprintf("Iteration: %d, Log-likelihood: %f", itr, l.current))

        # Test for convergence based on likelihood
        if (abs(l.current-l.prev) < tol){
            converged = TRUE
            if(verbose) message(sprintf("Converged with likelihood: %f", l.current))
            break
        } else if (lambda > lambda.max){
            warning(sprintf("NOT Converged with likelihood: %f", l.current))
            break
        }
    }
    if(!converged) 
      warning(sprintf("Iteration limit reached. NOT Converged with likelihood: %f", l.current))
    
    # Consistent naming
    J = inf.function(X, beta, phi)
    mu.fit = c(inv.link(X%*%beta))
    
    if(!is.null(colnames(X))) {
        colnames(J) <- colnames(X)
        row.names(J) <- colnames(X)
        names(beta) <- colnames(X)
    }
    
    if(!is.null(row.names(X))){
      names(mu.fit) = row.names(X)
    }
    
    return(list(coefficients=list(mean=beta), fitted.values=mu.fit,
            converged=converged, phi=phi, loglik=l.current, iter=itr, vcov=J))
}