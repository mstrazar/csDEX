\name{nbreg.fit}
\alias{nbreg.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Negative Binomial regression model.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nbreg.fit(X, y, phi, beta.init = NULL, offset = 0, max.iter = 100, tol = 1e-06, lambda.0 = 1e-04, lambda.max = 100, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{beta.init}{
%%     ~~Describe \code{beta.init} here~~
}
  \item{offset}{
%%     ~~Describe \code{offset} here~~
}
  \item{max.iter}{
%%     ~~Describe \code{max.iter} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{lambda.0}{
%%     ~~Describe \code{lambda.0} here~~
}
  \item{lambda.max}{
%%     ~~Describe \code{lambda.max} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, y, phi, beta.init = NULL, offset = 0, max.iter = 100, 
    tol = 1e-06, lambda.0 = 1e-04, lambda.max = 100, verbose = FALSE) 
{
    nb.log.likelihood <- function(y, mu, phi) {
        l = sum(log(dnbinom(y, mu = mu, size = phi)))
        return(l)
    }
    b.der <- function(mu, phi) {
        return(phi/(mu * (mu + phi)))
    }
    b.der.2 <- function(mu, phi) {
        return(-(2 * mu * phi + phi^2)/(mu^2 * (mu + phi)^2))
    }
    c.der <- function(mu, phi) {
        return(-phi/(mu + phi))
    }
    c.der.2 <- function(mu, phi) {
        return(phi/(mu + phi)^2)
    }
    var.nb <- function(mu, phi) {
        v = (b.der.2(mu, phi) * c.der(mu, phi) - c.der.2(mu, 
            phi) * b.der(mu, phi))/b.der(mu, phi)^3
        v[is.na(v) | is.infinite(v)] = max(v[!(is.na(v) | is.infinite(v))])
        stopifnot(v >= 0)
        return(v)
    }
    U.nb <- function(X, y, beta, phi, inv.link = exp, inv.link.der = exp) {
        p = ncol(X)
        mu = inv.link(X \%*\% beta)
        v = var.nb(mu, phi)
        U = t(X) \%*\% ((y - mu)/v * inv.link.der(X \%*\% beta))
        return(U)
    }
    J.nb <- function(X, beta, phi, inv.link = exp, inv.link.der = exp) {
        p = ncol(X)
        mu = inv.link(X \%*\% beta)
        v = var.nb(mu, phi)
        D = diag(as.vector(1/v * inv.link.der(X \%*\% beta)^2))
        J = (t(X) \%*\% D) \%*\% X
        return(J)
    }
    n = nrow(X)
    p = ncol(X)
    link = log
    inv.link = exp
    log.likelihood = nb.log.likelihood
    score.function = U.nb
    inf.function = J.nb
    beta = beta.init
    if (is.null(beta.init)) {
        N = inv.link(offset)
        rate = y/N
        w = 1 * N/(1 + phi * N)
        weighted.avg.y = log(sum(w * rate)/sum(w))
        beta = qr.coef(qr(X), rep(weighted.avg.y, length = n))
    }
    mu = drop(inv.link(X \%*\% beta + offset))
    jnxs = 1:p
    converged = FALSE
    lambda = lambda.0
    for (itr in 1:max.iter) {
        beta.prev = beta
        beta.trial = rep(0, p)
        ymax = max(y)
        mu.prev = inv.link(X \%*\% beta.prev)
        l.prev = log.likelihood(y, mu.prev, phi)
        if (verbose) 
            message(sprintf("Likelihood (start): \%f", l.prev))
        U = score.function(X, y, beta, phi)
        J = inf.function(X, beta, phi)
        repeat {
            G = chol(J + lambda * diag(p), pivot = TRUE)
            while (is.null(attr(G, "rank")) || attr(G, "rank") < 
                p) {
                lambda = lambda * 10
                G = chol(J + lambda * diag(p), pivot = TRUE)
            }
            jnxs = attr(G, "pivot")
            b = backsolve(G, backsolve(G, U[jnxs], transpose = TRUE))
            beta.trial[jnxs] = beta.prev[jnxs] + b
            mu.trial = inv.link(X \%*\% beta.trial)
            l.trial = log.likelihood(y, mu.trial, phi)
            if (verbose) 
                message(sprintf("Likelihood (lambda=\%f): \%f", 
                  lambda, l.trial))
            if (l.trial > l.prev) {
                if (verbose) 
                  message(sprintf("Likelihood improvement achieved with lambda \%f", 
                    lambda))
                beta = beta.trial
                mu = mu.trial
                lambda = lambda/10
                break
            }
            if (lambda > lambda.max) {
                if (verbose) 
                  message(sprintf("Lambda too large: \%f", lambda))
                beta = beta.prev
                break
            }
            else {
                lambda = lambda * 10
                if (verbose) 
                  message(sprintf("Repeating with lambda \%f", 
                    lambda))
            }
        }
        l.current = log.likelihood(y, mu, phi)
        if (verbose) 
            message(sprintf("Iteration: \%d, Log-likelihood: \%f", 
                itr, l.current))
        if (abs(l.current - l.prev) < tol) {
            converged = TRUE
            if (verbose) 
                message(sprintf("Converged with likelihood: \%f", 
                  l.current))
            break
        }
        else if (lambda > lambda.max) {
            warning(sprintf("NOT Converged with likelihood: \%f", 
                l.current))
            break
        }
    }
    if (!converged) 
        warning(sprintf("NOT Converged with likelihood: \%f", 
            l.current))
    J = inf.function(X, beta, phi)
    mu.fit = c(inv.link(X \%*\% beta))
    if (!is.null(colnames(X))) {
        colnames(J) <- colnames(X)
        row.names(J) <- colnames(X)
        names(beta) <- colnames(X)
    }
    if (!is.null(row.names(X))) {
        names(mu.fit) = row.names(X)
    }
    return(list(coefficients = list(mean = beta), fitted.values = mu.fit, 
        converged = converged, phi = phi, loglik = l.current, 
        iter = itr, vcov = J))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
