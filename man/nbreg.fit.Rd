\name{nbreg.fit}
\alias{nbreg.fit}
\title{
    Negative Binomial regression.
}
\description{
    Fit a Negative Binomial regression model. Called internally by the geneModel function.
}
\usage{
    nbreg.fit(X, y, phi, beta.init = NULL, offset = 0, max.iter = 100, tol = 1e-06,
    lambda.0 = 1e-04, lambda.max = 100, verbose = FALSE)
}
\arguments{
  \item{X}{
    The design matrix with n rows and p columns.}

  \item{y}{
    A vector of n target values. 
}
  \item{phi}{
    Precision hyperparameter. Either a scalar or a vector of n values. 
}

  \item{beta.init}{
    Optional. Initial estimate of parameter values. 
}
  \item{offset}{
    Optional. Known offset of target values from a negative binomially distributed data.

}
  \item{max.iter}{
    Number of optimization iteration.
}

  \item{tol}{
    Relative change in log-likelihood to declare convergence.
}

  \item{lambda.0}{
    Initial Levenberg-Marquardt damping factor.
}

  \item{lambda.max}{
    Maximal Levenberg-Marquardt damping factor.
}
  \item{verbose}{
    Print diagnostic messages.
}
}
\details{
    Custom implementation of the Negative Binomial regression model fitting
with Levenberg-Marquardt update, enabling custom initialization values.
}

\value{
    \item{coefficients }{Parameter values for the mean of each fitted value.}
    \item{fitted.values }{Means of n fitted values.} 
    \item{converged }{Binary convergence factor. If FALSE, use verbose=TRUE to determine possible reasons.} 
    \item{phi }{Precision.} 
    \item{loglik }{Log-likelihood of model parameters given data.} 
    \item{iter }{Number of executed gradient steps.} 
    \item{vcov }{Parameter variance-covariance matrix.} 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Martin Stražar
}
\note{
}

\seealso{

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
