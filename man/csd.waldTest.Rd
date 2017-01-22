\name{csd.waldTest}
\alias{csd.waldTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wald test to save coefficients.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
csd.waldTest(mm0, model0, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mm0}{
%%     ~~Describe \code{mm0} here~~
}
  \item{model0}{
%%     ~~Describe \code{model0} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (mm0, model0, alpha = 0.05) 
{
    beta0 = model0$coefficients$mean
    scores = diag(model0$vcov)[1:length(beta0)]
    wald = rep(0, length(scores))
    names(wald) = names(scores)
    for (i in 1:length(wald)) {
        w = wald.test(b = beta0[i], scores[i], T = 1)
        wald[i] = w$result$chi2["P"]
    }
    sig = names(beta0)[wald < alpha]
    not.sig = names(beta0)[wald >= alpha]
    if (length(not.sig) > 1) {
        mu.prior = mm0[, not.sig] \%*\% beta0[not.sig]
        colnames(mu.prior) = "mu.prior"
        mm.prior = cbind(mu.prior, mm0[, sig])
        coef.prior = c(1, beta0[sig])
    }
    else {
        mm.prior = mm0
        coef.prior = beta0
    }
    return(list(sig = sig, not.sig = not.sig, mm.prior = mm.prior, 
        coef.prior = coef.prior))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
