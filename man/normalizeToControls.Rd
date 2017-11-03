\name{normalizeToControls}
\alias{normalizeToControls}

\title{
  Compare sample values to controls. 
}
\description{
  A typical experimenetal design includes a subset of control samples.
}
\usage{
  normalizeToControls(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
    A csDEXdataSet object.
  }
}
\details{
  Currently, only normalization for PSI-based data is implemented. 
  The control valus are subtracted to original values to produce \emph{delta PSI}.
}
\value{
  \item{obj }{A csDEXdataSet object with normalized epxression values. Control samples are
  removed from the data matrix.}
}
\references{
}
\author{
Martin Stražar
}
\note{
See the vignette for examples.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ normalization }% use one of  RShowDoc("KEYWORDS")
