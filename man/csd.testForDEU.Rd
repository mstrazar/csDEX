\name{csd.testForDEU}
\alias{csd.testForDEU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for differential usage of exons.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
csd.testForDEU(cdx, workers = 1, tmp.dir = NULL, min.cpm = NULL, alpha.wald = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdx}{
%%     ~~Describe \code{cdx} here~~
}
  \item{workers}{
%%     ~~Describe \code{workers} here~~
}
  \item{tmp.dir}{
%%     ~~Describe \code{tmp.dir} here~~
}
  \item{min.cpm}{
%%     ~~Describe \code{min.cpm} here~~
}
  \item{alpha.wald}{
%%     ~~Describe \code{alpha.wald} here~~
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
function (cdx, workers = 1, tmp.dir = NULL, min.cpm = NULL, alpha.wald = NULL) 
{
    gene.set = unique(rowData(cdx)$groupID)
    inputs = list()
    dist = dataType(cdx)
    if (dist == "count") {
        stopifnot(!is.null(colData(cdx)$size.factor))
        stopifnot(!is.null(rowData(cdx)$dispersion))
    }
    gene.cpm = NULL
    message("Constructing inputs ... \n")
    for (g in gene.set) {
        inxs = which(rowData(cdx)$groupID == g)
        gene.expr = exprData(cdx, normalized = (dist == "count"))[inxs, 
            ]
        gene.rowdata = rowData(cdx)[inxs, ]
        if (!is.null(min.cpm)) 
            gene.cpm = cpmData(cdx)[g, ]
        inputs[[length(inputs) + 1]] = list(expr = gene.expr, 
            rowdata = gene.rowdata, gene.cpm = gene.cpm)
    }
    message(sprintf("Dispatching on \%d workers, cache directory \%s. \n", 
        workers, tmp.dir))
    jobs = mclapply(inputs, csd.geneModel, min.cpm, tmp.dir, 
        dist, alpha.wald, mc.preschedule = FALSE, mc.cores = workers, 
        mc.silent = TRUE)
    results = data.frame()
    for (i in 1:length(jobs)) {
        if (!class(jobs[[i]]) == "try-error") 
            results = rbind(results, jobs[[i]])
    }
    results = results[order(results$pvalue), ]
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
