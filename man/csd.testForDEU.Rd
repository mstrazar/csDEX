\name{csd.testForDEU}
\alias{csd.testForDEU}
\alias{csd.geneModel}
\alias{csd.waldTest}

\title{
Test for differential usage of exons using the PSI or the count models. 
}
\description{
Difference of deviance tests for interactions between each feature and condition for an instance of a gene.
Produces a ranked list of significant candidate interactions. 
The function csd.geneModel invoke by csd.testForDEU and possibly dispatched to multiple cores. 
The csd.waldTest is invoked by csd.geneModel if model approximation is used.
}
\usage{
csd.testForDEU(cdx, workers = 1, tmp.dir = NULL, min.cpm = NULL, alpha.wald = NULL)
csd.geneModel(input, min.cpm = NULL, tmp.dir = NULL, dist="count", alpha.wald = NULL)
csd.waldTest(mm0, model0, alpha=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdx}{
    A csDEXdataSet object.
}
  \item{workers}{
    Number of cores.
}
  \item{tmp.dir}{
    Temporary directory. If set, results for individual genes are stored immediately as they are computed.
  }
  \item{min.cpm}{
    Minimal counts per million reads (CPM) for a gene to be considered. No positions within a gene with a CPM below the set threshold will be tested. An output data frame is still produced, with appropriate message in the msg field.
}
  \item{alpha.wald}{
  Significance threshold for model approximation based on the Wald test for parameter significance. If set,
  any parameters non-significantly different from zero in the null model are merged to a single "error" parameter.
  Can significantly reduce computation time for large and sparse datasets.
}
}
\details{
The core function of the package. Evaluate significance of all possible pairs of interactions between features and conditions by means of difference in deviance model selection. Generalized linear models (GLMs) based on the Negative Binomial distribution are assumed for the quantification based on read counts. The Beta distribution GLM is used to model data based on  percentage-spliced (PSI) quantification.  

The returned data frame containing interactions will be sorted according to statistical significance (determined by the p-value and padj). The residual value can be used to test the absolute deviation of the null model from the expected value. In the event of an optimization error (numerical problem in matrix inversion or non-convergence) in null or alternative model fitting, the feature-condition pair will be declared non-testable, with corresponding error message stored in the message field.
 
The read count model requires additional hyperparameters, to be stored in the csDEXDataSet object, which need to be precomputed by package methods estimateSizeFactors, estimateDispersions, estimateGeneCPM, in that order. The PSI model can be run on a constructed csDEXDataSet object directly.

The time complexity of GLM fitting is cubic in the number of parameters, which equal the sum of number of features and conditions. A model approximation based on Wald test of parameter significance can be used by setting the alpha.wald argument to a significance threshold in (0, 1), e.g. 0.05. The internal method csd.waldTest is invoked for model approximation. This leads to faster execution, without affecting the probability of Type I error (i.e. false-positive rate), but affecting the Type II error (i.e. false-negative rate). The returned values nrow, ncol and time can be used to compare the savings with respect to the full model. See references for a detailed description of the approximation.

}
\value{
  \item{featureID }{Feature identifier.}
  \item{groupID }{Gene identifier.}
  \item{condition }{Condition identifier.}
  \item{y}{The expression value of feature in the given experimental condition. 
    A number in open interval (0, 1) for the PSI model or a non-negative integer for the count model.}
  \item{cpm }{Gene count per million reads (CPM), used by the count model.}
  \item{pvalue }{Statistical significance of the interaction (probability of type 1 error).}
  \item{padj }{Bonferroni-adjusted p-value on a gene-level.}
  \item{residual }{Model residual. The difference between actual value y and the fitted value of the null model.}
  \item{testable }{Whether the interaction complies to used defined threshold and model optimization converged without error.}
  \item{loglik }{Alternative model log-likelihood.}
  \item{LR }{Difference of deviance.}
  \item{time }{Time for alternative model fitting.}
  \item{nrow }{Number of rows in the model design matrix.}
  \item{ncol }{Number of columns in the model design matrix.}
  \item{msg }{Error/warning message, if any.}
}
\references{
%% TODO: add preprints.
%% ~put references to the literature/web site here ~
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
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
