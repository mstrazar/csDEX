\name{estimateSizeFactors-methods}
\docType{methods}
\alias{estimateSizeFactors-methods}
\alias{estimateSizeFactors,csDEXdataSet-method}
\title{ ~~ Methods for Function \code{estimateSizeFactors}  ~~}
\description{
 ~~ Methods for function \code{estimateSizeFactors}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(obj = "csDEXdataSet")}}{
    Estimate the size factors hyperparameter for each feature by using the external DESeq2::estimateSizeFactorsForMatrix function. Required for normalization of samples with respect to sequencing depth and mandatory when fitting count models.
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
