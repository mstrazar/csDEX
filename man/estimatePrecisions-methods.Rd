\name{estimatePrecisions-methods}
\docType{methods}
\alias{estimatePrecisions-methods}
\alias{estimatePrecisions,csDEXdataSet-method}
\title{ ~~ Methods for Function \code{estimatePrecisions}  ~~}
\description{
 ~~ Methods for function \code{estimatePrecisions}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(obj = "csDEXdataSet")}}{
  Estimate the precisions hyperparameter for each feature by using the external edgeR::estimateDisp function. Required for fitting read-count based models. The function estimateSizeFactors should be called on the cdx object prior to estimatePrecisions call.
}
}}
\keyword{methods}
\keyword{models}
