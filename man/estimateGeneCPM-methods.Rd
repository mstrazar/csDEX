\name{estimateGeneCPM-methods}
\docType{methods}
\alias{estimateGeneCPM-methods}
\alias{estimateGeneCPM,csDEXdataSet-method}
\title{ ~~ Methods for Function \code{estimateGeneCPM}  ~~}
\description{
 ~~ Methods for function \code{estimateGeneCPM}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(obj = "csDEXdataSet")}}{
  Estimate counts per million reads (CPM) associated to each gene (groupID). An integer field named input.read.count is required in the design file, specifying the sequencing depth measured by the number of sequenced reads / read pairs. Can be used by both PSI and count models to skip testing for low-expressed genes and avoid Type I errors that are dut to changes in gene expression (and not differential splicing).
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
