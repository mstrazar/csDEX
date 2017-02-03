\name{csDEXdataSet-class}
\Rdversion{1.1}
\docType{class}
\alias{csDEXdataSet-class}
\alias{colData,csDEXdataSet-method}
\alias{colData<-,csDEXdataSet-method}
\alias{cpmData,csDEXdataSet-method}
\alias{cpmData<-,csDEXdataSet-method}
\alias{dataType,csDEXdataSet-method}
\alias{estimateDispersions,csDEXdataSet-method}
\alias{estimateGeneCPM,csDEXdataSet-method}
\alias{estimateSizeFactors,csDEXdataSet-method}
\alias{exprData,csDEXdataSet-method}
\alias{exprData<-,csDEXdataSet-method}
\alias{rowData,csDEXdataSet-method}
\alias{rowData<-,csDEXdataSet-method}

\title{Class \code{"csDEXdataSet"}}
\description{
A data structure with the features-by-condition matrix, storing the expression
values for percentage-spliced in or read counts, along with metadata for both
rows and columns. Optionally, gene expression is stored as counts-per-million
reads (CPM).
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("csDEXdataSet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{dataType}:}{Object of class \code{"character"} Either "PSI" (percentage-spliced in) or "count" (read counts), determining the underlying data distribution.}
    \item{\code{exprData}:}{Object of class \code{"matrix"} Feature-by-condition matrix containing expression data. }
    \item{\code{rowData}:}{Object of class \code{"data.frame"} Row (feature) metadata. }
    \item{\code{colData}:}{Object of class \code{"data.frame"} Column (condition) metadata, the experiment design. }
    \item{\code{cpmData}:}{Object of class \code{"matrix"} Optional. Counts-per-million reads (CPM) associated to a group (gene).}
  }
}
\section{Methods}{
  \describe{

    \item{dataType}{\code{signature(obj = "csDEXdataSet")}: The accessor
(getter) of the \code{csDEXdataSet} object data type (PSI or count). }

    \item{colData}{\code{signature(obj = "csDEXdataSet")}: The accessor
(getter) of column (condition) metadata provided by the experiment design file.}
    \item{colData<-}{\code{signature(obj = "csDEXdataSet")}: The accessor
(setter) of column (condition) metadata provided by the experiment design file.}
    
    \item{cpmData}{\code{signature(obj = "csDEXdataSet")}: The accessor
(getter) of the counts per million reads (CPM) associated to each gene
(\code{groupID}). Computable if the \code{input.read.count} field is set and provided in the experiment design file. }

    \item{cpmData<-}{\code{signature(obj = "csDEXdataSet")}: The accessor (setter) of the counts per million reads (CPM) associated to each gene (\code{groupID}).}
    
    \item{exprData}{\code{signature(obj = "csDEXdataSet")}: The accessor (getter) of the expression data matrix constructed by merging the experimental data in samples. }

    \item{exprData<-}{\code{signature(obj = "csDEXdataSet")}: The accessor (setter) of the expression data matrix constructed by merging the experimental data in samples.  }
    
    \item{rowData}{\code{signature(obj = "csDEXdataSet")}: The accessor
(getter) of row (feature) metadata provided by the experiment data files. }
    \item{rowData<-}{\code{signature(obj = "csDEXdataSet")}: The accessor
(setter) of row (feature) metadata provided by the experiment data files. }
    
    
    \item{estimatePrecisions}{\code{signature(obj = "csDEXdataSet")}: Estimate
the precisions hyperparameter for each feature by using the external
edgeR::estimateDisp function. Required for fitting read-count based models. The
function estimateSizeFactors should be called on the cdx object prior to
estimatePrecisions call. }
    
    \item{estimateGeneCPM}{\code{signature(obj = "csDEXdataSet")}:   Estimate
counts per million reads (CPM) associated to each gene (groupID). An integer
field named input.read.count is required in the design file, specifying the
sequencing depth measured by the number of sequenced reads / read pairs. Can be
used by both PSI and count models to skip testing for low-expressed genes and
avoid Type I errors that are dut to changes in gene expression (and not
differential splicing). }
    
    \item{estimateSizeFactors}{\code{signature(obj = "csDEXdataSet")}: Estimate
the size factors hyperparameter for each feature by using the external
DESeq2::estimateSizeFactorsForMatrix function. Required for normalization of
samples with respect to sequencing depth and mandatory when fitting count
models. }
    
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("csDEXdataSet")
}
\keyword{classes}
