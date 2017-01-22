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
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("csDEXdataSet", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{dataType}:}{Object of class \code{"character"} ~~ }
    \item{\code{exprData}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{rowData}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{colData}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{cpmData}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{colData}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{colData<-}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{cpmData}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{cpmData<-}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{dataType}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{estimateDispersions}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{estimateGeneCPM}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{estimateSizeFactors}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{exprData}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{exprData<-}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{rowData}{\code{signature(obj = "csDEXdataSet")}: ... }
    \item{rowData<-}{\code{signature(obj = "csDEXdataSet")}: ... }
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
