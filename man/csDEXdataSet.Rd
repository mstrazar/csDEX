\name{csDEXdataSet}
\alias{csDEXdataSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Object to store a csDEX dataset.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
csDEXdataSet(data.dir, design.file, type = "count", col.condition = "Experiment.target", col.replicate = "File.accession", data.file.ext = "txt", aggregation = NULL, min.bin.count = NULL, min.gene.count = NULL, zero.out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.dir}{
%%     ~~Describe \code{data.dir} here~~
}
  \item{design.file}{
%%     ~~Describe \code{design.file} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{col.condition}{
%%     ~~Describe \code{col.condition} here~~
}
  \item{col.replicate}{
%%     ~~Describe \code{col.replicate} here~~
}
  \item{data.file.ext}{
%%     ~~Describe \code{data.file.ext} here~~
}
  \item{aggregation}{
%%     ~~Describe \code{aggregation} here~~
}
  \item{min.bin.count}{
%%     ~~Describe \code{min.bin.count} here~~
}
  \item{min.gene.count}{
%%     ~~Describe \code{min.gene.count} here~~
}
  \item{zero.out}{
%%     ~~Describe \code{zero.out} here~~
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
function (data.dir, design.file, type = "count", col.condition = "Experiment.target", 
    col.replicate = "File.accession", data.file.ext = "txt", 
    aggregation = NULL, min.bin.count = NULL, min.gene.count = NULL, 
    zero.out = NULL) 
{
    if (!(type \%in\% c("count", "PSI"))) 
        stop("type must be one of ('count', 'PSI')")
    if (is.null(min.bin.count)) 
        min.bin.count = 0
    if (type == "PSI") {
        if (is.null(aggregation)) 
            aggregation = mean
        if (is.null(zero.out)) 
            zero.out = TRUE
        if (is.null(min.gene.count)) 
            min.gene.count = 1
    }
    else if (type == "count") {
        if (is.null(aggregation)) 
            aggregation = sum
        if (is.null(zero.out)) 
            zero.out = FALSE
        if (is.null(min.gene.count)) 
            min.gene.count = 0
    }
    zeroOut <- function(repData, min.count = 0) {
        genes = unlist(lapply(row.names(repData), function(x) strsplit(x, 
            ":")[[1]][1]))
        for (g in unique(genes)) {
            inxs = genes == g
            zeros = which(apply(repData[inxs, ], 2, sum) < min.count)
            repData[inxs, zeros] = NA
        }
        return(repData)
    }
    design = read.csv(design.file, sep = "\t", header = TRUE)
    stopifnot(col.condition \%in\% colnames(design))
    stopifnot(col.replicate \%in\% colnames(design))
    conditions = sort(unique(design[, col.condition]))
    n.con = length(conditions)
    message("Processing expression data")
    exprData = NULL
    lib.sizes = NULL
    for (i in 1:length(conditions)) {
        cond = conditions[i]
        message(sprintf("Condition \%s", cond))
        replicates = design[design[, col.condition] == cond, 
            col.replicate]
        cond.lib.size = NULL
        if (!is.null(design$input.read.count)) {
            cond.lib.size = aggregation(design[design[, col.condition] == 
                cond, "input.read.count"])
            lib.sizes = c(lib.sizes, cond.lib.size)
        }
        repData = NULL
        n.rep = length(replicates)
        for (j in 1:length(replicates)) {
            rep = replicates[j]
            rep.path = file.path(data.dir, paste(rep, data.file.ext, 
                sep = "."))
            y = read.table(rep.path, header = FALSE, comment.char = "_")
            y$V1 = as.character(y$V1)
            n.row = nrow(y)
            message(sprintf("     Replicate \%s, num. rows: \%d", 
                rep.path, n.row))
            if (is.null(repData)) {
                repData = matrix(0, ncol = n.rep, nrow = n.row)
                row.names(repData) = y$V1
            }
            repData[y$V1, j] = y$V2
        }
        repData[repData < min.bin.count] = 0
        if (zero.out || min.gene.count > 0) 
            repData = zeroOut(repData, min.gene.count)
        repVec = suppressWarnings(apply(repData, 1, aggregation, 
            na.rm = TRUE))
        repVec[is.infinite(repVec)] = 0
        repVec[is.na(repVec)] = 0
        if (is.null(exprData)) {
            exprData = matrix(0, ncol = n.con, nrow = n.row)
            row.names(exprData) = row.names(repData)
            colnames(exprData) = conditions
        }
        exprData[, i] = repVec
    }
    rowData = data.frame(featureID = row.names(exprData), groupID = unlist(lapply(row.names(repData), 
        function(x) strsplit(x, ":")[[1]][1])), binID = unlist(lapply(row.names(repData), 
        function(x) strsplit(x, ":")[[1]][2])))
    rowData$featureID = as.character(rowData$featureID)
    rowData$groupID = as.character(rowData$groupID)
    rowData$binID = as.character(rowData$binID)
    row.names(rowData) = rowData$featureID
    colData = data.frame(condition = colnames(exprData))
    if (!is.null(lib.sizes)) 
        colData$lib.size = lib.sizes
    new("csDEXdataSet", exprData = exprData, rowData = rowData, 
        colData = colData, dataType = type)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
