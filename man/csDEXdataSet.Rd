\name{csDEXdataSet}
\alias{csDEXdataSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Object to store a csDEX dataset.
%%  ~~function to do ... ~~
}
\description{
Initialize and instance of the csDEXdataSet object, by constructing a
feature-by-condition data matrix with associated metadata.
}
\usage{
csDEXdataSet(data.dir, design.file, type = "count", col.condition = "Experiment.target", col.replicate = "File.accession", col.testable = "testable", col.read.count = "input.read.count", col.additional = c(), data.file.ext = "txt", aggregation = NULL, min.bin.count = NULL, min.gene.count = NULL, zero.out= NULL))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.dir}{
    The directory containing the replicate (data) files.
}
  \item{design.file}{
    The design file with one line per replicate. Needs to contain at least two
columns, denoting file name (default: "File.accession") and
condition name (default: "Experiment.target"). An optional column can be
included to denote total read counts associated to a replicate (default:
"input.read.count"). Also, an optional binary column (default: "testable") 
can be used to test only selected conditions.
Other columns can be stored, but will not be used by
csDEX.
}
  \item{type}{
    Data type, either "PSI" (percent-spliced in) or "count" (read counts).
}
  \item{col.condition}{
    Name of the column in the design file denoting the experimental condition.
}
  \item{col.replicate}{
    Name of the column in the design file denoting unique file identifiers
(without extension).
}
  \item{col.read.count}{
  Column denotig original input read counts.  
}
  \item{col.additional}{
    Additional columns to import (for use in e.g. complex model designs).
}
  \item{col.testable}{
    Optional column with binary (TRUE/FALSE) information for which conditions to test.
  }
  \item{data.file.ext}{
    Replicate files extension (default: "txt").
}
  \item{aggregation}{
    Replicate aggregation function. A pointer to a function for merging values
of replicates into one per each sample (default: mean).
    
}
  \item{min.bin.count}{
    Minimum value at a feature to be considered non-zero. If defined, values
below this threshold will be set to zero.
}
  \item{min.gene.count}{
    Minimum sum of expression values accross a gene to be considered expressed.
If defined, all the features associated to the gene will be set to zero.
}

\item{zero.out}{
    Use for PSI data to prevent PSI being not equal to either 0 or 1 for constitutive features. 
    If feature value is less than min.count, it is set to NA and ignored by the aggregation function.
}

}
\details{
    The main task of this function is to construct a csDEXdataSet object by
aggregating the replicate data into one column per sample. Arbitrary
aggregation functions can be passed, the most common ones being mean, sum and
max. The options min.bin.count, min.gene.count and zero.out can be used for filtering.

All metadata provided by the design file is stored in the object colData. The
two mandatory columns are specified by col.condition and col.replicate, as
defined above. The total experiment read count can be provided along with each
replicate, to enable usage of estiamteGeneCPM function in the downstream
analysis.

The row metadata is parsed as featureID (individual regions of interest) and groupID (gene).

}
\value{
An initialized csDEXdataSet object.
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Martin Stražar.
}
\note{
%%  ~~further notes~~
}


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
