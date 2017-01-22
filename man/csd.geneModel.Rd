\name{csd.geneModel}
\alias{csd.geneModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs deviance ratio test for each interaction in a gene.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
csd.geneModel(input, min.cpm = NULL, tmp.dir = NULL, dist = "count", alpha.wald = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
}
  \item{min.cpm}{
%%     ~~Describe \code{min.cpm} here~~
}
  \item{tmp.dir}{
%%     ~~Describe \code{tmp.dir} here~~
}
  \item{dist}{
%%     ~~Describe \code{dist} here~~
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
function (input, min.cpm = NULL, tmp.dir = NULL, dist = "count", 
    alpha.wald = NULL) 
{
    suppressPackageStartupMessages({
    })
    write.gene.file <- function(tmp.dir, results, gene) {
        if (!is.null(tmp.dir)) {
            fname = file.path(tmp.dir, sprintf("\%s.csv", gene))
            write.csv(results, fname, quote = TRUE, row.names = FALSE)
        }
    }
    expr = input$expr
    rowdata = input$rowdata
    gene = rowdata$groupID[1]
    results = melt(expr)
    colnames(results) = c("featureID", "condition", "y")
    results = merge(results, rowdata, by = "featureID")
    results$featureID = as.factor(results$featureID)
    results$testable = TRUE
    results[, c("cpm", "pvalue", "padj", "residual", "loglik", 
        "LR", "time", "nrow", "ncol", "msg")] = NA
    if (var(results$y) == 0) {
        results$testable = FALSE
        results$msg = "zero gene variance"
        write.gene.file(tmp.dir, results, gene)
        return(results)
    }
    if (!is.null(min.cpm)) {
        gene.cpm = input$gene.cpm
        results$testable = gene.cpm[as.character(results$condition)] >= 
            min.cpm
        results$cpm = gene.cpm[as.character(results$condition)]
        results[results$cpm < min.cpm, "msg"] = "low cpm"
    }
    N = nrow(results)
    if (dist == "PSI") 
        results$y = ((N - 1) * results$y + 0.5)/N
    if (dist == "count") 
        results$y = round(results$y)
    frm0 = y ~ featureID + condition
    mm0 = model.matrix(frm0, results)
    if (dist == "count") {
        model0 = tryCatch(nbreg.fit(X = mm0, y = as.vector(results$y), 
            phi = results$dispersion, tol = 1e-04), warning = function(w) w, 
            error = function(e) e)
    }
    else if (dist == "PSI") {
        model0 = tryCatch(betareg.fit(x = mm0, y = as.vector(results$y)), 
            warning = function(w) w, error = function(e) e)
    }
    if (class(model0) != "list") {
        results$msg[is.na(results$msg)] = sprintf("null model fit error: \%s", 
            gsub("[\r\n\t]", "", model0))
        results$testable = FALSE
        write.gene.file(tmp.dir, results, gene)
        return(results)
    }
    else {
        start.coefs = model0$coefficients$mean
    }
    results$residual = results$y - model0$fitted.values
    results$interaction = 0
    if (is.null(alpha.wald)) {
        frm1 = y ~ featureID + condition + interaction
        mm1 = model.matrix(frm1, results)
    }
    else {
        wt.stats = csd.waldTest(mm0, model0, alpha = alpha.wald)
        start.coefs = wt.stats$coef.prior
        mm1 = cbind(wt.stats$mm.prior, as.vector(results$interaction))
        colnames(mm1)[ncol(mm1)] = "interaction"
    }
    results$nrow = nrow(mm1)
    results$ncol = ncol(mm1)
    for (i in 1:nrow(results)) {
        if (!results[i, "testable"]) 
            next
        mm1[, "interaction"] = 0
        mm1[i, "interaction"] = 1
        start.time = Sys.time()
        if (dist == "count") {
            coefs.nbinom.init = c(start.coefs, 0)
            model1 = tryCatch(nbreg.fit(X = mm1, y = as.vector(results$y), 
                phi = results$dispersion, tol = 1e-04, beta.init = coefs.nbinom.init, 
                verbose = FALSE), warning = function(w) w, error = function(e) e)
        }
        else if (dist == "PSI") {
            coefs.beta.init = model0$coefficients
            coefs.beta.init$mean = c(start.coefs, 0)
            model1 = tryCatch(betareg.fit(x = mm1, y = as.vector(results$y), 
                control = betareg.control(start = coefs.beta.init)), 
                warning = function(w) w, error = function(e) e)
        }
        if (class(model1) != "list") {
            results$msg[i] = sprintf("fit error: \%s", gsub("[\r\n\t]", 
                "", model1))
            results$testable[i] = FALSE
            next
        }
        results$time[i] = Sys.time() - start.time
        results$loglik[i] = model1$loglik
        LR = 2 * (model1$loglik - model0$loglik)
        results$LR[i] = LR
        if (LR > 0) {
            results$pvalue[i] = 1 - pchisq(LR, 1)
        }
        else {
            results$pvalue[i] = 1
        }
    }
    results$interaction = NULL
    results$padj = p.adjust(results$pvalue, method = "bonferroni")
    results = results[order(results$pvalue), ]
    write.gene.file(tmp.dir, results, gene)
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
