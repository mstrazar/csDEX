### Class definition
setClass("csDEXdataSet",
    slots=c(
        dataType="character",   # Type of the dataset (count/PSI).
        exprData="matrix",      # alias for either counts or PSI.
        rowData="data.frame",   # bin names, precision, etc.
        colData="data.frame",   # design files on conditions
        cpmData="matrix"        # computed cpm per group (gene)
    ))


### Getters
setGeneric("exprData", function(obj, normalized=FALSE) standardGeneric("exprData"))
setMethod("exprData", "csDEXdataSet",
    function(obj, normalized=FALSE){
        if(normalized){
            if(is.null(csDEX::colData(obj)$size.factor)){
                message("Size factors not provided yet. Run estimateSizeFactors first.\n")
                NULL
            } else {
                sweep(obj@exprData, 2, colData(obj)$size.factor, "/")
            }
        } else {
            obj@exprData
        }})

setGeneric("rowData", function(obj) standardGeneric("rowData"))
setMethod("rowData", "csDEXdataSet", function(obj) obj@rowData)

setGeneric("dataType", function(obj) standardGeneric("dataType"))
setMethod("dataType", "csDEXdataSet", function(obj) obj@dataType)

setGeneric("colData", function(obj) standardGeneric("colData"))
setMethod("colData", "csDEXdataSet", function(obj) obj@colData)

setGeneric("cpmData", function(obj) standardGeneric("cpmData"))
setMethod("cpmData", "csDEXdataSet", function(obj) obj@cpmData)

### Setters
setGeneric("exprData<-", function(obj, value) standardGeneric("exprData<-"))
setReplaceMethod("exprData", "csDEXdataSet",
    function(obj, value) {obj@exprData <- value; validObject(obj); obj})

setGeneric("rowData<-", function(obj, value) standardGeneric("rowData<-"))
setReplaceMethod("rowData", "csDEXdataSet",
    function(obj, value) {obj@rowData <- value; validObject(obj); obj})

setGeneric("colData<-", function(obj, value) standardGeneric("colData<-"))
setReplaceMethod("colData", "csDEXdataSet",
    function(obj, value) {obj@colData <- value; validObject(obj); obj})

setGeneric("cpmData<-", function(obj, value) standardGeneric("cpmData<-"))
setReplaceMethod("cpmData", "csDEXdataSet",
    function(obj, value) {obj@cpmData <- value; validObject(obj); obj})


### Data generating functions

# Return a csDEXdataSet with estimated size Factors
setGeneric("estimateSizeFactors",
    function(obj) standardGeneric("estimateSizeFactors"))
setMethod("estimateSizeFactors", "csDEXdataSet",
    function(obj) {
        csDEX::colData(obj)$size.factor = 1.0 / c(edgeR::calcNormFactors(obj@exprData))
        obj
        })

# Return a csDEXdataSet with estimated precisions
setGeneric("estimatePrecisions",
    function(obj) standardGeneric("estimatePrecisions"))
setMethod("estimatePrecisions", "csDEXdataSet",
    function(obj) {        
        expr = csDEX::exprData(obj, normalized = TRUE)
        csDEX::rowData(obj)$precision = 1.0 / edgeR::estimateDisp(expr)$tagwise.dispersion
        obj
        })

# Estimate *maximum* CPM in bins inside a gene
# Required for filtering untestable genes in particular conditions
# Implemented as a separate function to avoid unnecesary increase of memory
setGeneric("estimateGeneCPM",
    function(obj) standardGeneric("estimateGeneCPM"))
setMethod("estimateGeneCPM", "csDEXdataSet",
    function(obj) {
         if(is.null(colData(obj)$lib.size)){
            message("Design file must provide input.read.count (library size) per sample.\n")
            NULL
        } else {
            expr = csDEX::exprData(obj)
            lib.size = csDEX::colData(obj)$lib.size
            cpms = edgeR::cpm(expr, lib.size=lib.size)

            groups = unique(csDEX::rowData(obj)$groupID)
            cpm.max = matrix(0, nrow=length(groups), ncol=ncol(expr))
            colnames(cpm.max) = colnames(expr)
            row.names(cpm.max) = groups

            for (i in 1:length(groups))
                cpm.max[i, ] = apply(cpms[csDEX::rowData(obj)$groupID == groups[i],], 2, max)

            csDEX::cpmData(obj) = cpm.max
            obj}
        })


# Class should define way to handle replicates, etc.
# Enable passing a custom aggregation function, such as geom. mean, simple mean, sum, etc.
csDEXdataSet <- function(data.dir, design.file, type="count",
        col.condition="Experiment.target",
        col.replicate="File.accession",
        col.testable="testable",
        input.read.count="input.read.count",
        data.file.ext="txt",
        aggregation=NULL,
        min.bin.count=NULL,
        min.gene.count=NULL,
        zero.out=NULL){

    ### Input checks ###
    if(!(type %in% c("count", "PSI")))
      stop("type must be one of ('count', 'PSI')")
    
    ### Set recommended parameters
    if(is.null(min.bin.count)) min.bin.count = 0
    if(type == "PSI"){
      if(is.null(aggregation)) aggregation = mean
      if(is.null(zero.out)) zero.out = TRUE
      if(is.null(min.gene.count)) min.gene.count = 1
    } else if (type == "count"){
      if(is.null(aggregation)) aggregation = sum
      if(is.null(zero.out)) zero.out = FALSE
      if(is.null(min.gene.count)) min.gene.count = 0
    }
    
    ### Private functions ###
    zeroOut <- function(repData, min.count=0){
        # Zero-out non-expressed genes in individual replicates
        genes = unlist(lapply(row.names(repData),
                function(x) strsplit(x, ":")[[1]][1]))

        for(g in unique(genes)){
            inxs = genes == g
            zeros = which(apply(repData[inxs, ], 2, sum) < min.count)
            repData[inxs, zeros] = NA
        }
            return(repData)
    }

    # Parse design file and ensure columns are present
    design = read.csv(design.file, sep="\t", header=TRUE)
    stopifnot(col.condition %in% colnames(design))
    stopifnot(col.replicate %in% colnames(design))

    # Parse condition but retain order of rowss
    conditions = sort(unique(design[,col.condition]))
    n.con = length(conditions)

    # Parse expression data
    message("Processing expression data")

    exprData = NULL
    lib.sizes = NULL
    testable = rep(TRUE, n.con)
    
    for (i in 1:length(conditions)){
        cond = conditions[i]
        message(sprintf("Condition %s", cond))

        replicates = design[design[,col.condition] == cond, col.replicate]

        # Merge data on read counts if available
        # Aggregated in the same manner than counts
        cond.lib.size = NULL
        if (input.read.count %in% colnames(design)){
            cond.lib.size = aggregation(design[design[,col.condition] == cond, input.read.count])
            lib.sizes = c(lib.sizes, cond.lib.size)
        }
        
        # Information on testability if available.
        if (col.testable %in% colnames(design)){
          testable[i] = as.logical(sum(design[design[,col.condition] == cond, col.testable]))
        }
        
        # Merge data on replicates
        repData = NULL
        n.rep = length(replicates)

        for (j in 1:length(replicates)){
            rep = replicates[j]
            rep.path = file.path(data.dir, paste(rep, data.file.ext, sep="."))
            y = read.table(rep.path, header=FALSE, comment.char="_")
            y$V1 = as.character(y$V1)
            n.row = nrow(y)
            message(sprintf("     Replicate %s, num. rows: %d", rep.path, n.row))

            if(is.null(repData)){
                repData = matrix(0, ncol=n.rep, nrow=n.row)
                row.names(repData) = y$V1
            }
            repData[y$V1, j] = y$V2
        }

        # Zero out on individual bin level
        repData[repData < min.bin.count] = 0

        # Zero out on a gene level and merge replicates
        if(zero.out || min.gene.count > 0)
            repData = zeroOut(repData, min.gene.count)

        repVec = suppressWarnings(apply(repData, 1, aggregation, na.rm=TRUE))
        repVec[is.infinite(repVec)] = 0
        repVec[is.na(repVec)] = 0

        if(is.null(exprData)){
            exprData = matrix(0, ncol=n.con, nrow=n.row)
            row.names(exprData) = row.names(repData)
            colnames(exprData) = conditions
        }
        exprData[,i] = repVec
    }

    # Order of passed data files must match order of rows in design file
    rowData = data.frame(
            featureID = row.names(exprData),                # gene: exon
            groupID = unlist(lapply(row.names(repData),     # gene
                function(x) strsplit(x, ":")[[1]][1])),
            binID = unlist(lapply(row.names(repData),       # exon
                function(x) strsplit(x, ":")[[1]][2]))
            )
    rowData$featureID = as.character(rowData$featureID)
    rowData$groupID = as.character(rowData$groupID)
    rowData$binID = as.character(rowData$binID)
    row.names(rowData) = rowData$featureID

    colData = data.frame(condition=colnames(exprData))
    row.names(colData) <- colData$condition
    if(!is.null(lib.sizes)) colData$lib.size = lib.sizes
    colData$testable = testable

    new("csDEXdataSet", exprData=exprData, rowData=rowData, colData=colData, dataType=type)
}

waldTest <- function(mm0, model0, alpha=0.05){
  # Wald test for Beta model coefficients
  # Reduce non-significant coefficients to zero for efficient computation
  beta0 = model0$coefficients$mean
  
  # Chi2 test for coefficients greater than zero
  scores = diag(model0$vcov)[1:length(beta0)]
  wald = rep(0, length(scores))
  names(wald) = names(scores)
  
  # Test statistic
  for (i in 1:length(wald)){
    w = aod::wald.test(b=beta0[i], scores[i], T=1)
    wald[i] = w$result$chi2["P"]
  }
  
  # (Non-)significant coefficients
  sig = names(beta0)[wald < alpha]
  not.sig = names(beta0)[wald >= alpha]
  
  # Merge prior coefficients for correct computation of likelihood
  # New model matrix yields same likelihood
  if(length(not.sig) > 1){
    mu.prior = mm0[,not.sig] %*% beta0[not.sig] 
    colnames(mu.prior) = "mu.prior" 
    mm.prior = cbind(mu.prior, mm0[,sig])  
    coef.prior = c(1, beta0[sig])
  } else {
    mm.prior = mm0
    coef.prior = beta0
  } 
  return(list(sig=sig, not.sig=not.sig, mm.prior=mm.prior, coef.prior=coef.prior))
}


geneModel <- function(input, min.cpm=NULL, tmp.dir=NULL, dist="count", alpha.wald=NULL){
  
    # Writing to file
    write.gene.file <- function(tmp.dir, results, gene){
        if(!is.null(tmp.dir)){
          fname = file.path(tmp.dir, sprintf("%s.csv", gene))
          write.csv(results, fname, quote=TRUE, row.names=FALSE)}
    }

    # Extract input values
    expr = input$expr
    rowdata = input$rowdata
    coldata = input$coldata
    gene = rowdata$groupID[1]

    # Prepare a results data frame
    results = reshape::melt(expr)
    colnames(results) = c("featureID", "condition", "y")
    results = merge(results, rowdata, by="featureID")
    results$featureID = as.factor(results$featureID)
    results$testable = coldata[results$condition, "testable"]
    results[, c("cpm", "pvalue", "fdr", "fdr.all", "padj", "residual", "loglik", "ddev", "time", "nrow", "ncol", "msg")] = NA

    # Check variance
    if(var(results$y) == 0){
      results$testable = FALSE
      results$msg = "zero gene variance"
      write.gene.file(tmp.dir, results, gene)
      return(results)}
    
    # Define testability based on min.cpm
    if(!is.null(min.cpm)){
        gene.cpm = input$gene.cpm
        results$testable = results$testable & (gene.cpm[as.character(results$condition)] >= min.cpm)
        results$cpm = gene.cpm[as.character(results$condition)]
        results[results$cpm < min.cpm, "msg"] = "low cpm"
    }

    # Correct values to (0, 1) for the beta model
    N = nrow(results)
    if(dist == "PSI")
        results$y = ((N - 1) * results$y + 0.5) / N
    if(dist == "count")
        results$y = round(results$y)

    # Null hypothesis
    frm0 = y ~ featureID + condition
    mm0  = model.matrix(frm0, results)

    # Fit null model and store coefficients for better initialization
    if(dist == "count"){
        model0 = tryCatch(csDEX::nbreg.fit(X=mm0, 
                                   y=as.vector(results$y),
                                   phi=results$precision, 
                                   tol=0.0001),
                  warning = function(w) w,
                  error = function(e) e)
    } else if(dist == "PSI"){
        model0 = tryCatch(
                betareg::betareg.fit(x=mm0, y=as.vector(results$y)),
                warning = function(w) w,
                error = function(e) e)        
    }
    if(class(model0) != "list") {
      results$msg[is.na(results$msg)] = sprintf("null model fit error: %s", 
                                                gsub("[\r\n\t]", "", model0))
      results$testable = FALSE
      write.gene.file(tmp.dir, results, gene)
      return(results)
    } else {
      start.coefs = model0$coefficients$mean
    }
    
    # Filter testable values according to fit residual
    results$residual = results$y - model0$fitted.values
    
    # Compute Wald test statistics if alpha.wald
    results$interaction = 0
    if(is.null(alpha.wald)){
      # Full alternative model for gene
      frm1 = y ~ featureID + condition + interaction
      mm1 = model.matrix(frm1, results)
    } else {
      # Zero coefficients merged into one column
      wt.stats = csDEX::waldTest(mm0, model0, alpha=alpha.wald)
      start.coefs = wt.stats$coef.prior
      mm1 = cbind(wt.stats$mm.prior, as.vector(results$interaction))
      colnames(mm1)[ncol(mm1)] = "interaction"
    }
    
    # Store design matrix dimensions
    results$nrow = nrow(mm1)
    results$ncol = ncol(mm1)
    
    # Compute deviance and results
    for (i in 1:nrow(results)){
        if(!results[i, "testable"]) next
        mm1[, "interaction"] = 0
        mm1[i, "interaction"] = 1

        start.time = Sys.time()
        if(dist == "count"){
            coefs.nbinom.init = c(start.coefs, 0)
            model1 = tryCatch(csDEX::nbreg.fit(X = mm1,
                                       y=as.vector(results$y),
                                       phi = results$precision, 
                                       tol=0.0001,
                                       beta.init = coefs.nbinom.init,
                                       verbose=FALSE),
                warning = function(w) w,
                error = function(e) e)
        } else if(dist == "PSI"){
            coefs.beta.init = model0$coefficients
            coefs.beta.init$mean = c(start.coefs, 0)
            model1 = tryCatch(
                betareg::betareg.fit(x=mm1, y=as.vector(results$y),
                    control=betareg::betareg.control(start=coefs.beta.init)),
                warning = function(w) w,
                error = function(e) e)
        }
        if(class(model1) != "list") {
          results$msg[i] = sprintf("fit error: %s", 
                                   gsub("[\r\n\t]", "", model1))
          results$testable[i] = FALSE
          next
        }
          
        # Compute difference of deviance (likelihood ratio) test statistic
        results$time[i] = Sys.time() - start.time
        results$loglik[i] = model1$loglik

        # Likelihood ratio-test pvalue
        ddev = 2 * (model1$loglik - model0$loglik)
        results$ddev[i] = ddev
        if (ddev > 0){
            results$pvalue[i] = 1.0 - pchisq(ddev, 1)
        } else {
            results$pvalue[i] = 1
        }
    }

    # Correct for multiple testing
    results$interaction = NULL
    results$padj = p.adjust(results$pvalue, method="bonferroni")
    results$fdr = p.adjust(results$pvalue, method="BH")
    results = results[order(results$pvalue),]

    # Write intermediate results
    write.gene.file(tmp.dir, results, gene)
    return(results)
}


testForDEU <- function(cdx, workers=1, tmp.dir=NULL, min.cpm=NULL, alpha.wald=NULL){
    # Test for differential exon usage given a csDEX data object file,
    # with calculated precisions and size factors. Only test genes
    #   if they have any reads mapped onto them

    gene.set = unique(csDEX::rowData(cdx)$groupID)
    inputs = list()
    dist = csDEX::dataType(cdx)
    
    # Check input for NBINOM
    if(dist == "count"){
        stopifnot(!is.null(csDEX::colData(cdx)$size.factor))
        stopifnot(!is.null(csDEX::rowData(cdx)$precision))
    }

    # Add to queue if gene has not equal count at all cells
    gene.cpm = NULL
    message("Constructing inputs ... \n")
    for (g in gene.set){
        inxs = which(csDEX::rowData(cdx)$groupID == g)
        gene.expr = csDEX::exprData(cdx, normalized=(dist=="count"))[inxs,]
        gene.rowdata = csDEX::rowData(cdx)[inxs,]
        if(!is.null(min.cpm)) gene.cpm = csDEX::cpmData(cdx)[g,]

        inputs[[length(inputs)+1]] = list(expr=gene.expr, 
                                          rowdata=gene.rowdata, 
                                          coldata=csDEX::colData(cdx),
                                          gene.cpm=gene.cpm)
    }
    
    # If workers == 1, run within the same process
    if(workers > 1){
      message(sprintf("Dispatching on %d workers, cache directory %s. \n", workers, tmp.dir))
      jobs = mclapply(inputs, csDEX::geneModel, min.cpm, tmp.dir, dist, alpha.wald,
          mc.preschedule=FALSE, mc.cores=workers, mc.silent=TRUE)
    } else {
      jobs = list()
      for(inp in inputs) 
        jobs[[length(jobs)+1]] = csDEX::geneModel(inp, min.cpm, tmp.dir, dist, alpha.wald)
    }

    results = data.frame()
    for (i in 1:length(jobs)){
        if (!class(jobs[[i]]) == "try-error")
            results = rbind(results, jobs[[i]])
    }
    
    results$fdr.all = p.adjust(results$pvalue, "BH")
    results=results[order(results$pvalue),]
    return(results)
}