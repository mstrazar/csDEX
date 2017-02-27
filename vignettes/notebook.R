require(pROC)
require(DEXSeq)
require(csDEX)

data.dir = "/Users/martin/Desktop/csdex-temp/"

# Run csDEX
run.csdex <- function(obj, alpha.wald=NULL){
  cdx = obj$data
  cdx = csDEX::estimateSizeFactors(cdx)
  cdx = csDEX::estimatePrecisions(cdx)
  results = csDEX::testForDEU(cdx, workers=3, alpha.wald=alpha.wald)
  row.names(results) = sprintf("%s:%s", results$featureID, results$condition)
  results
}

# Run DEXSeq ; repeat analysis for each case condition
run.dexseq <- function(obj){
  results = data.frame()
  ctrl = "cond_001"
  for(case in unique(obj$design$condition)){
    if (case == ctrl) next;
    message(sprintf("Testing for condition %s", case))
    inxs = obj$design$condition == ctrl 
    inxs = inxs | obj$design$condition == case
    sampleData = obj$design[inxs,]
    countfiles = file.path(obj$data.dir, "data", paste0(sampleData$File.accession, ".txt"))
    
    dex = DEXSeqDataSetFromHTSeq(
      countfiles = countfiles,
      sampleData = sampleData)
    dex = DEXSeq::estimateSizeFactors(dex)
    dex = DEXSeq::estimateDispersions(dex)
    dex = DEXSeq::testForDEU(dex)
    results.dex = DEXSeqResults(dex)
    results.dex$condition = case
    row.names(results.dex) <- sprintf("%s:%s:%s", 
                                      results.dex$groupID,
                                      gsub("E", "", results.dex$featureID), 
                                      case)
    results = rbind(results, results.dex)
  }
  results
}

# Score AUC for given object and results file
score.AUC <- function(obj, results){
  truth = obj$coefficients$interacting[row.names(results), "interaction"]
  score.auc = auc(truth != 0, -log(results$pvalue))
  c(score.auc)
}


# Produce a sample output
repeats = 5
conds = c(3, 5, 10, 20, 50)
results = data.frame()
for (r in 1:repeats){
  for(nc in conds){
    obj = generate(exons=20, conditions=nc, 
                   interacting=20, replicates=2, genes=3, 
                   type="count", data.dir=data.dir)  
    results.dex = run.dexseq(obj)
    results.cdx = run.csdex(obj)
    auc.dex = score.AUC(obj, results.dex)
    auc.cdx = score.AUC(obj, results.cdx)
    
    message(sprintf("r:%d, c:%d, csDEX:%0.3f, DEXSeq: %0.3f", r, nc, auc.cdx, auc.dex))
    df = data.frame(rep=r, conditions=nc, AUC=auc.cdx, method="csDEX")
    results = rbind(results, df)    
    df = data.frame(rep=r, conditions=nc, AUC=auc.dex, method="DEXSeq")
    results = rbind(results, df)    
  }
}


