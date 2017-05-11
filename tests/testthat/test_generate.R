context("Test csDEXdataSet generate")

# TODO: solve problem in this correlated 
# variance computation
tmp.dir = "./temp"
dir.create(tmp.dir, recursive=TRUE)
for (typ in c("count", "PSI")){
    n.exons = 16
    data = generate(exons=n.exons, conditions=20,
                    interacting=0, replicates=2, genes=2, 
                    type=typ, data.dir=tmp.dir, 
                    seed=NULL, dispersions=NULL)
    
    # Exon factors are independent for this purpose
    cdx = data$data
    if(typ == "count"){
      cdx = estimateSizeFactors(cdx)
      cdx = estimatePrecisions(cdx)
    }
    results = testForDEU(cdx)
    expr = exprData(cdx)
    vars = apply(expr, 1, var)
    num.sig = sum(results$pvalue < 0.05) / nrow(results)
    message(sprintf("Number of cases p<0.05: %.4f", num.sig))
    
    # Compute minimum p-value in relation to perceived variance
    min.pvalue = aggregate(results$pvalue, by=list(featureID=results$featureID),  FUN=min)
    y = min.pvalue$x
    x = vars[min.pvalue$featureID]
    cr = cor(x, y, method="spearman")
    message(sprintf("Correlation between variance and min. p-value: %.4f", cr))
}
unlink(tmp.dir, recursive = TRUE)