# Generate a synthetic csDEX dataset
# TODO: add unit tests
generate <- function(exons=16, conditions=20, interacting=2, replicates=2, genes=2, 
  type="count", data.dir=NULL, seed=NULL, dispersions=NULL){
    
  # Internal functions to manipulate the Beta distribution
  muvarToAlpha <- function(mu, var) {
    # Mean/variance to alpha
    ((1 - mu) / var - 1 / mu) * mu ^ 2
  }
  
  muvarToBeta <- function(mu, var){
    # Mean/variance to beta
    alpha = muvarToAlpha(mu, var)
    alpha * (1 / mu - 1)
  }
  
  muphiToVar <- function(mu, phi){
    # Mean/precision to variance
    mu * (1 - mu) / (1 + phi)
  }
  
  muphiToAlpha <- function(mu, phi){
    # Mean/precision to alpha
    mu * phi
  }
  muphiToBeta <- function(mu, phi){
    # Mean/precision to beta
    phi * (1 - mu)
  }
  
  inv.logit <- function(nu){
    # Inverse logit function
    as.numeric(exp(nu) / (1 + exp(nu)))
  }
  
  # Input check
  stopifnot(type %in% c("PSI", "count"))
  stopifnot(conditions >= 2)
  set.seed(seed)
  
  # Latent parameters
  beta_ex      = runif(exons*genes, -1, 2)
  beta_con     = runif(conditions,  -1, 2)
  
  if(is.null(dispersions)) 
    dispersions  = rgamma(exons*genes, 1, 2)
  precisions   = 1.0 / dispersions
  size.factors = rgamma(conditions, 1, 2)
  
  Mu = matrix(0, exons*genes, conditions)     # Means
  I  = matrix(0, exons*genes, conditions)     # Interactions

  
  # Inferred parameters - PSI
  alphas = matrix(0, exons*genes, conditions)    # Beta shape 1; Alphas
  betas  = matrix(0, exons*genes, conditions)    # Beta shape 2; Betas
  variances = matrix(0, exons*genes, conditions) # Variance is a transformed parameter
  
  # Randomly input interacting pairs ; strong interaction effect
  ix = sample(2:(exons*genes), interacting, replace=TRUE)
  iy = sample(2:conditions, interacting, replace=TRUE)
  for (k in 1:interacting)
    I[ix[k], iy[k]] = sample(c(runif(1, -8, -6), runif(1, 1, 2)), 1)
  
  # Observable parameters in multiple replicates
  # max is used to allow generating one replicate
  X = array(0, c(max(2, replicates), exons * genes, conditions))  # Observables
  for (r in 1:replicates){
    for (i in 1:(exons * genes)){
      for (j in 1:conditions){
        if(type == "count"){
          Mu[i, j]    = exp(size.factors[j] + beta_ex[i] + beta_con[j] + I[i, j])
          X[r, i, j]  = rnbinom(1, mu = Mu[i, j], size=precisions[i])
          
        } else if(type == "PSI"){
          # inv.logit is a squashing function
          phi = precisions[i]
          Mu[i, j] = inv.logit(beta_ex[i] + beta_con[j] + I[i, j])
          variances[i, j] = muphiToVar(Mu[i, j], phi)
          alphas[i, j] = muvarToAlpha(Mu[i, j], variances[i, j])
          betas[i, j] = muvarToBeta(Mu[i, j], variances[i, j])
          X[r, i, j] = rbeta(1, shape1=alphas[i, j], shape2=betas[i, j])
        }
      }}}
  
  # Create metadata and store to disk
  cdx = NULL
  if(!is.null(data.dir)){
    count.dir = file.path(data.dir, "data")
    pars.dir = file.path(data.dir, "parameters")
    
    if(file.exists(data.dir)){
      message(sprintf("Destination %s exists, overwriting", data.dir))
      unlink(data.dir, recursive=TRUE)
    }
    dir.create(data.dir, recursive=TRUE)
    dir.create(count.dir, recursive=TRUE)
    dir.create(pars.dir, recursive=TRUE)
    
    message("Storing data")
    condition.ids = sprintf("cond_%03d", 1:conditions)
    conditions_replicated = c(t(replicate(replicates, condition.ids)))
    replicate.ids = sprintf(paste0(t(replicate(replicates, 
                        paste0(sprintf("cond_%03d_repl_", 1:conditions)))), "%02d"), 1:replicates)
    metadata = data.frame(File.accession=replicate.ids,
                            Experiment.target=conditions_replicated,
                            condition=conditions_replicated)
    write.table(metadata, file.path(data.dir, "metadata.tsv"), 
                sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)
    
    exons_genes = sprintf(paste0(t(replicate(exons, paste0(sprintf("GENE%03d:", 1:genes)))), "%03d"), 1:exons)
    for (r in 1:replicates){
      for (j in 1:conditions){
        cframe = data.frame(groupID=exons_genes, count=X[r,,j])
        fname = sprintf("cond_%03d_repl_%02d.txt", j, r)
        write.table(cframe, file.path(count.dir, fname), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
    }}
    
    message("Storing parameters")
    row.names(I) = exons_genes
    colnames(I) = condition.ids
    ints = melt(I)
    ints = ints[order(ints$value),]
    colnames(ints) = c("featureID", "condition", "interaction")
    row.names(ints) = sprintf("%s:%s", ints$featureID, ints$condition)
    names(dispersions) = exons_genes
    names(precisions) = exons_genes
    
    write.table(as.data.frame(beta_ex), 
                file.path(pars.dir, "beta_ex.tab"), sep="\t", row.names=exons_genes)
    write.table(as.data.frame(beta_con),  
                file.path(pars.dir, "beta_con.tab"), sep="\t", row.names=condition.ids)
    write.table(as.data.frame(dispersions), 
                file.path(pars.dir, "dispersions.tab"), sep="\t", row.names=exons_genes)
    write.table(as.data.frame(precisions), 
                file.path(pars.dir, "precisions.tab"), sep="\t", row.names=exons_genes)
    write.table(as.data.frame(size.factors), 
                file.path(pars.dir, "size_factors.tab"), sep="\t", row.names=condition.ids)
    write.table(ints, file.path(data.dir, "interactions.tab"), 
                sep="\t", row.names=FALSE)
    
    # Construct dataset
    cdx = csDEXdataSet(data.dir=file.path(data.dir, "data"), 
                       design.file=file.path(data.dir, "metadata.tsv"), 
                       type=type)
  }
  
  # Return a constructed synthetic dataset  
  list(data=cdx,
       coefficients=list(
         exons=beta_ex,
         conditions=beta_con,
         interacting=ints),
       dispersions=dispersions,
       precisions=precisions,
       size.factors=size.factors,
       design=metadata,
       data.dir=data.dir)
}