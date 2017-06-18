
# The csDEX package: Quick Start

  
csDEX (Condition-specific Differential Exon Expression) is a family of general
linear models (GLMs), used to model differential splicing given a set of
sequencing experiments based on RNA-seq. The methods assume tens to hundreds
experimental conditions to be compared, and the result are
changes in usage of splice sites that are specific to an
experimental condition. 

Capabilities:
* handling of sequence data in form of read counts or percentages (e.g. percent-spliced in, PSI),
* modelling all type of splicing events,
* discovery and ranking of condition-specific changes,
* parallel processing of multiple genes (groups),
* time-efficient model approximation with no change in false positive rate.

## Vignettes
* The vignete on general usage of csDEX - this page 
[(Rmd)](vignettes/csdex_vignette.Rmd).
* Comparison of csDEX and DEXSeq models on generated data [(Rmd)](vignettes/prediction.Rmd)[(html)](http://htmlpreview.github.com/?https://github.com/mstrazar/csDEX/blob/master/vignettes/prediction.nb.html).
* Comparison of different dispersion estimation methods [(Rmd)](vignettes/dispersion.Rmd)[(html)](http://htmlpreview.github.com/?https://github.com/mstrazar/csDEX/blob/master/vignettes/dispersion.nb.html).


## Installation

Install the latest version hosted on GitHub:

```r
require(devtools)
install_github("mstrazar/csDEX")
```

Test the installation:

```r
require(csDEX)
```

## Data preparation

The package assumes a directory with replicate data and a design (decoder) file
to link experimental conditions to appropriate files.


### Replicate data
 
Replicate data is stored as `.txt` files containing a table with two columns in the
following format:
```
    groupID:featureID    value
    groupID:featureID    value
```

The replicate files are tab- or space-delimited and have no header line. The
fields encode the following information:
* `groupID`: group of features belonging to the same transcriptional unit,
    most likely a gene.    
* `featureID`: alternative transcriptional unit, such as an exon, exonic
    part, alternative 5'/3' end, etc. Most commonly referring to a
    non-overlapping exonic bin, to which reads are mapped.
* `value`: expression quantification, either in form of read counts mapping
    to the feature (non-negative integers) or percentage spliced-in (PSI,
    real numbers 0.0-1.0, inclusive). Choice of quantification shall be
    consistent within and between all files.

Example for PSI:

```
    ENSG00000198912:003 0.13
    ENSG00000198912:004 1.00
    ENSG00000198912:005 0.02
    ENSG00000198912:006 1.00
    ENSG00000198912:007 0.97
    ENSG00000198912:008 0.13
```


An arbitrary number of replicates can be used for each of the experimental conditions.


#### Obtaining replicate files from BAM/SAM files. 

Given the number of available tools out there, we do not reinvent the wheel but
point users to existing tools. These will typically assume a genome annotation
listed in a `.gtf/.gff` format.

Read count data:
- QoRTs
- DEXSeq

PSI:
- cuffLinks
- rMATS
- manually from ENCODE pre-made transcript quantifications


### The design file.

The experimental design is stored as a tab-delimited table and links replicate filenames to
experimental conditions. It assumes (at least) two columns denoting
* `File.accession`: replicate file name within the input directory, without extension.
* `Experiment.target`: experimental condition identifier
* `input.read.count`: (optional) Number of reads from which the replicate was generated. Required to compute gene counts-per-million reads (CPM) if required.

```
    File.accession Experiment.target input.read.count
    ENCFF120ICS       DDX27-human          1380110
    ENCFF237IIP       DDX27-human          1593308
    ENCFF202WCG       NELFE-human          1472191
    ENCFF759XIJ       NELFE-human          2889128
    ENCFF781OAE        PUM2-human          2021596
    ENCFF803RGM        PUM2-human          2362944
```

Files sharing the same `Experiment.target` are assumed to be technical replicates
of the same experimental condition. The column names are not fixed and can be defined at
runtime. Other columns can be stored, but will be ignored by csDEX. The format
is consistent with the ENCODE metadata format.


Two exemplary small datasets are provided with the csDEX package, to be
used for testing purposes only.



## Using csDEX

The example dataset is accessed as follows.

```r
design.file = system.file("extdata/metadata.tsv", 
    package="csDEX", mustWork=TRUE)
data.dir.psi = system.file("extdata/psi", 
    package="csDEX", mustWork=TRUE)
data.dir.count = system.file("extdata/count", 
    package="csDEX", mustWork=TRUE)
```

### The PSI model

The percentage-spliced in model is based on the Beta regression GLM. It is somewhat
more streamlined and therefore presented first.  The data is loaded into a
`csDEXdataSet` object:


```r
cdx = csDEXdataSet(data.dir=data.dir.psi, 
                   design.file=design.file, 
                   type="PSI",
                   aggregation=mean)
```

The `aggregation` argument is a function determining how to merge replicate
measurements associated to a feature.  Other functions, e.g. `sum`, `max`,
`min`, ..., are possible, as well as custom ones. The metadata associated to
rows and columns can be viewed or modified using `rowData(cdx)` and
`colData(cdx)` extractor functions. See the methods documentation for a
complete list of preprocessing and filtering options.


The differential feature expression is run with:


```r
result = csd.testForDEU(cdx, 
    workers=1, 
    tmp.dir=NULL)
```


```r
head(result[c("featureID", "condition", "y", "pvalue", "padj")])
```

```
##               featureID   condition         y       pvalue         padj
## 145 ENSG00000069812:003 PTBP1-human 0.6927059 1.189677e-08 1.011226e-06
## 237 ENSG00000180758:005 DDX27-human 0.4706000 5.813390e-05 2.906695e-03
## 110 ENSG00000158109:001 DDX27-human 0.9640000 9.958324e-05 2.987497e-03
## 161 ENSG00000158109:004 DDX27-human 0.9640000 9.958324e-05 2.987497e-03
## 574 ENSG00000178585:012 DDX55-human 0.4507692 1.082377e-04 7.035453e-03
## 331 ENSG00000198912:007 NELFE-human 0.9485000 2.939775e-04 1.175910e-02
```

The result is a table of features and conditions ranked by statistical
significance, where the columns `pvalue`, `padj` (Bonferroni correction) and
`fdr` (Benjamini-Hochberg correction) are of interest. The argument `tmp.dir`
can be set to a valid system path where intermediary results will be stored
(one file per gene), which is useful for lengthy computations. The `workers`
parameter can be used to increase the number of processing cores, if available.


### The read count model

The read count model is based on Negative binomial regression. It requires 2-3
extra preprocessing steps to infer model hyperparameters. 


Data is loaded similarly as for the PSI model. Note the `type` argument.


```r
cdx = csDEXdataSet(data.dir=data.dir.count, 
    design.file=design.file, 
    type="count")
```

The following methods act on a `csDEXdataSet` object, by adding new fields to `rowData` or `colData`.

Due to possible unequal sequencing, size factors are computed.


```r
cdx = estimateSizeFactors(cdx)
```

Normalized count values can be accessed by setting the `normalized` parameter
to `exprData`.


```r
exprData(cdx, normalized=TRUE)
```

The NB model assumes the dispersion parameter (sometimes termed precision) to be known.
csDEX uses the `edgeR` package method `estimateDisp(...)$tagwise.dispersion`, but custom
dispersion models can be used.


```r
cdx = estimateDispersions(cdx)
```

```
## Design matrix not provided. Switch to the classic mode.
```

Finally, gene-wise read coverage, expressed as counts-per-million reads (CPM)
can be calculated, provided that a column name `input.read.count` is present
and valid in the design file. This step is optional, but can be used to 
filter out low coverage genes.


```r
cdx = estimateGeneCPM(cdx)
cpmData(cdx)
```

```
##                 DDX27-human DDX55-human NELFE-human PTBP1-human PUM2-human
## ENSG00000186891     0.00000     0.00000    22.92884    24.14742   45.61482
## ENSG00000158109    67.26265    54.52247    45.85769    48.29484   45.61482
## ENSG00000198912    67.26265    54.38616    45.85769    48.29484   45.38674
## ENSG00000236423    67.26265    54.52247    45.85769    48.29484   45.61482
## ENSG00000116237    65.91739    52.88679    44.48196    47.81189   43.79022
## ENSG00000069812    67.26265    40.89185    45.85769    34.53081   45.61482
## ENSG00000173662    33.63132    40.89185    22.92884    24.14742   45.61482
## ENSG00000180758    67.26265    27.26123    45.85769    48.29484   45.61482
## ENSG00000178585    66.92633    54.38616    45.85769    48.05337   45.61482
```

Having computed all the hyperparameters, testing for differential
usage is run. Note the `min.cpm` argument is set to filter out genes
based on CPM, as described above.


```r
result = csd.testForDEU(cdx, workers=1, tmp.dir=NULL, min.cpm=1)
```


```r
head(result[c("featureID", "condition", "y", "pvalue", "padj")])
```

```
##                featureID   condition y     pvalue      padj
## 319  ENSG00000173662:007 DDX27-human 0 0.01305160 0.8483537
## 365  ENSG00000173662:008 DDX27-human 0 0.01305160 0.8483537
## 1110 ENSG00000173662:003 DDX27-human 0 0.01409734 0.9163272
## 166  ENSG00000173662:004 DDX27-human 0 0.01409734 0.9163272
## 514  ENSG00000173662:011 DDX27-human 0 0.01409734 0.9163272
## 286  ENSG00000173662:006 NELFE-human 0 0.02140260 1.0000000
```

#### Session information

```r
sessionInfo()
```

```
## R version 3.3.2 (2016-10-31)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X El Capitan 10.11.6
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] csDEX_1.0                  reshape_0.8.6             
##  [3] betareg_3.1-0              aod_1.3                   
##  [5] edgeR_3.14.0               limma_3.28.21             
##  [7] DESeq2_1.12.4              SummarizedExperiment_1.2.3
##  [9] Biobase_2.32.0             GenomicRanges_1.24.3      
## [11] GenomeInfoDb_1.8.7         IRanges_2.6.1             
## [13] S4Vectors_0.10.3           BiocGenerics_0.18.0       
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.9          locfit_1.5-9.1       lattice_0.20-34     
##  [4] zoo_1.7-14           assertthat_0.1       digest_0.6.11       
##  [7] lmtest_0.9-34        plyr_1.8.4           backports_1.0.5     
## [10] acepack_1.4.1        RSQLite_1.1-2        evaluate_0.10       
## [13] ggplot2_2.2.1        zlibbioc_1.18.0      lazyeval_0.2.0      
## [16] data.table_1.10.0    annotate_1.50.1      rpart_4.1-10        
## [19] Matrix_1.2-7.1       checkmate_1.8.2      splines_3.3.2       
## [22] BiocParallel_1.6.6   geneplotter_1.50.0   stringr_1.1.0       
## [25] foreign_0.8-67       RCurl_1.95-4.8       munsell_0.4.3       
## [28] base64enc_0.1-3      htmltools_0.3.5      nnet_7.3-12         
## [31] tibble_1.2           gridExtra_2.2.1      htmlTable_1.8       
## [34] Hmisc_4.0-2          XML_3.98-1.5         bitops_1.0-6        
## [37] grid_3.3.2           xtable_1.8-2         gtable_0.2.0        
## [40] DBI_0.5-1            magrittr_1.5         scales_0.4.1        
## [43] stringi_1.1.2        XVector_0.12.1       genefilter_1.54.2   
## [46] flexmix_2.3-13       latticeExtra_0.6-28  sandwich_2.3-4      
## [49] Formula_1.2-1        RColorBrewer_1.1-2   tools_3.3.2         
## [52] survival_2.40-1      AnnotationDbi_1.34.4 colorspace_1.3-2    
## [55] cluster_2.0.5        memoise_1.0.0        knitr_1.15.1        
## [58] modeltools_0.2-21
```
