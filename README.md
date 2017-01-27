---
  title: "csDEX vignette"
date: "2017-01-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{csDEX vignette}
%\VignetteEngine{knitr::rmarkdown}
\usepackage[utf8]{inputenc}
---


# The csDEX package: Quick Start

  
Build with R CMD Sweave vignettes/csdex_vignette.Rmd
 
csDEX (Condition-specific Differential Exon Expression) is a family of general
linear models (GLMs), used to model differential splicing given a set of
sequencing experiments based on RNA-seq. The methods assume tens to hundreds
experimental conditions to be compared, and the result are
changes in usage of splice sites that are specific (particular) to an
experimental condition. 

Main features include:
* handling of sequence data in form of read counts or percentages (e.g. percent-spliced in, PSI),
* modelling all type of splicing events,
* discovery and ranking of condition-specific changes,
* parallel processing of multiple genes (groups),
* time-efficient model approximation with no change in false positive rate.

## Installation

The easiest installation method is via GitHub:

```r
# install.packages("")
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
* groupID: group of features belonging to the same transcriptional unit,
    most likely a gene.    
* featureID: alternative transcriptional unit, such as an exon, exonic
    part, alternative 5'/3' end, etc. Most commonly referring to a
    non-overlapping exonic bin, to which reads are mapped.
* value: expression quantification, either in form of read counts mapping
    to the feature (non-negative integers) or percentage spliced-in (PSI,
    real numbers 0.0-1.0, inclusive). Choice of quantification shall be
    consistent within and between all files.

Example for PSI:
```
    ENSG00000012345:ENSE00000067890 0.42
    ENSG00000012345:ENSE00000067891 0.71
    ENSG00000012346:ENSE00000067892 0
    ENSG00000012346:ENSE00000067893 1.0
```

An arbitrary number of replicates can be used for each experimental conditions.


#### Obtaining replicate files from BAM/SAM files. 

Given the number of available tools out there, we do not reinvent the wheel but
point users to existing tools. These will typically assume a genome annotation
listed in a .gtf/.gff format.
* Read count data:
** QoRTs
** DEXSeq
* PSI:
** cuffLinks
** rMATS
** manually from ENCODE pre-made transcript quantifications


### The design file.

The experimental design is stored as a tab-delimited table and links replicate filenames to
experimental conditions. It assumes (at least) two columns with arbitrary
names, denoting
* `File.accession`: replicate file name within the input directory, without extension.
* `Experiment.target`: experimental condition identifier
* `input.read.count`: (optional) Number of reads from which the replicate was generated. Required to compute gene counts-per-million reads (CPM) if required.

Files sharing the same `Experiment.target` are assumed to be technical replicates
of the same experimental condition. The column names are not fixed and can be defined at
runtime. Other columns can be stored, but will be ignored by csDEX. The format
is consistent with the ENCODE metadata format.

Two exemplary small datasets are provided with the csDEX package, to be
used for testing purposes only.



## Running csDEX

Once csDEX is installed, the setup can be tested using the provided data.

The example files are accessed as follows.

```r
require(csDEX)

design.file = system.file("extdata/metadata.tsv", 
    package="csDEX", mustWork=TRUE)
data.dir.psi = system.file("extdata/psi", 
    package="csDEX", mustWork=TRUE)
data.dir.count = system.file("extdata/count", 
    package="csDEX", mustWork=TRUE)
```

### The PSI model

The percentage-spliced in model is based on Beta regression GLM. Is is somewhat
more streamlined and therefore presented first.  The data is loaded into a
csDEXDataSet object:


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

The result is table of features and conditions ranked by statistical
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
