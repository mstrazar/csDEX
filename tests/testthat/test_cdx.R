context("Test csDEXdataSet")

design.file = system.file("extdata/metadata.tsv", package="csDEX", mustWork=TRUE)
data.dir.psi = system.file("extdata/psi", package="csDEX", mustWork=TRUE)
data.dir.count = system.file("extdata/count", package="csDEX", mustWork=TRUE)

### Test loading
test_that("test csDEX dataset - input error", {
  expect_error(
    csDEXdataSet(data.dir=NULL, design.file=NULL, type="")
  )
})

### 1 worker workflows

test_that("test csDEX dataset - PSI (standard workflow)", {
  cdx = csDEXdataSet(data.dir=data.dir.psi, design.file=design.file, type="PSI")
  expect_s4_class(cdx, "csDEXdataSet")
  result = testForDEU(cdx, workers=1, tmp.dir=NULL)
})

test_that("test csDEX dataset - PSI (reduced workflow with Wald test)", {
  cdx = csDEXdataSet(data.dir=data.dir.psi, design.file=design.file, type="PSI")
  expect_s4_class(cdx, "csDEXdataSet")
  result = testForDEU(cdx, workers=1, tmp.dir=NULL, alpha.wald=0.05)
})

test_that("test csDEX dataset - PSI (data with offset)", {
  # Simulate deltaPSI in (-1, 1)
  cdx = csDEXdataSet(data.dir=data.dir.psi, design.file=design.file, type="PSI")
  exprData(cdx) = 2 * exprData(cdx) - 1
  expect_s4_class(cdx, "csDEXdataSet")
  resultDP = testForDEU(cdx, workers=1, tmp.dir=NULL, scale=2, offset=-1)
})

test_that("test csDEX dataset - PSI (data with offset and controls)", {
  # Simulate deltaPSI in (-1, 1)
  cdx = csDEXdataSet(data.dir=data.dir.psi, 
                     design.file=design.file, 
                     type="PSI", 
                     col.condition = "Experiment.accession",
                     col.additional =  c("Controlled.by"))
  cdxn = normalizeToControls(cdx)
  expect_s4_class(cdxn, "csDEXdataSet")
  resultDP = testForDEU(cdxn, workers=1, tmp.dir=NULL, scale=2, offset=-1)
})

test_that("test csDEX dataset - count (standard workflow)", {
  cdx = csDEXdataSet(data.dir=data.dir.count, design.file=design.file, type="count")
  expect_s4_class(cdx, "csDEXdataSet")
  cdx = estimateSizeFactors(cdx)
  
  # Test accessor functions
  expr = exprData(cdx, normalized=FALSE)
  expr = exprData(cdx, normalized=TRUE)
  rowdata = rowData(cdx)
  coldata = colData(cdx)
  rowData(cdx) <- rowdata
  colData(cdx) <- coldata
  
  cdx = estimatePrecisions(cdx)
  cdx = estimateGeneCPM(cdx)
  cpm = cpmData(cdx)  
  cpmData(cdx) <- cpm
  
  result = testForDEU(cdx, workers=1, tmp.dir=NULL, min.cpm=1)
})

test_that("test csDEX dataset - PSI (reduced workflow with Wald test)", {
  cdx = csDEXdataSet(data.dir=data.dir.psi, design.file=design.file, type="PSI")
  expect_s4_class(cdx, "csDEXdataSet")
  result = testForDEU(cdx, workers=1, tmp.dir=NULL, alpha.wald=0.05)
})

### 2 worker workflows
test_that("test csDEX dataset - PSI (standard workflow)", {
  cdx = csDEXdataSet(data.dir=data.dir.psi, design.file=design.file, type="PSI")
  expect_s4_class(cdx, "csDEXdataSet")
  result = testForDEU(cdx, workers=2, tmp.dir=NULL)
})

test_that("test csDEX dataset - PSI (data with offset)", {
  # Simulate deltaPSI in (-1, 1)
  cdx = csDEXdataSet(data.dir=data.dir.psi, design.file=design.file, type="PSI")
  exprData(cdx) = 2 * exprData(cdx) - 1
  expect_s4_class(cdx, "csDEXdataSet")
  resultDP = testForDEU(cdx, workers=2, tmp.dir=NULL, scale=2, offset=-1)
})


### Additional fields
test_that("test csDEX dataset - PSI (additional column)", {
  cdx = csDEXdataSet(data.dir=data.dir.psi, design.file=design.file, type="PSI",
                     col.condition="Experiment.accession",
                     col.additional=c("Experiment.target", "Controlled.by"))
  expect_s4_class(cdx, "csDEXdataSet")
  result = testForDEU(cdx, workers=1, tmp.dir=NULL, 
                      formula = y ~ featureID + Experiment.target + Controlled.by)
})

### Additional fields - invalid
test_that("test csDEX dataset - PSI (wrong additional column)", {
  expect_error(
    csDEXdataSet(data.dir=data.dir.psi, design.file=design.file, type="PSI",
                     col.condition="Experiment.target",
                     col.additional=c("Derived.from"))
  )
})