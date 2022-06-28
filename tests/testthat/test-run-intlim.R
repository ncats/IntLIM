# If a value other than a data frame is input, an error should be thrown.
testthat::test_that("Inputting the wrong class causes early termination.", {
  inputData <- "Hello World"
  testthat::expect_error(IntLIM::RunIntLim(inputData, stype = "Feat1"), "Input must be an IntLimData object", 
               ignore.case = TRUE)
})

# If the input has more than 2 classes and is not continuous, an error should be thrown.
testthat::test_that("Inputting more than 2 classes causes early termination.", {
  
  # Create toy data.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0),
                      "Level"=c("Low", "Medium", "High", "Medium"))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  geneData <- data.frame("P1"=c(0,0,0), "P2"=c(0,0,0), "P3"=c(0,0,0),
                         "P4"=c(0,0,0))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(0,0,0), "P2"=c(0,0,0), "P3"=c(0,0,0),
                          "P4"=c(0,0,0))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)
  
  testthat::expect_error(IntLIM::RunIntLim(dat, stype = "Level"), 
               paste("IntLim currently requires only two categories.  Make sure the column",
                     "Level only has two unique values. Did you mean to set continuous to TRUE?"),
               ignore.case = TRUE)
})

# Check invalid parameters for input or outcome type.
testthat::test_that("Invalid parameters cause an error.", {
  
  # Create toy data.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0),
                      "Level"=c("Low", "Medium", "Low", "Medium"))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  geneData <- data.frame("P1"=c(0,0,0), "P2"=c(0,0,0), "P3"=c(0,0,0),
                         "P4"=c(0,0,0))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(0,0,0), "P2"=c(0,0,0), "P3"=c(0,0,0),
                          "P4"=c(0,0,0))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)
  
  testthat::expect_error(IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = "Scooby Dooby Doo",
                                 outcome = "Where are you?"), paste(
                                   "Error! independent.var.type and outcome.type must",
                                   "both be either 1 or 2 in RunIntLim."), ignore.case = TRUE)
})

# If the data is missing, an error should be thrown.
testthat::test_that("Inputting missing data types causes early termination.", {
  # Create toy data.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0),
                      "Level"=c("Low", "Medium", "Low", "Medium"))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  geneData <- data.frame("P1"=c(0,0,0), "P2"=c(0,0,0), "P3"=c(0,0,0),
                         "P4"=c(0,0,0))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(0,0,0), "P2"=c(0,0,0), "P3"=c(0,0,0),
                          "P4"=c(0,0,0))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat1 <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                       analyteType2=as.matrix(geneData),
                       analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                       analyteType2MetaData = geneMetaData,
                       sampleMetaData = pData)
  dat2 <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                       analyteType2=matrix(, nrow = 0, ncol = 0),
                       analyteType1MetaData = metabMetaData,
                       analyteType2MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                       sampleMetaData = pData)

  # Check for missing data.
  testthat::expect_error(IntLIM::RunIntLim(dat1, stype = "Level", independent.var.type = 1,
                                 outcome = 2),
               "One type of analyte data is missing. Cannot run.",
               ignore.case = TRUE)
  testthat::expect_error(IntLIM::RunIntLim(dat2, stype = "Level", independent.var.type = 2,
                                 outcome = 1),
               "One type of analyte data is missing. Cannot run.",
               ignore.case = TRUE)

  # Check for missing data when analyte type is same.
  testthat::expect_error(IntLIM::RunIntLim(dat1, stype = "Level", independent.var.type = 1,
                                 outcome = 1), "Analyte type 1 is missing. Cannot run.",
               ignore.case = TRUE)
  testthat::expect_error(IntLIM::RunIntLim(dat2, stype = "Level", independent.var.type = 2,
                                 outcome = 2), "Analyte type 2 is missing. Cannot run.",
               ignore.case = TRUE)
})

# Check that all fields are populated correctly for discrete data.
testthat::test_that("Fields are populated appropriately (discrete).", {

  # Create toy data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8),
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c("Low", "Medium", "Low", "Medium", "Medium", "Low",
                                "Low", "Medium"))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  geneData <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                          "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                          "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                          "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)

  # Run with different analytes and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 2,
                               suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)

  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)), sort(c("interaction.adj.pvalues", 
                                                    "interaction.coefficients",
                                                    "model.rsquared", 
                                                    "covariate.pvalues", 
                                                    "covariate.coefficients", 
                                                    "corr", "filt.results",
                                                    "warnings", "stype", 
                                                    "outcome", 
                                                    "independent.var.type", 
                                                    "covar", "interaction.pvalues",
                                                    "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)

  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)

  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 1, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 2, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)
})

# Check that all fields are populated correctly with plus signs in the feature data.
testthat::test_that("Fields are populated appropriately when feature data contains plus signs.", {
  
  # Create toy data.
  pData <- data.frame("Feat1"=c("one", "two+three", "one", "z", "one",
                                "z", "two+three", "two+three"),
                      "Level"=c("Low", "Medium", "Low", "Medium", "Medium", "Low",
                                "Low", "Medium"))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  geneData <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                          "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                          "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                          "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)
  
  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1"),
                               independent.var.type = 2, outcome = 1, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)), sort(c("interaction.adj.pvalues", 
                                                    "interaction.coefficients",
                                                    "model.rsquared", "covariate.pvalues", 
                                                    "covariate.coefficients", "corr", 
                                                    "filt.results",
                                                    "warnings", "stype", "outcome", 
                                                    "independent.var.type", "covar", 
                                                    "interaction.pvalues",
                                                    "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)
  
  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1"),
                               independent.var.type = 1, outcome = 1, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1"),
                               independent.var.type = 2, outcome = 2, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 0)
})

# Check that all fields are populated correctly for continuous data.
testthat::test_that("Fields are populated appropriately (continuous).", {

  # Create toy data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8),
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  geneData <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                          "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                          "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                          "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)

  # Run with different analytes and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 2,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues","continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)), sort(c("interaction.adj.pvalues", 
                                                    "interaction.coefficients",
                                                    "model.rsquared", 
                                                    "covariate.pvalues", 
                                                    "covariate.coefficients", 
                                                    "corr", "filt.results",
                                                    "warnings", "stype", "outcome", 
                                                    "independent.var.type", 
                                                    "covar", "interaction.pvalues",
                                                    "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 1, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 2, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)
})

# Check that all fields are populated correctly when covariates are factors.
testthat::test_that("Fields are populated appropriately with factor covariates.", {

  # Create toy data.
  pData <- data.frame("Feat1"=c("one", "two", "three", "one", "two", "three",
                                "one", "two"),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  geneData <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                          "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                          "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                          "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)

  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1"),
                               independent.var.type = 2, outcome = 1, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)), sort(c("interaction.adj.pvalues", 
                                                    "interaction.coefficients",
                                                    "model.rsquared", 
                                                    "covariate.pvalues", 
                                                    "covariate.coefficients", 
                                                    "corr", "filt.results",
                                                    "warnings", "stype", 
                                                    "outcome", 
                                                    "independent.var.type", 
                                                    "covar", "interaction.pvalues",
                                                    "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1"),
                               independent.var.type = 1, outcome = 1, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1"),
                               independent.var.type = 2, outcome = 2, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)
})

# Check that all fields are populated correctly when some covariates are factors.
testthat::test_that("Fields are populated appropriately with a mix of numeric and factor covariates.", {

  # Create toy data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8),
                      "Feat2"=c("one", "one", "two", "two", "one", "one", "two",
                                "two"),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  geneData <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                          "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                          "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                          "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)

  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2"),
                               independent.var.type = 2, outcome = 1, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)), sort(c("interaction.adj.pvalues", 
                                                    "interaction.coefficients",
                                                    "model.rsquared", 
                                                    "covariate.pvalues", 
                                                    "covariate.coefficients", 
                                                    "corr", "filt.results",
                                                    "warnings", "stype", 
                                                    "outcome", 
                                                    "independent.var.type", 
                                                    "covar", 
                                                    "interaction.pvalues",
                                                    "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2"),
                               independent.var.type = 1, outcome = 1, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2"),
                               independent.var.type = 2, outcome = 2, suppressWarnings = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results","warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)
})

# Check remove.duplicates.
testthat::test_that("Duplicates are removed when appropriate (continuous).", {

  # Create toy data (discrete).
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8),
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  geneData <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                          "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                          "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                          "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)

  # Run with different analytes and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1,
                               outcome = 2, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_warning(IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1,
                                   outcome = 2, remove.duplicates = TRUE,
                                   continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2,
                               outcome = 1, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_warning(IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2,
                                   outcome = 1, remove.duplicates = TRUE, continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_warning(IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                                   independent.var.type = 1, outcome = 2, remove.duplicates = TRUE, continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_warning(IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                                   independent.var.type = 2, outcome = 1, remove.duplicates = TRUE, continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1,
                               outcome = 1, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  testthat::expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  testthat::expect_equal(length(which(is.na(results@model.rsquared))), 6)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2,
                               outcome = 2, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  testthat::expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  testthat::expect_equal(length(which(is.na(results@model.rsquared))), 6)
  testthat::expect_equal(results@continuous, 1)

  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               remove.duplicates = TRUE, independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  testthat::expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  testthat::expect_equal(length(which(is.na(results@model.rsquared))), 6)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               remove.duplicates = TRUE, independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", 
                          "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", 
                          "covar", "interaction.pvalues", "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  testthat::expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  testthat::expect_equal(length(which(is.na(results@model.rsquared))), 6)
  testthat::expect_equal(results@continuous, 1)
})

testthat::test_that("Duplicates are removed when appropriate (continuous).", {

  # Create toy data (discrete).
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8),
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  geneData <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                          "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                          "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                          "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)

  # Run with different analytes and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1,
                               outcome = 2, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  testthat::expect_warning(IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1,
                                   outcome = 2, remove.duplicates = TRUE,
                                   continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2,
                               outcome = 1, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  testthat::expect_warning(IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2,
                                   outcome = 1, remove.duplicates = TRUE,
                                   continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  testthat::expect_warning(IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                                   independent.var.type = 1, outcome = 2, remove.duplicates = TRUE,
                                   continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  testthat::expect_warning(IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                                   independent.var.type = 2, outcome = 1, remove.duplicates = TRUE,
                                   continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(results@continuous, 1)

  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1,
                               outcome = 1, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  testthat::expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  testthat::expect_equal(length(which(is.na(results@model.rsquared))), 6)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2,
                               outcome = 2, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  testthat::expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  testthat::expect_equal(length(which(is.na(results@model.rsquared))), 6)
  testthat::expect_equal(results@continuous, 1)

  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               remove.duplicates = TRUE, independent.var.type = 1, outcome = 1,
                               continuous = TRUE, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  testthat::expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  testthat::expect_equal(length(which(is.na(results@model.rsquared))), 6)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               remove.duplicates = TRUE, independent.var.type = 2, outcome = 2,
                               continuous = TRUE, suppressWarnings = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(length(results@covariate.pvalues), 0)
  testthat::expect_equal(length(results@covariate.coefficients), 0)
  testthat::expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  testthat::expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  testthat::expect_equal(length(which(is.na(results@model.rsquared))), 6)
  testthat::expect_equal(results@continuous, 1)
})

# Check saving covariate pvals.
testthat::test_that("Covariate p-values are saved (discrete).", {

  # Create toy data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8),
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c("Low", "Medium", "Low", "Medium", "Medium", "Low",
                                "Low", "Medium"))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  geneData <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                          "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                          "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                          "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)

  # Run with different analytes and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 2,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  testthat::expect_equal(results@continuous, 0)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  testthat::expect_equal(results@continuous, 0)

  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE,
                               save.covar.pvals = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  testthat::expect_equal(results@continuous, 0)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  testthat::expect_equal(results@continuous, 0)

  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  testthat::expect_equal(results@continuous, 0)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  testthat::expect_equal(results@continuous, 0)

  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  testthat::expect_equal(results@continuous, 0)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  testthat::expect_equal(results@continuous, 0)
})

# Check saving covariate pvals.
testthat::test_that("Covariate p-values are saved (continuous).", {

  # Create toy data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8),
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  geneData <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                          "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                          "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                          "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)

  # Run with different analytes and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 2,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  testthat::expect_equal(results@continuous, 1)

  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE,
                               save.covar.pvals = TRUE, continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  testthat::expect_equal(results@continuous, 1)

  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  testthat::expect_equal(results@continuous, 1)

  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  testthat::expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  testthat::expect_equal(results@continuous, 1)

  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  testthat::expect_identical(sort(slotNames(results)),
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients",
                          "corr", "filt.results", "warnings", "stype", "outcome",
                          "independent.var.type", "covar", "interaction.pvalues",
                          "continuous")))
  testthat::expect_equal(length(results@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(results@interaction.coefficients), 9)
  testthat::expect_equal(length(results@model.rsquared), 9)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  testthat::expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  testthat::expect_equal(results@continuous, 1)
})

# Check for correlated variables.
testthat::test_that("Check for correlated variables.", {
  # Create toy data.
  pData1 <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8),
                       "Feat2"=c(seq(1:8)),
                       "Feat3"=c(seq(1:8)),
                       "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData1) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                        "P7", "P8")
  pData2 <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8),
                       "Feat2"=c(seq(1:8)),
                       "Feat3"=c(rev(seq(1:8))),
                       "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData2) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                        "P7", "P8")
  geneData <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                          "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                          "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                          "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat1 <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                       analyteType2=as.matrix(geneData),
                       analyteType1MetaData = metabMetaData,
                       analyteType2MetaData = geneMetaData,
                       sampleMetaData = pData1)
  dat2 <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                       analyteType2=as.matrix(geneData),
                       analyteType1MetaData = metabMetaData,
                       analyteType2MetaData = geneMetaData,
                       sampleMetaData = pData2)

  # Run
  run1 <- IntLIM::RunIntLim(dat1, stype = "Level", independent.var.type = 1,
                            outcome = 2, covar = c("Feat1", "Feat2", "Feat3"),
                            continuous = TRUE)
  testthat::expect_identical(run1@warnings[[1]], "Using pseudoinverse for Metab1")
  testthat::expect_identical(run1@warnings[[2]], paste0("The following covariates have correlation > 0.9 : ",
                                              "(Feat2, Feat3)"))
  run2 <- IntLIM::RunIntLim(dat2, stype = "Level", independent.var.type = 1,
                            outcome = 2, covar = c("Feat1", "Feat2", "Feat3"),
                            continuous = TRUE)
  testthat::expect_identical(run2@warnings[[1]], "Using pseudoinverse for Metab1")
  testthat::expect_identical(run2@warnings[[3]], paste0("The following covariates have correlation < -0.9 : ",
                                              "(Feat2, Feat3)"))
})

# Removes genes and metabolites with a standard deviation of zero
testthat::test_that("If standard deviation is zero, analytes are removed.", {
  # With and without standard deviation of zero
  pData1 <- data.frame("Feat1"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8),
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8))
  rownames(pData1) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                        "P7", "P8")
  pData2 <- data.frame("Feat1"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8),
                       "Feat2"=c(rep(0,8)),
                       "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8))
  rownames(pData2) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                        "P7", "P8")
  pData3 <- data.frame("Feat1"=c(rep(0,8)),
                       "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                       "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8))
  rownames(pData3) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                        "P7", "P8")

  # Genes with and without standard deviation of zero.
  geneData1 <- data.frame("P1"=c(46.1,20.2,59.3), "P2"=c(11.1,34.2,19.3),
                         "P3"=c(28.1,71.2,94.3), "P4"=c(51.1,91.2,32.3),
                         "P5"=c(73.1,26.2,40.3), "P6"=c(91.1,99.2,12.3),
                         "P7"=c(38.1,44.2,60.3), "P8"=c(91.1,93.2,63.3))
  rownames(geneData1) <- c("Gene1", "Gene2", "Gene3")
  geneData2 <- data.frame("P1"=c(0,20.2,59.3), "P2"=c(0,34.2,19.3),
                          "P3"=c(0,71.2,94.3), "P4"=c(0,91.2,32.3),
                          "P5"=c(0,26.2,40.3), "P6"=c(0,99.2,12.3),
                          "P7"=c(0,44.2,60.3), "P8"=c(0,93.2,63.3))
  rownames(geneData2) <- c("Gene1", "Gene2", "Gene3")

  # Metabolites with and without standard deviation of zero.
  metabData1 <- data.frame("P1"=c(0,32.2,81.3), "P2"=c(0,58.2,45.3),
                           "P3"=c(0,61.2,67.3), "P4"=c(0,7.2,79.3),
                           "P5"=c(0,87.2,91.3), "P6"=c(0,87.2,91.3),
                           "P7"=c(0,10.2,85.3), "P8"=c(0,14.2,76.3))
  rownames(metabData1) <- c("Metab1", "Metab2", "Metab3")
  metabData2 <- data.frame("P1"=c(60.1,32.2,81.3), "P2"=c(68.1,58.2,45.3),
                           "P3"=c(30.1,61.2,67.3), "P4"=c(36.1,7.2,79.3),
                           "P5"=c(5.1,87.2,91.3), "P6"=c(5.1,87.2,91.3),
                           "P7"=c(99.1,10.2,85.3), "P8"=c(51.1,14.2,76.3))
  rownames(metabData2) <- c("Metab1", "Metab2", "Metab3")

  # Metadata
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))

  # Input data
  dat1 <- methods::new("IntLimData", analyteType1=as.matrix(metabData1),
                       analyteType2=as.matrix(geneData1),
                       analyteType1MetaData = metabMetaData,
                       analyteType2MetaData = geneMetaData,
                       sampleMetaData = pData1)
  dat2 <- methods::new("IntLimData", analyteType1=as.matrix(metabData2),
                       analyteType2=as.matrix(geneData2),
                       analyteType1MetaData = metabMetaData,
                       analyteType2MetaData = geneMetaData,
                       sampleMetaData = pData1)
  dat3 <- methods::new("IntLimData", analyteType1=as.matrix(metabData2),
                       analyteType2=as.matrix(geneData1),
                       analyteType1MetaData = metabMetaData,
                       analyteType2MetaData = geneMetaData,
                       sampleMetaData = pData2)
  dat4 <- methods::new("IntLimData", analyteType1=as.matrix(metabData2),
                       analyteType2=as.matrix(geneData1),
                       analyteType1MetaData = metabMetaData,
                       analyteType2MetaData = geneMetaData,
                       sampleMetaData = pData3)

  # Metabolite with standard deviation of 0
  run1 <- IntLIM::RunIntLim(dat1, stype = "Feat1", independent.var.type = 1,
                            outcome = 2, continuous = TRUE, covar = c("Feat2", "Feat3"))
  testthat::expect_identical(run1@warnings[[1]],
                   "Removed 1 analytes of type 1 that had a standard deviation of 0: Metab1")
  testthat::expect_equal(length(run1@interaction.adj.pvalues), 6)
  testthat::expect_equal(length(run1@interaction.coefficients), 6)
  testthat::expect_equal(length(run1@model.rsquared), 6)
  testthat::expect_equal(length(run1@continuous), 1)

  # Gene with standard deviation of 0
  run2 <- IntLIM::RunIntLim(dat2, stype = "Feat1", independent.var.type = 1,
                            outcome = 2, continuous = TRUE, covar = c("Feat2", "Feat3"))
  testthat::expect_identical(run2@warnings[[1]],
                   "Removed 1 analytes of type 2 that had a standard deviation of 0: Gene1")
  testthat::expect_equal(length(run2@interaction.adj.pvalues), 6)
  testthat::expect_equal(length(run2@interaction.coefficients), 6)
  testthat::expect_equal(length(run2@model.rsquared), 6)
  testthat::expect_equal(length(run2@continuous), 1)

  # Covariates with standard deviation of 0
  run3 <- IntLIM::RunIntLim(dat3, stype = "Feat1", independent.var.type = 1,
                            outcome = 2, continuous = TRUE, covar = c("Feat2", "Feat3"))
  testthat::expect_identical(run3@warnings[[1]],
                   "Removed 1 covariates that had a standard deviation of 0: Feat2")
  testthat::expect_equal(length(run3@interaction.adj.pvalues), 9)
  testthat::expect_equal(length(run3@interaction.coefficients), 9)
  testthat::expect_equal(length(run3@model.rsquared), 9)
  testthat::expect_equal(length(run3@covar), 1)
  testthat::expect_equal(length(run3@continuous), 1)

  # Phenotype with standard deviation of 0
  testthat::expect_error(IntLIM::RunIntLim(dat4, stype = "Feat1", independent.var.type = 1,
                                 outcome = 2, continuous = TRUE, covar = c("Feat2", "Feat3")),
               "stype variable has a standard deviation of zero. Cannot run.")
})