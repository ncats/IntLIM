# If a value other than an IntLimResults object is input, an error should be thrown.
testthat::test_that("Inputting the wrong class causes early termination.", {
  
  # Generate input data.
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
  inputDataGood <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                                analyteType2=as.matrix(geneData),
                                analyteType1MetaData = metabMetaData,
                                analyteType2MetaData = geneMetaData,
                                sampleMetaData = pData)
  inputDataBad <- "Hello World"
  
  # Generate result data.
  pvals <- matrix(rep(0.2, 9), nrow = 3, ncol = 3)
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  coef <- matrix(rep(8, 9), nrow = 3, ncol = 3)
  rsq <- coef <- matrix(rep(0.6, 9), nrow = 3, ncol = 3)
  inputResultsGood <- methods::new("IntLimResults", 
                                   interaction.pvalues=pvals, 
                                   interaction.adj.pvalues=adj_pvals,
                                   interaction.coefficients=coef, 
                                   model.rsquared = rsq, 
                                   covariate.pvalues = data.frame(matrix(, nrow = 0, ncol = 0)),
                                   covariate.coefficients = data.frame(matrix(, nrow = 0, ncol = 0)), 
                                   corr=data.frame(matrix(, nrow = 0, ncol = 0)),
                                   filt.results=data.frame(matrix(, nrow = 0, ncol = 0)),
                                   warnings=list(),
                                   stype="Level", 
                                   outcome=1,
                                   independent.var.type=2,
                                   covar="",
                                   continuous=1)
  inputResultsBad <- "Hello World"
  
  
  testthat::expect_error(IntLIM::ProcessResults(inputResultsBad, inputDataGood), 
               paste("Results must be an IntLIMResults object"), 
               ignore.case = TRUE)
  testthat::expect_error(IntLIM::ProcessResults(inputResultsGood, inputDataBad), 
               paste("Data must be an IntLIMData object"), 
               ignore.case = TRUE)
})

# Types for both outcome and independent analyte must be either 1 or 2.
testthat::test_that("Check that outcome and independent analyte types are appropriate.", {
  
  # Generate input data.
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
  inputData <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                            analyteType2=as.matrix(geneData),
                            analyteType1MetaData = metabMetaData,
                            analyteType2MetaData = geneMetaData,
                            sampleMetaData = pData)
  
  # Generate result data.
  pvals <- matrix(rep(0.2, 9), nrow = 3, ncol = 3)
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  coef <- matrix(rep(8, 9), nrow = 3, ncol = 3)
  rsq <- coef <- matrix(rep(0.6, 9), nrow = 3, ncol = 3)
  inputResults1 <- methods::new("IntLimResults",
                                interaction.pvalues=pvals,
                                interaction.adj.pvalues=adj_pvals,
                                interaction.coefficients=coef,
                                model.rsquared = rsq,
                                covariate.pvalues = data.frame(matrix(, nrow = 0, ncol = 0)),
                                covariate.coefficients = data.frame(matrix(, nrow = 0, ncol = 0)),
                                corr=data.frame(matrix(, nrow = 0, ncol = 0)),
                                filt.results=data.frame(matrix(, nrow = 0, ncol = 0)),
                                warnings=list(),
                                stype="Level",
                                outcome=50,
                                independent.var.type=2,
                                covar="",
                                continuous=1)
  inputResults2 <- methods::new("IntLimResults",
                                interaction.pvalues=pvals,
                                interaction.adj.pvalues=adj_pvals,
                                interaction.coefficients=coef,
                                model.rsquared = rsq,
                                covariate.pvalues = data.frame(matrix(, nrow = 0, ncol = 0)),
                                covariate.coefficients = data.frame(matrix(, nrow = 0, ncol = 0)),
                                corr=data.frame(matrix(, nrow = 0, ncol = 0)),
                                filt.results=data.frame(matrix(, nrow = 0, ncol = 0)),
                                warnings=list(),
                                stype="Level",
                                outcome=1,
                                independent.var.type=-3,
                                covar="",
                                continuous=1)
  
  
  testthat::expect_error(IntLIM::ProcessResults(inputResults1, inputData), 
               paste("Independent variable and outcome must both",
                     "be either 1 or 2."),
               ignore.case = TRUE)
  testthat::expect_error(IntLIM::ProcessResults(inputResults2, inputData), 
               paste("Independent variable and outcome must both",
                     "be either 1 or 2."),
               ignore.case = TRUE)
})


# The outcome and independent analytes must both be found in the InputData.
testthat::test_that("An error is thrown if an analyte type is used that is not present in the data.", {
  
  # Generate input data.
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
  inputData1 <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                             analyteType2=matrix(, nrow = 0, ncol = 0),
                             analyteType1MetaData = metabMetaData,
                             analyteType2MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                             sampleMetaData = pData)
  inputData2 <- methods::new("IntLimData", analyteType1= matrix(, nrow = 0, ncol = 0),
                             analyteType2=as.matrix(geneData),
                             analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                             analyteType2MetaData = geneMetaData,
                             sampleMetaData = pData)
  
  # Generate result data.
  pvals <- matrix(rep(0.2, 9), nrow = 3, ncol = 3)
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  coef <- matrix(rep(8, 9), nrow = 3, ncol = 3)
  rsq <- coef <- matrix(rep(0.6, 9), nrow = 3, ncol = 3)
  inputResults <- methods::new("IntLimResults",
                               interaction.pvalues=pvals,
                               interaction.adj.pvalues=adj_pvals,
                               interaction.coefficients=coef,
                               model.rsquared = rsq,
                               covariate.pvalues = data.frame(matrix(, nrow = 0, ncol = 0)),
                               covariate.coefficients = data.frame(matrix(, nrow = 0, ncol = 0)),
                               corr=data.frame(matrix(, nrow = 0, ncol = 0)),
                               filt.results=data.frame(matrix(, nrow = 0, ncol = 0)),
                               warnings=list(),stype="Level",
                               outcome=1,
                               independent.var.type=2,
                               covar="",
                               continuous=1)
  
  testthat::expect_error(IntLIM::ProcessResults(inputResults, inputData2), 
               paste("Outcome type is not present in original data"),
               ignore.case = TRUE)
  testthat::expect_error(IntLIM::ProcessResults(inputResults, inputData1), 
               paste("Independent data type is not present in",
                     "original data"),
               ignore.case = TRUE)
})

# The type variable must have two levels or be continuous.
testthat::test_that("More than two levels causes an error.", {
  
  # Generate input data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c("Low","Medium","Low","Medium","Low","Medium",
                                "Low","High"))
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
  inputData <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                            analyteType2=as.matrix(geneData),
                            analyteType1MetaData = metabMetaData,
                            analyteType2MetaData = geneMetaData,
                            sampleMetaData = pData)
  
  # Generate result data.
  pvals <- matrix(rep(0.2, 9), nrow = 3, ncol = 3)
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  coef <- matrix(rep(8, 9), nrow = 3, ncol = 3)
  rsq <- coef <- matrix(rep(0.6, 9), nrow = 3, ncol = 3)
  inputResults <- methods::new("IntLimResults",
                               interaction.pvalues=pvals,
                               interaction.adj.pvalues=adj_pvals,
                               interaction.coefficients=coef,
                               model.rsquared = rsq,
                               covariate.pvalues = data.frame(matrix(, nrow = 0, ncol = 0)),
                               covariate.coefficients = data.frame(matrix(, nrow = 0, ncol = 0)),
                               corr=data.frame(matrix(, nrow = 0, ncol = 0)),
                               filt.results=data.frame(matrix(, nrow = 0, ncol = 0)),
                               warnings=list(),
                               stype="Level",
                               outcome=1,
                               independent.var.type=2,
                               covar="",
                               continuous=0)
  
  testthat::expect_error(IntLIM::ProcessResults(inputResults, inputData), paste(
    "IntLim requires two categories only for correlation analysis. Make sure the column",
    "Level only has two unique values or is continuous"), ignore.case = TRUE)
})

# Check that out-of-bounds values for interaction coefficient, r-squared,
# and p-value are not allowed.
testthat::test_that("Out of bounds values are not allowed.", {
  
  # Generate input data.
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
  inputData <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                            analyteType2=as.matrix(geneData),
                            analyteType1MetaData = metabMetaData,
                            analyteType2MetaData = geneMetaData,
                            sampleMetaData = pData)
  
  # Generate result data.
  pvals <- matrix(rep(0.2, 9), nrow = 3, ncol = 3)
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  coef <- matrix(rep(8, 9), nrow = 3, ncol = 3)
  rsq <- coef <- matrix(rep(0.6, 9), nrow = 3, ncol = 3)
  inputResults <- methods::new("IntLimResults",
                               interaction.pvalues=pvals,
                               interaction.adj.pvalues=adj_pvals,
                               interaction.coefficients=coef,
                               model.rsquared = rsq,
                               covariate.pvalues = data.frame(matrix(, nrow = 0, ncol = 0)),
                               covariate.coefficients = data.frame(matrix(, nrow = 0, ncol = 0)),
                               corr=data.frame(matrix(, nrow = 0, ncol = 0)),
                               filt.results=data.frame(matrix(, nrow = 0, ncol = 0)),
                               warnings=list(),
                               stype="Level",
                               outcome=1,
                               independent.var.type=2,
                               covar="",
                               continuous=0)
  
  # Check boundaries
  testthat::expect_error(IntLIM::ProcessResults(inputResults, inputData, pvalcutoff = -2), paste(
    "P-value must be between 0 and 1"), ignore.case = TRUE)
  testthat::expect_error(IntLIM::ProcessResults(inputResults, inputData, pvalcutoff = 2), paste(
    "P-value must be between 0 and 1"), ignore.case = TRUE)
  testthat::expect_error(IntLIM::ProcessResults(inputResults, inputData, rsquaredCutoff = -2), paste(
    "R-squared value must be between 0 and 1"), ignore.case = TRUE)
  testthat::expect_error(IntLIM::ProcessResults(inputResults, inputData, rsquaredCutoff = 2), paste(
    "R-squared value must be between 0 and 1"), ignore.case = TRUE)
  testthat::expect_error(IntLIM::ProcessResults(inputResults, inputData, interactionCoeffPercentile = -2), paste(
    "Interaction coefficient percentile must be between 0 and 1"), ignore.case = TRUE)
  testthat::expect_error(IntLIM::ProcessResults(inputResults, inputData, interactionCoeffPercentile = 2), paste(
    "Interaction coefficient percentile must be between 0 and 1"), ignore.case = TRUE)
})

# Check that everything is returned when there is no filtering.
testthat::test_that("Data is returned appropriately with no filtering.", {

  # Generate input data (discrete).
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
  inputData <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                            analyteType2=as.matrix(geneData),
                            analyteType1MetaData = metabMetaData,
                            analyteType2MetaData = geneMetaData,
                            sampleMetaData = pData)

  # Generate input data (continuous).
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  inputDataCont <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                            analyteType2=as.matrix(geneData),
                            analyteType1MetaData = metabMetaData,
                            analyteType2MetaData = geneMetaData,
                            sampleMetaData = pData)

  # Generate result data.
  pvals <- matrix(rep(0.2, 9), nrow = 3, ncol = 3)
  rownames(pvals) <- c("Gene1", "Gene2", "Gene3")
  colnames(pvals) <- c("Metab1", "Metab2", "Metab3")
  adj_pvals <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  rownames(adj_pvals) <- c("Gene1", "Gene2", "Gene3")
  colnames(adj_pvals) <- c("Metab1", "Metab2", "Metab3")
  coef <- matrix(rep(8, 9), nrow = 3, ncol = 3)
  rownames(coef) <- c("Gene1", "Gene2", "Gene3")
  colnames(coef) <- c("Metab1", "Metab2", "Metab3")
  rsq <- coef <- matrix(rep(0.6, 9), nrow = 3, ncol = 3)
  rownames(rsq) <- c("Gene1", "Gene2", "Gene3")
  colnames(rsq) <- c("Metab1", "Metab2", "Metab3")
  inputResults <- methods::new("IntLimResults",
                               interaction.pvalues=pvals,
                               interaction.adj.pvalues=adj_pvals,
                               interaction.coefficients=coef,
                               model.rsquared = rsq,
                               covariate.pvalues = data.frame(matrix(, nrow = 0, ncol = 0)),
                               covariate.coefficients = data.frame(matrix(, nrow = 0, ncol = 0)),
                               corr=data.frame(matrix(, nrow = 0, ncol = 0)),
                               filt.results=data.frame(matrix(, nrow = 0, ncol = 0)),
                               warnings=list(),
                               stype="Level",
                               outcome=1,
                               independent.var.type=2,
                               covar="",
                               continuous=0)

  # Check results
  results <- IntLIM::ProcessResults(inputResults, inputData, pvalcutoff = 1,
		interactionCoeffPercentile = 0, rsquaredCutoff = 0)
  testthat::expect_identical(sort(unlist(inputResults@interaction.adj.pvalues)), sort(results$FDRadjPval))
  testthat::expect_identical(sort(unlist(inputResults@interaction.pvalues)), sort(results$Pval))
  testthat::expect_identical(sort(unlist(inputResults@model.rsquared)), sort(results$rsquared))
  testthat::expect_identical(sort(unlist(inputResults@interaction.coefficients)), sort(results$interaction_coeff))
  inputResults@continuous <- 1
  results <- IntLIM::ProcessResults(inputResults, inputDataCont, pvalcutoff = 1,
                                    interactionCoeffPercentile = 0, rsquaredCutoff = 0)
  testthat::expect_identical(sort(unlist(inputResults@interaction.adj.pvalues)), sort(results$FDRadjPval))
  testthat::expect_identical(sort(unlist(inputResults@interaction.pvalues)), sort(results$Pval))
  testthat::expect_identical(sort(unlist(inputResults@model.rsquared)), sort(results$rsquared))
  testthat::expect_identical(sort(unlist(inputResults@interaction.coefficients)), sort(results$interaction_coeff))
})

# Check coefficient filtering.
testthat::test_that("Check that coefficients are filtered as expected.", {

  # Generate input data (discrete).
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
  inputData <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                            analyteType2=as.matrix(geneData),
                            analyteType1MetaData = metabMetaData,
                            analyteType2MetaData = geneMetaData,
                            sampleMetaData = pData)

  # Generate input data (continuous).
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  inputDataC <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                                analyteType2=as.matrix(geneData),
                                analyteType1MetaData = metabMetaData,
                                analyteType2MetaData = geneMetaData,
                                sampleMetaData = pData)

  # Generate result data.
  pvals <- matrix(cbind(c(0.0005, 0.001, 0.002), c(0.003, 0.004, 0.005), c(0.006, 0.007, 0.008)),
                  nrow = 3, ncol = 3)
  rownames(pvals) <- c("Gene1", "Gene2", "Gene3")
  colnames(pvals) <- c("Metab1", "Metab2", "Metab3")
  adj_pvals <- matrix(cbind(c(0.05, 0.1, 0.2), c(0.3, 0.4, 0.5), c(0.6, 0.7, 0.8)),
                      nrow = 3, ncol = 3)
  rownames(adj_pvals) <- c("Gene1", "Gene2", "Gene3")
  colnames(adj_pvals) <- c("Metab1", "Metab2", "Metab3")
  coef <- matrix(cbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)),
                 nrow = 3, ncol = 3)
  rownames(coef) <- c("Gene1", "Gene2", "Gene3")
  colnames(coef) <- c("Metab1", "Metab2", "Metab3")
  rsq <- matrix(cbind(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6), c(0.7, 0.8, 0.9)),
                nrow = 3, ncol = 3)
  rownames(rsq) <- c("Gene1", "Gene2", "Gene3")
  colnames(rsq) <- c("Metab1", "Metab2", "Metab3")
  inputResults <- methods::new("IntLimResults",
                               interaction.pvalues=pvals,
                               interaction.adj.pvalues=adj_pvals,
                               interaction.coefficients=coef,
                               model.rsquared = rsq,
                               covariate.pvalues = data.frame(matrix(, nrow = 0, ncol = 0)),
                               covariate.coefficients = data.frame(matrix(, nrow = 0, ncol = 0)),
                               corr=data.frame(matrix(, nrow = 0, ncol = 0)),
                               filt.results=data.frame(matrix(, nrow = 0, ncol = 0)),
                               warnings=list(),
                               stype="Level",
                               outcome=1,
                               independent.var.type=2,covar="",
                               continuous=0)

  # Check results (discrete).
  results <- IntLIM::ProcessResults(inputResults, inputData, pvalcutoff = 1,
		interactionCoeffPercentile = 0.7, rsquaredCutoff = 0)
  testthat::expect_equal(max(results$FDRadjPval), 0.8)
  testthat::expect_equal(length(results$FDRadjPval), 3)
  testthat::expect_equal(max(results$Pval), 0.008)
  testthat::expect_equal(length(results$Pval), 3)
  testthat::expect_equal(max(results$interaction_coeff), 9)
  testthat::expect_equal(length(results$interaction_coeff), 3)
  testthat::expect_equal(max(results$rsquared), 0.9)
  testthat::expect_equal(length(results$rsquared), 3)

  # Check results (continuous).
  inputResults@continuous <- 1
  results <- IntLIM::ProcessResults(inputResults, inputDataC, pvalcutoff = 1,
		interactionCoeffPercentile = 0.7, rsquaredCutoff = 0)
  testthat::expect_equal(max(results$FDRadjPval), 0.8)
  testthat::expect_equal(length(results$FDRadjPval), 3)
  testthat::expect_equal(max(results$Pval), 0.008)
  testthat::expect_equal(length(results$Pval), 3)
  testthat::expect_equal(max(results$interaction_coeff), 9)
  testthat::expect_equal(length(results$interaction_coeff), 3)
  testthat::expect_equal(max(results$rsquared), 0.9)
  testthat::expect_equal(length(results$rsquared), 3)
})

# Check p-value filtering.
testthat::test_that("Check that p-values are filtered as expected.", {

  # Generate input data (discrete).
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
  inputData <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                            analyteType2=as.matrix(geneData),
                            analyteType1MetaData = metabMetaData,
                            analyteType2MetaData = geneMetaData,
                            sampleMetaData = pData)
  
  # Generate input data (continuous).
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  inputDataC <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                             analyteType2=as.matrix(geneData),
                             analyteType1MetaData = metabMetaData,
                             analyteType2MetaData = geneMetaData,
                             sampleMetaData = pData)

  # Generate result data.
  pvals <- matrix(cbind(c(0.0005, 0.001, 0.002), c(0.003, 0.004, 0.005), c(0.006, 0.007, 0.008)),
                  nrow = 3, ncol = 3)
  rownames(pvals) <- c("Gene1", "Gene2", "Gene3")
  colnames(pvals) <- c("Metab1", "Metab2", "Metab3")
  adj_pvals <- matrix(cbind(c(0.05, 0.1, 0.2), c(0.3, 0.4, 0.5), c(0.6, 0.7, 0.8)),
                      nrow = 3, ncol = 3)
  rownames(adj_pvals) <- c("Gene1", "Gene2", "Gene3")
  colnames(adj_pvals) <- c("Metab1", "Metab2", "Metab3")
  coef <- matrix(cbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)),
                 nrow = 3, ncol = 3)
  rownames(coef) <- c("Gene1", "Gene2", "Gene3")
  colnames(coef) <- c("Metab1", "Metab2", "Metab3")
  rsq <- matrix(cbind(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6), c(0.7, 0.8, 0.9)),
                nrow = 3, ncol = 3)
  rownames(rsq) <- c("Gene1", "Gene2", "Gene3")
  colnames(rsq) <- c("Metab1", "Metab2", "Metab3")
  inputResults <- methods::new("IntLimResults",
                               interaction.pvalues=pvals,
                               interaction.adj.pvalues=adj_pvals,
                               interaction.coefficients=coef,
                               model.rsquared = rsq,
                               covariate.pvalues = data.frame(matrix(, nrow = 0, ncol = 0)),
                               covariate.coefficients = data.frame(matrix(, nrow = 0, ncol = 0)),
                               corr=data.frame(matrix(, nrow = 0, ncol = 0)),
                               filt.results=data.frame(matrix(, nrow = 0, ncol = 0)),
                               warnings=list(),
                               stype="Level",
                               outcome=1,
                               independent.var.type=2,
                               covar="",
                               continuous=0)

  # Check results (discrete).
  results <- IntLIM::ProcessResults(inputResults, inputData, pvalcutoff = 0.3, 
		interactionCoeffPercentile = 0, rsquaredCutoff = 0)
  testthat::expect_equal(max(results$FDRadjPval), 0.3)
  testthat::expect_equal(length(results$FDRadjPval), 4)
  testthat::expect_equal(max(results$Pval), 0.003)
  testthat::expect_equal(length(results$Pval), 4)
  testthat::expect_equal(max(results$interaction_coeff), 4)
  testthat::expect_equal(length(results$interaction_coeff), 4)
  testthat::expect_equal(max(results$rsquared), 0.4)
  testthat::expect_equal(length(results$rsquared), 4)

  # Check results (continuous).
  inputResults@continuous <- 1
  results <- IntLIM::ProcessResults(inputResults, inputDataC, pvalcutoff = 0.3,
		interactionCoeffPercentile = 0, rsquaredCutoff = 0)
  testthat::expect_equal(max(results$FDRadjPval), 0.3)
  testthat::expect_equal(length(results$FDRadjPval), 4)
  testthat::expect_equal(max(results$Pval), 0.003)
  testthat::expect_equal(length(results$Pval), 4)
  testthat::expect_equal(max(results$interaction_coeff), 4)
  testthat::expect_equal(length(results$interaction_coeff), 4)
  testthat::expect_equal(max(results$rsquared), 0.4)
  testthat::expect_equal(length(results$rsquared), 4)
})

# Check r-squared filtering
testthat::test_that("Check that R-squared values are filtered as expected.", {

  # Generate input data (discrete).
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
  inputData <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                            analyteType2=as.matrix(geneData),
                            analyteType1MetaData = metabMetaData,
                            analyteType2MetaData = geneMetaData,
                            sampleMetaData = pData)
  
  # Generate input data (continuous).
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
  inputDataC <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                             analyteType2=as.matrix(geneData),
                             analyteType1MetaData = metabMetaData,
                             analyteType2MetaData = geneMetaData,
                             sampleMetaData = pData)
  
  # Generate result data.
  pvals <- matrix(cbind(c(0.0005, 0.001, 0.002), c(0.003, 0.004, 0.005), c(0.006, 0.007, 0.008)),
                  nrow = 3, ncol = 3)
  rownames(pvals) <- c("Gene1", "Gene2", "Gene3")
  colnames(pvals) <- c("Metab1", "Metab2", "Metab3")
  adj_pvals <- matrix(cbind(c(0.05, 0.1, 0.2), c(0.3, 0.4, 0.5), c(0.6, 0.7, 0.8)),
                      nrow = 3, ncol = 3)
  rownames(adj_pvals) <- c("Gene1", "Gene2", "Gene3")
  colnames(adj_pvals) <- c("Metab1", "Metab2", "Metab3")
  coef <- matrix(cbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)),
                 nrow = 3, ncol = 3)
  rownames(coef) <- c("Gene1", "Gene2", "Gene3")
  colnames(coef) <- c("Metab1", "Metab2", "Metab3")
  rsq <- matrix(cbind(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6), c(0.7, 0.8, 0.9)),
                nrow = 3, ncol = 3)
  rownames(rsq) <- c("Gene1", "Gene2", "Gene3")
  colnames(rsq) <- c("Metab1", "Metab2", "Metab3")
  inputResults <- methods::new("IntLimResults", 
                               interaction.pvalues=pvals, 
                               interaction.adj.pvalues=adj_pvals,
                               interaction.coefficients=coef, 
                               model.rsquared = rsq, 
                               covariate.pvalues = data.frame(matrix(, nrow = 0, ncol = 0)),
                               covariate.coefficients = data.frame(matrix(, nrow = 0, ncol = 0)), 
                               corr=data.frame(matrix(, nrow = 0, ncol = 0)),
                               filt.results=data.frame(matrix(, nrow = 0, ncol = 0)),
                               warnings=list(),
                               stype="Level", 
                               outcome=1,
                               independent.var.type=2,
                               covar="",
                               continuous=1)

  # Check results (discrete).
  results <- IntLIM::ProcessResults(inputResults, inputData, pvalcutoff = 1,
		interactionCoeffPercentile = 0, rsquaredCutoff = 0.4)
  testthat::expect_equal(max(results$FDRadjPval), 0.8)
  testthat::expect_equal(length(results$FDRadjPval), 6)
  testthat::expect_equal(max(results$Pval), 0.008)
  testthat::expect_equal(length(results$Pval), 6)
  testthat::expect_equal(max(results$interaction_coeff), 9)
  testthat::expect_equal(length(results$interaction_coeff), 6)
  testthat::expect_equal(max(results$rsquared), 0.9)
  testthat::expect_equal(length(results$rsquared), 6)

  # Check results (continuous).
  results <- IntLIM::ProcessResults(inputResults, inputDataC, pvalcutoff = 1,
		interactionCoeffPercentile = 0, rsquaredCutoff = 0.4)
  testthat::expect_equal(max(results$FDRadjPval), 0.8)
  testthat::expect_equal(length(results$FDRadjPval), 6)
  testthat::expect_equal(max(results$Pval), 0.008)
  testthat::expect_equal(length(results$Pval), 6)
  testthat::expect_equal(max(results$interaction_coeff), 9)
  testthat::expect_equal(length(results$interaction_coeff), 6)
  testthat::expect_equal(max(results$rsquared), 0.9)
  testthat::expect_equal(length(results$rsquared), 6)
})