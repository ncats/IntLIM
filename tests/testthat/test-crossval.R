# When input data is not valid, should terminate early.
testthat::test_that("Errors on wrong input type",{
  testthat::expect_error(RunCrossValidation(inputData = "?!%", folds = 4, suppressWarnings=TRUE), 
               "input data is not a IntLimData class", ignore.case = TRUE)
})

# When the number of folds is not sufficient, should terminate early.
testthat::test_that("Insufficient fold count causes an error",{
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
  testthat::expect_error(RunCrossValidation(inputData = dat, folds = 1, suppressWarnings=TRUE), 
               "At least 2 folds are required.", ignore.case = TRUE)
})

# When there are more folds than samples, building data folds should terminate early.
testthat::test_that("Inputting too many folds causes early termination.", {
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
  testthat::expect_error(RunCrossValidation(inputData = dat, folds = 9, suppressWarnings=TRUE), 
               "The number of folds is greater than the number of samples!", ignore.case = TRUE)
})

# Check that the multi-omic function works with meta-data.
testthat::test_that("Function returns all folds correctly in multi-omic case.", {
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
  
  res <- RunCrossValidation(inputData = dat, folds = 4, suppressWarnings=TRUE,
                            independent.var.type = 2, outcome = 1, pvalcutoff = 1,
                            interactionCoeffPercentile = 0, rsquaredCutoff = 0.4,
                            stype = "Level")
  
  testthat::expect_equal(length(res$folds), 4)
})

# Check that we are still able to run without the metadata.
testthat::test_that("Function still returns all folds when metadata is missing.", {

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
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                      analyteType2=as.matrix(geneData),
                      analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                      analyteType2MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                      sampleMetaData = pData)
  
  res <- RunCrossValidation(inputData = dat, folds = 4, suppressWarnings=TRUE,
                            independent.var.type = 2, outcome = 1, pvalcutoff = 1, 
                            interactionCoeffPercentile = 0, rsquaredCutoff = 0.4,
                            stype = "Level")
  
  testthat::expect_equal(length(res$folds), 4)
})

# Check that we are able to build folds with both single-omic and multi-omic data.
testthat::test_that("Single-omic data gives expected results.", {
  # Create toy data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c("Low", "Medium", "Low", "Medium", "Medium", "Low",
                                "Low", "Medium"))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P5", "P6",
                       "P7", "P8")
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
                      analyteType2=matrix(, nrow = 0, ncol = 0),
                      analyteType1MetaData = metabMetaData,
                      analyteType2MetaData = geneMetaData,
                      sampleMetaData = pData)
  
  res <- RunCrossValidation(inputData = dat, folds = 4, suppressWarnings=TRUE,
                            independent.var.type = 1, outcome = 1, pvalcutoff = 1, 
                            interactionCoeffPercentile = 0, rsquaredCutoff = 0.4,
                            stype = "Level")
  
  testthat::expect_equal(length(res$folds), 4)
})

# If data is not divisible by fold number
testthat::test_that("When number of samples is not divisible by fold count, make sure the samples
          are evenly divided", {
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
    
    res <- RunCrossValidation(inputData = dat, folds = 3, suppressWarnings=TRUE,
                              independent.var.type = 2, outcome = 1, pvalcutoff = 1, 
                              interactionCoeffPercentile = 0, rsquaredCutoff = 0.4,
                              stype = "Level")
    
    testthat::expect_equal(length(res$folds), 3)
    testthat::expect_lte(nrow(res$folds[[1]]$training@sampleMetaData), 6)
    testthat::expect_lte(nrow(res$folds[[2]]$training@sampleMetaData), 6)
    testthat::expect_lte(nrow(res$folds[[3]]$training@sampleMetaData), 6)
    testthat::expect_lte(nrow(res$folds[[1]]$testing@sampleMetaData), 3)
    testthat::expect_lte(nrow(res$folds[[2]]$testing@sampleMetaData), 3)
    testthat::expect_lte(nrow(res$folds[[3]]$testing@sampleMetaData), 3)
})