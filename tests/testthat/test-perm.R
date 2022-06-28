# When input data is not valid, should terminate early.
testthat::test_that("Errors on wrong input type",{
  testthat::expect_error(PermuteIntLIM(data = "?!%", stype="something", outcome = 1,
                             independent.var.type = 2, num.permutations = 4), 
               "The data must be an IntLimData object", ignore.case = TRUE)
})

# When the number of permutations is not sufficient, should terminate early.
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
  testthat::expect_error(PermuteIntLIM(data = dat, stype="Level", outcome = 1,
                             independent.var.type = 2, num.permutations = 0), 
               "The number of permutations must be greater than or equal to 1", 
               ignore.case = TRUE)
})

# Check that the multi-omic function works with meta-data.
testthat::test_that("Function works correctly in multi-omic case.", {
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
  
  # Discrete 
  res <- PermuteIntLIM(data = dat, stype="Level", outcome = 1,
                       independent.var.type = 2, num.permutations = 4)
  testthat::expect_equal(nrow(res[[1]]), 4)
  testthat::expect_equal(ncol(res[[1]]), 2)
  testthat::expect_equal(length(res[[2]]), 4)
  
  # Continuous
  res <- PermuteIntLIM(data = dat, stype="Feat1", outcome = 1,
                       independent.var.type = 2, num.permutations = 4,
                       continuous = TRUE)
  testthat::expect_equal(nrow(res[[1]]), 4)
  testthat::expect_equal(ncol(res[[1]]), 2)
  testthat::expect_equal(length(res[[2]]), 4)
})

# Check that we are still able to run without the metadata.
testthat::test_that("Function still works when metadata is missing.", {

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
  
  # Discrete
  res <- PermuteIntLIM(data = dat, stype="Level", outcome = 1,
                       independent.var.type = 2, num.permutations = 4)
  testthat::expect_equal(nrow(res[[1]]), 4)
  testthat::expect_equal(ncol(res[[1]]), 2)
  testthat::expect_equal(length(res[[2]]), 4)
  
  # Continuous
  res <- PermuteIntLIM(data = dat, stype="Feat1", outcome = 1,
                       independent.var.type = 2, num.permutations = 4,
                       continuous = TRUE)
  testthat::expect_equal(nrow(res[[1]]), 4)
  testthat::expect_equal(ncol(res[[1]]), 2)
  testthat::expect_equal(length(res[[2]]), 4)
})

# Check that we are able to run with both single-omic and multi-omic data.
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
  
  # Discrete
  res <- PermuteIntLIM(data = dat, stype="Level", outcome = 1,
                       independent.var.type = 1, num.permutations = 4)
  testthat::expect_equal(nrow(res[[1]]), 4)
  testthat::expect_equal(ncol(res[[1]]), 2)
  testthat::expect_equal(length(res[[2]]), 4)
  
  # Continuous
  res <- PermuteIntLIM(data = dat, stype="Feat1", outcome = 1,
                       independent.var.type = 1, num.permutations = 4,
                       continuous = TRUE)
  testthat::expect_equal(nrow(res[[1]]), 4)
  testthat::expect_equal(ncol(res[[1]]), 2)
  testthat::expect_equal(length(res[[2]]), 4)
})