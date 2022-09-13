# If a data set is input that is not of the IntLimData type, the program
# should terminate.
testthat::test_that("Inputting data other than a IntLimData causes early termination.", {
  testthat::expect_error(IntLIM::FilterData("?!%"), "input data is not a IntLimData class", ignore.case = TRUE)
})
# If a data set does not contain either expression data or metabolite data,
# the program should terminate.
testthat::test_that("Inputting data set without expression or metabolite data causes early termination.", {
  # sample metadata
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  
  message <- paste("input data must contain assayData of at least one type of analyte.
	     Try reading in the data with the ReadData function")
  dat <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                               analyteType2=matrix(, nrow = 0, ncol = 0),
                               analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               analyteType2MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               sampleMetaData = pData)
  testthat::expect_error(IntLIM::FilterData(dat, message, ignore.case = TRUE))
})

# If parameters are out of bounds, the program should terminate.
testthat::test_that("Inputting out-of-bounds parameters causes early termination.", {

  # Create toy data set.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
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
							   
  # Check for errors.
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = -1, analyteType2perc = 0, analyteMiss = 0), 
               "analyteType1perc parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 100, analyteType2perc = 0, analyteMiss = 0), 
               "analyteType1perc parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = -1, analyteMiss = 0), 
               "analyteType2perc parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 100, analyteMiss = 0), 
               "analyteType2perc parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, analyteMiss = -1), 
               "analyteMiss parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, analyteMiss = 100), 
               "analyteMiss parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_warning(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, analyteMiss = 0), 
               "No filtering parameters were set so the data remains unfiltered", 
               ignore.case = TRUE)
			   
  # Perform the same checks for single-analyte data.
  dat <- methods::new("IntLimData", analyteType1=as.matrix(geneData),
                               analyteType2=matrix(, nrow = 0, ncol = 0),
                               analyteType1MetaData = geneMetaData,
                               analyteType2MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               sampleMetaData = pData) 
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = -1, analyteType2perc = 0, analyteMiss = 0), 
               "analyteType1perc parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 100, analyteType2perc = 0, analyteMiss = 0), 
               "analyteType1perc parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_warning(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, analyteMiss = 0), 
               "No filtering parameters were set so the data remains unfiltered", 
               ignore.case = TRUE)
			   
  # Perform the same checks for single-analyte data.
  dat <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                               analyteType2=as.matrix(metabData),
                               analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               analyteType2MetaData = metabMetaData,
                               sampleMetaData = pData) 
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = -1, analyteMiss = 0), 
               "analyteType2perc parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 100, analyteMiss = 0), 
               "analyteType2perc parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, analyteMiss = -1), 
               "analyteMiss parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, analyteMiss = 100), 
               "analyteMiss parameter must be between 0 and 1", ignore.case = TRUE)
  testthat::expect_warning(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, analyteMiss = 0), 
               "No filtering parameters were set so the data remains unfiltered", 
               ignore.case = TRUE)
})

# If all data is in the correct format and no filtering is done, check that
# all fields are populated appropriately.
testthat::test_that("Check that output is correct without filtering.", {

  # Create toy data set.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
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

  # Check for warning.
  testthat::expect_warning(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, analyteMiss = 0),
               "No filtering parameters were set so the data remains unfiltered",
               ignore.case = TRUE)
  testthat::expect_identical(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                      analyteMiss = 0, suppressWarnings = TRUE),dat)

  # Create single-omic toy data and check for warning.
  dat <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                               analyteType2=as.matrix(geneData),
                               analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               analyteType2MetaData = geneMetaData,
                               sampleMetaData = pData)
  testthat::expect_warning(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, analyteMiss = 0),
               "No filtering parameters were set so the data remains unfiltered",
               ignore.case = TRUE)
  filtdata <- IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                 analyteMiss = 0, suppressWarnings = TRUE)
  testthat::expect_identical(filtdata,dat)

  # Create single-omic toy data and check for warning.
  dat <- methods::new("IntLimData", analyteType1=as.matrix(metabData),
                               analyteType2=matrix(, nrow = 0, ncol = 0),
                               analyteType1MetaData = metabData,
                               analyteType2MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               sampleMetaData = pData)
  testthat::expect_warning(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, analyteMiss = 0),
               "No filtering parameters were set so the data remains unfiltered",
               ignore.case = TRUE)
  filtdata <- IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                 analyteMiss = 0, suppressWarnings = TRUE)
  testthat::expect_identical(filtdata,dat)
})

# Check that missing metabolites are filtered properly.
testthat::test_that("Check that output is correct when we remove metabolites with NA.", {

  # Create toy data set.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  geneData <- data.frame("P1"=c(0,0,0), "P2"=c(0,0,0), "P3"=c(0,0,0), 
                         "P4"=c(0,0,0))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(0,0,0), "P2"=c(0,NA,0), "P3"=c(NA,0,0), 
                          "P4"=c(0,0,0))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(geneData),
                               analyteType2=as.matrix(metabData),
                               analyteType1MetaData = geneMetaData,
                               analyteType2MetaData = metabMetaData,
                               sampleMetaData = pData)
  # Check filtering.
  filtdata <- IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                 analyteMiss = 0.25, suppressWarnings = TRUE)
  testthat::expect_identical(slotNames(filtdata), slotNames(dat))
  testthat::expect_identical(filtdata@analyteType1, as.matrix(geneData))
  testthat::expect_identical(rownames(filtdata@analyteType2), "Metab3")
  testthat::expect_identical(filtdata@sampleMetaData, pData)

  # Create toy data set with only analyte 2 and check filtering.
  dat <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                               analyteType2=as.matrix(metabData),
                               analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               analyteType2MetaData = metabMetaData,
                               sampleMetaData = pData)
  filtdata <- IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                 analyteMiss = 0.25, suppressWarnings = TRUE)
  testthat::expect_identical(slotNames(filtdata), slotNames(dat))
  testthat::expect_identical(filtdata@analyteType1, matrix(, nrow = 0, ncol = 0))
  testthat::expect_identical(rownames(filtdata@analyteType2), "Metab3")
  testthat::expect_identical(filtdata@sampleMetaData, pData)
})

# Check coefficient of variation filtering.
testthat::test_that("Check that metabolites are filtered properly when cov filtering is added.", {
  
  # Create toy data set.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  geneData <- data.frame("P1"=c(0,0,0), "P2"=c(0,0,0), "P3"=c(0,0,0), 
                         "P4"=c(0,0,0))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(0,0,0), "P2"=c(0,NA,0), "P3"=c(NA,0,0), 
                          "P4"=c(0,0,0))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(geneData),
                      analyteType2=as.matrix(metabData),
                      analyteType1MetaData = geneMetaData,
                      analyteType2MetaData = metabMetaData,
                      sampleMetaData = pData)
  # Check filtering.
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                 analyteMiss = 0.25, cov.cutoff = 0.3, suppressWarnings = TRUE),
               "All analytes have been removed from your type 2 data! Change your filtering criteria.")
  
  # Create toy data set with only analyte 2 and check filtering.
  dat <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                      analyteType2=as.matrix(metabData),
                      analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                      analyteType2MetaData = metabMetaData,
                      sampleMetaData = pData)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                 analyteMiss = 0.25, cov.cutoff = 0.3,
                                 suppressWarnings = TRUE),
               "All analytes have been removed from your type 2 data! Change your filtering criteria.")
})

# Check that warnings are thrown when coefficient of variation is applied to log-scaled data.
testthat::test_that("Check for warning when data is log-scaled.", {
  
  # Create toy data set.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  geneData <- data.frame("P1"=c(-1,1,2), "P2"=c(1,2,3), "P3"=c(1,2,3), 
                         "P4"=c(1,2,3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(1,2,3), "P2"=c(1,2,3), "P3"=c(1,2,3), 
                          "P4"=c(1,2,-3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(geneData),
                      analyteType2=as.matrix(metabData),
                      analyteType1MetaData = geneMetaData,
                      analyteType2MetaData = metabMetaData,
                      sampleMetaData = pData)
  # Check filtering.
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                    analyteMiss = 0, cov.cutoff = -1, suppressWarnings = TRUE),
               "cov.cutoff parameter must be between 0 and 1")
  testthat::expect_warning(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                  analyteMiss = 0, cov.cutoff = 0.3),
                 paste("Coefficient of variation filtering will not be applied",
                 "to analyte type 1 because data is log-scaled"))
  
  # Create toy data set with only analyte 2 and check filtering.
  dat <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                      analyteType2=as.matrix(metabData),
                      analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                      analyteType2MetaData = metabMetaData,
                      sampleMetaData = pData)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                  analyteMiss = 0, cov.cutoff = 2, suppressWarnings = TRUE),
               "cov.cutoff parameter must be between 0 and 1")
  testthat::expect_warning(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                  analyteMiss = 0, cov.cutoff = 0.3),
               paste("Coefficient of variation filtering will not be applied",
                     "to analyte type 2 because data is log-scaled"))
  
  # Create toy data set with only analyte 1 and check filtering.
  dat <- methods::new("IntLimData", analyteType2=matrix(, nrow = 0, ncol = 0),
                      analyteType1=as.matrix(geneData),
                      analyteType2MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                      analyteType1MetaData = geneMetaData,
                      sampleMetaData = pData)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                  analyteMiss = 0, cov.cutoff = 2, suppressWarnings = TRUE),
               "cov.cutoff parameter must be between 0 and 1")
  testthat::expect_warning(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0, 
                                    analyteMiss = 0, cov.cutoff = 0.3),
                 paste("Coefficient of variation filtering will not be applied",
                       "to analyte type 1 because data is log-scaled"))
})

# Check that lowest percentile analytes are filtered.
testthat::test_that("Check that output is correct when we filter the lowest percentile analytes.", {

  # Create toy data set.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  geneData <- data.frame("P1"=c(1,2,3), "P2"=c(1,2,3), "P3"=c(1,2,3), 
                         "P4"=c(1,2,3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(1,2,3), "P2"=c(1,2,3), "P3"=c(1,2,3), 
                          "P4"=c(1,2,3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(geneData),
                               analyteType2=as.matrix(metabData),
                               analyteType1MetaData = geneMetaData,
                               analyteType2MetaData = metabMetaData,
                               sampleMetaData = pData)

  # Check filtering.
  filtdata <- IntLIM::FilterData(dat, analyteType1perc = 0.4, analyteType2perc = 0.4, 
                                 analyteMiss = 0, suppressWarnings = TRUE)
  testthat::expect_identical(slotNames(filtdata), slotNames(dat))
  testthat::expect_identical(rownames(filtdata@analyteType1), c("Gene2", "Gene3"))
  testthat::expect_identical(rownames(filtdata@analyteType2), c("Metab2", "Metab3"))
  testthat::expect_identical(filtdata@sampleMetaData, pData)

  # Create single-omic toy dataset and check filtering.
  dat <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                               analyteType2=as.matrix(metabData),
                               analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               analyteType2MetaData = metabMetaData,
                               sampleMetaData = pData)
  filtdata <- IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0.4, 
                                 analyteMiss = 0, suppressWarnings = TRUE)
  testthat::expect_identical(slotNames(filtdata), slotNames(dat))
  testthat::expect_identical(filtdata@analyteType1, matrix(, nrow = 0, ncol = 0))
  testthat::expect_identical(rownames(filtdata@analyteType2), c("Metab2", "Metab3"))
  testthat::expect_identical(filtdata@sampleMetaData, pData)

  # Create single-omic toy dataset and check filtering.
  dat <- methods::new("IntLimData", analyteType1= as.matrix(geneData),
                               analyteType2=matrix(, nrow = 0, ncol = 0),
                               analyteType1MetaData = geneMetaData,
                               analyteType2MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               sampleMetaData = pData)
  filtdata <- IntLIM::FilterData(dat, analyteType1perc = 0.4, analyteType2perc = 0, 
                                 analyteMiss = 0, suppressWarnings = TRUE)
  testthat::expect_identical(slotNames(filtdata), slotNames(dat))
  testthat::expect_identical(filtdata@analyteType2, matrix(, nrow = 0, ncol = 0))
  testthat::expect_identical(rownames(filtdata@analyteType1), c("Gene2", "Gene3"))
  testthat::expect_identical(filtdata@sampleMetaData, pData)
})

# Check that lowest percentile analytes are filtered and NA is removed.
testthat::test_that("Check that output is correct when we filter the lowest percentile analytes and remove NA.", {

  # Create toy data set.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  geneData <- data.frame("P1"=c(1,2,3), "P2"=c(1,2,3), "P3"=c(1,2,3), 
                         "P4"=c(1,2,3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(1,2,3), "P2"=c(1,2,3), "P3"=c(1,2,3), 
                          "P4"=c(1,NA,3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(geneData),
                               analyteType2=as.matrix(metabData),
                               analyteType1MetaData = geneMetaData,
                               analyteType2MetaData = metabMetaData,
                               sampleMetaData = pData)

  # Check filtering.
  filtdata <- IntLIM::FilterData(dat, analyteType1perc = 0.4, analyteType2perc = 0.4, 
                                 analyteMiss = 0.25, suppressWarnings = TRUE)
  testthat::expect_identical(slotNames(filtdata), slotNames(dat))
  testthat::expect_identical(rownames(filtdata@analyteType1), c("Gene2", "Gene3"))
  testthat::expect_identical(rownames(filtdata@analyteType2), "Metab3")
  testthat::expect_identical(filtdata@sampleMetaData, pData)

  # Create single-omic toy dataset and check filtering.
  dat <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                               analyteType2=as.matrix(metabData),
                               analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               analyteType2MetaData = metabMetaData,
                               sampleMetaData = pData)
  filtdata <- IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0.4, 
                                 analyteMiss = 0.25, suppressWarnings = TRUE)
  testthat::expect_identical(slotNames(filtdata), slotNames(dat))
  testthat::expect_identical(filtdata@analyteType1, matrix(, nrow = 0, ncol = 0))
  testthat::expect_identical(rownames(filtdata@analyteType2), "Metab3")
  testthat::expect_identical(filtdata@sampleMetaData, pData)
})

# Check that an error is given if all data is filtered out.
testthat::test_that("Check that an error is given if all data is filtered out.", {

  # Create toy data set.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  geneData <- data.frame("P1"=c(1,2,3), "P2"=c(1,2,3), "P3"=c(1,2,3), 
                         "P4"=c(1,2,3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("P1"=c(1,2,3), "P2"=c(1,2,NA), "P3"=c(1,2,3), 
                          "P4"=c(1,NA,3))
  rownames(metabData) <- c("Metab1", "Metab2", "Metab3")
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  dat <- methods::new("IntLimData", analyteType1=as.matrix(geneData),
                               analyteType2=as.matrix(metabData),
                               analyteType1MetaData = geneMetaData,
                               analyteType2MetaData = metabMetaData,
                               sampleMetaData = pData)

  # Check filtering.
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0.4, 
                                    analyteMiss = 0.25, suppressWarnings = TRUE),
                 "All analytes have been removed from your type 2 data! Change your filtering criteria.")

  # Create single-omic toy dataset and check filtering.
  dat <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                               analyteType2=as.matrix(metabData),
                               analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               analyteType2MetaData = metabMetaData,
                               sampleMetaData = pData)
  testthat::expect_error(IntLIM::FilterData(dat, analyteType1perc = 0, analyteType2perc = 0.4, 
                                  analyteMiss = 0.25, suppressWarnings = TRUE),
	"All analytes have been removed from your type 2 data! Change your filtering criteria.")
})
