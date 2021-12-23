# If a value other than a data frame is input, an error should be thrown.
test_that("Inputting the wrong class causes early termination.", {
  inputData <- "Hello World"
  expect_error(IntLIM::RunIntLim(inputData, stype = "Feat1"), "Input must be an IntLimData object", 
               ignore.case = TRUE)
})

# If the input has more than 2 classes and is not continuous, an error should be thrown.
test_that("Inputting more than 2 classes causes early termination.", {
  
  # Create toy data.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0),
                      "Level"=c("Low", "Medium", "High", "Medium"))
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
  geneData <- data.frame("Fred"=c(0,0,0), "Wilma"=c(0,0,0), "Pebbles"=c(0,0,0),
                         "Bambam"=c(0,0,0))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("Fred"=c(0,0,0), "Wilma"=c(0,0,0), "Pebbles"=c(0,0,0),
                          "Bambam"=c(0,0,0))
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
  
  expect_error(IntLIM::RunIntLim(dat, stype = "Level"), 
               paste("IntLim currently requires only two categories.  Make sure the column",
                     "Level only has two unique values. Did you mean to set continuous to TRUE?"),
               ignore.case = TRUE)
})

# Check invalid parameters for input or outcome type.
test_that("Invalid parameters cause an error.", {
  
  # Create toy data.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0),
                      "Level"=c("Low", "Medium", "Low", "Medium"))
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
  geneData <- data.frame("Fred"=c(0,0,0), "Wilma"=c(0,0,0), "Pebbles"=c(0,0,0),
                         "Bambam"=c(0,0,0))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("Fred"=c(0,0,0), "Wilma"=c(0,0,0), "Pebbles"=c(0,0,0),
                          "Bambam"=c(0,0,0))
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
  
  expect_error(IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = "Scooby Dooby Doo",
                                 outcome = "Where are you?"), paste("Error! independent.var.type and outcome.type must",
                                                                    "both be either 1 or 2 in RunIntLim."), ignore.case = TRUE)
})

# If the data is missing, an error should be thrown.
test_that("Inputting missing data types causes early termination.", {
  # Create toy data.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0),
                      "Level"=c("Low", "Medium", "Low", "Medium"))
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
  geneData <- data.frame("Fred"=c(0,0,0), "Wilma"=c(0,0,0), "Pebbles"=c(0,0,0),
                         "Bambam"=c(0,0,0))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("Fred"=c(0,0,0), "Wilma"=c(0,0,0), "Pebbles"=c(0,0,0),
                          "Bambam"=c(0,0,0))
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
  expect_error(IntLIM::RunIntLim(dat1, stype = "Level", independent.var.type = 1,
                                 outcome = 2), 
               "One type of analyte data is missing. Cannot run.",
               ignore.case = TRUE)
  expect_error(IntLIM::RunIntLim(dat2, stype = "Level", independent.var.type = 2,
                                 outcome = 1), 
               "One type of analyte data is missing. Cannot run.",
               ignore.case = TRUE)
  
  # Check for missing data when analyte type is same.
  expect_error(IntLIM::RunIntLim(dat1, stype = "Level", independent.var.type = 1,
                                 outcome = 1), "Analyte type 1 is missing. Cannot run.",
               ignore.case = TRUE)
  expect_error(IntLIM::RunIntLim(dat2, stype = "Level", independent.var.type = 2,
                                 outcome = 2), "Analyte type 2 is missing. Cannot run.",
               ignore.case = TRUE)
})

# Check that all fields are populated correctly for discrete data.
test_that("Fields are populated appropriately (discrete).", {
  
  # Create toy data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c("Low", "Medium", "Low", "Medium", "Medium", "Low",
                                "Low", "Medium"))
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                       "Dino", "Hoppy")
  geneData <- data.frame("Fred"=c(46.1,20.2,59.3), "Wilma"=c(11.1,34.2,19.3), 
                         "Pebbles"=c(28.1,71.2,94.3), "Bambam"=c(51.1,91.2,32.3),
                         "Betty"=c(73.1,26.2,40.3), "Barney"=c(91.1,99.2,12.3),
                         "Dino"=c(38.1,44.2,60.3), "Hoppy"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("Fred"=c(60.1,32.2,81.3), "Wilma"=c(68.1,58.2,45.3), 
                          "Pebbles"=c(30.1,61.2,67.3), "Bambam"=c(36.1,7.2,79.3),
                          "Betty"=c(5.1,87.2,91.3), "Barney"=c(5.1,87.2,91.3),
                          "Dino"=c(99.1,10.2,85.3), "Hoppy"=c(51.1,14.2,76.3))
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
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1, suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), sort(c("interaction.adj.pvalues", "interaction.coefficients",
                                                    "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                                                    "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 1, suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 2, suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
})

# Check that all fields are populated correctly for continuous data.
test_that("Fields are populated appropriately (continuous).", {
  
  # Create toy data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                       "Dino", "Hoppy")
  geneData <- data.frame("Fred"=c(46.1,20.2,59.3), "Wilma"=c(11.1,34.2,19.3), 
                         "Pebbles"=c(28.1,71.2,94.3), "Bambam"=c(51.1,91.2,32.3),
                         "Betty"=c(73.1,26.2,40.3), "Barney"=c(91.1,99.2,12.3),
                         "Dino"=c(38.1,44.2,60.3), "Hoppy"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("Fred"=c(60.1,32.2,81.3), "Wilma"=c(68.1,58.2,45.3), 
                          "Pebbles"=c(30.1,61.2,67.3), "Bambam"=c(36.1,7.2,79.3),
                          "Betty"=c(5.1,87.2,91.3), "Barney"=c(5.1,87.2,91.3),
                          "Dino"=c(99.1,10.2,85.3), "Hoppy"=c(51.1,14.2,76.3))
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
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE, continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE, 
                               continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1, suppressWarnings = TRUE,
                               continuous = TRUE)
  expect_identical(sort(slotNames(results)), sort(c("interaction.adj.pvalues", "interaction.coefficients",
                                                    "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                                                    "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 1, suppressWarnings = TRUE,
                               continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 2, suppressWarnings = TRUE,
                               continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results","warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
})

# Check remove.duplicates.
test_that("Duplicates are removed when appropriate (continuous).", {
  
  # Create toy data (discrete).
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                       "Dino", "Hoppy")
  geneData <- data.frame("Fred"=c(46.1,20.2,59.3), "Wilma"=c(11.1,34.2,19.3), 
                         "Pebbles"=c(28.1,71.2,94.3), "Bambam"=c(51.1,91.2,32.3),
                         "Betty"=c(73.1,26.2,40.3), "Barney"=c(91.1,99.2,12.3),
                         "Dino"=c(38.1,44.2,60.3), "Hoppy"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("Fred"=c(60.1,32.2,81.3), "Wilma"=c(68.1,58.2,45.3), 
                          "Pebbles"=c(30.1,61.2,67.3), "Bambam"=c(36.1,7.2,79.3),
                          "Betty"=c(5.1,87.2,91.3), "Barney"=c(5.1,87.2,91.3),
                          "Dino"=c(99.1,10.2,85.3), "Hoppy"=c(51.1,14.2,76.3))
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
  expect_warning(IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, 
                                   outcome = 2, remove.duplicates = TRUE,
                                   continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, 
                               outcome = 1, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  expect_warning(IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, 
                                   outcome = 1, remove.duplicates = TRUE, continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  expect_warning(IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                                   independent.var.type = 1, outcome = 2, remove.duplicates = TRUE, continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  expect_warning(IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                                   independent.var.type = 2, outcome = 1, remove.duplicates = TRUE, continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, 
                               outcome = 1, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  expect_equal(length(which(is.na(results@model.rsquared))), 6)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, 
                               outcome = 2, remove.duplicates = TRUE,
                               suppressWarnings = TRUE, continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  expect_equal(length(which(is.na(results@model.rsquared))), 6)
  
  
  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"), 
                               remove.duplicates = TRUE, independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  expect_equal(length(which(is.na(results@model.rsquared))), 6)
  
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"), 
                               remove.duplicates = TRUE, independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", "corr", "filt.results",
                          "warnings", "stype", "outcome", "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  expect_equal(length(which(is.na(results@model.rsquared))), 6)
})

test_that("Duplicates are removed when appropriate (continuous).", {
  
  # Create toy data (discrete).
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                       "Dino", "Hoppy")
  geneData <- data.frame("Fred"=c(46.1,20.2,59.3), "Wilma"=c(11.1,34.2,19.3), 
                         "Pebbles"=c(28.1,71.2,94.3), "Bambam"=c(51.1,91.2,32.3),
                         "Betty"=c(73.1,26.2,40.3), "Barney"=c(91.1,99.2,12.3),
                         "Dino"=c(38.1,44.2,60.3), "Hoppy"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("Fred"=c(60.1,32.2,81.3), "Wilma"=c(68.1,58.2,45.3), 
                          "Pebbles"=c(30.1,61.2,67.3), "Bambam"=c(36.1,7.2,79.3),
                          "Betty"=c(5.1,87.2,91.3), "Barney"=c(5.1,87.2,91.3),
                          "Dino"=c(99.1,10.2,85.3), "Hoppy"=c(51.1,14.2,76.3))
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
  expect_warning(IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, 
                                   outcome = 2, remove.duplicates = TRUE,
                                   continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, 
                               outcome = 1, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  expect_warning(IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, 
                                   outcome = 1, remove.duplicates = TRUE,
                                   continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  expect_warning(IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                                   independent.var.type = 1, outcome = 2, remove.duplicates = TRUE,
                                   continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  expect_warning(IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                                   independent.var.type = 2, outcome = 1, remove.duplicates = TRUE,
                                   continuous = TRUE),
                 paste("remove.duplicates only applies if the independent variable",
                       "and outcome are of the same analyte type. Duplicates will not be removed."))
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  
  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, 
                               outcome = 1, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  expect_equal(length(which(is.na(results@model.rsquared))), 6)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, 
                               outcome = 2, remove.duplicates = TRUE,
                               continuous = TRUE, suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  expect_equal(length(which(is.na(results@model.rsquared))), 6)
  
  
  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"), 
                               remove.duplicates = TRUE, independent.var.type = 1, outcome = 1,
                               continuous = TRUE, suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  expect_equal(length(which(is.na(results@model.rsquared))), 6)
  
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"), 
                               remove.duplicates = TRUE, independent.var.type = 2, outcome = 2,
                               continuous = TRUE, suppressWarnings = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(length(results@covariate.pvalues), 0)
  expect_equal(length(results@covariate.coefficients), 0)
  expect_equal(length(which(is.na(results@interaction.adj.pvalues))), 6)
  expect_equal(length(which(is.na(results@interaction.coefficients))), 6)
  expect_equal(length(which(is.na(results@model.rsquared))), 6)
})

# Check saving covariate pvals.
test_that("Covariate p-values are saved (discrete).", {
  
  # Create toy data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c("Low", "Medium", "Low", "Medium", "Medium", "Low",
                                "Low", "Medium"))
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                       "Dino", "Hoppy")
  geneData <- data.frame("Fred"=c(46.1,20.2,59.3), "Wilma"=c(11.1,34.2,19.3), 
                         "Pebbles"=c(28.1,71.2,94.3), "Bambam"=c(51.1,91.2,32.3),
                         "Betty"=c(73.1,26.2,40.3), "Barney"=c(91.1,99.2,12.3),
                         "Dino"=c(38.1,44.2,60.3), "Hoppy"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("Fred"=c(60.1,32.2,81.3), "Wilma"=c(68.1,58.2,45.3), 
                          "Pebbles"=c(30.1,61.2,67.3), "Bambam"=c(36.1,7.2,79.3),
                          "Betty"=c(5.1,87.2,91.3), "Barney"=c(5.1,87.2,91.3),
                          "Dino"=c(99.1,10.2,85.3), "Hoppy"=c(51.1,14.2,76.3))
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
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  
  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE,
                               save.covar.pvals = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  
  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  
  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
})

# Check saving covariate pvals.
test_that("Covariate p-values are saved (continuous).", {
  
  # Create toy data.
  pData <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8),
                      "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                       "Dino", "Hoppy")
  geneData <- data.frame("Fred"=c(46.1,20.2,59.3), "Wilma"=c(11.1,34.2,19.3), 
                         "Pebbles"=c(28.1,71.2,94.3), "Bambam"=c(51.1,91.2,32.3),
                         "Betty"=c(73.1,26.2,40.3), "Barney"=c(91.1,99.2,12.3),
                         "Dino"=c(38.1,44.2,60.3), "Hoppy"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("Fred"=c(60.1,32.2,81.3), "Wilma"=c(68.1,58.2,45.3), 
                          "Pebbles"=c(30.1,61.2,67.3), "Bambam"=c(36.1,7.2,79.3),
                          "Betty"=c(5.1,87.2,91.3), "Barney"=c(5.1,87.2,91.3),
                          "Dino"=c(99.1,10.2,85.3), "Hoppy"=c(51.1,14.2,76.3))
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
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  
  # Run with different analytes and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 2, suppressWarnings = TRUE,
                               save.covar.pvals = TRUE, continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  
  # Run with the same analyte type and no covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 36)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 36)
  
  # Run with the same analyte type and covariates.
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 1, outcome = 1,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  expect_equal(ncol(results@covariate.coefficients) * nrow(results@covariate.coefficients), 63)
  
  results <- IntLIM::RunIntLim(dat, stype = "Level", covar = c("Feat1", "Feat2", "Feat3"),
                               independent.var.type = 2, outcome = 2,
                               suppressWarnings = TRUE, save.covar.pvals = TRUE,
                               continuous = TRUE)
  expect_identical(sort(slotNames(results)), 
                   sort(c("interaction.adj.pvalues", "interaction.coefficients",
                          "model.rsquared", "covariate.pvalues", "covariate.coefficients", 
                          "corr", "filt.results", "warnings", "stype", "outcome", 
                          "independent.var.type", "covar", "interaction.pvalues")))
  expect_equal(length(results@interaction.adj.pvalues), 9)
  expect_equal(length(results@interaction.coefficients), 9)
  expect_equal(length(results@model.rsquared), 9)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
  expect_equal(ncol(results@covariate.pvalues) * nrow(results@covariate.pvalues), 63)
})

# TODO: Check whether a pseudoinverse is calculated for a singular matrix.
# XtX must be singular, where X is the model matrix for an individual analyte.
# The model matrix is formed using the formula and the clinical data. The formula is Y ~ a + type + a:type + covar.
# The clinical data includes the outcome type, covariate matrix, and the matrix of analytes used as independent
# variables.

# Check for correlated variables.
test_that("Check for correlated variables.", {
  # Create toy data.
  pData1 <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                       "Feat2"=c(seq(1:8)),
                       "Feat3"=c(seq(1:8)),
                       "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData1) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                        "Dino", "Hoppy")
  pData2 <- data.frame("Feat1"=c(47.1,26.2,84.3,98.4,43.5,82.6,13.7,87.8), 
                       "Feat2"=c(seq(1:8)),
                       "Feat3"=c(rev(seq(1:8))),
                       "Level"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8))
  rownames(pData2) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                        "Dino", "Hoppy")
  geneData <- data.frame("Fred"=c(46.1,20.2,59.3), "Wilma"=c(11.1,34.2,19.3), 
                         "Pebbles"=c(28.1,71.2,94.3), "Bambam"=c(51.1,91.2,32.3),
                         "Betty"=c(73.1,26.2,40.3), "Barney"=c(91.1,99.2,12.3),
                         "Dino"=c(38.1,44.2,60.3), "Hoppy"=c(91.1,93.2,63.3))
  rownames(geneData) <- c("Gene1", "Gene2", "Gene3")
  metabData <- data.frame("Fred"=c(60.1,32.2,81.3), "Wilma"=c(68.1,58.2,45.3), 
                          "Pebbles"=c(30.1,61.2,67.3), "Bambam"=c(36.1,7.2,79.3),
                          "Betty"=c(5.1,87.2,91.3), "Barney"=c(5.1,87.2,91.3),
                          "Dino"=c(99.1,10.2,85.3), "Hoppy"=c(51.1,14.2,76.3))
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
  expect_identical(run1@warnings[[1]], "Using pseudoinverse for Metab1")
  expect_identical(run1@warnings[[2]], paste0("The following covariates have correlation > 0.9 : ",
                                              "(Feat2, Feat3)"))
  run2 <- IntLIM::RunIntLim(dat2, stype = "Level", independent.var.type = 1, 
                            outcome = 2, covar = c("Feat1", "Feat2", "Feat3"),
                            continuous = TRUE)
  expect_identical(run2@warnings[[1]], "Using pseudoinverse for Metab1")
  expect_identical(run2@warnings[[3]], paste0("The following covariates have correlation < -0.9 : ",
                                              "(Feat2, Feat3)"))
})

# Removes genes and metabolites with a standard deviation of zero
test_that("If standard deviation is zero, analytes are removed.", {
  # With and without standard deviation of zero
  pData1 <- data.frame("Feat1"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8), 
                      "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                      "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8))
  rownames(pData1) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                        "Dino", "Hoppy")
  pData2 <- data.frame("Feat1"=c(62.1,44.2,42.3,14.4,58.5,95.6,91.7,1.8), 
                       "Feat2"=c(rep(0,8)),
                       "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8))
  rownames(pData2) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                        "Dino", "Hoppy")
  pData3 <- data.frame("Feat1"=c(rep(0,8)), 
                       "Feat2"=c(37.1,40.2,80.3,83.4,6.5,12.6,43.7,75.8),
                       "Feat3"=c(14.1,74.2,11.3,19.4,73.5,55.6,18.7,91.8))
  rownames(pData3) <- c("Fred", "Wilma", "Pebbles", "Bambam", "Betty", "Barney",
                        "Dino", "Hoppy")
  
  # Genes with and without standard deviation of zero.
  geneData1 <- data.frame("Fred"=c(46.1,20.2,59.3), "Wilma"=c(11.1,34.2,19.3), 
                         "Pebbles"=c(28.1,71.2,94.3), "Bambam"=c(51.1,91.2,32.3),
                         "Betty"=c(73.1,26.2,40.3), "Barney"=c(91.1,99.2,12.3),
                         "Dino"=c(38.1,44.2,60.3), "Hoppy"=c(91.1,93.2,63.3))
  rownames(geneData1) <- c("Gene1", "Gene2", "Gene3")
  geneData2 <- data.frame("Fred"=c(0,20.2,59.3), "Wilma"=c(0,34.2,19.3), 
                          "Pebbles"=c(0,71.2,94.3), "Bambam"=c(0,91.2,32.3),
                          "Betty"=c(0,26.2,40.3), "Barney"=c(0,99.2,12.3),
                          "Dino"=c(0,44.2,60.3), "Hoppy"=c(0,93.2,63.3))
  rownames(geneData2) <- c("Gene1", "Gene2", "Gene3")
  
  # Metabolites with and without standard deviation of zero.
  metabData1 <- data.frame("Fred"=c(0,32.2,81.3), "Wilma"=c(0,58.2,45.3), 
                           "Pebbles"=c(0,61.2,67.3), "Bambam"=c(0,7.2,79.3),
                           "Betty"=c(0,87.2,91.3), "Barney"=c(0,87.2,91.3),
                           "Dino"=c(0,10.2,85.3), "Hoppy"=c(0,14.2,76.3))
  rownames(metabData1) <- c("Metab1", "Metab2", "Metab3")
  metabData2 <- data.frame("Fred"=c(60.1,32.2,81.3), "Wilma"=c(68.1,58.2,45.3), 
                           "Pebbles"=c(30.1,61.2,67.3), "Bambam"=c(36.1,7.2,79.3),
                           "Betty"=c(5.1,87.2,91.3), "Barney"=c(5.1,87.2,91.3),
                           "Dino"=c(99.1,10.2,85.3), "Hoppy"=c(51.1,14.2,76.3))
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
  expect_identical(run1@warnings[[1]], 
                   "Removed 1 analytes of type 1 that had a standard deviation of 0: Metab1")
  expect_equal(length(run1@interaction.adj.pvalues), 6)
  expect_equal(length(run1@interaction.coefficients), 6)
  expect_equal(length(run1@model.rsquared), 6)
  
  # Gene with standard deviation of 0
  run2 <- IntLIM::RunIntLim(dat2, stype = "Feat1", independent.var.type = 1, 
                            outcome = 2, continuous = TRUE, covar = c("Feat2", "Feat3"))
  expect_identical(run2@warnings[[1]], 
                   "Removed 1 analytes of type 2 that had a standard deviation of 0: Gene1")
  expect_equal(length(run2@interaction.adj.pvalues), 6)
  expect_equal(length(run2@interaction.coefficients), 6)
  expect_equal(length(run2@model.rsquared), 6)
  
  # Covariates with standard deviation of 0
  run3 <- IntLIM::RunIntLim(dat3, stype = "Feat1", independent.var.type = 1, 
                            outcome = 2, continuous = TRUE, covar = c("Feat2", "Feat3"))
  expect_identical(run3@warnings[[1]], 
                   "Removed 1 covariates that had a standard deviation of 0: Feat2")
  expect_equal(length(run3@interaction.adj.pvalues), 9)
  expect_equal(length(run3@interaction.coefficients), 9)
  expect_equal(length(run3@model.rsquared), 9)
  expect_equal(length(run3@covar$covariate), 1)
  
  # Phenotype with standard deviation of 0
  expect_error(IntLIM::RunIntLim(dat4, stype = "Feat1", independent.var.type = 1, 
                                 outcome = 2, continuous = TRUE, covar = c("Feat2", "Feat3")),
               "stype variable has a standard deviation of zero. Cannot run.")
})