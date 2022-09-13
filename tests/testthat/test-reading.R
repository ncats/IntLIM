# If a file is passed as input that does not exist in the system, an error
# should be thrown.
testthat::test_that("Nonexistent files cause early termination.", {
  testthat::expect_error(IntLIM::ReadData("?!%"), "CSV input file does not exist", ignore.case = TRUE)
})

# If a file with inappropriate column names is passed as input, an error
# should be thrown.
testthat::test_that("Incorrect column names cause early termination.", {

  # Create a data frame with inappropriate column names.
  incorrect_colname_df = data.frame("P2"=c(1,2,3), "P1"=c(0,0,0))

  # Save the data frame as a file.
  fname <- paste(getwd(), "incorrect_colname_file.csv", sep = "/")
  write.csv(incorrect_colname_df, file=fname, quote=FALSE, row.names = FALSE)

  # Check for error.
  testthat::expect_error(IntLIM::ReadData(fname),
               "Check column names of input files.  'type' and 'filenames' are required",
               ignore.case = TRUE)
  file.remove(fname)
})

# If one or more data type headings is missing, an error should be thrown.
testthat::test_that("Missing data types cause early termination.", {

  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  # Check when each type is missing.
  for(heading in expected){

    # Create a file with one missing heading.
    missing_heading_df = data.frame("filenames"=c(0,0,0,0))
    rownames(missing_heading_df) <- setdiff(expected, heading)
    fname <- paste(getwd(), "missing_heading_file.csv", sep = "/")
    write.csv(missing_heading_df, file=fname, quote=FALSE)

    # Check that error occurs.
    message <- paste("The column 'type' contains non-allowed entries (See Description). The",
                     "CSV input file must contain 6 rows (if optional meta data files for analytes",
                     "are not to be input, have the corresponding filenames be blanks.")
    testthat::expect_error(IntLIM::ReadData(fname), message, fixed = TRUE)
    file.remove(fname)
  }
})

# If no data is included, an error should be thrown.
testthat::test_that("Absence of data causes early termination.", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  # Create a file with both gene and metabolite data missing.
  missing_data_df = data.frame("filenames"=c("","",0,0,0))
  rownames(missing_data_df) <- expected
  fname <- paste(getwd(), "missing_data_file.csv", sep = "/")
  write.csv(missing_data_df, file=fname, quote=FALSE)

  # Check that error occurs.
  testthat::expect_error(IntLIM::ReadData(fname), "No data provided.", ignore.case = TRUE)
  file.remove(fname)
})

# If only one data type is missing, a warning should be given.
# However, fields should still be populated appropriately.
testthat::test_that("Absence of one omic type results in a warning.", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  #Create metabolite data file.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  # Create metabolite metadata file.
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  # Create gene data file.
  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  # Create gene metadata file.
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                                c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  # Create patient data file.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  # Create a file with only metabolite data missing.
  missing_metab_df = data.frame("filenames"=c("","gene_file.csv",
                                              "metab_metadata_file.csv",
                                              "gene_metadata_file.csv",
                                              "pdata_file.csv"))
  rownames(missing_metab_df) <- expected
  fname <- paste(getwd(), "missing_metabs_file.csv", sep = "/")
  write.csv(missing_metab_df, file = fname, quote=FALSE)

  # Check that warning occurs.
  message <- paste("No data provided for Analyte Type 1. This means you cannot run",
                   "analyses involving this analyte type.")
  testthat::expect_warning(IntLIM::ReadData(fname, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric")),
                 message, ignore.case = TRUE)

  # Check that fields are still populated appropriately.
  data_nometab <- IntLIM::ReadData(fname, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"),
                                 suppressWarnings = TRUE)
  expected_names <- c("analyteType1", "analyteType2", "analyteType1MetaData",
                      "analyteType2MetaData", "sampleMetaData")
  testthat::expect_identical(slotNames(data_nometab), expected_names)
  testthat::expect_identical(colnames(data_nometab@sampleMetaData),
                   c("Feat1", "Feat2", "Feat3"))
  testthat::expect_identical(rownames(data_nometab@analyteType2),
                   c("Gene1", "Gene2","Gene3"))
  testthat::expect_identical(rownames(data_nometab@analyteType2MetaData),
                   c("Gene1", "Gene2","Gene3"))
  testthat::expect_equal(length(data_nometab@analyteType2), 12)
  testthat::expect_equal(length(data_nometab@analyteType1),0)

  # Create a file with only gene data missing.
  missing_gene_df = data.frame("filenames"=c("metab_file.csv","",
                                             "metab_metadata_file.csv",
                                             "gene_metadata_file.csv",
                                             "pdata_file.csv"))
  rownames(missing_gene_df) <- expected
  fname <- paste(getwd(), "missing_genes_file.csv", sep = "/")
  write.csv(missing_gene_df, file = fname, quote=FALSE)

  # Check that warning occurs.
  message <- paste("No data provided for Analyte Type 2. This means you cannot run",
                   "analyses involving this analyte type.")
  testthat::expect_warning(IntLIM::ReadData(fname, class.feat = list(Feat1 = "numeric",
                                                         Feat2 = "numeric",
                                                         Feat3 = "numeric")),
                 message, ignore.case = TRUE)

  # Check that fields are still populated appropriately.
  data_nogene <- IntLIM::ReadData(fname, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"),
                                  suppressWarnings = TRUE)
  expected_names <- c("analyteType1", "analyteType2", "analyteType1MetaData",
                      "analyteType2MetaData", "sampleMetaData")
  testthat::expect_identical(slotNames(data_nogene), expected_names)
  testthat::expect_identical(colnames(data_nogene@sampleMetaData),
                   c("Feat1", "Feat2", "Feat3"))
  testthat::expect_identical(rownames(data_nogene@analyteType1),
                   c("Metab1", "Metab2","Metab3"))
  testthat::expect_identical(rownames(data_nogene@analyteType1MetaData),
                   c("Metab1", "Metab2","Metab3"))
  testthat::expect_equal(length(data_nogene@analyteType1),12)
  testthat::expect_equal(length(data_nogene@analyteType2),0)

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(fname)
})

# If there are two analytes with the same name, an error should be thrown.
testthat::test_that("Duplicate analytes cause early termination", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                            "analyteType2MetaData","sampleMetaData")

  #Create metabolite data file.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  # Create metabolite data file with duplicates.
  metabData <- data.frame("1"=c("Metab1",0,0,0,0), "2"=c("Metab2",0,0,0,0),
                          "3"=c("Metab2",0,0,0,0))
  rownames(metabData) <- c("","P1", "P2", "P3", "P4")
  fname_metab_dup <- paste(getwd(), "metab_file_dup.csv", sep = "/")
  write.table(t(metabData), file = fname_metab_dup, quote=FALSE, row.names = FALSE,
              sep = ",")

  # Create metabolite metadata file.
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  # Create metabolite metadata file with duplicates.
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab2"), "metabname"=
                                c("Metab1", "Metab2", "Metab2"))
  fname_metab_meta_dup <- paste(getwd(), "metab_metadata_file_dup.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta_dup, quote=FALSE, row.names = FALSE)

  # Create gene data file.
  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  # Create gene data file with duplicates.
  geneData <- data.frame("1"=c("Gene1",0,0,0,0), "2"=c("Gene2",0,0,0,0),
                         "3"=c("Gene2",0,0,0,0))
  rownames(geneData) <- c("","P1", "P2", "P3", "P4")
  fname_gene_dup <- paste(getwd(), "gene_file_dup.csv", sep = "/")
  write.table(t(geneData), file = fname_gene_dup, quote=FALSE, row.names = FALSE,
              sep = ",")

  # Create gene metadata file.
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  # Create gene metadata file with duplicates.
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene2"), "genename"=
                               c("Gene1", "Gene2", "Gene2"))
  fname_gene_meta_dup <- paste(getwd(), "gene_metadata_file_dup.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta_dup, quote=FALSE, row.names = FALSE)

  # Create patient data file.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  # Create file with duplicated metabolite reference.
  metab_dup_df = data.frame("filenames"=c("metab_file_dup.csv","gene_file.csv",
                                              "metab_metadata_file_dup.csv"
                                              ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(metab_dup_df) <- expected
  fname_metab_dup_ref <- paste(getwd(), "metab_file_dup_ref.csv", sep = "/")
  write.csv(metab_dup_df, file=fname_metab_dup_ref, quote=FALSE)

  # Create file with duplicated gene reference.
  gene_dup_df = data.frame("filenames"=c("metab_file.csv","gene_file_dup.csv",
                                         "metab_metadata_file.csv"
                                         ,"gene_metadata_file_dup.csv","pdata_file.csv"))
  rownames(gene_dup_df) <- expected
  fname_gene_dup_ref <- paste(getwd(), "gene_file_dup_ref.csv", sep = "/")
  write.csv(gene_dup_df, file=fname_gene_dup_ref, quote=FALSE)

  # Check that error is thrown.
  fname_metab_dup <- paste(getwd(), "metab_file_dup.csv", sep = "/")
  testthat::expect_error(IntLIM::ReadData(fname_metab_dup_ref),
               paste("Error: your input file",fname_metab_dup,"has duplicate",
                     "entries in column 1. Please make sure you have one row per",
                     "analyte"), fixed = TRUE)

  # Check that error is thrown.
  fname_gene_dup <- paste(getwd(), "gene_file_dup.csv", sep = "/")
  testthat::expect_error(IntLIM::ReadData(fname_gene_dup_ref),
                paste("Error: your input file",fname_gene_dup,"has duplicate",
                     "entries in column 1. Please make sure you have one row per analyte"),
               fixed = TRUE)

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_dup)
  file.remove(fname_metab_meta)
  file.remove(fname_metab_meta_dup)
  file.remove(fname_gene)
  file.remove(fname_gene_dup)
  file.remove(fname_gene_meta)
  file.remove(fname_gene_meta_dup)
  file.remove(fname_pdata)
  file.remove(fname_metab_dup_ref)
  file.remove(fname_gene_dup_ref)
})

# If a data file is not accessible on the system, an error should be thrown.
testthat::test_that("Inaccessible data files cause early termination.", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  # Create metabolite file.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  # Create metabolite metadata file.
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE)

  # Create gene file.
  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  # Create gene metadata file.
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE)

  # Create patient data file.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  # Create reference file that points to nonexistent metabolite file.
  missing_metab_df = data.frame("filenames"=c("?!%","gene_file.csv","metab_metadata_file.csv"
                                              ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(missing_metab_df) <- expected
  fname_missing_metab <- paste(getwd(), "missing_metabdata_file.csv", sep = "/")
  write.csv(missing_metab_df, file = fname_missing_metab, quote=FALSE)

  # Check that an error is thrown.
  testthat::expect_error(IntLIM::ReadData("missing_metabdata_file.csv"),
                 paste("File", paste0(base::dirname("missing_metabdata_file.csv"), "/?!%"),
                       "does not exist"), fixed = TRUE)

  # Create reference file that points to nonexistent gene file.
  missing_gene_df = data.frame("filenames"=c("metab_file.csv","?!%", "metabolite_metadata_file.csv"
                                             ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(missing_gene_df) <- expected
  fname_missing_gene <- paste(getwd(), "missing_gene_file.csv", sep = "/")
  write.csv(missing_gene_df, file = fname_missing_gene, quote=FALSE)

  # Check that an error is thrown.
  testthat::expect_error(IntLIM::ReadData("missing_gene_file.csv"),
                 paste("File", paste0(base::dirname("missing_gene_file.csv"), "/?!%"),
                       "does not exist"), fixed = TRUE)

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(fname_missing_metab)
  file.remove(fname_missing_gene)
})

# If analyte metadata files are missing, this should throw a warning but not lead
# to early termination.
# Fields should still be populated.
testthat::test_that("Missing metadata causes a warning.", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                            "analyteType2MetaData","sampleMetaData")

  # Save metabolite file.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  # Save gene file.
  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  # Save patient data file.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  # Save metadata files.
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                                c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  # Check when only metabolite data is missing.
  missing_metab_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv", "",
                                              "gene_metadata_file.csv","pdata_file.csv"))
  rownames(missing_metab_df) <- expected
  fname_metab_missing <- paste(getwd(), "missing_metab_metadata_file.csv", sep = "/")
  write.csv(missing_metab_df, file = fname_metab_missing, quote=FALSE)
  testthat::expect_warning(IntLIM::ReadData(fname_metab_missing, class.feat = list(Feat1 = "numeric",
                                                                       Feat2 = "numeric",
                                                                       Feat3 = "numeric")),
                 "No metadata provided for Analyte Type 1", fixed = TRUE)
  # Check that fields are still populated appropriately.
  data_metab <- IntLIM::ReadData(fname_metab_missing, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"),
                                 suppressWarnings = TRUE)
  expected_names <- c("analyteType1", "analyteType2", "analyteType1MetaData",
                      "analyteType2MetaData", "sampleMetaData")
  testthat::expect_identical(slotNames(data_metab), expected_names)
  testthat::expect_identical(colnames(data_metab@sampleMetaData),
                   c("Feat1", "Feat2", "Feat3"))
  testthat::expect_identical(rownames(data_metab@analyteType1),
                   c("Metab1", "Metab2","Metab3"))
  testthat::expect_identical(rownames(data_metab@analyteType2),
                   c("Gene1", "Gene2","Gene3"))
  testthat::expect_equal(length(data_metab@analyteType1MetaData), 0)
  testthat::expect_identical(rownames(data_metab@analyteType2MetaData),
				   c("Gene1", "Gene2","Gene3"))

  # Check when only gene data is missing.
  missing_gene_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv",
                                             "metab_metadata_file.csv","",
                                             "pdata_file.csv"))
  rownames(missing_gene_df) <- expected
  fname_gene_missing <- paste(getwd(), "missing_gene_metadata_file.csv", sep = "/")
  write.csv(missing_gene_df, file = fname_gene_missing, quote=FALSE)
  testthat::expect_warning(IntLIM::ReadData(fname_gene_missing, class.feat = list(Feat1 = "numeric",
                                                                        Feat2 = "numeric",
                                                                        Feat3 = "numeric")),
                 "No metadata provided for Analyte Type 2", fixed = TRUE)

  # Check that fields are still populated appropriately.
  data_gene <- IntLIM::ReadData(fname_gene_missing, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"),
                                suppressWarnings = TRUE)
  expected_names <- c("analyteType1", "analyteType2", "analyteType1MetaData",
                      "analyteType2MetaData", "sampleMetaData")
  testthat::expect_identical(slotNames(data_gene), expected_names)
  testthat::expect_identical(colnames(data_gene@sampleMetaData),
                   c("Feat1", "Feat2", "Feat3"))
  testthat::expect_identical(rownames(data_gene@analyteType1),
                   c("Metab1", "Metab2","Metab3"))
  testthat::expect_identical(rownames(data_gene@analyteType2),
                   c("Gene1", "Gene2","Gene3"))
  testthat::expect_equal(length(data_gene@analyteType2MetaData), 0)
  testthat::expect_identical(rownames(data_gene@analyteType1MetaData),
				   c("Metab1", "Metab2","Metab3"))

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(fname_metab_missing)
  file.remove(fname_gene_missing)
})

# If patient data file is missing, this should cause early termination of the
# program.
testthat::test_that("Missing patient data leads to early termination", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                            "analyteType2MetaData","sampleMetaData")

  # Save files.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), file = "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene_meta, quote=FALSE)

  # Check when patient data is missing.
  missing_pdata_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                              ,"gene_metadata_file.csv","?%!"))
  rownames(missing_pdata_df) <- expected
  fname_pdata_missing <- paste(getwd(), "missing_pdata_file.csv", sep = "/")
  write.csv(missing_pdata_df, file = fname_pdata_missing, quote=FALSE)
  testthat::expect_error(IntLIM::ReadData(fname_pdata_missing, class.feat = list(Feat1 = "numeric",
                                                                       Feat2 = "numeric",
                                                                       Feat3 = "numeric")),
                 paste("File", paste(getwd(), "?%!", sep = "/"),
                       "does not exist"), fixed = TRUE)

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(fname_pdata_missing)
})

# If the id column is missing from the metadata, this should lead to early
# termination of the program.
testthat::test_that("Missing 'id' column in metadata leads to early termination.",{
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  # Save files.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  metabMetaDataMissing <- data.frame("metabname" = c("Metab1","Metab2","Metab3"))
  fname_metab_meta_missing <- paste(getwd(), "metab_metadata_missing_file.csv", sep = "/")
  write.csv(metabMetaDataMissing, file = fname_metab_meta_missing, quote=FALSE,
            row.names = FALSE)

  geneMetaDataMissing <- data.frame("genename" = c("Gene1","Gene2","Gene3"))
  fname_gene_meta_missing <- paste(getwd(), "gene_metadata_missing_file.csv", sep = "/")
  write.csv(geneMetaDataMissing, file = fname_gene_meta_missing, quote=FALSE,
            row.names = FALSE)

  # Check that error is thrown.
  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_missing_file.csv"
                                              ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  fname_metab_meta_noid <- paste(getwd(), "metab_metadata_noid_file.csv", sep = "/")
  write.csv(all_df, file = fname_metab_meta_noid, quote=FALSE)
  testthat::expect_error(IntLIM::ReadData(fname_metab_meta_noid, class.feat = list(Feat1 = "numeric",
                                                                         Feat2 = "numeric",
                                                                         Feat3 = "numeric")),
               paste("analyteType1id provided id does not exist in",
                     "Analyte Type 1 meta data file"),
                 fixed = TRUE)

  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                    ,"gene_metadata_missing_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  fname_gene_meta_noid <- paste(getwd(), "gene_metadata_noid_file.csv", sep = "/")
  write.csv(all_df, file = fname_gene_meta_noid, quote=FALSE)
  testthat::expect_error(IntLIM::ReadData(fname_gene_meta_noid, class.feat = list(Feat1 = "numeric",
                                                                        Feat2 = "numeric",
                                                                        Feat3 = "numeric")),
               paste("analyteType2id provided id does not exist in",
                     "Analyte Type 2 meta data file"),
               fixed = TRUE)

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(fname_metab_meta_missing)
  file.remove(fname_gene_meta_missing)
  file.remove(fname_metab_meta_noid)
  file.remove(fname_gene_meta_noid)
})

# If the ID's in the metadata file do not match the ID's in the analyte
# file, program should terminate.
testthat::test_that("Discrepancy between metabolite names and metabolite metadata
          leads to early termination.",{
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                            "analyteType2MetaData","sampleMetaData")

  # Save files.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  metabMetaDataWrong <- data.frame("id"=c("Metab1","FruityP3","Metab3"),
                              "metabname" = c("Metab1","Metab2","Metab3"))
  fname_metab_meta_wrong <- paste(getwd(), "metab_metadata_file_wrong.csv", sep = "/")
  write.csv(metabMetaDataWrong, file = fname_metab_meta_wrong, quote=FALSE, row.names = FALSE)

  geneMetaDataWrong <- data.frame("id"=c("Gene1","CocoaP3","Gene3"),
                                   "genename" = c("Gene1","Gene2","Gene3"))
  fname_gene_meta_wrong <- paste(getwd(), "gene_metadata_file_wrong.csv", sep = "/")
  write.csv(geneMetaDataWrong, file = fname_gene_meta_wrong, quote=FALSE, row.names = FALSE)

  # Check when only metabolite data is missing.
  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file_wrong.csv"
                                    ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  fname_metab_meta_wrong_ref <- paste(getwd(), "gene_metadata_wrong_file.csv", sep = "/")
  write.csv(all_df, file = fname_metab_meta_wrong_ref, quote=FALSE)
  testthat::expect_error(IntLIM::ReadData(fname_metab_meta_wrong_ref, class.feat = list(Feat1 = "numeric",
                                                                              Feat2 = "numeric",
                                                                              Feat3 = "numeric")),
               "Analytes in Type 1 data file and meta data files are not equal",
                 fixed = TRUE)

  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                    ,"gene_metadata_file_wrong.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  fname_gene_meta_wrong_ref <- paste(getwd(), "gene_genedata_wrong_file.csv", sep = "/")
  write.csv(all_df, file = fname_gene_meta_wrong_ref, quote=FALSE)
  testthat::expect_error(IntLIM::ReadData(fname_gene_meta_wrong_ref, class.feat = list(Feat1 = "numeric",
                                                                             Feat2 = "numeric",
                                                                             Feat3 = "numeric")),
               "Analytes in Type 2 data file and meta data files are not equal",
               fixed = TRUE)

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(fname_gene_meta_wrong)
  file.remove(fname_metab_meta_wrong)
  file.remove(fname_metab_meta_wrong_ref)
  file.remove(fname_gene_meta_wrong_ref)
})

# When all fields are present, an IntLimData object should be created.
# All fields should be populated appropriately.
testthat::test_that("All fields are present", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  # Save files.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                    ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  ref_file <- paste(getwd(), "all_ref.csv", sep = "/")
  write.csv(all_df, file = ref_file, quote=FALSE)

  # Check the characteristics of the data set (both gene and metabolite).
  dataset <- IntLIM::ReadData(ref_file, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"))
  expected_names <- c("analyteType1", "analyteType2", "analyteType1MetaData",
                      "analyteType2MetaData", "sampleMetaData")
  testthat::expect_identical(slotNames(dataset), expected_names)
  testthat::expect_identical(colnames(dataset@sampleMetaData), c("Feat1", "Feat2", "Feat3"))
  testthat::expect_identical(rownames(dataset@sampleMetaData), c("P1","P2", "P3",
                                                "P4"))
  testthat::expect_identical(rownames(dataset@analyteType1), c("Metab1", "Metab2","Metab3"))
  testthat::expect_identical(colnames(dataset@analyteType1), c("P1","P2", "P3",
                                                       "P4"))
  testthat::expect_identical(rownames(dataset@analyteType2), c("Gene1","Gene2","Gene3"))
  testthat::expect_identical(colnames(dataset@analyteType2), c("P1","P2", "P3",
                                                     "P4"))

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(ref_file)
})

# When all fields are present, an IntLimData object should be created.
# This should work even when a single covariate is included.
testthat::test_that("All fields are present", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")
  
  # Save files.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)
  
  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)
  
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)
  
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)
  
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)
  
  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                    ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  ref_file <- paste(getwd(), "all_ref.csv", sep = "/")
  write.csv(all_df, file = ref_file, quote=FALSE)
  
  # Check the characteristics of the data set (both gene and metabolite).
  dataset <- IntLIM::ReadData(ref_file, class.feat = list(Feat1 = "numeric"))
  expected_names <- c("analyteType1", "analyteType2", "analyteType1MetaData",
                      "analyteType2MetaData", "sampleMetaData")
  testthat::expect_identical(slotNames(dataset), expected_names)
  testthat::expect_identical(colnames(dataset@sampleMetaData), c("Feat1"))
  testthat::expect_identical(rownames(dataset@sampleMetaData), c("P1","P2", "P3",
                                                       "P4"))
  testthat::expect_identical(rownames(dataset@analyteType1), c("Metab1", "Metab2","Metab3"))
  testthat::expect_identical(colnames(dataset@analyteType1), c("P1","P2", "P3",
                                                     "P4"))
  testthat::expect_identical(rownames(dataset@analyteType2), c("Gene1","Gene2","Gene3"))
  testthat::expect_identical(colnames(dataset@analyteType2), c("P1","P2", "P3",
                                                     "P4"))
  
  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(ref_file)
})

# When no covariates are specified, all data in sampleMetaData should
# be read.
testthat::test_that("All fields are present", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")
  
  # Save files.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)
  
  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)
  
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)
  
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)
  
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)
  
  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                    ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  ref_file <- paste(getwd(), "all_ref.csv", sep = "/")
  write.csv(all_df, file = ref_file, quote=FALSE)
  
  # Check the characteristics of the data set (both gene and metabolite).
  dataset <- IntLIM::ReadData(ref_file)
  expected_names <- c("analyteType1", "analyteType2", "analyteType1MetaData",
                      "analyteType2MetaData", "sampleMetaData")
  testthat::expect_identical(slotNames(dataset), expected_names)
  testthat::expect_identical(colnames(dataset@sampleMetaData), c("Feat1", "Feat2", "Feat3"))
  testthat::expect_identical(rownames(dataset@sampleMetaData), c("P1","P2", "P3",
                                                       "P4"))
  testthat::expect_identical(rownames(dataset@analyteType1), c("Metab1", "Metab2","Metab3"))
  testthat::expect_identical(colnames(dataset@analyteType1), c("P1","P2", "P3",
                                                     "P4"))
  testthat::expect_identical(rownames(dataset@analyteType2), c("Gene1","Gene2","Gene3"))
  testthat::expect_identical(colnames(dataset@analyteType2), c("P1","P2", "P3",
                                                     "P4"))
  
  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(ref_file)
})

# When log scaling is requested for positive valued data, data should be
# appropriately log-scaled.
testthat::test_that("Data is log-scaled", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  # Save files.
  metabData <- data.frame("Metab1"=c(0,1,2,4), "Metab2"=c(0,1,2,4), "Metab3"=c(0,1,2,4))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  geneData <- data.frame("Gene1"=c(0,1,2,4), "Gene2"=c(0,1,2,4), "Gene3"=c(0,1,2,4))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  pData <- data.frame("Feat1"=c(0,1,2,4), "Feat2"=c(0,1,2,4), "Feat3"=c(0,1,2,4))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                    ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  ref_file <- paste(getwd(), "all_ref.csv", sep = "/")
  write.csv(all_df, file = ref_file, quote=FALSE)

  # Check the characteristics of the data set (both gene and metabolite).
  dataset <- IntLIM::ReadData(ref_file, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"),
                              logAnalyteType1 = TRUE, logAnalyteType2 = TRUE)
  cutoff <- 0.0000001
  testthat::expect_lt(dataset@analyteType1[1,1],0)
  testthat::expect_lt(dataset@analyteType1[2,1],0)
  testthat::expect_lt(dataset@analyteType1[3,1],0)
  testthat::expect_equal(dataset@analyteType1[1,2],0, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType1[2,2],0, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType1[3,2],0, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType1[1,3],1, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType1[2,3],1, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType1[3,3],1, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType1[1,4],2, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType1[2,4],2, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType1[3,4],2, tolerance=cutoff)

  testthat::expect_lt(dataset@analyteType2[1,1],0)
  testthat::expect_lt(dataset@analyteType2[2,1],0)
  testthat::expect_lt(dataset@analyteType2[3,1],0)
  testthat::expect_equal(dataset@analyteType2[1,2],0, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType2[2,2],0, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType2[3,2],0, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType2[1,3],1, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType2[2,3],1, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType2[3,3],1, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType2[1,4],2, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType2[2,4],2, tolerance=cutoff)
  testthat::expect_equal(dataset@analyteType2[3,4],2, tolerance=cutoff)

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(ref_file)
})

# Log-scaling for negative data should result in a warning.
testthat::test_that("Negative data is not log-scaled", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  # Save files.
  metabData <- data.frame("Metab1"=c(-1,1,2,4), "Metab2"=c(0,1,2,4), "Metab3"=c(0,1,2,4))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  geneData <- data.frame("Gene1"=c(0,1,2,4), "Gene2"=c(-1,1,2,4), "Gene3"=c(0,1,2,4))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  pData <- data.frame("Feat1"=c(0,1,2,4), "Feat2"=c(0,1,2,4), "Feat3"=c(0,1,2,4))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                    ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  ref_file <- paste(getwd(), "all_ref.csv", sep = "/")
  write.csv(all_df, file = ref_file, quote=FALSE)

  # Check the characteristics of the data set (both gene and metabolite).
  dataset <- IntLIM::ReadData(ref_file, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"),
														  suppressWarnings = TRUE, logAnalyteType1 = TRUE, 
														  logAnalyteType2 = TRUE)
  testthat::expect_warning(IntLIM::ReadData(ref_file, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"),
                                  logAnalyteType1 = TRUE, logAnalyteType2 = TRUE),
														  "Analyte Type 1 data has negative values. Continuing without log-scaling.", fixed = TRUE)
  testthat::expect_warning(IntLIM::ReadData(ref_file, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"),
                                  logAnalyteType1 = TRUE, logAnalyteType2 = TRUE),
														  "Analyte Type 2 data has negative values. Continuing without log-scaling.", fixed = TRUE)
  testthat::expect_equal(dataset@analyteType1,t(metabData))
  testthat::expect_equal(dataset@analyteType2,t(geneData))

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(ref_file)
})

# Test that output is appropriate when a different ID is used.
  # Data types expected
testthat::test_that("Other ID's also work", {
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  # Save files.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  metabMetaData <- data.frame("metabIdentifier"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  geneMetaData <- data.frame("geneIdentifier"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                    ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  ref_file <- paste(getwd(), "all_ref.csv", sep = "/")
  write.csv(all_df, file = ref_file, quote=FALSE)

  # Check the characteristics of the data set (both gene and metabolite).
  dataset <- IntLIM::ReadData(ref_file, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"),
														  analyteType1id = "metabIdentifier",
														  analyteType2id = "geneIdentifier")
  expected_names <- c("analyteType1", "analyteType2", "analyteType1MetaData",
                      "analyteType2MetaData", "sampleMetaData")
  testthat::expect_identical(slotNames(dataset), expected_names)
  testthat::expect_identical(colnames(dataset@sampleMetaData), c("Feat1", "Feat2", "Feat3"))
  testthat::expect_identical(rownames(dataset@sampleMetaData), c("P1","P2", "P3",
                                                "P4"))
  testthat::expect_identical(rownames(dataset@analyteType1), c("Metab1", "Metab2","Metab3"))
  testthat::expect_identical(colnames(dataset@analyteType1), c("P1","P2", "P3",
                                                       "P4"))
  testthat::expect_identical(rownames(dataset@analyteType2), c("Gene1","Gene2","Gene3"))
  testthat::expect_identical(colnames(dataset@analyteType2), c("P1","P2", "P3",
                                                     "P4"))

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(ref_file)
 })

# Test that samples not shared between patient and analyte files are removed.
testthat::test_that("Samples not shared are removed", {

  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  # Save files.
  metabData <- data.frame("Metab1"=c(0,0,0,0,0), "Metab2"=c(0,0,0,0,0), "Metab3"=c(0,0,0,0,0))
  rownames(metabData) <- c("P1", "P2", "P3", "P4", "P7")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  geneData <- data.frame("Gene1"=c(0,0,0,0,0), "Gene2"=c(0,0,0,0,0), "Gene3"=c(0,0,0,0,0))
  rownames(geneData) <- c("P1", "P2", "P3", "P4", "P5")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  pData <- data.frame("Feat1"=c(0,0,0,0,0), "Feat2"=c(0,0,0,0,0), "Feat3"=c(0,0,0,0,0))
  rownames(pData) <- c("P1", "P2", "P3", "P4", "P6")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                    ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  ref_file <- paste(getwd(), "all_ref.csv", sep = "/")
  write.csv(all_df, file = ref_file, quote=FALSE)

  # Check for the warning.
  testthat::expect_warning(IntLIM::ReadData(ref_file, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric")),
				paste("The following samples were not shared in all data types and were removed:",
				"P6, P7, P5"),
				fixed = TRUE)

  # Check the characteristics of the data set (both gene and metabolite).
  dataset <- IntLIM::ReadData(ref_file, class.feat = list(Feat1 = "numeric",
                                                          Feat2 = "numeric",
                                                          Feat3 = "numeric"),
														  suppressWarnings = TRUE)
  expected_names <- c("analyteType1", "analyteType2", "analyteType1MetaData",
                      "analyteType2MetaData", "sampleMetaData")
  testthat::expect_identical(slotNames(dataset), expected_names)
  testthat::expect_identical(colnames(dataset@sampleMetaData), c("Feat1", "Feat2", "Feat3"))
  testthat::expect_identical(rownames(dataset@sampleMetaData), c("P1","P2", "P3",
                                                "P4"))
  testthat::expect_identical(rownames(dataset@analyteType1), c("Metab1", "Metab2","Metab3"))
  testthat::expect_identical(colnames(dataset@analyteType1), c("P1","P2", "P3",
                                                       "P4"))
  testthat::expect_identical(rownames(dataset@analyteType2), c("Gene1","Gene2","Gene3"))
  testthat::expect_identical(colnames(dataset@analyteType2), c("P1","P2", "P3",
                                                     "P4"))

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(ref_file)
 })

# Test that names are converted appropriately.
testthat::test_that("Names with special symbols are converted appropriately", {
  # Data types expected
  expected <- c("analyteType1","analyteType2","analyteType1MetaData",
                "analyteType2MetaData","sampleMetaData")

  # Save files.
  metabData <- data.frame("Metab-1"=c(0,0,0,0), "Met-ab2"=c(0,0,0,0), "Metab&3"=c(0,0,0,0))
  rownames(metabData) <- c("P1%%", "4P2", "*P3-", "Bam*bam")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  geneData <- data.frame("Gene*1"=c(0,0,0,0), "Gene&&2"=c(0,0,0,0), "Gene()3"=c(0,0,0,0))
  rownames(geneData) <- c("P1%%", "4P2", "*P3-", "Bam*bam")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  pData <- data.frame("Feat-1"=c(0,0,0,0), "Feat~2"=c(0,0,0,0), "Feat`3"=c(0,0,0,0))
  rownames(pData) <- c("P1%%", "4P2", "*P3-", "Bam*bam")
  fname_pdata <- paste(getwd(), "pdata_file.csv", sep = "/")
  write.csv(pData, file = fname_pdata, quote=FALSE)

  metabMetaData <- data.frame("id"=c("Metab-1", "Met-ab2", "Metab&3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)

  geneMetaData <- data.frame("id"=c("Gene*1", "Gene&&2", "Gene()3"), "genename"=
                               c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE, row.names = FALSE)

  all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
                                    ,"gene_metadata_file.csv","pdata_file.csv"))
  rownames(all_df) <- expected
  ref_file <- paste(getwd(), "all_ref.csv", sep = "/")
  write.csv(all_df, file = ref_file, quote=FALSE)

  # Check the characteristics of the data set (both gene and metabolite).
  dataset <- IntLIM::ReadData(ref_file, class.feat = list(Feat.1 = "numeric",
                                                          Feat.2 = "numeric",
                                                          Feat.3 = "numeric"))
  expected_names <- c("analyteType1", "analyteType2", "analyteType1MetaData",
                      "analyteType2MetaData", "sampleMetaData")
  testthat::expect_identical(slotNames(dataset), expected_names)
  testthat::expect_identical(colnames(dataset@sampleMetaData), c("Feat.1", "Feat.2", "Feat.3"))
  testthat::expect_identical(rownames(dataset@sampleMetaData), c("P1..", "X4P2", "X.P3.", "Bam.bam"))
  testthat::expect_identical(rownames(dataset@analyteType1), c("Metab.1", "Met.ab2","Metab.3"))
  testthat::expect_identical(colnames(dataset@analyteType1), c("P1..", "X4P2", "X.P3.", "Bam.bam"))
  testthat::expect_identical(rownames(dataset@analyteType2), c("Gene.1","Gene..2","Gene..3"))
  testthat::expect_identical(colnames(dataset@analyteType2), c("P1..", "X4P2", "X.P3.", "Bam.bam"))

  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(ref_file)
})