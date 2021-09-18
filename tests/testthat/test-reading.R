# If a file is passed as input that does not exist in the system, an error
# should be thrown.
test_that("Nonexistent files cause early termination.", {
  expect_error(IntLIM::ReadData("?!%"), "CSV input file does not exist", ignore.case = TRUE)
})

# If a file with inappropriate column names is passed as input, an error
# should be thrown.
test_that("Incorrect column names cause early termination.", {

  # Create a data frame with inappropriate column names.
  incorrect_colname_df = data.frame("Wilma"=c(1,2,3), "Fred"=c(0,0,0))

  # Save the data frame as a file.
  fname <- paste(getwd(), "incorrect_colname_file.csv", sep = "/")
  write.csv(incorrect_colname_df, file=fname, quote=FALSE, row.names = FALSE)
  
  # Check for error.
  expect_error(IntLIM::ReadData(fname),
               "Check column names of input files.  'type' and 'filenames' are required", 
               ignore.case = TRUE)
  file.remove(fname)
})

# If one or more data type headings is missing, an error should be thrown.
test_that("Missing data types cause early termination.", {

  # Data types expected
  expected <- c("metabData","geneData","metabMetaData","geneMetaData","sampleMetaData")

  # Check when each type is missing.
  for(heading in expected){

    # Create a file with one missing heading.
    missing_heading_df = data.frame("filenames"=c(0,0,0,0))
    rownames(missing_heading_df) <- setdiff(expected, heading)
    fname <- paste(getwd(), "missing_heading_file.csv", sep = "/")
    write.csv(missing_heading_df, file=fname, quote=FALSE)

    # Check that error occurs.
    message <- paste("The column 'type' contains non-allowed entries (See Description). The",
                     "CSV input file must contain 6 rows (if optional meta data files for metabolites",
                     "and genes are not to be input, have the corresponding filenames be blanks.")
    expect_error(IntLIM::ReadData(fname), message, fixed = TRUE)
    file.remove(fname)
  }
})

# If no data is included, an error should be thrown.
test_that("Absence of data causes early termination.", {
  # Data types expected
  expected <- c("metabData","geneData","metabMetaData","geneMetaData",
                "sampleMetaData")

  # Create a file with both gene and metabolite data missing.
  missing_data_df = data.frame("filenames"=c("","",0,0,0))
  rownames(missing_data_df) <- expected
  fname <- paste(getwd(), "missing_data_file.csv", sep = "/")
  write.csv(missing_data_df, file=fname, quote=FALSE)

  # Check that error occurs.
  expect_error(IntLIM::ReadData(fname), "No data provided.", ignore.case = TRUE)
  file.remove(fname)
})

# If only one data type is missing, a warning should be given.
test_that("Absence of one omic type results in a warning.", {
  # Data types expected
  expected <- c("metabData","geneData","metabMetaData","geneMetaData",
                "sampleMetaData")
  
  #Create metabolite data file.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)
  
  # Create metabolite metadata file.
  metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
                                c("Metab1", "Metab2", "Metab3"))
  fname_metab_meta <- paste(getwd(), "metab_metadata_file.csv", sep = "/")
  write.csv(metabMetaData, file = fname_metab_meta, quote=FALSE, row.names = FALSE)
    
  # Create gene data file.
  geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
  rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)
  
  # Create gene metadata file.
  geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
                                c("Gene1", "Gene2", "Gene3"))
  fname_gene_meta <- paste(getwd(), "gene_metadata_file.csv", sep = "/")
  write.csv(geneMetaData, file = fname_gene_meta, quote=FALSE)

  # Create patient data file.
  pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
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
  message <- paste("No data provided for metabolites. This means you cannot run",
                "metabolite-metabolite or gene-metabolite analyses.\n")
  expect_warning(IntLIM::ReadData(fname), message, ignore.case = TRUE)
  file.remove(fname)

  # Create a file with only gene data missing.
  missing_gene_df = data.frame("filenames"=c("metab_file.csv","",
                                             "metab_metadata_file.csv",
                                             "gene_metadata_file.csv",
                                             "pdata_file.csv"))
  rownames(missing_gene_df) <- expected
  fname <- paste(getwd(), "missing_genes_file.csv", sep = "/")
  write.csv(missing_gene_df, file = fname, quote=FALSE)

  # Check that warning occurs.
  message <- paste("No data provided for genes. This means you cannot run",
                   "gene-gene or gene-metabolite analyses.\n")
  expect_warning(IntLIM::ReadData(fname), message, ignore.case = TRUE)
  
  # Remove files.
  file.remove(fname_metab)
  file.remove(fname_metab_meta)
  file.remove(fname_gene)
  file.remove(fname_gene_meta)
  file.remove(fname_pdata)
  file.remove(fname)
})

# If there are two analytes with the same name, an error should be thrown.
test_that("Duplicate analytes cause early termination", {
  # Data types expected
  expected <- c("metabData","geneData","metabMetaData","geneMetaData",
                "sampleMetaData")

  #Create metabolite data file.
  metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
  rownames(metabData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
  fname_metab <- paste(getwd(), "metab_file.csv", sep = "/")
  write.csv(t(metabData), file = fname_metab, quote=FALSE)

  # Create metabolite data file with duplicates.
  metabData <- data.frame("1"=c("Metab1",0,0,0,0), "2"=c("Metab2",0,0,0,0), 
                          "3"=c("Metab2",0,0,0,0))
  rownames(metabData) <- c("","Fred", "Wilma", "Pebbles", "Bambam")
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
  rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
  fname_gene <- paste(getwd(), "gene_file.csv", sep = "/")
  write.csv(t(geneData), file = fname_gene, quote=FALSE)

  # Create gene data file with duplicates.
  geneData <- data.frame("1"=c("Gene1",0,0,0,0), "2"=c("Gene2",0,0,0,0), 
                         "3"=c("Gene2",0,0,0,0))
  rownames(geneData) <- c("","Fred", "Wilma", "Pebbles", "Bambam")
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
  rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
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
  expect_error(IntLIM::ReadData(fname_metab_dup_ref),
                 paste("Error: your input file",fname_metab_dup, "has duplicate",
                       "entries in column 1. Please make sure you have one row per metabolite"),
                 ignore.case = TRUE)

  # Check that error is thrown.
  fname_gene_dup <- paste(getwd(), "gene_file_dup.csv", sep = "/")
  expect_error(IntLIM::ReadData(fname_gene_dup_ref),
                 paste("Error: your input file",fname_gene_dup, "has duplicate",
                       "entries in column 1. Please make sure you have one row per gene"),
                 ignore.case = TRUE)
  
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
# 
# # If a data file is not accessible on the system, an error should be thrown.
# test_that("Inaccessible data files cause early termination.", {
#   # Data types expected
#   expected <- c("metabData","geneData","metabMetaData","geneMetaData",
#                 "sampleMetaData")
#   
#   # Create metabolite file.
#   metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
#   rownames(metabData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("metab_file.csv", t(metabData), quote=FALSE)
#   
#   # Create metabolite metadata file.
#   metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
#                                 c("Metab1", "Metab2", "Metab3"))
#   write.csv("metab_metadata_file.csv", metabMetaData, quote=FALSE)
#   
#   # Create gene file.
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   
#   # Create gene metadata file.
#   geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
#                                c("Gene1", "Gene2", "Gene3"))
#   write.csv("gene_metadata_file.csv", geneMetaData, quote=FALSE)
#   
#   # Create patient data file.
#   pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
#   rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("pdata_file.csv", pData, quote=FALSE)
#   
#   # Create reference file that points to nonexistent metabolite file.
#   missing_metab_df = data.frame("filenames"=c("?!%","gene_file.csv","metabolite_metadata_file.csv"
#                                               ,"gene_metadata_file.csv","pdata_file.csv"))
#   rownames(missing_metab_df) <- expected
#   write.csv("missing_metabdata_file.csv", missing_metab_df, quote=FALSE)
#   
#   # Check that an error is thrown.
#   expect_error(IntLIM::ReadData("missing_metabdata_file.csv"), 
#                  paste("File", base::dirname("missing_metabdata_file.csv"), "/?!%",
#                        "does not exist"), ignore.case = TRUE)
#   
#   # Create reference file that points to nonexistent gene file.
#   missing_gene_df = data.frame("filenames"=c("metabolite_file.csv","?!%", "metabolite_metadata_file.csv"
#                                              ,"gene_metadata_file.csv","pdata_file.csv"))
#   rownames(missing_gene_df) <- expected
#   write.csv("missing_gene_file.csv", missing_gene_df, quote=FALSE)
#   
#   # Check that an error is thrown.
#   expect_error(IntLIM::ReadData("missing_gene_file.csv"), 
#                  paste("File", base::dirname("missing_gene_file.csv"), "/?!%",
#                        "does not exist"), ignore.case = TRUE)
# })
# 
# # If analyte metadata files are missing, this should throw a warning but not lead
# # to early termination.
# test_that("Missing metadata causes a warning.", {
#   # Data types expected
#   expected <- c("metabData","geneData","metabMetaData","geneMetaData",
#                 "sampleMetaData")
#   
#   # Save files.
#   metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
#   rownames(metabData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("metab_file.csv", t(metabData), quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
#   rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("pdata_file.csv", pData, quote=FALSE)
#   metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
#                                 c("Metab1", "Metab2", "Metab3"))
#   write.csv("metab_metadata_file.csv", metabMetaData, quote=FALSE)
#   geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
#                                 c("Gene1", "Gene2", "Gene3"))
#   write.csv("gene_metadata_file.csv", geneMetaData, quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   
#   # Check when only metabolite data is missing.
#   missing_metab_df = data.frame("filenames"=c("metab_file","gene_file.csv","?!%"
#                                               ,"gene_metadata_file.csv","pdata_file.csv"))
#   rownames(missing_metab_df) <- expected
#   write.csv("missing_metab_metadata_file.csv", missing_metab_df, quote=FALSE)
#   expect_warning(IntLIM::ReadData("missing_metab_metadata_file.csv"), 
#                  paste("File", base::dirname("missing_metab_metadata_file.csv"), "/?!%",
#                        "does not exist"), ignore.case = TRUE)
#   
#   # Check when only gene data is missing.
#   missing_gene_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
#                                               ,"?!%","pdata_file.csv"))
#   rownames(missing_gene_df) <- expected
#   write.csv("missing_gene_metadata_file.csv", missing_gene_df, quote=FALSE)
#   expect_warning(IntLIM::ReadData("missing_gene_metadata_file.csv"), 
#                  paste("File", base::dirname("missing_gene_metadata_file.csv"), "/?!%",
#                        "does not exist"), ignore.case = TRUE)
# })
# 
# # If patient data file is missing, this should cause early termination of the
# # program.
# test_that("Missing patient data leads to early termination", {
#   # Data types expected
#   expected <- c("metabData","geneData","metabMetaData","geneMetaData",
#                 "sampleMetaData")
#   
#   # Save files.
#   metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
#   rownames(metabData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("metab_file.csv", t(metabData), quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
#   rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("pdata_file.csv", pData, quote=FALSE)
#   metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
#                                 c("Metab1", "Metab2", "Metab3"))
#   write.csv("metab_metadata_file.csv", metabMetaData, quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   
#   # Check when patient data is missing.
#   missing_pdata_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
#                                               ,"gene_metadata_file.csv","?%!"))
#   rownames(missing_pdata_df) <- expected
#   write.csv("missing_pdata_file.csv", missing_pdata_df, quote=FALSE)
#   expect_error(IntLIM::ReadData("missing_pdata_file.csv"), 
#                  paste("File", base::dirname("missing_pdata_file.csv"), "/?!%",
#                        "does not exist"), ignore.case = TRUE)
# })
# 
# # If the id column is missing from the metabolite metadata, this should lead to early
# # termination of the program.
# test_that("Missing 'id' column in metabolite metadata leads to early termination.",{
#   # Data types expected
#   expected <- c("metabData","geneData","metabMetaData","geneMetaData",
#                 "sampleMetaData")
#   
#   # Save files.
#   metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
#   rownames(metabData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("metab_file.csv", t(metabData), quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
#   rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("pdata_file.csv", pData, quote=FALSE)
#   metabMetaData <- data.frame("metabname" = c("Metab1","Metab2","Metab3"))
#   write.csv("metab_file.csv", t(metabData), quote=FALSE, rownames=FALSE)
#   write.csv("gene_metadata_file.csv", geneMetaData, quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   
#   # Remove the ID column.
#   metabMetaData <- data.frame("metabname"=
#                                 c("Metab1", "Metab2", "Metab3"))
#   write.csv("metab_metadata_file.csv", metabMetaData, quote=FALSE)
#   geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
#                                 c("Gene1", "Gene2", "Gene3"))
#   
#   # Check that error is thrown.
#   all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
#                                               ,"gene_metadata_file.csv","pdata_file.csv"))
#   rownames(all_df) <- expected
#   write.csv("all.csv", all_df, quote=FALSE)
#   expect_error(IntLIM::ReadData("all.csv"), 
#                  stop(paste("metabid provided id",
#                                   "does not exist in metabolite meta data file")), 
#                  ignore.case = TRUE)
# })
# 
# # If the id column is missing from the gene metadata, this should lead to early
# # termination of the program.
# test_that("Missing 'id' column in gene metadata leads to early termination.",{
#   # Data types expected
#   expected <- c("metabData","geneData","metabMetaData","geneMetaData",
#                 "sampleMetaData")
#   
#   # Save files.
#   metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
#   rownames(metabData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("metab_file.csv", t(metabData), quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
#   rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("pdata_file.csv", pData, quote=FALSE)
#   metabMetaData <- data.frame("metabname" = c("Metab1","Metab2","Metab3"))
#   write.csv("metab_file.csv", t(metabData), quote=FALSE, rownames=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
#                                 c("Metab1", "Metab2", "Metab3"))
#   write.csv("metab_metadata_file.csv", metabMetaData, quote=FALSE)
#   
#   # Remove the ID column.
# 
#   geneMetaData <- data.frame("genename"=
#                                c("Gene1", "Gene2", "Gene3"))
#   write.csv("gene_metadata_file.csv", geneMetaData, quote=FALSE)
#   
#   # Check that error is thrown.
#   all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
#                                     ,"gene_metadata_file.csv","pdata_file.csv"))
#   rownames(all_df) <- expected
#   write.csv("all.csv", all_df, quote=FALSE)
#   expect_error(IntLIM::ReadData("all.csv"), 
#                stop(paste("metabid provided id",
#                           "does not exist in metabolite meta data file")), 
#                ignore.case = TRUE)
# })
# 
# # If the ID's in the metadata file do not match the ID's in the metabolite
# # file, program should terminate.
# test_that("Discrepancy between metabolite names and metabolite metadata
#           leads to early termination.",{
#   # Data types expected
#   expected <- c("metabData","geneData","metabMetaData","geneMetaData",
#                 "sampleMetaData")
#   
#   # Save files.
#   metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
#   rownames(metabData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("metab_file.csv", t(metabData), quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
#   rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("pdata_file.csv", pData, quote=FALSE)
#   metabMetaData <- data.frame("id"=c("Metab1","FruityPebbles","Metab3"),
#                               "metabname" = c("Metab1","Metab2","Metab3"))
#   write.csv("metab_file.csv", t(metabData), quote=FALSE, rownames=FALSE)
#   geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
#                                c("Gene1", "Gene2", "Gene3"))
#   write.csv("gene_metadata_file.csv", geneMetaData, quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   
#   # Check when only metabolite data is missing.
#   all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
#                                     ,"gene_metadata_file.csv","pdata_file.csv"))
#   rownames(all_df) <- expected
#   write.csv("all.csv", all_df, quote=FALSE)
#   expect_warning(IntLIM::ReadData("all.csv"), 
#                  stop("Metabolites in abundance data file and metabolite meta 
#                  data file are not equal"), 
#                  ignore.case = TRUE)
# })
# 
# # If the patient ID does not match the ID's in the analyte data, the program should
# # terminate early.
# test_that("Discrepancy between samples in patient data and analyte data files
#           leads to early termination.",{
#   # Data types expected
#   expected <- c("metabData","geneData","metabMetaData","geneMetaData",
#                 "sampleMetaData")
#   
#   # Save files.
#   metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
#   rownames(metabData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("metab_file.csv", t(metabData), quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
#   rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Dino")
#   write.csv("pdata_file.csv", pData, quote=FALSE)
#   metabMetaData <- data.frame("id"=c("Metab1","Metab2","Metab3"),
#                               "metabname" = c("Metab1","Metab2","Metab3"))
#   write.csv("metab_file.csv", t(metabData), quote=FALSE, rownames=FALSE)
#   geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
#                                c("Gene1", "Gene2", "Gene3"))
#   write.csv("gene_metadata_file.csv", geneMetaData, quote=FALSE)
#   
#   all_df = data.frame("filenames"=c("metab_file.csv","gene_file.csv","metab_metadata_file.csv"
#                                     ,"gene_metadata_file.csv","pdata_file.csv"))
#   rownames(all_df) <- expected
#   write.csv("all.csv", all_df, quote=FALSE)
#   expect_warning(IntLIM::ReadData("all.csv"), 
#                  stop("All samples in abundance data file must be in metabolite 
#                       meta data file"), 
#                  ignore.case = TRUE)
# })
# 
# # When all fields are present, a MultiDataSet object should be created.
# # All fields should be populated appropriately.
# # NEEDS TO BE FILLED IN
# test_that("All fields are present", {
#   # Data types expected
#   expected <- c("metabData","geneData","metabMetaData","geneMetaData",
#                 "sampleMetaData")
#   
#   # Save files.
#   metabData <- data.frame("Metab1"=c(0,0,0,0), "Metab2"=c(0,0,0,0), "Metab3"=c(0,0,0,0))
#   rownames(metabData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("metabolite_file.csv", t(metabData), quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   pData <- data.frame("Feat1"=c(0,0,0,0), "Feat2"=c(0,0,0,0), "Feat3"=c(0,0,0,0))
#   rownames(pData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("pdata_file.csv", pData, quote=FALSE)
#   correct_df = data.frame("filenames"=c("metab_file","gene_file.csv","metab_metadata_file"
#                                               ,"gene_metadata_file.csv","pdata_file.csv"))
#   rownames(correct_df) <- expected
#   write.csv("correct_data.csv", correct_df, quote=FALSE)
#   metabMetaData <- data.frame("id"=c("Metab1", "Metab2", "Metab3"), "metabname"=
#                                 c("Metab1", "Metab2", "Metab3"))
#   write.csv("metab_metadata_file.csv", metabMetaData, quote=FALSE)
#   geneMetaData <- data.frame("id"=c("Gene1", "Gene2", "Gene3"), "genename"=
#                                 c("Gene1", "Gene2", "Gene3"))
#   write.csv("gene_metadata_file.csv", geneMetaData, quote=FALSE)
#   geneData <- data.frame("Gene1"=c(0,0,0,0), "Gene2"=c(0,0,0,0), "Gene3"=c(0,0,0,0))
#   rownames(geneData) <- c("Fred", "Wilma", "Pebbles", "Bambam")
#   write.csv("gene_file.csv", t(geneData), quote=FALSE)
#   
#   # Build the data set that should be built.
#   multi <- MultiDataSet::createMultiDataSet()
#   metab.set <- methods::new("MetaboliteSet",metabData = metabdata,
#                          phenoData = metabphenoData, featureData = metabfeatureData)
#   gene.set <- Biobase::ExpressionSet(assayData=as.matrix(genedata))
#   multi <- add_metabolite(multi, metab.set)
#   multi <- MultiDataSet::add_genexp(multi, gene.set)
#   
#   # Check that data format is correct.
#   expect_identical(IntLIM::ReadData("missing_pdata_file.csv"), 
#                  paste("File", base::dirname("missing_pdata_file.csv"), "/?!%",
#                        "does not exist"), ignore.case = TRUE)
# })
