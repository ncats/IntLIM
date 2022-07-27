#' Output data into individual CSV files.  All data will be zipped into one file with all data.
#'
#' @param inputData data output from ReadData() or FilterData() function
#' @param filename name of file to be output (default: 'tempdir/output.zip') 
#'
#' @return the filename of the CSV file with results named with cohort
#' @export
OutputData <- function (inputData,filename="~/output.zip"){
  type1 <- inputData@analyteType1
  type2 <- inputData@analyteType2
	phenoData <- inputData@sampleMetaData
	utils::write.csv(type1,path.expand("~/AnalyteType1Data.csv"),quote=FALSE)
  utils::write.csv(type2,path.expand("~/AnalyteType2Data.csv"),quote=FALSE)
  utils::write.csv(phenoData,path.expand("~/MetaSamples.csv"),quote=FALSE)
	utils::zip(zipfile=path.expand(filename),
	           flags="-b ~",
	           files=c(path.expand("~/AnalyteType1Data.csv"),
	                   path.expand("~/AnalyteType2Data.csv"),
	                   path.expand("~/MetaSamples.csv")))
	file.remove(path.expand("~/AnalyteType1Data.csv"))
	file.remove(path.expand("~/AnalyteType2Data.csv"))
	file.remove(path.expand("~/MetaSamples.csv"))
}

#' Output results into a zipped CSV file.  Results include gene and metabolite pairs, along with model interaction p-values, and correlations in each group being evaluated.
#'
#' @param inputResults IntLimResults object with model results (output of ProcessResults())
#' @param filename name of file to be output (default: 'tempdir/results.csv')
#'
#' @return the filename of the CSV file with results named with cohort
#' @export
OutputResults <- function (inputResults,filename="results.csv"){
	utils::write.csv(inputResults,filename,quote=TRUE,row.names=FALSE)
}


