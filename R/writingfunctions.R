#' Output data into individual CSV files.  All data will be zipped into one file with all data.
#'
#' @param inputData data output from ReadData() or FilterData() function
#' @param filename name of file to be output (default: '~/output.zip') 
#'
#' @return the filename of the CSV file with results named with cohort
#' @export

OutputData <- function (inputData=NULL,filename="~/output.zip"){
  mygenes <- inputData$gene
  mymetab <- inputData$metab
	phenoData <- inputData$covar_matrix
	utils::write.csv(mygenes,path.expand("~/GeneData.csv"),quote=F)
  utils::write.csv(mymetab,path.expand("~/MetabData.csv"),quote=F)
  utils::write.csv(phenoData,path.expand("~/MetaSamples.csv"),quote=F)
	utils::zip(zipfile=path.expand(filename),
	           flags="-b ~",
	           files=c(path.expand("~/GeneData.csv"),
	                   path.expand("~/MetabData.csv"),
	                   path.expand("~/MetaSamples.csv")))
	file.remove(path.expand("~/GeneData.csv"))
	file.remove(path.expand("~/MetabData.csv"))
	file.remove(path.expand("~/MetaSamples.csv"))
}

#' Output results into a zipped CSV file.  Results include gene and metabolite pairs, along with model interaction p-values, and correlations in each group being evaluated.
#'
#' @param inputResults IntLimResults object with model results (output of ProcessResults())
#' @param filename name of file to be output (default: '~/results.csv')
#'
#' @return the filename of the CSV file with results named with cohort
#' @export

OutputResults <- function (inputResults=NULL,filename="results.csv"){
	if(is.null(inputResults)) {stop("Input results from ProcessResults()")}
	utils::write.csv(inputResults,filename,quote=T,row.names=F)
}


