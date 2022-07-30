#' Output data into individual CSV files.  All data will be zipped into one file with all data.
#'
#' @param inputData data output from ReadData() or FilterData() function
#' @param filename name of file to be output (default: 'tempdir/output.zip') 
#'
#' @return the filename of the CSV file with results named with cohort
#' @export
OutputData <- function (inputData,filename=""){
  
  # Set up temporary directory if filename is not specified.
  if(filename == ""){
    tmp<- tempdir()
    if(!file.exists(tmp)){
      dir.create(tmp)
    }
    filename<-paste0(tmp, "/output.zip")
  }
  
  # Write data.
  type1 <- inputData@analyteType1
  type2 <- inputData@analyteType2
	phenoData <- inputData@sampleMetaData
	tmp2<- tempdir()
	if(!file.exists(tmp2)){
	  dir.create(tmp2)
	}
	utils::write.csv(type1,paste0(tmp2, "/AnalyteType1Data.csv"),quote=FALSE)
  utils::write.csv(type2,paste0(tmp2, "/AnalyteType2Data.csv"),quote=FALSE)
	utils::write.csv(phenoData,paste0(tmp2, "/MetaSamples.csv"),quote=FALSE)
	
	# Move data to zipped folder.
	utils::zip(zipfile=path.expand(filename),
	           flags=paste("-b", tmp2),
	                 files=c(paste0(tmp2, "/AnalyteType1Data.csv"),
	                   paste0(tmp2, "/AnalyteType2Data.csv"),
	                         paste0(tmp2, "/MetaSamples.csv")))
	file.remove(paste0(tmp2, "/AnalyteType1Data.csv"))
	                 file.remove(paste0(tmp2, "/AnalyteType2Data.csv"))
	file.remove(paste0(tmp2, "/MetaSamples.csv"))
	message(paste("Wrote results to", filename))
}

#' Output results into a zipped CSV file.  Results include gene and metabolite pairs, along with model interaction p-values, and correlations in each group being evaluated.
#'
#' @param inputResults IntLimResults object with model results (output of ProcessResults())
#' @param filename name of file to be output (default: 'tempdir/results.csv')
#'
#' @return the filename of the CSV file with results named with cohort
#' @export
OutputResults <- function (inputResults,filename=""){
  
  # Set up temporary directory if filename is not specified.
  if(filename == ""){
    tmp<- tempdir()
    if(!file.exists(tmp)){
      dir.create(tmp)
    }
    filename<-paste0(tmp, "/results.csv")
  }
  
  # Write data.
  utils::write.csv(inputResults,filename,quote=TRUE,row.names=FALSE)
  message(paste("Wrote results to", filename))
}