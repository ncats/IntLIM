#' Read in CSV file
#'
#' The metadata associated with data files to be analyzed in IntLim is supplied
#' as a CSV file with two columns and 6 rows: 
#'    type,filenames
#'    analyteType1,myfilename
#'    analyteType2,myfilename
#'    analyteType1MetaData,myfilename (optional)
#'    analyteType2MetaData,myfilename (optional)
#'    sampleMetaData,myfilename
#'
#' Note that all files supplied in the CSV file, and the CSV file itself should 
#' be placed in the same folder.  The software assumes will automatically retrieve 
#' the file path of
#' the input files (based on location of CSV files).  
#' Note also that the input data files should be in a specific format:
#'	analyteType1/2: rows are analytes, columns are samples
#'	analyteType1/2MetaData: rows are analytes, features are columns
#'	sampleMetaData: rows are samples, features are columns
#' In addition, the first column of the sampleMetaData file is assumed to be the sample id, 
#' and those sample ids should match the columns of analyteType1/2 (e.g. it is required
#' that all sample ids in the analyteType1/2 are also in the sampleMetaData).
#'
#' @include internalfunctions.R
#'
#' @param inputFile input file in CSV format (see Despcription)
#' @param analyteType1id name of column from Analyte Type 1 meta data to be used as id
#'      (required if an Analyte Type 1 meta data file is present, 
#'      must match Analyte Type 1 data)
#' @param analyteType2id name of column from Analyte Type 2 meta data to be used as id
#'      (required if an Analyte Type 2 meta data file is present, 
#'      must match Analyte Type 2 data)
#' @param logAnalyteType1 whether or not to log values for Analyte Type 1(T/F)
#' @param logAnalyteType2 whether or not to log values for Analyte Type 2(T/F)
#' @param class.feat class ("factor" or "numeric") for each covariate. The following format
#' is required: list(covar1="numeric", covar2="factor")
#' @param suppressWarnings whether or not to suppress warnings
#' @return IntLimData object with input data
#'
#' @export
ReadData <- function(inputFile,analyteType1id="id",analyteType2id="id", 
                     logAnalyteType1=FALSE,logAnalyteType2=FALSE, class.feat = NULL,
                     suppressWarnings = FALSE){
      
  # Initialize output.
  intlimData <- NULL
  
  # Check that file exists
  if (!file.exists(inputFile)) {
    stop("CSV input file does not exist")
  }
  # Make into df to make access easier
  csvfile <- as.data.frame(utils::read.csv(inputFile, header=TRUE,row.names=1))
  type1Data <- NULL
  type2Data <- NULL
  
  # Check column names are correct
  if (colnames(csvfile)!="filenames") {
    stop("Check column names of input files.  'type' and 'filenames' are required")
  }

  # Check that all types required are present
  mytypes <- c("analyteType1","analyteType2","analyteType1MetaData",
               "analyteType2MetaData",
               "sampleMetaData")
  mymatches <- as.numeric(lapply(mytypes,function(x) 
    length(which(rownames(csvfile)==x))))
  if(sum(mymatches)!=5) {
    stop(paste("The column 'type' contains non-allowed entries (See Description). The",
               "CSV input file must contain 6 rows (if optional meta data files for analytes",
               "are not to be input, have the corresponding filenames be blanks."))
  }
  
  mydir <- base::dirname(inputFile)
  
  if((is.na(as.character(csvfile['analyteType1',])) && is.na(as.character(csvfile['analyteType2',])))
     || ((as.character(csvfile['analyteType1',]) == "") && as.character(csvfile['analyteType2',]) == "")){
    stop("No data provided.")
  }
  else{
    # Read the data for Analyte Type 2.
    if(is.na(csvfile['analyteType2',]) || as.character(csvfile['analyteType2',])=="") {
      if(suppressWarnings == FALSE){
        warning(paste("No data provided for Analyte Type 2. This means you cannot run",
                      "analyses involving this analyte type.\n"))
        ;type2Data<-NULL;
      }
      else{
        type2Data<-NULL
      }
    }else{
      temp <- paste0(mydir,"/",as.character(csvfile['analyteType2',]))
      if(!file.exists(temp)) {
        stop(paste("File", temp, "does not exist"))
      } 
      else {
        ids <- utils::read.csv(temp,check.names=F)[,1]
        if(length(ids) != length(unique(ids))) {
          stop(paste("Error: your input file",temp,"has duplicate", 
                     "entries in column 1. Please make sure you have one row per", 
                     "analyte"))
        } 
        else {
          type2Data<-utils::read.csv(temp,row.names = 1,check.names=F)
        }
      }
    }
    # Read the data for Analyte Type 1.
    if(is.na(csvfile['analyteType1',]) || as.character(csvfile['analyteType1',])=="") {
      if(suppressWarnings == FALSE){
        warning(paste("No data provided for Analyte Type 1. This means you cannot run",
                      "analyses involving this analyte type.\n"))
        ;type1Data<-NULL;
      }
      else{
        type1Data<-NULL
      }
    }else {
      temp <- paste0(mydir,"/",as.character(csvfile['analyteType1',]))
      if(!file.exists(temp)) {
        stop(paste("File", temp, "does not exist"))
      } else {
        ids <- utils::read.csv(temp,check.names=F)[,1]
        if(length(ids) != length(unique(ids))) {
          stop(paste("Error: your input file",temp,"has duplicate", 
                     "entries in column 1. Please make sure you have one row per analyte"))
        } else {
          type1Data<-utils::read.csv(temp,row.names = 1,check.names=F)
        }
      }
    }

    # Read the Analyte Type 2 metadata.
    temp <- paste0(mydir,"/",as.character(csvfile['analyteType2MetaData',]))
    if(as.character(csvfile['analyteType2MetaData',])=="" ||
       as.character(csvfile['analyteType2MetaData',])==mydir) {
      if(suppressWarnings == FALSE){
        warning("No metadata provided for Analyte Type 2");type2MetaData<-NULL;analyteType2id=NULL; 
      }
      else{
        type2MetaData<-NULL
        analyteType2id=NULL
      }
    }else if(!file.exists(temp)) {
      stop(paste("File", temp, "does not exist"))
    } else {
      type2MetaData<-utils::read.csv(temp)
      colnames(type2MetaData)[which(colnames(type2MetaData)==analyteType2id)]="id"
    }

    # Read the Analyte Type 1 metadata.
    temp <- paste0(mydir,"/",as.character(csvfile['analyteType1MetaData',]))
    if(as.character(csvfile['analyteType1MetaData',])==""|| 
       as.character(csvfile['analyteType1MetaData',])==mydir) {
      if(suppressWarnings == FALSE){
        warning("No metadata provided for Analyte Type 1");type1MetaData<-NULL;
        analyteType1id=NULL; 
      }
      else{
        type1MetaData<-NULL
        analyteType1id=NULL
      }
    }else if(!file.exists(temp)) {
      stop(paste("File", temp, "does not exist"))
    } else {
      type1MetaData<-utils::read.csv(temp)
      colnames(type1MetaData)[which(colnames(type1MetaData)==analyteType1id)]="id"
    }

    # Read the sample data.
    temp <- paste0(mydir,"/",as.character(csvfile['sampleMetaData',]))
    if(!file.exists(temp)) {
      stop(paste("File", temp, "does not exist"))
    } else {
      pData<-utils::read.csv(temp,row.names = 1)
    }
    
    # Check that Analyte Type 1 ID's match metadata.
    if(!is.null(type1Data)){
      
      # Check that data and metadata match.
      if (!is.null(type1MetaData)) {
        if(length(which(colnames(type1MetaData)=='id'))!=1) {
          stop(paste("analyteType1id provided",analyteType1id,"does not exist in",
                     "Analyte Type 1 meta data file"))
        } else if (length(intersect(rownames(type1Data),as.character(type1MetaData[,"id"])))
                   <nrow(type1Data)){
          stop("Analytes in Type 1 data file and meta data files are not equal")
        }
        rownames(type1MetaData)=as.character(type1MetaData[,'id'])
      }
      
      # Make sure order of feature data is the same as the data matrix:
      if(is.null(type1MetaData)) {
        type1MetaData <- data.frame(id = rownames(type1Data),stringsAsFactors=FALSE)
        rownames(type1MetaData) <- type1MetaData[,1]
      }
      myind=as.numeric(lapply(rownames(type1Data),function(x) which(type1MetaData[,"id"]==x)))
      type1MetaData <- data.frame(type1MetaData[myind,],stringsAsFactors=FALSE)
      rownames(type1MetaData) <- type1MetaData[,1]
    }
    
    # Check that Analyte Type 2 ID's match metadata.
    if(!is.null(type2Data)){
      
      # Check that data and metadata match.
      if (!is.null(type2MetaData)) {
        if(length(which(colnames(type2MetaData)=='id'))!=1) {
          stop(paste("analyteType2id provided",analyteType2id,"does not exist in",
                     "Analyte Type 2 meta data file"))
        } else if (length(intersect(rownames(type2Data),as.character(type2MetaData[,"id"])))
                   <nrow(type2Data)){
          stop("Analytes in Type 2 data file and meta data files are not equal")
        }
        rownames(type2MetaData)=as.character(type2MetaData[,'id'])
      }
      
      # Make sure order of feature data is the same as the data matrix:
      if(is.null(type2MetaData)) {
        type2MetaData <- data.frame(id = rownames(type2Data),stringsAsFactors=FALSE)
        rownames(type2MetaData) <- type2MetaData[,1]
      }
      myind=as.numeric(lapply(rownames(type2Data),function(x) which(type2MetaData[,"id"]==x)))
      type2MetaData <- data.frame(type2MetaData[myind,],stringsAsFactors=FALSE)
      rownames(type2MetaData) <- type2MetaData[,1]
    }
    
    # Log data if applicable.
    cutoff <- 0.0000001
    if (logAnalyteType1 == TRUE && !is.null(type1Data)){
      if(min(type1Data) < 0){
        warning("Analyte Type 1 data has negative values. Continuing without log-scaling.")
      }else{
        type1Data[which(type1Data == 0)] <- cutoff
        type1Data <- log2(type1Data)
      }
    }
    if (logAnalyteType2 == TRUE && !is.null(type2Data)){
      if(min(type2Data) < 0){
        warning("Analyte Type 2 data has negative values. Continuing without log-scaling.")
      }else{
        type2Data[which(type2Data == 0)] <- cutoff
        type2Data <- log2(type2Data)
      }
    }
    
    # Get common samples.
    pDataOld <- pData
    type1DataOld <- type1Data
    type2DataOld <- type2Data
    myind <- rownames(pData)
    if(!is.null(type1Data)){
      type1Samps <- colnames(type1Data)
      myind <- intersect(myind, type1Samps)
    }
    if(!is.null(type2Data)){
      type2Samps <- colnames(type2Data)
      myind <- intersect(myind, type2Samps)
    }
    pData<-pData[myind,]
    pDataSpecific <- setdiff(rownames(pDataOld), myind)
    if(length(pDataSpecific) > 0){
      warning(paste("The following samples were only included in the sample data",
                    "and were removed:", pDataSpecific))
    }
    if(!is.null(type1Data)){
      type1Data <- type1Data[,myind]
      type1Specific <- setdiff(colnames(type1DataOld), myind)
      if(length(type1Specific) > 0){
        warning(paste("The following samples were only included in the Analyte 1 data",
                      "and were removed:", type1Specific))
      }
    }
    if(!is.null(type2Data)){
      type2Data <- type2Data[,myind]
      type2Specific <- setdiff(colnames(type2DataOld), myind)
      if(length(type2Specific) > 0){
        warning(paste("The following samples were only included in the Analyte 2 data",
                      "and were removed:", type2Specific))
      }
    }
    
    # Fix names.
    if(!is.null(type2Data)){
      colnames(type2Data) <- make.names(colnames(type2Data))
      rownames(type2Data) <- make.names(rownames(type2Data))
    }
    if(!is.null(type1Data)){
      colnames(type1Data) <- make.names(colnames(type1Data))
      rownames(type1Data) <- make.names(rownames(type1Data))
    }
    if(!is.null(type1MetaData)){
      colnames(type1MetaData) <- make.names(colnames(type1MetaData))
      rownames(type1MetaData) <- make.names(rownames(type1MetaData))
    }
    if(!is.null(type2MetaData)){
      colnames(type2MetaData) <- make.names(colnames(type2MetaData))
      rownames(type2MetaData) <- make.names(rownames(type2MetaData))
    }
    colnames(pData) <- make.names(colnames(pData))
    rownames(pData) <- make.names(rownames(pData))
    
    # Extract sampleMetaData and coerce to numeric or factor.
    # Coerce sampleMetaData classes.
    names(class.feat) <- make.names(names(class.feat))
    covarMatrix <- pData[,names(class.feat)]
    for(covar in names(class.feat)){
      if(class.feat[covar] == "numeric"){
        covarMatrix[,covar] <- as.numeric(as.character(covarMatrix[,covar]))
      }else if(class.feat[covar] == "factor"){
        covarMatrix[,covar] <- as.factor(as.character(covarMatrix[,covar]))
      }else{
        stop(paste(class.feat[covar], "is not a valid class for covariate", covar))
      }
    }
    
    # Convert data to matrix.
    if(!is.null(type1Data)){
      type1Data <- as.matrix(type1Data)
    }else{
      type1Data <- matrix(, nrow = 0, ncol = 0)
    }
    if(!is.null(type2Data)){
      type2Data <- as.matrix(type2Data)
    }else{
      type2Data <- matrix(, nrow = 0, ncol = 0)
    }
    
    # Set null metadata to data frame.
    if(is.null(type1MetaData)){
      type1MetaData <- as.data.frame(matrix(, nrow = 0, ncol = 0))
    }
    if(is.null(type2MetaData)){
      type2MetaData <- as.data.frame(matrix(, nrow = 0, ncol = 0))
    }
    
    intlimData <- methods::new("IntLimData",analyteType1=type1Data,
                               analyteType2=type2Data,
                               analyteType1MetaData = type1MetaData,
                               analyteType2MetaData = type2MetaData,
                               sampleMetaData = covarMatrix)
    
    print("IntLIMData created")
  }
  return(intlimData)
}