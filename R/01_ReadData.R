#' Read in CSV file
#'
#' The metadata associated with data files to be analyzed in IntLim is supplied
#' as a CSV file with two columns and 6 rows: 
#'    type,filenames
#'    metabData,myfilename
#'    geneData,myfilename
#'    metabMetaData,myfilename (optional)
#'    geneMetaData,myfilename (optional)
#'    sampleMetaData,myfilename
#'
#' Note that all files supplied in the CSV file, and the CSV file itself should be placed in the same folder.  The software assumes will automatically retrieve the file path of
#' the input files (based on location of CSV files).  
#' Note also that the input data files should be in a specific format:
#'	metabData: rows are metabolites, columns are samples
#'	geneData: rows are genes, columns are samples
#'	metabMetaData: rows are metabolites, features are columns
#'	geneMetaData: rows are genes, features are columns
#'	sampleMetaData: rows are samples, features are columns
#' In addition, the first column of the sampleMetaData file is assumed to be the sample id, 
#' and those sample ids should match the columns of metabData and geneData (e.g. it is required
#' that all sample ids in the metabData and geneData are also in the sampleMetaDatafile).
#'
#' @include MetaboliteSet_addMetabolite.R
#' @include internalfunctions.R
#'
#' @param inputFile input file in CSV format (see Despcription)
#' @param metabid name of column from metabolite meta data to be used as id
#'      (required if a metabolite meta dadta file is present, must match metabolite abundances data)
#' @param geneid name of column from gene meta data to be used as id
#'	(required if a gene meta data file is present, must match gene expression data)
#' @param logmetab whether or not to log metabolite values (T/F)
#' @param loggene whether or not to log gene values (T/F)
#' @return MultiDataSet object with input data
#'
#' @examples
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(inputFile = csvfile,metabid='id',geneid='id')
#' @export
ReadData <- function(inputFile,metabid=NULL,geneid=NULL, logmetab=FALSE,loggene=FALSE){
      
    #Create Multi
    pieces <- CreateIntLimObjectPieces(inputFile,metabid,geneid, logmetab,loggene)
    GMdata <- CreateIntLimObject(genefdata=pieces[["GmetaData"]], 
                                 metabfdata=pieces[["MmetaData"]],
                                 metabid=pieces[["metabid"]],
                                 geneid=pieces[["geneid"]],
                                 pdata=pieces[["pData"]],
                                 metabdata=pieces[["MData"]],
                                 genedata=pieces[["GData"]],
                                 logmetab=pieces[["logmetab"]],
                                 loggene=pieces[["loggene"]])
    

    print("CreateMultiDataSet created")
    return(GMdata)
}

#' This is the helper function for ReadData. It is also used in CreateCrossValFolds.
#'
#' @include MetaboliteSet_addMetabolite.R
#' @include internalfunctions.R
#'
#' @param inputFile input file in CSV format (see Despcription)
#' @param metabid name of column from metabolite meta data to be used as id
#'      (required if a metabolite meta dadta file is present, must match metabolite abundances data)
#' @param geneid name of column from gene meta data to be used as id
#'	(required if a gene meta data file is present, must match gene expression data)
#' @param logmetab whether or not to log metabolite values (T/F)
#' @param loggene whether or not to log gene values (T/F)
#' @param suppressWarnings whether to suppress warnings
#' @return Named list of components for MultiDataSet
CreateIntLimObjectPieces <- function(inputFile,metabid=NULL,geneid=NULL, logmetab=FALSE,loggene=FALSE,
                                     suppressWarnings=FALSE){
  # Check that file exists
  if (!file.exists(inputFile)) {
    stop("CSV input file does not exist")
  }
  # Make into df to make access easier
  csvfile <- as.data.frame(utils::read.csv(inputFile, header=TRUE,row.names=1))
  GMData <- NULL
  GData <- NULL
  MData <- NULL
  
  # Check column names are correct
  if (colnames(csvfile)!="filenames") {
    stop("Check column names of input files.  'type' and 'filenames' are required")
  }
  
  # Check that all types required are present
  mytypes <- c("metabData","geneData","metabMetaData","geneMetaData",
               "sampleMetaData")
  mymatches <- as.numeric(lapply(mytypes,function(x) 
    length(which(rownames(csvfile)==x))))
  if(sum(mymatches)!=5) {
    stop(paste("The column 'type' contains non-allowed entries (See Description). The",
               "CSV input file must contain 6 rows (if optional meta data files for metabolites",
               "and genes are not to be input, have the corresponding filenames be blanks."))
  }
  
  mydir <- base::dirname(inputFile)
  
  if((is.na(as.character(csvfile['geneData',])) && is.na(as.character(csvfile['metabData',])))
     || ((as.character(csvfile['geneData',]) == "") && as.character(csvfile['metabData',]) == "")){
    stop("No data provided.")
  }
  else{
    # Check that files exist then read them in one by one
    if(is.na(csvfile['metabData',]) || as.character(csvfile['metabData',])=="") {
      if(suppressWarnings == FALSE){
        warning(paste("No data provided for metabolites. This means you cannot run",
                      "metabolite-metabolite or gene-metabolite analyses.\n"))
        ;MData<-NULL;
      }
      else{
        MData<-NULL
      }
    }else{
      temp <- paste0(mydir,"/",as.character(csvfile['metabData',]))
      if(!file.exists(temp)) {
        stop(paste("File", temp, "does not exist"))
      } 
      else {
        ids <- utils::read.csv(temp,check.names=F)[,1]
        if(length(ids) != length(unique(ids))) {
          stop(paste("Error: your input file",temp,"has duplicate", 
                       "entries in column 1. Please make sure you have one row per", 
                       "metabolite"))
        } 
        else {
          MData<-utils::read.csv(temp,row.names = 1,check.names=F)
        }
      }
    }
    if(is.na(csvfile['geneData',]) || as.character(csvfile['geneData',])=="") {
      if(suppressWarnings == FALSE){
        warning(paste("No data provided for genes. This means you cannot run",
                      "gene-gene or gene-metabolite analyses.\n"))
        ;GData<-NULL;
      }
      else{
        GData<-NULL
      }
    }else {
      temp <- paste0(mydir,"/",as.character(csvfile['geneData',]))
      if(!file.exists(temp)) {
        stop(paste("File", temp, "does not exist"))
      } else {
        ids <- utils::read.csv(temp,check.names=F)[,1]
        if(length(ids) != length(unique(ids))) {
          stop(paste("Error: your input file",temp,"has duplicate", 
                         "entries in column 1. Please make sure you have one row per gene"))
        } else {
          GData<-utils::read.csv(temp,row.names = 1,check.names=F)
        }
      }
    }
    
    temp <- paste0(mydir,"/",as.character(csvfile['metabMetaData',]))
    if(as.character(csvfile['metabMetaData',])=="" || as.character(csvfile['metabMetaData',])==mydir) {
      if(suppressWarnings == FALSE){
        warning("No metadata provided for metabolites");MmetaData<-NULL;metabid=NULL; 
      }
      else{
        MmetaData<-NULL
        metabid=NULL
      }
    }else if(!file.exists(temp)) {
      stop(paste("File", temp, "does not exist"))
    } else {
        MmetaData<-utils::read.csv(temp)
        colnames(MmetaData)[which(colnames(MmetaData)==metabid)]="id"
    }
    
    temp <- paste0(mydir,"/",as.character(csvfile['geneMetaData',]))
    if(as.character(csvfile['geneMetaData',])==""|| as.character(csvfile['metabMetaData',])==mydir) {
      if(suppressWarnings == FALSE){
        warning("No metadata provided for genes");GmetaData<-NULL;geneid=NULL; 
      }
      else{
        GmetaData<-NULL
        metabid=NULL
      }
    }else if(!file.exists(temp)) {
      stop(paste("File", temp, "does not exist"))
    } else {
        GmetaData<-utils::read.csv(temp)
        colnames(GmetaData)[which(colnames(GmetaData)==geneid)]="id"
    }
    temp <- paste0(mydir,"/",as.character(csvfile['sampleMetaData',]))
    if(!file.exists(temp)) {
      stop(paste("File", temp, "does not exist"))
    } else {
        pData<-utils::read.csv(temp,row.names = 1)
    }
    
    # Return data.
    pieces <- list("GmetaData" = GmetaData, "MmetaData" = MmetaData, 
                   "metabid" = metabid, "geneid" = geneid, "pData" = pData,
                   "MData" = MData, "GData" = GData, "logmetab" = logmetab,
                   "loggene" = loggene)
  }
  return(pieces)
}

#' Filter input data by abundance values (gene and metabolite data) and number of missing values (metabolite data only).
#'
#' Filter data by abundance (with user-input percentile cutoff) of missing values (with user-input percent cutoff). Missing values are commonly found in metabolomics data so the parameter currently only applies to metabolomics data.
#'
#' @param inputFile input file in CSV format (see Despcription)
#' @param metabid name of column from metabolite meta data to be used as id
#'      (required if a metabolite meta dadta file is present, must match metabolite abundances data)
#' @param geneid name of column from gene meta data to be used as id
#'	(required if a gene meta data file is present, must match gene expression data)
#' @param logmetab whether or not to log metabolite values (T/F)
#' @param loggene whether or not to log gene values (T/F)
#' @param folds number of folds to create
#' @param suppressWarnings whether to suppress warnings
#' @return A set of MultiDataSet training and testing sets, of the following format:
#' list(list("train" = MultiDataSet, "test" = MultiDataSet), ... list("train" = MultiDataSet,
#' "test" = MultiDataSet))
#' @export
CreateCrossValFolds <- function(inputFile,metabid=NULL,geneid=NULL, 
                                logmetab=FALSE,loggene=FALSE,
                                folds, suppressWarnings=FALSE) {
  
  # Create the components of the input.
  pieces <- CreateIntLimObjectPieces(inputFile,metabid,geneid, logmetab,loggene,
                                     suppressWarnings)
  
  # Extract all samples.
  samps <- rownames(pieces[["pData"]])
  
  # Stop if the number of folds is greater than the number of samples.
  if(folds > length(samps)){
    stop("ERROR: The number of folds is greater than the number of samples!")
  }

  # Permute samples and divide into folds.
  sets_of <- ceiling(length(samps) / folds)
  perm_samps <- sample(samps, length(samps), replace = FALSE)
  fold_samps <- split(perm_samps, ceiling(seq_along(perm_samps)/sets_of))

  # For each fold, extract the samples from the data set.
  trainTestObjects <- lapply(fold_samps, function(fold){
    # Initialize.
    GmetaData <- pieces$GmetaData
    MmetaData <- pieces$MmetaData
    metabid <- pieces$metabid
    geneid <- pieces$geneid
    pData_train <- pieces$pData[setdiff(samps, fold),]
    pData_test <- pieces$pData[fold,]
    MData_train <- NULL
    MData_test <- NULL
    GData_train <- NULL
    GData_test <- NULL
    logmetab_train <- NULL
    logmetab_test <- NULL
    loggene_train <- NULL
    loggene_test <- NULL
    
    # Include all but the current fold in the training data.
    if(!is.null(pieces$MData)){
      MData_train <- pieces$MData[,setdiff(samps, fold)]
      MData_test <- as.data.frame(pieces$MData[,fold])
      logmetab_train <- pieces$logmetab
      logmetab_test <- pieces$logmetab
      rownames(MData_test) <- rownames(pieces$MData)
      colnames(MData_test) <- rownames(pData_test)
    }
    if(!is.null(pieces$GData)){
      GData_train <- pieces$GData[,setdiff(samps, fold)]
      GData_test <- as.data.frame(pieces$GData[,fold])
      loggene_train <- pieces$loggene
      loggene_test <- pieces$loggene
      rownames(GData_test) <- rownames(pieces$GData)
      colnames(GData_test) <- rownames(pData_test)
    }
    
    training <- CreateIntLimObject(genefdata=GmetaData, 
                                 metabfdata=MmetaData,
                                 metabid=metabid,
                                 geneid=geneid,
                                 pdata=pData_train,
                                 metabdata=MData_train,
                                 genedata=GData_train,
                                 logmetab=logmetab,
                                 loggene=loggene)
    

    # Include the current fold in the testing data.
    testing <- CreateIntLimObject(genefdata=GmetaData, 
                                   metabfdata=MmetaData,
                                   metabid=metabid,
                                   geneid=geneid,
                                   pdata=pData_test,
                                   metabdata=MData_test,
                                   genedata=GData_test,
                                   logmetab=logmetab,
                                   loggene=loggene)
    
    return(list("training"=training, "testing"=testing))
  })
  return(trainTestObjects)
}

