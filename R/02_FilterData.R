#' Filter input data by abundance values (analyte data) and number of missing 
#' values.
#'
#' Filter data by abundance (with user-input percentile cutoff) of missing values 
#' (with user-input percent cutoff). Missing values are commonly found in 
#' metabolomics data.
#' @param csvfile The path to the CSV file that contains the full dataset.
#' @param inputDataFolds List of IntLimData objects (output of ReadData()) with analylte levels and
#'  associated meta-data
#' @param logAnalyteType1 whether or not to log values for Analyte Type 1(T/F).
#' Must coincide with parameter used during CreateCrossValFolds.
#' @param logAnalyteType2 whether or not to log values for Analyte Type 2(T/F).
#' Must coincide with parameter used during CreateCrossValFolds.
#' @param class.feat class ("factor" or "numeric") for each covariate. The following format
#' is required: list(covar1="numeric", covar2="factor"). 
#' Must coincide with parameter used during CreateCrossValFolds.
#' @param analyteType1id name of column from Analyte Type 1 meta data to be used as id
#'      (required if an Analyte Type 1 meta data file is present, 
#'      must match Analyte Type 1 data), Must coincide with parameter used during 
#'      CreateCrossValFolds.
#' @param analyteType2id name of column from Analyte Type 2 meta data to be used as id
#'      (required if an Analyte Type 2 meta data file is present, 
#'      must match Analyte Type 2 data). Must coincide with parameter used during 
#'      CreateCrossValFolds.
#' @param analyteType1perc percentile cutoff (0-1) for filtering analyte type 1 (e.g. 
#' remove analytes with mean values < 'analyteType1perc' percentile) (default: 0)
#' @param analyteType2perc percentile cutoff (0-1) for filtering analyte type 2 
#' (default: no filtering of analytes) (default:0)
#' @param analyteType2miss missing value percent cutoff (0-1) for filtering analyte type 2 
#' (analytes with > 80\% missing values will be removed) (default:0)
#' @return filtData IntLimData object with input data after filtering
#' @export
FilterDataFolds <- function(csvfile,inputDataFolds,analyteType1perc=0,
                            analyteType2perc=0, analyteType2miss=0, analyteType1id="id",
                            analyteType2id="id", logAnalyteType1 = F,
                            logAnalyteType2 = F, class.feat) {
	# First, read the original data set.
	inputData <- ReadData(csvfile,analyteType1id,analyteType2id, logAnalyteType1,
	                                logAnalyteType2, class.feat)
	filtdata <- FilterData(inputData,analyteType1perc=analyteType1perc,
	                       analyteType2perc=analyteType2perc, 
	                       analyteType2miss=analyteType2miss)
	
	# Filter the original data set.
	type1keepers = NULL
	type2keepers = NULL
	if(length(filtdata@analyteType1)>0){
	  type1Data <- filtdata@analyteType1
	  type1keepers <- rownames(type1Data)
	}
	if(length(filtdata@analyteType2)>0){
	  type2Data <- filtdata@analyteType2
	  type2keepers <- rownames(type2Data)
	}

	# Now, filter each fold so that the analytes in the set match the analytes.
	for(i in 1:length(inputDataFolds)){

		# Filter training and testing data.
		inputDataFolds[[i]]$training <- FilterFromKeepers(inputData = inputDataFolds[[i]]$training,
		                                                  type2keepers = type2keepers,
		                                                  type1keepers = type1keepers)
		inputDataFolds[[i]]$testing <- FilterFromKeepers(inputData = inputDataFolds[[i]]$testing,
		                                                 type2keepers = type2keepers,
		                                                 type1keepers = type1keepers)
	}
	return(inputDataFolds)
}

#' Filter a data set given the analytes to keep.
#'
#' @param inputData IntLimData object (output of ReadData()) with analyte leves
#'  and associated meta-data
#' @param stype phenotype or outcome of interest.
#' @param type2keepers Type 2 analytes to keep.
#' @param type1keepers Type 1 analytes to keep.
#' @param stype.type Either "factor" or "numeric". The outcome type.
#' @param class.covar class ("numeric" or "character") for each covariate. The following format
#' is required: list(covar1="numeric", covar2="character")
#' @return filtData IntLimData object with input data after filtering
FilterFromKeepers <- function(inputData,stype,type1keepers, type2keepers, stype.type,
                              class.covar = NULL) {
  # Initialize
  type1Data <- NULL
  type2Data <- NULL
  
  # Filter analyte type 1
  if(length(inputData@analyteType1)>0){
    # Filter the assay data.
    type1Data <- inputData@analyteType1
    samps <- colnames(type1Data)
    type1Data <- as.matrix(type1Data[type1keepers,])
    colnames(type1Data) <- samps
  }

  # Filter analyte type 2
  if(length(inputData@analyteType2)>0){
    # Filter the assay data.
    type2Data <- inputData@analyteType2
    samps <- colnames(type2Data)
    type2Data <- as.matrix(type2Data[type2keepers,])
    colnames(type2Data) <- samps
  }

  # Create new data set.
  filteredData <- methods::new("IntLimData", analyteType1=type1Data,
               analyteType2=type2Data,
               analyteType1MetaData = inputData@analyteType1MetaData,
               analyteType2MetaData = inputData@analyteType2MetaData,
               sampleMetaData = inputData@sampleMetaData)
  return(filteredData)
}

#' Filter input data by abundance values and number of missing values.
#'
#' Filter data by abundance (with user-input percentile cutoff) of missing values 
#' (with user-input percent cutoff). Missing values are commonly found in metabolomics 
#' data.
#'
#' @param inputData IntLimData object (output of ReadData()) with analylte levels and
#'  associated meta-data
#' @param analyteType1perc percentile cutoff (0-1) for filtering analyte type 1 (e.g. 
#' remove analytes with mean values < 'analyteType1perc' percentile) (default: 0)
#' @param analyteType2perc percentile cutoff (0-1) for filtering analyte type 2 
#' (default: no filtering of analytes) (default:0)
#' @param analyteType2miss missing value percent cutoff (0-1) for filtering analyte type 2 
#' (analytes with > 80\% missing values will be removed) (default:0)
#' @return filtData IntLimData object with input data after filtering
#' @export
FilterData <- function(inputData,analyteType1perc=0,analyteType2perc=0, analyteType2miss=0) {

  # Check that input is a IntLimData
  if (class(inputData) != "IntLimData") {
	   stop("input data is not a IntLimData class")
  }
  filtdata <- NULL
  if(length(inputData@analyteType1)==0 && length(inputData@analyteType2)==0) {
	   stop("input data must contain assayData of at least one type of analyte.
	     Try reading in the data with the ReadData function")
  }	else if(length(inputData@analyteType1)>0 && length(inputData@analyteType2)>0){
    if(!is.null(analyteType1perc) && analyteType1perc > 1) {
      stop("analyteType1perc parameter must be between 0 and 1")
    }
    if(!is.null(analyteType2perc) && analyteType2perc > 1) {
      stop("analyteType2perc parameter must be between 0 and 1")
    }
    if(!is.null(analyteType2miss) && analyteType2miss > 1) {
      stop("analyteType2miss parameter must be between 0 and 1")
    }
  
    # Check that at least one parameter is not null
    len <- length(c(analyteType1perc,analyteType2perc,analyteType2miss))
    if ((analyteType1perc+analyteType2perc+analyteType2miss) ==0) {
      warning("All filtering parameters are NULL so the data remains unfiltered")
      filtdata <- inputData
    } else {
      type1Data <- inputData@analyteType1
      type2Data <- inputData@analyteType2
      if(analyteType1perc > 0) {
        if(analyteType1perc>1) {analyteType1perc=analyteType1perc}
        mymean <- as.numeric(apply(type1Data,1, function(x) mean(x,na.rm=T)))
        keepers <- which(mymean > stats::quantile(mymean,analyteType1perc))
        type1Data <- type1Data[keepers,]
      } else {
        print("No filtering by percentile is applied for analyte type 1")
      }
      if(analyteType2perc > 0) {
        if(analyteType2perc>1) {analyteType2perc=analyteType2perc}
          mymean <- as.numeric(apply(type2Data,1, function(x) mean(x,na.rm=T)))
          keepers <- which(mymean > stats::quantile(mymean,analyteType2perc))
          type2Data <- type2Data[keepers,]
      } else {
        print("No filtering by percentile is applied for analyte type 2")
      }
      if(analyteType2miss > 0) {
        missnum <- as.numeric(apply(type2Data,1,function(x) length(which(x==min(x,na.rm=T)))))-1
        mycut <- analyteType2miss * ncol(type2Data)
        keepers <- which(missnum < mycut)
        type2Data <- type2Data[keepers,]
      } else {
        print("No filtering by missing values is applied for analyte type 2")
      }
    
      # Now reconstruct the multidataset object
      filtdata <- methods::new("IntLimData", analyteType1=type1Data,
                                              analyteType2=type2Data,
                                              analyteType1MetaData = inputData@analyteType1MetaData,
                                              analyteType2MetaData = inputData@analyteType2MetaData,
                                              sampleMetaData = inputData@sampleMetaData)
    }
  }else if(length(inputData@analyteType1)>0){
    if(!is.null(analyteType1perc) && analyteType1perc > 1) {
      stop("analyteType1perc parameter must be between 0 and 1")
    }
    
    # Check that at least one parameter is not null
    len <- length(c(analyteType1perc,analyteType2perc,analyteType2miss))
    if (is.null(analyteType1perc)) {
      warning("All filtering parameters are NULL so the data remains unfiltered")
      filtdata <- inputData
    }
    else {
      type1Data <- inputData@analyteType1
      if(analyteType1perc > 0) {
        if(analyteType1perc>1) {analyteType1perc=analyteType1perc}
        mymean <- as.numeric(apply(type1Data,1, function(x)
          mean(x,na.rm=T)))
        keepers <- which(mymean > stats::quantile(mymean,analyteType1perc))
        type1Data <- type1Data[keepers,]
      } else {
        print("No percentile filtering is applied for analyte type 1")
      }
      
      # Now reconstruct the multidataset object
      filtdata <- methods::new("IntLimData", phenotype=inputData@phenotype,
                               analyteType1=type1Data,
                               analyteType2=NULL,
                               analyteType1MetaData = inputData@analyteType1MetaData,
                               analyteType2MetaData = NULL,
                               sampleMetaData = inputData@sampleMetaData,
                               stype=inputData@stype,
                               stype.type=inputData@stype.type)
    }
  } else if(length(inputData@analyteType2)>0){
    if(!is.null(analyteType2perc) && analyteType2perc > 1) {
      stop("analyteType2perc parameter must be between 0 and 1")
    }
    if(!is.null(analyteType2miss) && analyteType2miss > 1) {
      stop("analyteType2miss parameter must be between 0 and 1")
    }
    
    # Check that at least one parameter is not null
    len <- length(c(analyteType2perc,analyteType2miss))
    if ((analyteType2perc+analyteType2miss) ==0) {
      warning("All filtering parameters are NULL so the data remains unfiltered")
      filtdata <- inputData
    }
    else {
      type2Data <- inputData@analyteType2
      if(analyteType2perc > 0) {
        if(analyteType2perc>1) {analyteType2perc=analyteType2perc}
        mymean <- as.numeric(apply(type2Data,1, function(x)
          mean(x,na.rm=T)))
        keepers <- which(mymean > stats::quantile(mymean,analyteType2perc))
        type2Data <- type2Data[keepers,]
      } else {
        print("No filtering by percentile is applied for analyte type 2")
      }
      if(analyteType2miss > 0) {
        missnum <- as.numeric(apply(type2Data,1,function(x) length(which(x==min(x,na.rm=T)))))-1
        mycut <- analyteType2miss * ncol(type2Data)
        keepers <- which(missnum < mycut)
        type2Data <- type2Data[keepers,]
      } else {
        print("No filtering by missing values is applied for analyte type 2")
      }
      
      # Now reconstruct the multidataset object
      filtdata <- methods::new("IntLimData", phenotype=inputData@phenotype,
                               analyteType1=NULL,
                               analyteType2=type2Data,
                               analyteType1MetaData = NULL,
                               analyteType2MetaData = inputData@analyteType2MetaData,
                               sampleMetaData = inputData@sampleMetaData,
                               stype=inputData@stype,
                               stype.type=inputData@stype.type)
    }
  }
  
	return(filtdata)
}

