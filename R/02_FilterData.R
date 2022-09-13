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
#' @param analyteMiss missing value percent cutoff (0-1) for filtering both analyte types 
#' (analytes with > 80\% missing values will be removed) (default:0)
#' @param cov.cutoff percentile cutoff (0-1) for the covariances of the anaytes (default: 0.30)
#' @param suppressWarnings whether or not to print warnings. If TRUE, warnings will
#' not be printed.
#' @return filtData IntLimData object with input data after filtering
#' @export
FilterData <- function(inputData,analyteType1perc=0,analyteType2perc=0, analyteMiss=0,
                       suppressWarnings = FALSE, cov.cutoff=0) {
  # Check that input is a IntLimData
  if (!methods::is(inputData, "IntLimData")) {
    stop("input data is not a IntLimData class")
  }
  filtdata <- NULL
  if(length(inputData@analyteType1)==0 && length(inputData@analyteType2)==0) {
    stop("input data must contain assayData of at least one type of analyte.
	     Try reading in the data with the ReadData function")
  }	else if(length(inputData@analyteType1)>0 && length(inputData@analyteType2)>0){
    if(analyteType1perc < 0 || analyteType1perc > 1) {
      stop("analyteType1perc parameter must be between 0 and 1")
    }
    if(analyteType2perc < 0 || analyteType2perc > 1) {
      stop("analyteType2perc parameter must be between 0 and 1")
    }
    if(analyteMiss < 0 || analyteMiss > 1) {
      stop("analyteMiss parameter must be between 0 and 1")
    }
    if(cov.cutoff < 0 || cov.cutoff > 1) {
      stop("cov.cutoff parameter must be between 0 and 1")
    }

    # Check that at least one parameter is not null
    if ((analyteType1perc+analyteType2perc+analyteMiss+cov.cutoff) ==0) {
      if(suppressWarnings == FALSE){
        warning("No filtering parameters were set so the data remains unfiltered")
      }
      filtdata <- inputData
    } else {
      type1Data <- inputData@analyteType1
      type2Data <- inputData@analyteType2
      if(analyteType1perc > 0) {
        # Filter
        mymean <- as.numeric(apply(type1Data,1, function(x) mean(x,na.rm=TRUE)))
        keepers <- which(mymean > stats::quantile(mymean,analyteType1perc, na.rm = TRUE))
        oldData <- type1Data
        type1Data <- as.matrix(type1Data[keepers,])
        # Transpose if only one column is remaining
        if(ncol(type1Data)==1 && nrow(type1Data) == ncol(oldData)){
          type1Data <- t(type1Data)
        }
        # Set rownames
        rownames(type1Data) <- rownames(oldData)[keepers]
        # Create empty matrix if no columns are remaining
        if(nrow(type1Data)==0){
          type1Data <- matrix(, nrow = 0, ncol = 0)
        }
      } else {
        if(suppressWarnings == FALSE){
          warning("No filtering by percentile is applied for analyte type 1")
        }
      }
      if(analyteType2perc > 0) {
        # Filter
        mymean <- as.numeric(apply(type2Data,1, function(x) mean(x,na.rm=TRUE)))
        keepers <- which(mymean > stats::quantile(mymean,analyteType2perc, na.rm = TRUE))
        oldData <- type2Data
        type2Data <- as.matrix(type2Data[keepers,])
        # Transpose if only one column is remaining
        if(ncol(type2Data)==1 && nrow(type2Data) == ncol(oldData)){
          type2Data <- t(type2Data)
        }
        # Set rownames
        rownames(type2Data) <- rownames(oldData)[keepers]
        # Create empty matrix if no columns are remaining
        if(nrow(type2Data)==0){
          type2Data <- matrix(, nrow = 0, ncol = 0)
        }
      } else {
        if(suppressWarnings == FALSE){
          warning("No filtering by percentile is applied for analyte type 2")
        }
      }
      if(analyteMiss > 0) {
        # Filter
        missnum <- as.numeric(apply(type2Data,1,function(x) length(which(is.na(x)))))
        mycut <- analyteMiss * ncol(type2Data)
        keepers <- which(missnum < mycut)
        oldData <- type2Data
        type2Data <- as.matrix(type2Data[keepers,])
        # Transpose if only one column is remaining
        if(ncol(type2Data)==1 && nrow(type2Data) == ncol(oldData)){
          type2Data <- t(type2Data)
        }
        # Set rownames
        rownames(type2Data) <- rownames(oldData)[keepers]
        # Create empty matrix if no columns are remaining
        if(nrow(type2Data)==0){
          type2Data <- matrix(, nrow = 0, ncol = 0)
        }
      } else {
        if(suppressWarnings == FALSE){
          warning("No filtering by missing values is applied for analyte type 2")
        }
      }
      
      #Checking for Log Scaling and Filtering the Data
      if(length(which(!is.na(type1Data))) > 0 && min(type1Data) >= 0
         && cov.cutoff > 0) {
        analyte1mean <- apply(type1Data, MARGIN = 1, mean)
        analyte1sd <- apply(type1Data, MARGIN = 1, mean)
        analyte1cov <- analyte1sd / analyte1mean
        cov.quant <- stats::quantile(analyte1cov, probs = cov.cutoff, na.rm = TRUE)
        type1Data <- type1Data[which(analyte1cov <= cov.quant),]
      }else if(suppressWarnings == FALSE && min(type1Data) < 0 && cov.cutoff > 0){
        warning(paste("Coefficient of variation filtering will not be applied",
                "to analyte type 1 because data is log-scaled"))
      }
      
      if(length(which(!is.na(type2Data))) > 0 && min(type2Data) >= 0
         && cov.cutoff > 0) {
        analyte2mean <- apply(type2Data, MARGIN = 1, mean)
        analyte2sd <- apply(type2Data, MARGIN = 1, mean)
        analyte2cov <- analyte2sd / analyte2mean
        cov.quant <- stats::quantile(analyte2cov, probs = cov.cutoff, na.rm = TRUE)
        type2Data <- type2Data[which(analyte2cov <= cov.quant),]
      }else if(suppressWarnings == FALSE && min(type2Data) < 0 && cov.cutoff > 0){
        warning(paste("Coefficient of variation filtering will not be applied",
                "to analyte type 2 because data is log-scaled"))
      }
      
      # Now reconstruct the multidataset object
      filtdata <- methods::new("IntLimData", analyteType1=type1Data,
                               analyteType2=type2Data,
                               analyteType1MetaData = inputData@analyteType1MetaData,
                               analyteType2MetaData = inputData@analyteType2MetaData,
                               sampleMetaData = inputData@sampleMetaData)
    }
  }else if(length(inputData@analyteType1)>0){
    if(analyteType1perc < 0 || analyteType1perc > 1) {
      stop("analyteType1perc parameter must be between 0 and 1")
    }
    if(analyteMiss < 0 || analyteMiss > 1){
      stop("analyteMiss parameter must be between 0 and 1")
    }
    if(cov.cutoff < 0 || cov.cutoff > 1) {
      stop("cov.cutoff parameter must be between 0 and 1")
    }
    
    # Check that at least one parameter is not null
    if (analyteType1perc+cov.cutoff+analyteMiss ==0) {
      if(suppressWarnings == FALSE){
        warning("No filtering parameters were set so the data remains unfiltered")
      }
      filtdata <- inputData
    }
    else {
      type1Data <- inputData@analyteType1
      if(analyteType1perc > 0) {
        type1Data <- inputData@analyteType1
        mymean <- as.numeric(apply(type1Data,1, function(x)
          mean(x,na.rm=TRUE)))
        keepers <- which(mymean > stats::quantile(mymean,analyteType1perc, na.rm = TRUE))
        oldData <- type1Data
        type1Data <- as.matrix(type1Data[keepers,])
        # Transpose if only one column is remaining
        if(ncol(type1Data)==1 && nrow(type1Data) == ncol(oldData)){
          type1Data <- t(type1Data)
        }
        # Set rownames
        rownames(type1Data) <- rownames(oldData)[keepers]
        # Create empty matrix if no columns are remaining
        if(nrow(type1Data)==0){
          type1Data <- matrix(, nrow = 0, ncol = 0)
        }
      } else {
        if(suppressWarnings == FALSE){
          warning("No filtering by percentile is applied for analyte type 1")
        }
      }
      if(analyteMiss > 0) {
        # Filter
        missnum <- as.numeric(apply(type1Data,1,function(x) length(which(is.na(x)))))
        mycut <- analyteMiss * ncol(type1Data)
        keepers <- which(missnum < mycut)
        oldData <- type1Data
        type1Data <- as.matrix(type1Data[keepers,])
        # Transpose if only one column is remaining
        if(ncol(type1Data)==1 && nrow(type1Data) == ncol(oldData)){
          type1Data <- t(type1Data)
        }
        # Set rownames
        rownames(type1Data) <- rownames(oldData)[keepers]
        # Create empty matrix if no columns are remaining
        if(nrow(type1Data)==0){
          type1Data <- matrix(, nrow = 0, ncol = 0)
        }
      } else {
        if(suppressWarnings == FALSE){
          warning("No filtering by missing values is applied for analyte type 1")
        }
      }

      #Checking for Log Scaling and Filtering the Data
      if(length(which(!is.na(type1Data))) > 0 && min(type1Data) >= 0
         && cov.cutoff > 0) {
        analyte1mean <- apply(type1Data, MARGIN = 1, mean)
        analyte1sd <- apply(type1Data, MARGIN = 1, mean)
        analyte1cov <- analyte1sd / analyte1mean
        cov.quant <- stats::quantile(analyte1cov, probs = cov.cutoff, na.rm = TRUE)
        type1Data <- type1Data[which(analyte1cov <= cov.quant),]
      }else if(suppressWarnings == FALSE && min(type1Data) < 0 && cov.cutoff > 0){
        warning(paste("Coefficient of variation filtering will not be applied",
                "to analyte type 1 because data is log-scaled"))
      }
      
      # Now reconstruct the object
      filtdata <- methods::new("IntLimData", analyteType1=type1Data,
                               analyteType2=matrix(, nrow = 0, ncol = 0),
                               analyteType1MetaData = inputData@analyteType1MetaData,
                               analyteType2MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               sampleMetaData = inputData@sampleMetaData)
      
    }
  } else if(length(inputData@analyteType2)>0){
    if(analyteType2perc < 0 || analyteType2perc > 1) {
      stop("analyteType2perc parameter must be between 0 and 1")
    }
    if(analyteMiss < 0 || analyteMiss > 1) {
      stop("analyteMiss parameter must be between 0 and 1")
    }
    if(cov.cutoff < 0 || cov.cutoff > 1) {
      stop("cov.cutoff parameter must be between 0 and 1")
    }
    
    # Check that at least one parameter is not null
    if ((analyteType2perc+analyteMiss+cov.cutoff) ==0) {
      if(suppressWarnings == FALSE){
        warning("No filtering parameters were set so the data remains unfiltered")
      }
      filtdata <- inputData
    }
    else {
      type2Data <- inputData@analyteType2
      if(analyteType2perc > 0) {
        # Filter
        mymean <- as.numeric(apply(type2Data,1, function(x)
          mean(x,na.rm=TRUE)))
        keepers <- which(mymean > stats::quantile(mymean,analyteType2perc, na.rm = TRUE))
        oldData <- type2Data
        # Transpose if only one column is remaining
        type2Data <- as.matrix(type2Data[keepers,])
        if(ncol(type2Data)==1 && nrow(type2Data) == ncol(oldData)){
          type2Data <- t(type2Data)
        }
        # Set rownames
        rownames(type2Data) <- rownames(oldData)[keepers]
        # Create empty matrix if no columns are remaining
        if(nrow(type2Data)==0){
          type2Data <- matrix(, nrow = 0, ncol = 0)
        }
      } else {
        if(suppressWarnings == FALSE){
          warning("No filtering by percentile is applied for analyte type 2")
        }
      }
      if(analyteMiss > 0) {
        # Filter
        missnum <- as.numeric(apply(type2Data,1,function(x) length(which(is.na(x)))))
        mycut <- analyteMiss * ncol(type2Data)
        keepers <- which(missnum < mycut)
        oldData <- type2Data
        type2Data <- as.matrix(type2Data[keepers,])
        # Transpose if only one column is remaining
        if(ncol(type2Data)==1 && nrow(type2Data) == ncol(oldData)){
          type2Data <- t(type2Data)
        }
        # Set rownames
        rownames(type2Data) <- rownames(oldData)[keepers]
        # Create empty matrix if no columns are remaining
        if(nrow(type2Data)==0){
          type2Data <- matrix(, nrow = 0, ncol = 0)
        }
      } else {
        if(suppressWarnings == FALSE){
          warning("No filtering by missing values is applied for analyte type 2")
        }
      }
      
      if(analyteType1perc < 0 || analyteType1perc > 1) {
        stop("analyteType1perc parameter must be between 0 and 1")
      }
      if(analyteType2perc < 0 || analyteType2perc > 1) {
        stop("analyteType2perc parameter must be between 0 and 1")
      }
      if(analyteMiss < 0 || analyteMiss > 1) {
        stop("analyteMiss parameter must be between 0 and 1")
      }
      if(cov.cutoff < 0 || cov.cutoff > 1) {
        stop("cov.cutoff parameter must be between 0 and 1")
      }
      
      #Checking for Log Scaling and Filtering the Data
      if(length(which(!is.na(type2Data))) > 0 && min(type2Data) >= 0
         && cov.cutoff > 0) {
        analyte2mean <- apply(type2Data, MARGIN = 1, mean)
        analyte2sd <- apply(type2Data, MARGIN = 1, mean)
        analyte2cov <- analyte2sd / analyte2mean
        cov.quant <- stats::quantile(analyte2cov, probs = cov.cutoff, na.rm = TRUE)
        type2Data <- type2Data[which(analyte2cov <= cov.quant),]
      }else if(suppressWarnings == FALSE && min(type2Data) < 0 && cov.cutoff > 0){
        warning(paste("Coefficient of variation filtering will not be applied",
                "to analyte type 2 because data is log-scaled"))
      }
      
      # Now reconstruct the object
      filtdata <- methods::new("IntLimData", analyteType1=matrix(, nrow = 0, ncol = 0),
                               analyteType2=type2Data,
                               analyteType1MetaData = as.data.frame(matrix(, nrow = 0, ncol = 0)),
                               analyteType2MetaData = inputData@analyteType2MetaData,
                               sampleMetaData = inputData@sampleMetaData)
    }
  }
  if(length(filtdata@analyteType2) == 0 && length(inputData@analyteType2) > 0){
    stop(paste("All analytes have been removed from your type 2 data! Change your filtering criteria."))
  }
  if(length(filtdata@analyteType1) == 0 && length(inputData@analyteType1) > 0){
    stop(paste("All analytes have been removed from your type 1 data! Change your filtering criteria."))
  }
  return(filtdata)
}