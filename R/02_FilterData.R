#' Filter input data by abundance values (gene and metabolite data) and number of missing values (metabolite data only).
#'
#' Filter data by abundance (with user-input percentile cutoff) of missing values (with user-input percent cutoff). Missing values are commonly found in metabolomics data so the parameter currently only applies to metabolomics data.
#'
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression, 
#' metabolite abundances, and associated meta-data
#' @param geneperc percentile cutoff (0-1) for filtering genes (e.g. remove genes with mean values 
#' < 'geneperc' percentile) (default: 0)
#' @param metabperc percentile cutoff (0-1) for filtering metabolites (default: no filtering of metabolites) (default:0)
#' @param metabmiss missing value percent cutoff (0-1) for filtering metabolites (metabolites with > 80\% missing values will be removed) (default:0)
#' @return filtData MultiDataSet object with input data after filtering
#'
#' @examples
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' inputData <- ReadData(csvfile,metabid='id',geneid='id')
#' inputDatafilt <- FilterData(inputData,geneperc=0.5)
#' @export
FilterData <- function(inputData,analyteType1perc=0,analyteType2perc=0,analyteType2miss=0,cov.cutoff=0.30) {
  
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
      
      analyte1mean <- apply(type1Data, MARGIN = 1, mean)
      analyte1sd <- apply(type1Data, MARGIN = 1, mean)
      analyte1cov <- analyte1sd / analyte1mean
      cov.quant <- quantile(analyte1cov, probs = cov.cutoff)
      type1Data <- type1Data[which(analyte1cov <= cov.quant),]
      
      analyte2mean <- apply(type2Data, MARGIN = 1, mean)
      analyte2sd <- apply(type2Data, MARGIN = 1, mean)
      analyte2cov <- analyte2sd / analyte2mean
      cov.quant <- quantile(analyte2cov, probs = cov.cutoff)
      type2Data <- type2Data[which(analyte1cov <= cov.quant),]
      
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
