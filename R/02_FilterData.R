#' Filter input data by abundance values (gene and metabolite data) and number of missing values (metabolite data only).
#'
#' Filter data by abundance (with user-input percentile cutoff) of missing values (with user-input percent cutoff). Missing values are commonly found in metabolomics data so the parameter currently only applies to metabolomics data.
#'
#' @param csvfile The path to the CSV file that contains the full dataset.
#' @param stype phenotype or outcome of interest.
#' @param inputDataFolds list of MultiDataSet objects (output of CreateCrossValFolds()) with gene expression, 
#' metabolite abundances, and associated meta-data
#' @param geneperc percentile cutoff (0-1) for filtering genes (e.g. remove genes with mean values 
#' < 'geneperc' percentile) (default: 0)
#' @param metabperc percentile cutoff (0-1) for filtering metabolites (default: no filtering of metabolites) (default:0)
#' @param metabmiss missing value percent cutoff (0-1) for filtering metabolites (metabolites with > 80\% missing values will be removed) (default:0)
#' @param class.covar class ("numeric" or "character") for each covariate. The following format
#' is required: list(covar1="numeric", covar2="character")
#' @param stype.type Either "factor" or "numeric". The outcome type.
#' @return filtData list of MultiDataSet objects with input data after filtering
#' @export
FilterDataFolds <- function(csvfile, stype,inputDataFolds,geneperc=0,metabperc=0, 
                            metabmiss=0,class.covar = NULL,
                            stype.type) {
	# First, read the original data set.
	inputData <- ReadData(csvfile,metabid='id',geneid='id')
	filtdata <- FilterData(inputData,stype=stype, geneperc=geneperc, 
	                       metabperc=metabperc, metabmiss = metabmiss,
	                       class.covar = class.covar,
	                       stype.type)
	
	# Filter the original data set.
	genekeepers = NULL
	metabkeepers = NULL
	if("gene" %in% names(filtdata)){
	  genes <- filtdata$gene
	  genekeepers <- rownames(genes)
	}
	if("metab" %in% names(filtdata)){
	  metabolites <- filtdata$metab
	  metabkeepers <- rownames(metabolites)
	}

	# Now, filter each fold so that the metabolites and genes in the set match the metabolites
	# and genes.
	for(i in 1:length(inputDataFolds)){

		# Filter training and testing data.
		inputDataFolds[[i]]$training <- FilterFromKeepers(inputData = inputDataFolds[[i]]$training,
		                                                  stype = stype,
		                                                  metabkeepers = metabkeepers,
		                                                  genekeepers = genekeepers,
		                                                  stype.type = stype.type, 
		                                                  class.covar = class.covar)
		inputDataFolds[[i]]$testing <- FilterFromKeepers(inputData = inputDataFolds[[i]]$testing,
		                                                 stype = stype,
		                                                  metabkeepers = metabkeepers,
		                                                  genekeepers = genekeepers,
		                                                 stype.type = stype.type,
		                                                 class.covar = class.covar)
	}
	return(inputDataFolds)
}

#' Filter a data set given the metabolites and genes to keep.
#'
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param stype phenotype or outcome of interest.
#' @param metabkeepers Metabolites to keep.
#' @param genekeepers Genes to keep.
#' @param stype.type Either "factor" or "numeric". The outcome type.
#' @param class.covar class ("numeric" or "character") for each covariate. The following format
#' is required: list(covar1="numeric", covar2="character")
#' @return filtData MultiDataSet object with input data after filtering
FilterFromKeepers <- function(inputData,stype,genekeepers, metabkeepers, stype.type,
                              class.covar = NULL) {
  # Extract metabolite and gene data to keep.
  mytypes <- names(Biobase::assayData(inputData))

  # Create new data set.
  filteredData <- MultiDataSet::createMultiDataSet()

  # Filter genes.
  if(any(mytypes=="expression")){
    # Filter the assay data.
    mygenes <- Biobase::assayDataElement(inputData[["expression"]], 'exprs')
    samps <- colnames(mygenes)
    mygenes <- as.matrix(mygenes[genekeepers,])
    colnames(mygenes) <- samps

    # Create a new expression set from the filtered assay data.
    gene.set <- Biobase::ExpressionSet(assayData=mygenes)
    fgenes <- Biobase::fData(inputData[["expression"]])[genekeepers,]
    Biobase::fData(gene.set) <- fgenes
    Biobase::pData(gene.set) <- Biobase::pData(inputData[["expression"]])

    # Add gene expression to the filtered data set.
    filteredData <- MultiDataSet::add_genexp(filteredData, gene.set)
  }

  # Filter metabolites.
  if(any(mytypes=="metabolite")){
    # Filter the assay data.
    mymetab <- Biobase::assayDataElement(inputData[["metabolite"]], 'metabData')
    samps <- colnames(mymetab)
    mymetab <- as.matrix(mymetab[metabkeepers,])
    colnames(mymetab) <- samps
    rownames(mymetab) <- metabkeepers

    # Create a new
    fmetab <- Biobase::AnnotatedDataFrame(data = Biobase::fData(inputData[["metabolite"]]))[metabkeepers,]
    metab.set <- methods::new("MetaboliteSet",metabData = mymetab,
                              phenoData =  Biobase::AnnotatedDataFrame(data = Biobase::pData(inputData[["metabolite"]])),
                              featureData =  fmetab)
    filteredData <- add_metabolite(filteredData, metab.set)
  }

  # Get common analytes.
  mytypes <- names(Biobase::assayData(filteredData))
  incommon <- filteredData
  if(any(mytypes == "expression") && any(mytypes == "metabolite")){
    incommon <- getCommon(inputData=filteredData,stype=stype)
  }else if(any(mytypes == "metabolite")){
    incommon <- formatSingleOmicInput(filteredData, stype,type = "metabolite")
  }else if(any(mytypes == "expression")){
    incommon <- formatSingleOmicInput(filteredData, stype,type = "expression")
  }
  
  # Coerce stype.
  if(stype.type == "factor"){
    incommon$p = as.factor(incommon$p)
  }else if(stype.type == "numeric"){
    incommon$p = as.numeric(as.character(incommon$p))
  }
  
  # Coerce covariate classes.
  for(covar in names(class.covar)){
    if(class.covar[covar] == "numeric"){
      incommon$covar_matrix[,covar] <- as.numeric(as.character(incommon$covar_matrix[,covar]))
    }else if(class.covar[covar] == "factor"){
      incommon$covar_matrix[,covar] <- as.factor(as.character(incommon$covar_matrix[,covar]))
    }else{
      stop(paste(class.covar[covar], "is not a valid class for covariate", covar))
    }
  }

  return(incommon)
}

#' Filter input data by abundance values (gene and metabolite data) and number of missing values (metabolite data only).
#'
#' Filter data by abundance (with user-input percentile cutoff) of missing values (with user-input percent cutoff). Missing values are commonly found in metabolomics data so the parameter currently only applies to metabolomics data.
#'
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression, 
#' metabolite abundances, and associated meta-data
#' @param stype phenotype or outcome of interest.
#' @param geneperc percentile cutoff (0-1) for filtering genes (e.g. remove genes with mean values 
#' < 'geneperc' percentile) (default: 0)
#' @param metabperc percentile cutoff (0-1) for filtering metabolites (default: no filtering of metabolites) (default:0)
#' @param metabmiss missing value percent cutoff (0-1) for filtering metabolites (metabolites with > 80\% missing values will be removed) (default:0)
#' @param class.covar class ("numeric" or "character") for each covariate. The following format
#' is required: list(covar1="numeric", covar2="character")
#' @param stype.type Either "factor" or "numeric". The outcome type.
#' @return filtData MultiDataSet object with input data after filtering
#'
#' @examples
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' inputData <- ReadData(csvfile,metabid='id',geneid='id')
#' inputDatafilt <- FilterData(inputData = inputData,stype = "PBO_vs_Leukemia", 
#' geneperc=0.5, stype.type = "factor")
#' @export
FilterData <- function(inputData,stype,geneperc=0,metabperc=0, metabmiss=0,
# We assume that the data is appropriately normalized...in the ReadData documentation
                       class.covar = NULL,
                       stype.type) {

  # Check that input is a MultiDataSet
  if (class(inputData) != "MultiDataSet") {
	   stop("input data is not a MultiDataSet class")
  }
  mytypes <- names(Biobase::assayData(inputData))
  filtdata <- NULL
  if(!any(mytypes=="expression") && !any(mytypes=="metabolite")) {
	   stop("input data must contain assayData of type 'metabolite' or 'expression'.
	     Try reading in the data with the ReadData function")
  }	else if(any(mytypes == "expression") && any(mytypes == "metabolite")){
    if(!is.null(geneperc) && geneperc > 1) {stop("geneperc parameter must be between 0 and 1")}
    if(!is.null(metabperc) && metabperc > 1) {stop("metabperc parameter must be between 0 and 1")}
    if(!is.null(metabmiss) && metabmiss > 1) {stop("metabmiss parameter must be between 0 and 1")}
  
    # Check that at least one parameter is not null
    len <- length(c(geneperc,metabperc,metabmiss))
    if ((geneperc+metabperc+metabmiss) ==0) {
      warning("All filtering parameters are NULL so the data remains unfiltered")
      filtdata <- inputData
    } else {
      mygenes <- Biobase::assayDataElement(inputData[["expression"]], 'exprs')
      mymetab <- Biobase::assayDataElement(inputData[["metabolite"]], 'metabData')
      if(geneperc > 0) {
        if(geneperc>1) {geneperc=geneperc}
        mymean <- as.numeric(apply(mygenes,1, function(x) mean(x,na.rm=T)))
        keepers <- which(mymean > stats::quantile(mymean,geneperc))
        mygenes <- mygenes[keepers,]
        pgenes <- Biobase::pData(inputData[["expression"]])
        fgenes <- Biobase::fData(inputData[["expression"]])[keepers,]
      } else {
        print("No gene filtering is applied")
        mygenes <- mygenes
        fgenes <- Biobase::fData(inputData[["expression"]])
      }
      if(metabperc > 0) {
        if(metabperc>1) {metabperc=metabperc}
          mymean <- as.numeric(apply(mymetab,1, function(x) mean(x,na.rm=T)))
          keepers <- which(mymean > stats::quantile(mymean,metabperc))
          mymetab <- mymetab[keepers,]
          fmetab <- Biobase::AnnotatedDataFrame(data = Biobase::fData(inputData[["metabolite"]]))[keepers,]
      } else {
        print("No metabolite filtering by percentile is applied")
        mymetab <- mymetab
        fmetab <- Biobase::AnnotatedDataFrame(data = Biobase::fData(inputData[["metabolite"]]))
      }
      if(metabmiss > 0) {
        missnum <- as.numeric(apply(mymetab,1,function(x) length(which(x==min(x,na.rm=T)))))-1
        mycut <- metabmiss * ncol(mymetab)
        keepers <- which(missnum < mycut)
        mymetab <- mymetab[keepers,]
        fmetab <- fmetab[keepers,]
      } else {
        print("No metabolite filtering by missing values is applied")
        mymetab <- mymetab
        fmetab <- fmetab
      }
    
      # Now reconstruct the multidataset object
      gene.set <- Biobase::ExpressionSet(assayData=mygenes)
      Biobase::fData(gene.set) <- fgenes
      Biobase::pData(gene.set) <- Biobase::pData(inputData[["expression"]])
    
      metab.set <- methods::new("MetaboliteSet",metabData = mymetab,
                              phenoData =  Biobase::AnnotatedDataFrame(data = Biobase::pData(inputData[["metabolite"]])), 
                              featureData =  fmetab)
    
      multi <- MultiDataSet::createMultiDataSet()
      multi1 <- MultiDataSet::add_genexp(multi, gene.set)
      filtdata <- add_metabolite(multi1, metab.set)
    }
  }else if(any(mytypes == "expression")){
    if(!is.null(geneperc) && geneperc > 1) {stop("geneperc parameter must be between 0 and 1")}
    
    # Check that at least one parameter is not null
    len <- length(c(geneperc,metabperc,metabmiss))
    if (is.null(geneperc)) {
      warning("All filtering parameters are NULL so the data remains unfiltered")
      filtdata <- inputData
    }
    else {
      mygenes <- Biobase::assayDataElement(inputData[["expression"]], 'exprs')
      if(geneperc > 0) {
        if(geneperc>1) {geneperc=geneperc}
        mymean <- as.numeric(apply(mygenes,1, function(x)
          mean(x,na.rm=T)))
        keepers <- which(mymean > stats::quantile(mymean,geneperc))
        mygenes <- mygenes[keepers,]
        pgenes <- Biobase::pData(inputData[["expression"]])
        fgenes <- Biobase::fData(inputData[["expression"]])[keepers,]
      } else {
        print("No gene filtering is applied")
        mygenes <- mygenes
        fgenes <- Biobase::fData(inputData[["expression"]])
      }
      
      # Now reconstruct the multidataset object
      gene.set <- Biobase::ExpressionSet(assayData=mygenes)
      Biobase::fData(gene.set) <- fgenes
      Biobase::pData(gene.set) <- Biobase::pData(inputData[["expression"]])
      
      multi <- MultiDataSet::createMultiDataSet()
      filtdata <- MultiDataSet::add_genexp(multi, gene.set)
    }
  } else if(any(mytypes == "metabolite")){
    if(!is.null(metabperc) && metabperc > 1) {stop("metabperc parameter must be between 0 and 1")}
    if(!is.null(metabmiss) && metabmiss > 1) {stop("metabmiss parameter must be between 0 and 1")}
    
    # Check that at least one parameter is not null
    len <- length(c(metabperc,metabmiss))
    if ((metabperc+metabmiss) ==0) {
      warning("All filtering parameters are NULL so the data remains unfiltered")
      filtdata <- inputData
    }
    else {
      mymetab <- Biobase::assayDataElement(inputData[["metabolite"]], 'metabData')
      if(metabperc > 0) {
        if(metabperc>1) {metabperc=metabperc}
        mymean <- as.numeric(apply(mymetab,1, function(x)
          mean(x,na.rm=T)))
        keepers <- which(mymean > stats::quantile(mymean,metabperc))
        mymetab <- mymetab[keepers,]
        fmetab <- Biobase::AnnotatedDataFrame(data = Biobase::fData(inputData[["metabolite"]]))[keepers,]
      } else {
        print("No metabolite filtering by percentile is applied")
        mymetab <- mymetab
        fmetab <- Biobase::AnnotatedDataFrame(data = Biobase::fData(inputData[["metabolite"]]))
      }
      if(metabmiss > 0) {
        missnum <- as.numeric(apply(mymetab,1,function(x) length(which(x==min(x,na.rm=T)))))-1
        mycut <- metabmiss * ncol(mymetab)
        keepers <- which(missnum < mycut)
        mymetab <- mymetab[keepers,]
        fmetab <- fmetab[keepers,]
      } else {
        print("No metabolite filtering by missing values is applied")
        mymetab <- mymetab
        fmetab <- fmetab
      }
      
      # Now reconstruct the multidataset object
      metab.set <- methods::new("MetaboliteSet",metabData = mymetab,
                                phenoData =  Biobase::AnnotatedDataFrame(data = Biobase::pData(inputData[["metabolite"]])), 
                                featureData =  fmetab)
      
      multi <- MultiDataSet::createMultiDataSet()
      filtdata <- add_metabolite(multi, metab.set)
    }
  }
  
  # Get common analytes.
  mytypes <- names(Biobase::assayData(filtdata))
  incommon <- filtdata
  if(any(mytypes == "expression") && any(mytypes == "metabolite")){
    incommon <- getCommon(inputData=filtdata,stype=stype)
  }else if(any(mytypes == "metabolite")){
    incommon <- formatSingleOmicInput(filtdata, stype,type = "metabolite")
  }else if(any(mytypes == "expression")){
    incommon <- formatSingleOmicInput(filtdata, stype,type = "expression")
  }
  
  # Coerce stype.
  if(stype.type == "factor"){
    incommon$p = as.factor(incommon$p)
  }else if(stype.type == "numeric"){
    incommon$p = as.numeric(as.character(incommon$p))
  }
  
  # Coerce covariate classes.
  for(covar in names(class.covar)){
    if(class.covar[covar] == "numeric"){
      incommon$covar_matrix[,covar] <- as.numeric(as.character(incommon$covar_matrix[,covar]))
    }else if(class.covar[covar] == "factor"){
      incommon$covar_matrix[,covar] <- as.factor(as.character(incommon$covar_matrix[,covar]))
    }else{
      stop(paste(class.covar[covar], "is not a valid class for covariate", covar))
    }
  }
  
	return(incommon)
}

