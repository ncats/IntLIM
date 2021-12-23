#' Runs the cross-validation end-to-end using the following steps:
#' 1. Create multiple cross-validation folds from the data.
#' 2. Filter each fold using the filtering criteria applied to the entire dataset.
#' 3. Run IntLIM for all folds.
#' 4. Process the results for all folds.
#'
#' @param inputFile input file in CSV format (see Description)
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
#' @param folds number of folds to create
#' @param suppressWarnings whether to suppress warnings
#' @param logAnalyteType1 whether or not to log values for Analyte Type 1(T/F).
#' Must coincide with parameter used during CreateCrossValFolds.
#' @param logAnalyteType2 whether or not to log values for Analyte Type 2(T/F).
#' Must coincide with parameter used during CreateCrossValFolds.
#' @param class.feat class ("factor" or "numeric") for each covariate. The following format
#' is required: list(covar1="numeric", covar2="factor"). 
#' Must coincide with parameter used during CreateCrossValFolds.
#' @param folds Number of folds to run.
#' @param analyteType1perc percentile cutoff (0-1) for filtering analyte type 1 (e.g. 
#' remove analytes with mean values < 'analyteType1perc' percentile) (default: 0)
#' @param analyteType2perc percentile cutoff (0-1) for filtering analyte type 2 
#' (default: no filtering of analytes) (default:0)
#' @param analyteType2miss missing value percent cutoff (0-1) for filtering analyte type 2 
#' (analytes with > 80\% missing values will be removed) (default:0)
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param diffcorr cutoff of differences in correlations for filtering (default 0.5)
#' @param corrtype spearman or pearson or other parameters allowed by cor() function 
#' (default spearman)
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient 
#' (default bottom 10 percent (high negative coefficients) and top 10 percent 
#' (high positive coefficients))
#' @param treecuts user-selected number of clusters (of pairs) 
#' to cut the tree into
#' @param rsquaredCutoff cutoff for lowest r-squared value
#' @param stype column name that represents sample type (by default, it will be used
#' in the interaction term). Only 2 categories are currently supported.
#' @param outcome list of outcomes to run. '1' or '2' must be set as outcome/independent variable
#' (default is '1')
#' @param covar Additional variables from the phenotypic data that be integrated into linear model
#' @param continuous boolean to indicate whether the data is continuous or discrete
#' @param save.covar.pvals boolean to indicate whether or not to save the p-values of all covariates,
#' which can be analyzed later but will also lengthen computation time. The default is FALSE.
#' @param independent.var.type list of independent variable types to run. '1' or '2' 
#' must be set as independent variable
#' (default is '1')
#' @param remove.duplicates boolean to indicate whether or not to remove the 
#'  pair with the highest p-value across two duplicate models (e.g. m1~m2 and m2~m1)
#' @return List of IntResults object with model results (now includes correlations)
RunCrossValidation <- function(inputFile,
                               analyteType1id="id",
                               analyteType2id="id",
                               logAnalyteType1=FALSE,
                               logAnalyteType2=FALSE, 
                               class.feat = "factor",
                               folds, 
                               analyteType1perc=0, 
                               analyteType2perc=0,
                               analyteType2miss=0,
                               stype="",
                               outcome=c(1), 
                               covar=c(), 
                               continuous = FALSE, 
                               save.covar.pvals=FALSE, 
                               independent.var.type=c(1), 
                               remove.duplicates = FALSE,
                               pvalcutoff=0.05,
                               diffcorr=0.5,
                               corrtype="spearman",
                               interactionCoeffPercentile=0.5,
                               rsquaredCutoff = 0.0,
                               treecuts = 0,
                               suppressWarnings=FALSE) {
  # Create the folds.
  inputDataFolds <- CreateCrossValFolds(inputFile,analyteType1id,analyteType2id, 
                                  logAnalyteType1,logAnalyteType2, class.feat,
                                  folds, suppressWarnings)
  
  # Filter the folds.
  inputDataFilt <- FilterDataFolds(inputFile,inputDataFolds,analyteType1perc,
                                analyteType2perc, analyteType2miss, analyteType1id,
                                analyteType2id, logAnalyteType1,
                                logAnalyteType2, class.feat)
  
  # Run IntLIM with the types specified.
  inputResults <- lapply(1:length(outcome), function(i){
    result <- RunIntLimAllFolds(inputDataFilt,stype,outcome[i], covar, 
                                  continuous, save.covar.pvals, independent.var.type[i], 
                                  remove.duplicates)
    return(result)
  })
  
  # Process all results.
  sigResults <- ProcessResultsAllFolds(inputResults, inputDataFilt, pvalcutoff=0.05,
                                         diffcorr, corrtype,
                                         interactionCoeffPercentile,
                                         rsquaredCutoff, treecuts)
    
  # Return everything.
  return(list(folds = inputDataFolds, filtered = inputDataFilt, results = inputResults,
             processed = sigResults))
    
}

#' Creates multiple cross-validation folds from the data. Format is a list of
#' IntLIMData training and testing pairs. The "training" slot contains all data
#' except that in the given fold, and the "testing" contains all data in the fold.
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
#' @param folds number of folds to create
#' @param suppressWarnings whether to suppress warnings
#' @return A set of IntLimData training and testing sets, of the following format:
#' list(list("train" = IntLimData, "test" = IntLimData), ... list("train" = IntLimData,
#' "test" = IntLimData))
#' 
#' @return List of IntLimModel objects with model results
CreateCrossValFolds <- function(inputFile,analyteType1id="id",analyteType2id="id", 
                                logAnalyteType1=FALSE,logAnalyteType2=FALSE, class.feat = "factor",
                                folds, suppressWarnings=FALSE) {
  
  # Create the components of the input.
  pieces <- ReadData(inputFile,analyteType1id,analyteType2id, logAnalyteType1,
                     logAnalyteType2, class.feat, suppressWarnings)
  
  # Extract all samples.
  samps <- rownames(pieces@sampleMetaData)
  
  # Stop if the number of folds is greater than the number of samples.
  if(folds > length(samps)){
    stop("ERROR: The number of folds is greater than the number of samples!")
  }
  
  # Permute samples and divide into folds.
  sets_of <- floor(length(samps) / folds)
  perm_samps <- sample(samps, length(samps), replace = FALSE)
  group_assignment <- ceiling(seq_along(perm_samps)/sets_of)
  # If uneven, assign remaining.
  if(length(unique(group_assignment)) > folds){
    which_extra <- which(group_assignment > folds)
    new_grp <- 1
    for(i in which_extra){
      group_assignment[i] <- new_grp
      new_grp <- new_grp + 1
    }
  }
  fold_samps <- split(perm_samps, group_assignment)
  
  # For each fold, extract the samples from the data set.
  trainTestObjects <- lapply(fold_samps, function(fold){
    # Initialize.
    not_fold <- sort(setdiff(samps, fold))
    type1MetaData <- pieces@analyteType1MetaData
    type2MetaData <- pieces@analyteType2MetaData
    covar_train <- pieces@sampleMetaData[not_fold,]
    covar_test <- pieces@sampleMetaData[fold,]
    type1_train <- pieces@analyteType1
    type1_test <- pieces@analyteType1
    type2_train <- pieces@analyteType2
    type2_test <- pieces@analyteType2
    
    # Include all but the current fold in the training data.
    if(length(pieces@analyteType1)>0){
      type1_train <- pieces@analyteType1[,not_fold]
      type1_test <- as.matrix(as.data.frame(pieces@analyteType1[,fold]))
      rownames(type1_test) <- rownames(pieces@analyteType1)
      colnames(type1_test) <- rownames(covar_test)
    }
    if(length(pieces@analyteType2)>0){
      type2_train <- pieces@analyteType2[,setdiff(samps, fold)]
      type2_test <- as.matrix(as.data.frame(pieces@analyteType2[,fold]))
      rownames(type2_test) <- rownames(pieces@analyteType2)
      colnames(type2_test) <- rownames(covar_test)
    }
    training <- methods::new("IntLimData",analyteType1=type1_train,
                             analyteType2=type2_train,
                             analyteType1MetaData = type1MetaData,
                             analyteType2MetaData = type2MetaData,
                             sampleMetaData = covar_train)
    
    
    # Include the current fold in the testing data.
    testing <- methods::new("IntLimData",analyteType1=type1_test,
                            analyteType2=type2_test,
                            analyteType1MetaData = type1MetaData,
                            analyteType2MetaData = type2MetaData,
                            sampleMetaData = covar_test)
    
    return(list("training"=training, "testing"=testing))
  })
  return(trainTestObjects)
}

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
#' @param type2keepers Type 2 analytes to keep.
#' @param type1keepers Type 1 analytes to keep.
#' @return filtData IntLimData object with input data after filtering
FilterFromKeepers <- function(inputData,type1keepers, type2keepers) {
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

#' Run linear models for all data folds. This is a wrapper to RunIntLim.
#'
#' @include internalfunctions.R
#'
#' @param inputData MultiDataSet object (output of ReadData()) with analyte levels
#'  and associated meta-data
#' @param stype column name that represents sample type (by default, it will be used
#' in the interaction term). Only 2 categories are currently supported.
#' @param outcome '1' or '2' must be set as outcome/independent variable
#' (default is '1')
#' @param covar Additional variables from the phenotypic data that be integrated into linear model
#' @param continuous boolean to indicate whether the data is continuous or discrete
#' @param save.covar.pvals boolean to indicate whether or not to save the p-values of all covariates,
#' which can be analyzed later but will also lengthen computation time. The default is FALSE.
#' @param independent.var.type '1' or '2' must be set as independent variable
#' (default is '1')
#' @param remove.duplicates boolean to indicate whether or not to remove the 
#'  pair with the highest p-value across two duplicate models (e.g. m1~m2 and m2~m1)
#' 
#' @return List of IntLimModel objects with model results
RunIntLimAllFolds <- function(inputData,stype="",outcome=1, covar=c(), 
                              continuous = FALSE, save.covar.pvals=FALSE, independent.var.type=1, 
                              remove.duplicates = FALSE){
  myres <- list()
  for(i in 1:length(inputData)){
    myres[i] <- IntLIM::RunIntLim(inputData = inputData[[i]]$training,
                                  stype=stype, 
                                  save.covar.pvals = save.covar.pvals, 
                                  outcome = outcome, 
                                  independent.var.type = independent.var.type,
                                  continuous = continuous)
  }
  return(myres)
}

#' Retrieve significant pairs, based on adjusted p-values, interaction
#' coefficient percentile, and r-squared values. This is a wrapper for ProcessResults.
#'
#' @include internalfunctions.R
#'
#' @param inputResults List of IntLimResults object with model results (output of RunIntLimAllFolds())
#' @param inputData List of MultiDataSet objects (output of CreateCrossValFolds()) 
#' with analyte levels and associated meta-data
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param diffcorr cutoff of differences in correlations for filtering (default 0.5)
#' @param corrtype spearman or pearson or other parameters allowed by cor() function 
#' (default spearman)
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient 
#' (default bottom 10 percent (high negative coefficients) and top 10 percent 
#' (high positive coefficients))
#' @param treecuts user-selected number of clusters (of pairs) 
#' to cut the tree into
#' @param rsquaredCutoff cutoff for lowest r-squared value
#' @return List of IntResults object with model results (now includes correlations)
ProcessResultsAllFolds <- function(inputResults,
                                   inputData,
                                   pvalcutoff=0.05,
                                   diffcorr=0.5,
                                   corrtype="spearman",
                                   interactionCoeffPercentile=0.5,
                                   rsquaredCutoff = 0.0,
                                   treecuts = 0){
  sig <- lapply(1:length(inputResults), function(i){
    return(lapply(1:length(inputResults[[1]]), function(j){
      return(IntLIM::ProcessResults(inputResults = inputResults[[i]][[j]], 
                                              inputData = inputData[[j]]$training, 
                                              pvalcutoff = pvalcutoff, 
                                              interactionCoeffPercentile = interactionCoeffPercentile, 
                                              rsquaredCutoff = rsquaredCutoff, 
                                              diffcorr = diffcorr, 
                                              corrtype = corrtype))
    }))
  })
  return(sig)
}
