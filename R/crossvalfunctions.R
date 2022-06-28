#' Runs the cross-validation end-to-end using the following steps:
#' 1. Create multiple cross-validation folds from the data.
#' 2. Filter each fold using the filtering criteria applied to the entire dataset.
#' 3. Run IntLIM for all folds.
#' 4. Process the results for all folds.
#'
#' @param inputData IntLimData object (output of ReadData()) with analylte levels and
#'  associated meta-data
#' @param folds number of folds to create
#' @param suppressWarnings whether to suppress warnings
#' @param analyteType1perc percentile cutoff (0-1) for filtering analyte type 1 (e.g. 
#' remove analytes with mean values < 'analyteType1perc' percentile) (default: 0)
#' @param analyteType2perc percentile cutoff (0-1) for filtering analyte type 2 
#' (default: no filtering of analytes) (default:0)
#' @param analyteMiss missing value percent cutoff (0-1) for filtering analytes 
#' (analytes with > 80\% missing values will be removed) (default:0)
#' @param cov.cutoff percentile cutoff (0-1) for the covariances of the anaytes (default: 0.30)
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient
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
#' @export
RunCrossValidation <- function(inputData,
                               folds, 
                               analyteType1perc=0, 
                               analyteType2perc=0,
                               analyteMiss=0,
                               cov.cutoff=0,
                               stype="",
                               outcome=c(1), 
                               covar=c(), 
                               continuous = FALSE, 
                               save.covar.pvals=FALSE, 
                               independent.var.type=c(1), 
                               remove.duplicates = FALSE,
                               pvalcutoff=0.05,
                               interactionCoeffPercentile=0,
                               rsquaredCutoff = 0.0,
                               treecuts = 0,
                               suppressWarnings=FALSE) {

  # Create the folds.
  inputDataFolds <- CreateCrossValFolds(inputData=inputData,
                                        folds=folds)
  
  # Filter the folds.
  inputDataFilt <- FilterDataFolds(inputDataFolds=inputDataFolds,
                                   analyteType1perc=analyteType1perc,
                                   analyteType2perc=analyteType2perc, 
                                   analyteMiss=analyteMiss,
                                   cov.cutoff=cov.cutoff,
                                   suppressWarnings = suppressWarnings)
  
  # Run IntLIM with the types specified.
  inputResults <- RunIntLimAllFolds(inputData=inputDataFilt,
                                    stype=stype,
                                    outcome=outcome, 
                                    covar=covar,
                                    continuous=continuous, 
                                    save.covar.pvals=save.covar.pvals, 
                                    independent.var.type=independent.var.type,
                                    remove.duplicates=remove.duplicates,
                                    suppressWarnings = suppressWarnings)
  
  # Process all results.
  sigResults <- ProcessResultsAllFolds(inputResults=inputResults, 
                                       inputData=inputDataFilt, 
                                       pvalcutoff=pvalcutoff,
                                       interactionCoeffPercentile=interactionCoeffPercentile,
                                       rsquaredCutoff=rsquaredCutoff)
    
  # Return everything.
  return(list(folds = inputDataFolds, filtered = inputDataFilt, results = inputResults,
             processed = sigResults))
    
}

#' Creates multiple cross-validation folds from the data. Format is a list of
#' IntLIMData training and testing pairs. The "training" slot contains all data
#' except that in the given fold, and the "testing" contains all data in the fold.
#'
#' @param inputData IntLimData object (output of ReadData()) with analylte levels and
#'  associated meta-data
#' @param folds number of folds to create
#' @return A set of IntLimData training and testing sets, of the following format:
#' list(list("train" = IntLimData, "test" = IntLimData), ... list("train" = IntLimData,
#' "test" = IntLimData))
#' 
#' @return List of IntLimModel objects with model results
CreateCrossValFolds <- function(inputData,folds) {
  
  # Check that input is a IntLimData
  if (!methods::is(inputData, "IntLimData")) {
    stop("input data is not a IntLimData class")
  }
  
  # Extract all samples.
  samps <- rownames(inputData@sampleMetaData)
  
  # Stop if the number of folds is greater than the number of samples.
  if(folds > length(samps)){
    stop("The number of folds is greater than the number of samples!")
  }
  
  # Stop if at least two folds are not specified.
  if(folds < 2){
    stop("At least 2 folds are required.")
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
    type1MetaData <- inputData@analyteType1MetaData
    type2MetaData <- inputData@analyteType2MetaData
    covar_train <- inputData@sampleMetaData[not_fold,]
    covar_test <- inputData@sampleMetaData[fold,]
    # If there are no covariates other than phenotype,
    # reformat the covariate matrices appropriately.
    if(ncol(inputData@sampleMetaData == 1)){
      # Training set
      covar_train <- as.data.frame(covar_train)
      colnames(covar_train) <- colnames(inputData@sampleMetaData)
      rownames(covar_train) <- not_fold
      # Testing set
      covar_test <- as.data.frame(covar_test)
      colnames(covar_test) <- colnames(inputData@sampleMetaData)
      rownames(covar_test) <- fold
    }
    type1_train <- inputData@analyteType1
    type1_test <- inputData@analyteType1
    type2_train <- inputData@analyteType2
    type2_test <- inputData@analyteType2
    
    # Include all but the current fold in the training data.
    if(length(inputData@analyteType1)>0){
      type1_train <- inputData@analyteType1[,not_fold]
      type1_test <- as.matrix(as.data.frame(inputData@analyteType1[,fold]))
      rownames(type1_test) <- rownames(inputData@analyteType1)
      colnames(type1_test) <- rownames(covar_test)
    }
    if(length(inputData@analyteType2)>0){
      type2_train <- inputData@analyteType2[,setdiff(samps, fold)]
      type2_test <- as.matrix(as.data.frame(inputData@analyteType2[,fold]))
      rownames(type2_test) <- rownames(inputData@analyteType2)
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
#' @param inputDataFolds List of IntLimData objects (output of ReadData()) with analylte levels and
#'  associated meta-data
#' @param analyteType1perc percentile cutoff (0-1) for filtering analyte type 1 (e.g. 
#' remove analytes with mean values < 'analyteType1perc' percentile) (default: 0)
#' @param analyteType2perc percentile cutoff (0-1) for filtering analyte type 2 
#' (default: no filtering of analytes) (default:0)
#' @param analyteMiss missing value percent cutoff (0-1) for filtering analytes 
#' (analytes with > 80\% missing values will be removed) (default:0)
#' @param suppressWarnings whether to suppress warnings
#' @param cov.cutoff percentile cutoff (0-1) for the covariances of the anaytes (default: 0.30)
#' @return filtData IntLimData object with input data after filtering
FilterDataFolds <- function(inputDataFolds,analyteType1perc=0,
                            analyteType2perc=0, analyteMiss=0,
                            cov.cutoff=0,
                            suppressWarnings=FALSE) {
  
  # Now, filter each fold so that the analytes in the set match the analytes.
  for(i in 1:length(inputDataFolds)){
    
    # Filter training and testing data.
    inputDataFolds[[i]]$training <- FilterData(inputData=inputDataFolds[[i]]$training,
                                               analyteType1perc=analyteType1perc,
                                               analyteType2perc=analyteType2perc, 
                                               analyteMiss=analyteMiss,
                                               cov.cutoff=cov.cutoff,
                                               suppressWarnings=suppressWarnings)
    inputDataFolds[[i]]$testing <- inputDataFolds[[i]]$testing
  }
  return(inputDataFolds)
}

#' Run linear models for all data folds. This is a wrapper to RunIntLim.
#'
#' @include internalfunctions.R
#'
#' @param inputData IntLimData object (output of ReadData()) with analyte levels
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
#' @param suppressWarnings whether to suppress warnings
#' @return List of IntLimModel objects with model results
RunIntLimAllFolds <- function(inputData,stype="",outcome=1, covar=c(), 
                              continuous = FALSE, save.covar.pvals=FALSE, independent.var.type=1, 
                              remove.duplicates = FALSE, suppressWarnings = FALSE){
  myres <- lapply(1:length(inputData), function(i){
    res<-IntLIM::RunIntLim(inputData=inputData[[i]]$training,
                           stype=stype, 
                           save.covar.pvals = save.covar.pvals, 
                           covar = covar,
                           outcome = outcome, 
                           independent.var.type = independent.var.type,
                           continuous = continuous,
                           suppressWarnings = suppressWarnings)
    return(res)
  })
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
                                   interactionCoeffPercentile=0.5,
                                   rsquaredCutoff = 0.0,
                                   treecuts = 0){
  sig <- lapply(1:length(inputResults), function(i){
      return(IntLIM::ProcessResults(inputResults = inputResults[[i]], 
                                              inputData = inputData[[i]]$training, 
                                              pvalcutoff = pvalcutoff, 
                                              interactionCoeffPercentile = interactionCoeffPercentile, 
                                              rsquaredCutoff = rsquaredCutoff))
    })
  return(sig)
}
