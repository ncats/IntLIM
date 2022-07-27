#' Run permutations of the IntLIM code to search for random cross-omic associations in dataset
#'
#' This function allows users to test different permutations of the metadata with their analytes to ensure
#' that any pairs being deemed significant by IntIM are not being suggested due to random chance, as is 
#' sometimes a problem in correlative associations.
#'
#' @param data IntLimData object (output of ReadData()) with analylte levels and
#'  associated sample meta-data
#' @param stype column name that represents sample type (by default, it will be used
#' in the interaction term). Only 2 categories are currently supported.
#' @param outcome '1' or '2' must be set as outcome/independent variable
#' (default is '1')
#' @param covar Additional variables from the phenotypic data that be integrated into linear model
#' @param save.covar.pvals boolean to indicate whether or not to save the p-values of all covariates,
#' which can be analyzed later but will also lengthen computation time. The default is FALSE.
#' @param independent.var.type '1' or '2' must be set as independent variable
#' (default is '1')
#' @param continuous boolean to indicate whether the data is continuous or discrete
#' @param pvalcutoff FDR adjusted p-value cutoff for number of significant multi-omic 
#' pairs (default = 0.20)
#' @param interactionCoeffPercentile Interaction coefficient cutoff for the IntLIM linear model (default = 0.10)
#' @param rsquaredCutoff Cutoff for the R-squared values for the models as a quality control (default = 0.50)
#' @param num.permutations Number of permutations to be ran (default = 1)
#' @param seed set.seed paramter allowing for custom random number generation seeds 
#' @return List object with 1st slot populated with dataframe containing the R^2 values of the models, 
#' and number of significant pairs before and after p-value adjustment. The 2nd slot in the list contains 
#' a string vector of the IDs of the significant pairs.
#' @export
PermuteIntLIM <- function(data, 
                          stype="", 
                          outcome=1, 
                          independent.var.type=1,
                          covar = c(), 
                          save.covar.pvals = FALSE, 
                          continuous = FALSE,
                          pvalcutoff = 0.05,
                          interactionCoeffPercentile = 0,
                          rsquaredCutoff = 0,
                          num.permutations = 1,
                          seed = 1) {
  
  avg.r2.vals <- c()
  sig.pairs <- c()
  sig.list <- vector(mode = "list", length = num.permutations)
  
  set.seed(seed)
  
  if(!methods::is(data, "IntLimData")) {
    stop("The data must be an IntLimData object")
  }
  
  if(num.permutations < 1) {
    stop("The number of permutations must be greater than or equal to 1.")
  }
  
  if(length(data@sampleMetaData) == 0) {
    stop("Sample metadata must be included in the IntlimData object")
  }
  
  for (i in 1:num.permutations) {
    
    analyte1 <- data.frame(data@analyteType1, check.names = FALSE)
    analyte2 <- data.frame(data@analyteType2, check.names = FALSE)
    sample.meta <- data.frame(data@sampleMetaData, check.names = FALSE)
    
    #Randomizing Analyte1 and Analyte2 Data
    New_Samples <- sample(rownames(sample.meta), size = nrow(sample.meta), replace = FALSE)
    stype_perm <- sample.meta[New_Samples,stype]
    sample.meta[,stype] <- stype_perm
    
    #Recompiling the InputData
    data@analyteType1 <- as.matrix(analyte1)
    data@analyteType2 <- as.matrix(analyte2)
    data@sampleMetaData <- as.data.frame(sample.meta)

    #Running IntLIM
    IntLIMResults <- IntLIM::RunIntLim(inputData = data, stype=stype,
                                       covar = covar, outcome = outcome,
                                       independent.var.type = independent.var.type,
                                       save.covar.pvals = save.covar.pvals,
                                       continuous = continuous)

    ProcessedResults <- IntLIM::ProcessResults(inputData = data, inputResults = IntLIMResults,
                                               pvalcutoff = pvalcutoff,
                                               interactionCoeffPercentile = interactionCoeffPercentile,
                                               rsquaredCutoff = rsquaredCutoff)

    #Compiling Results
    r.squared <- mean(IntLIMResults@model.rsquared)

    model.p.values <- ProcessedResults$FDRadjPval
    num.unadj.sig <- length(which(model.p.values <= pvalcutoff))

    avg.r2.vals <- c(avg.r2.vals, r.squared)
    sig.pairs <- c(sig.pairs, num.unadj.sig)

    pair.index <- which(model.p.values <= pvalcutoff)

    sig.pairs.ids <- c()
    for(j in 1:length(pair.index)) {
      current.id <- paste0(ProcessedResults$Analyte1[pair.index[j]], "__V__",
                           ProcessedResults$Analyte2[pair.index[j]])
      sig.pairs.ids <- c(sig.pairs.ids, current.id)
    }

    sig.list[[i]] <- as.vector(sig.pairs.ids)
  }
  
  permuted.df <- data.frame(avg.r2.vals, sig.pairs)
  colnames(permuted.df) <- c("Avg_R_Squared_Values", "Num_Significant_Pairs")
  
  permutedResults <- list(numSigPairs=permuted.df, listOfSigPairs=sig.list)
  return(permutedResults)
}