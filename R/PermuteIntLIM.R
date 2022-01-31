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
#' @param analyte.metadata boolean value to indicate whether the user has included analyte 
#' metadata in the analysis
#' @param filter.cutoff percentile cutoff based on variance for both analytes to remove lowest variances in
#' dataset (default = 0.20)
#' @param pval.cutoff unadjusted p-value cutoff for number of significant multi-omic pairs (default = 0.05)
#' @param adj.pval.cutoff FDR adjusted p-value cutoff for number of significant multi-omic 
#' pairs (default = 0.20)
#' @param num.permutations Number of permutations to be ran (default = 1)
#' @param seed set.seed paramter allowing for custom random number generation seeds 
#' @return List object with 1st slot populated with dataframe containing the R^2 values of the models, 
#' and number of significant pairs before and after p-value adjustment. The 2nd slot in the list contains 
#' a string vector of the IDs of the significant pairs
#' 
#' @examples 
#' \dontrun{
#' permutedResults <- PermuteIntLIM(data = mydata, stype = "PBO_vs_Leukemia", outcome = 1, 
#'                                  independent.var.type = 2, analyte.metadata = TRUE,
#'                                  num.permutations = 100)
#' }
#' @export
PermuteIntLIM <- function(data, 
                          stype, 
                          outcome, 
                          independent.var.type,
                          analyte.metadata,
                          covar = NULL, 
                          save.covar.pvals = TRUE, 
                          continuous = TRUE,
                          filter.cutoff = 0.20,
                          pval.cutoff = 0.05,
                          adj.pval.cutoff = 0.20,
                          num.permutations = 1,
                          seed = 1) {
  
  avg.r2.vals <- c()
  sig.pairs.unadj <- c()
  sig.pairs.adj <- c()
  sig.list <- vector(mode = "list", length = num.permutations)
  
  set.seed(seed)
  
  if(class(data) != "IntLimData") {
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
    
    if(analyte.metadata == TRUE) {
      if(length(data@analyteType1MetaData) == 0) {
        stop("Analyte 1 metadata needs to be present if analyte.metadata is set to TRUE")
      }
      if(length(data@analyteType2MetaData) == 0) {
        stop("Analyte 2 metadata needs to be present if analyte.metadata is set to TRUE")
      }
      analyte1.meta <- data.frame(data@analyteType1MetaData, check.names = FALSE)
      analyte2.meta <- data.frame(data@analyteType2MetaData, check.names = FALSE)
    }
    
    #Randomizing Analyte1 and Analyte2 Data
    New_Analytes_1 <- sample(rownames(analyte1), size = nrow(analyte1), replace = FALSE)
    rownames(analyte1) <- New_Analytes_1
    
    New_Analytes_2 <- sample(rownames(analyte2), size = nrow(analyte2), replace = FALSE)
    rownames(analyte2) <- New_Analytes_2
    
    #Filtering for Low Variance
    analyte2.var <- apply(analyte2, MARGIN = 1, FUN = stats::var)
    analyte1.var <- apply(analyte1, MARGIN = 1, FUN = stats::var)
    
    analyte2.cutoff <- stats::quantile(analyte2.var, probs = filter.cutoff)
    analyte1.cutoff <- stats::quantile(analyte1.var, probs = filter.cutoff)
    
    analyte2.index <- which(analyte2.var <= analyte2.cutoff)
    analyte1.index <- which(analyte1.var <= analyte1.cutoff)
    
    analyte2.filt <- analyte2[-analyte2.index,]
    analyte1.filt <- analyte1[-analyte1.index,]
    
    if(analyte.metadata == TRUE) {
      analyte2.meta.filt <- analyte2.meta[-analyte2.index,]
      analyte1.meta.filt <- analyte1.meta[-analyte1.index,]
    }
    
    #Recompiling the InputData
    data@analyteType1 <- as.matrix(analyte1.filt)
    data@analyteType2 <- as.matrix(analyte2.filt)
    data@sampleMetaData <- as.data.frame(sample.meta)
    
    if(analyte.metadata == TRUE) {
      data@analyteType1MetaData <- as.data.frame(analyte1.meta.filt)
      data@analyteType2MetaData <- as.data.frame(analyte2.meta.filt)
    }
    
    #Running IntLIM
    IntLIMResults <- IntLIM::RunIntLim(inputData = data, stype=stype, 
                                       covar = covar, outcome = outcome, 
                                       independent.var.type = independent.var.type,
                                       save.covar.pvals = save.covar.pvals, 
                                       continuous = continuous)
    
    #Compiling Results
    r.squared <- mean(IntLIMResults@model.rsquared)
    
    model.p.values <- IntLIMResults@interaction.pvalues
    num.unadj.sig <- length(which(model.p.values <= pval.cutoff)) 
    
    model.adj.values <- IntLIMResults@interaction.adj.pvalues
    num.adj.sig <- length(which(model.adj.values <= adj.pval.cutoff))
    
    avg.r2.vals <- c(avg.r2.vals, r.squared)
    sig.pairs.unadj <- c(sig.pairs.unadj, num.unadj.sig)
    sig.pairs.adj <- c(sig.pairs.adj, num.adj.sig)
    
    sig.pairs.ids <- c()
    for(k in 1:ncol(model.p.values)) {
      for(j in 1:nrow(model.p.values)) {
        if(model.p.values[j,k] <= pval.cutoff) {
          current.id <- paste0(colnames(model.p.values)[k], "__V__", rownames(model.p.values)[j])
          sig.pairs.ids <- c(sig.pairs.ids, current.id)
        }
      }
    }
    
    sig.list[[i]] <- as.vector(sig.pairs.ids)
  }
  
  permuted.df <- data.frame(avg.r2.vals, sig.pairs.unadj, sig.pairs.adj)
  colnames(permuted.df) <- c("Avg_R_Squared_Values", "Significant_Pairs_Unadjusted", 
                             "Significant_Pairs_Adjusted")
  
  permutedResults <- list(permuted.df, sig.list)
  return(permutedResults)
}