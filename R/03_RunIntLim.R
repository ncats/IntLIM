#' Run linear models and retrieve relevant statistics
#'
#' @include internalfunctions.R
#'
#' @param inputData Named list (output of 
#' FilterData()) with analyte abundances, 
#' and associated meta-data
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
#' pair with the highest p-value across two duplicate models (e.g. m1~m2 and m2~m1)
#' 
#' @return IntLimModel object with model results
#' @export
RunIntLim <- function(inputData,stype=NULL,outcome=1, covar=NULL, 
                      continuous = FALSE, 
                      save.covar.pvals=FALSE, independent.var.type=1, remove.duplicates = FALSE){

    if(!continuous & length(unique(stats::na.omit(inputData@sampleMetaData[,stype]))) != 2) {
	    stop(paste("IntLim currently requires only two categories.  Make sure the 
	               column",stype,"only has two unique values. Did you mean to set",
	               "continuous to TRUE?"))
    }
    print("Running the analysis on")
    if(continuous == FALSE){
        print(table(inputData@sampleMetaData[,stype]))
    }else{
        print(range(inputData@sampleMetaData[,stype]))
    }

    ptm <- proc.time()

    myres <- NULL
    removeDupWarning <- "remove.duplicates only applies if the independent variable
                    and outcome are of the same analyte type. Duplicates will
                    not be removed."
    if(independent.var.type == 1 && outcome == 2){
        if(remove.duplicates == TRUE){
            warning(removeDupWarning)
        }
        if(length(inputData@analyteType1)>0 && length(inputData@analyteType2)>0) {
            myres <- RunLM(inputData,outcome=outcome,
                           independentVariable = independent.var.type,
                    type=inputData@sampleMetaData[,stype],covar=covar, continuous = continuous, 
                    save.covar.pvals=save.covar.pvals)
        }
        else{
            stop("One type of analyte data is missing. Cannot run.\n")
        }
    }else if(independent.var.type == 2 && outcome == 1){
        if(remove.duplicates == TRUE){
            warning(removeDupWarning)
        }
        if(length(inputData@analyteType1)>0 && length(inputData@analyteType2)>0) {
            myres <- RunLM(inputData,outcome=outcome,
                           independentVariable = independent.var.type,
                           type=inputData@sampleMetaData[,stype],covar=covar, 
                           continuous = continuous, 
                           save.covar.pvals=save.covar.pvals)
        }
        else{
            stop("One type of analyte data is missing. Cannot run.\n")
        }
    }else if(independent.var.type == 1 && outcome == 1){
        if(length(inputData@analyteType1)>0) {
            myres <- RunLM(inputData,outcome=outcome,
                           independentVariable = independent.var.type,
                           type=inputData@sampleMetaData[,stype],covar=covar, 
                           continuous = continuous, 
                           save.covar.pvals=save.covar.pvals, 
                           keep.highest.pval = remove.duplicates)
        }
        else{
            stop("Analyte Type 1 is missing. Cannot run.\n")
        }
    }else if(independent.var.type == 2 && outcome == 2){
        if(length(inputData@analyteType2)>0) {
            myres <- RunLM(inputData,outcome=outcome,
                                    independentVariable = independent.var.type,
                                    type=inputData@sampleMetaData[,stype],covar=covar, 
                                    continuous = continuous, 
                                    save.covar.pvals=save.covar.pvals, 
                                    keep.highest.pval = remove.duplicates)
        }
        else{
            stop("Analyte Type 2 is missing. Cannot run.\n")
        }
    }else{
        print(paste("Error! independent.var.type and outcome.type must both be either",
        "1 or 2 in RunIntLim."))
    }
    print(proc.time() - ptm)
    myres@stype=stype
    myres@outcome=outcome
    myres@independent.var.type=independent.var.type

    if(!is.null(covar)){
        covariate <- covar
        class.var <- c()
        for(i in 1:length(covar)){
            class.var[i] <- class(covar)
        }

        myres@covar <- data.frame(covariate,class.var)
    }
    return(myres)
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
#' @export
RunIntLimAllFolds <- function(inputData,stype=NULL,outcome=1, covar=NULL, 
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
