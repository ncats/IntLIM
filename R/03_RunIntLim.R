#' Run linear models and retrieve relevant statistics
#'
#' @include internalfunctions.R
#' @include AllClasses.R
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
#' @param suppressWarnings whether or not to print warnings. If TRUE, do not print.
#' 
#' @return IntLimResults object with model results
#' @export
RunIntLim <- function(inputData,stype="",outcome=1, covar=c(), 
                      continuous = FALSE, 
                      save.covar.pvals=FALSE, independent.var.type=1, 
                      remove.duplicates = FALSE, suppressWarnings = FALSE){
    
    # If the wrong data type, stop.
    if(!methods::is(inputData, "IntLimData")){
        stop("Input must be an IntLimData object")
    }

    if(!continuous & length(unique(stats::na.omit(inputData@sampleMetaData[,stype]))) != 2) {
	    stop(paste("IntLim currently requires only two categories.  Make sure the column",
	               stype,"only has two unique values. Did you mean to set",
	               "continuous to TRUE?"))
    }
    if(continuous == FALSE){
        message(do.call(paste, list("Running the analysis on discrete data with group sizes", table(inputData@sampleMetaData[,stype])[1],
                                    table(inputData@sampleMetaData[,stype])[2])))
    }else{
        message(do.call(paste, list("Running the analysis on continuous data with range", range(inputData@sampleMetaData[,stype])[1],
                                    range(inputData@sampleMetaData[,stype])[2])))
    }

    ptm <- proc.time()

    myres <- NULL
    removeDupWarning <- paste("remove.duplicates only applies if the independent variable",
                    "and outcome are of the same analyte type. Duplicates will",
                    "not be removed.")
    if(independent.var.type == 1 && outcome == 2){
        if(remove.duplicates == TRUE && suppressWarnings == FALSE){
            warning(removeDupWarning)
        }
        if(length(inputData@analyteType1)>0 && length(inputData@analyteType2)>0) {
            myres <- RunLM(inputData,outcome=outcome,
                           independentVariable = independent.var.type,
                    type=inputData@sampleMetaData[,stype],covar=covar, continuous = continuous, 
                    save.covar.pvals=save.covar.pvals,
                    suppressWarnings = suppressWarnings)
        }
        else{
            stop("One type of analyte data is missing. Cannot run.\n")
        }
    }else if(independent.var.type == 2 && outcome == 1){
        if(remove.duplicates == TRUE && suppressWarnings == FALSE){
            warning(removeDupWarning)
        }
        if(length(inputData@analyteType1)>0 && length(inputData@analyteType2)>0) {
            myres <- RunLM(inputData,outcome=outcome,
                           independentVariable = independent.var.type,
                           type=inputData@sampleMetaData[,stype],covar=covar, 
                           continuous = continuous, 
                           save.covar.pvals=save.covar.pvals,
                           suppressWarnings = suppressWarnings)
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
                           keep.highest.pval = remove.duplicates,
                           suppressWarnings = suppressWarnings)
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
                                    keep.highest.pval = remove.duplicates,
                           suppressWarnings = suppressWarnings)
        }
        else{
            stop("Analyte Type 2 is missing. Cannot run.\n")
        }
    }else{
        stop(paste("Error! independent.var.type and outcome.type must both be either",
        "1 or 2 in RunIntLim."))
    }
    myres@stype=stype
    myres@outcome=outcome
    myres@independent.var.type=independent.var.type
    myres@continuous=0
    if(continuous == TRUE){myres@continuous=1}

    return(myres)
}