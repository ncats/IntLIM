#' Run linear models and retrieve relevant statistics
#'
#' @include internalfunctions.R
#'
#' @param inputData Named list (output of 
#' FilterData()) with gene expression, metabolite abundances, 
#' and associated meta-data
#' @param stype column name that represents sample type (by default, it will be used
#' in the interaction term). Only 2 categories are currently supported.
#' @param outcome 'metabolite' or 'gene' must be set as outcome/independent variable
#' (default is 'metabolite')
#' @param covar Additional variables from the phenotypic data that be integrated into linear model
#' @param continuous boolean to indicate whether the data is continuous or discrete
#' @param metabolite.pairs boolean to indicate whether to return metabolite-metabolite pairs (TRUE) or gene-metabolite pairs (FALSE)
#' @param save.covar.pvals boolean to indicate whether or not to save the p-values of all covariates,
#' which can be analyzed later but will also lengthen computation time. The default is FALSE.
#' @param independent.var.type 'metabolite' or 'gene' must be set as independent variable
#' (default is 'metabolite')
#' @param remove.duplicates boolean to indicate whether or not to remove the metabolite-metabolite
#' or gene-gene pair with the highest p-value across two duplicate models (e.g. m1~m2 and m2~m1)
#' 
#' @return IntLimModel object with model results
#' @export
RunIntLim <- function(inputData,stype=NULL,outcome="metabolite", covar=NULL, 
                      continuous = FALSE, metabolite.pairs=FALSE,
                      save.covar.pvals=FALSE, independent.var.type="gene", remove.duplicates = FALSE){

    if(!continuous & length(unique(stats::na.omit(inputData$p))) != 2) {
	    stop(paste("IntLim currently requires only two categories.  Make sure the 
	               column",stype,"only has two unique values"))
    }
    print("Running the analysis on")
    if(continuous == FALSE){
        print(table(inputData$p))
    }else{
        print(range(inputData$p))
    }

    ptm <- proc.time()

    myres <- NULL
    if(independent.var.type == "gene" && outcome == "metabolite"){
        if("gene" %in% names(inputData) && "metab" %in% names(inputData)) {
            myres <- RunLM(inputData,outcome=outcome,
                    type=inputData$p,covar=covar, continuous = continuous, 
                    save.covar.pvals=save.covar.pvals)
        }
        else{
            stop("Either the gene data or the metabolite data is missing. Cannot run.\n")
        }
    }else if(independent.var.type == "metabolite" && outcome == "gene"){
        if("gene" %in% names(inputData) && "metab" %in% names(inputData)) {
            myres <- RunLM(inputData,outcome=outcome,type=inputData$p,covar=covar, 
                           continuous = continuous, save.covar.pvals=save.covar.pvals)
        }
        else{
            stop("Either the gene data or the metabolite data is missing. Cannot run.\n")
        }
    }else if(independent.var.type == "metabolite" && outcome == "metabolite"){
        if("metab" %in% names(inputData)) {
            myres <- RunLMMetabolitePairs(inputData,type=inputData$p,covar=covar, 
                                          continuous = continuous, 
                                          save.covar.pvals=save.covar.pvals, 
                                          keep.highest.pval = remove.duplicates)
        }
        else{
            stop("The metabolite data is missing. Cannot run.\n")
        }
    }else if(independent.var.type == "gene" && outcome == "gene"){
        if("gene" %in% names(inputData)) {
            myres <- RunLMGenePairs(inputData,type=inputData$p,covar=covar, 
                                    continuous = continuous, 
                                    save.covar.pvals = save.covar.pvals, 
                                    keep.highest.pval = remove.duplicates)
        }
        else{
            stop("The gene data is missing. Cannot run.\n")
        }
    }else{
        print(paste("Error! independent.var.type and outcome.type must both be either",
        "gene or metabolite in RunIntLim."))
    }
    print(proc.time() - ptm)
    myres@stype=stype
    myres@outcome=outcome
    myres@independent.var.type=independent.var.type

    if(!is.null(covar)){
        covariate <- covar
        class.var <- c()
        for(i in 1:length(covar)){
            class.var[i] <- class(inputData$covar_matrix[,i])
        }

        myres@covar <- data.frame(covariate,class.var)
    }
    return(myres)
}

#' Run linear models for all data folds. This is a wrapper to RunIntLim.
#'
#' @include internalfunctions.R
#'
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param stype column name that represents sample type (by default, it will be used
#' in the interaction term). Only 2 categories are currently supported.
#' @param outcome 'metabolite' or 'gene' must be set as outcome/independent variable
#' (default is 'metabolite')
#' @param covar Additional variables from the phenotypic data that be integrated into linear model
#' @param continuous boolean to indicate whether the data is continuous or discrete
#' @param metabolite.pairs boolean to indicate whether to return metabolite-metabolite pairs (TRUE) or gene-metabolite pairs (FALSE)
#' @param save.covar.pvals boolean to indicate whether or not to save the p-values of all covariates,
#' which can be analyzed later but will also lengthen computation time. The default is FALSE.
#' @param independent.var.type 'metabolite' or 'gene' must be set as independent variable
#' (default is 'metabolite')
#' @param remove.duplicates boolean to indicate whether or not to remove the metabolite-metabolite
#' or gene-gene pair with the highest p-value across two duplicate models (e.g. m1~m2 and m2~m1)
#' 
#' @return List of IntLimModel objects with model results
#' @export
RunIntLimAllFolds <- function(inputData,stype=NULL,outcome="metabolite", covar=NULL, 
                      continuous = FALSE, metabolite.pairs=FALSE,
                      save.covar.pvals=FALSE, independent.var.type="gene", remove.duplicates = FALSE){
    myres <- list()
    for(i in 1:length(inputData)){
        # Gene-metabolite
        myres[i] <- IntLIM::RunIntLim(inputData = inputData[[i]]$training,
                                         stype=stype, 
                                         save.covar.pvals = save.covar.pvals, 
                                         outcome = outcome, 
                                         independent.var.type = independent.var.type,
                                         continuous = continuous)
    }
    return(myres)
}
