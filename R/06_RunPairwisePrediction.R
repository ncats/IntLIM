#' Given each significant pairwise model and the input data, predict the phenotype
#' for each sample. Recall that IntLIM models take the following form:
#' m ~ beta0 + beta1(g) + beta2(phenotype) + beta3(g:phenotype) + beta4...n(covariates)
#' Therefore, to predict phenotype given the betas learned by IntLIM, we use the
#' following model:
#' p ~ (m - (beta0 + beta1(g) + beta4...n(covariates)) / (beta2 + beta3(g))
#' @param inputResults The data frame of filtered results from 
#'  an IntLimResults object. Each object must include
#'  model and processing results (output of ProcessResults()). All results must
#'  include learned covariate weights (i.e. must be run with save.covar.pvals = TRUE)
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param stype The phenotype (outcome) to predict. This can be either a categorical
#' or numeric outcome.
#' @param covar The clinical covariates to include in the model. These should be the same
#' covariates that were included when running the IntLIM linear models.
#' @param independent.var.type 'metabolite' or 'gene' must be set as independent variable
#' (default is 'metabolite')
#' @param outcome 'metabolite' or 'gene' must be set as outcome/independent variable
#' (default is 'metabolite')
#' @export
RunPairwisePrediction <- function(inputResults, inputData, stype=NULL, covar=NULL,
                                  independent.var.type="gene", outcome="metabolite"){

  # Extract data needed for model.
  covariates <- inputData$covar_matrix
  independent.vars <- NULL
  dependent.vars <- NULL
  if(independent.var.type == "metabolite"){
    independent.vars <- as.data.frame(inputData$metab)
  }else{
    independent.vars <- as.data.frame(inputData$gene)
  }
  if(outcome == "metabolite"){
    dependent.vars <- as.data.frame(inputData$metab)
  }else{
    dependent.vars <- as.data.frame(inputData$gene)
  }
  
  # Extract coefficients.
  which_start <- which(colnames(inputResults) == "rsquared")[1]
  if(is.na(which_start)){
    which_start <- which(colnames(inputResults) == "FDRadjPval")[1]
  }
  which_start <- which_start + 1    
  coefficients <- inputResults[,c(1:2, which_start:(dim(inputResults)[2]))]
  coefficients$Analyte1 <- as.character(coefficients$Analyte1)
  coefficients$Analyte2 <- as.character(coefficients$Analyte2)
  which_interact <- which(grepl(":", colnames(coefficients)) == TRUE)
  
  # Construct matrix of coefficients.
  intercept <- matrix(coefficients[,"(Intercept)"], 
                      nrow=length(coefficients[,"(Intercept)"]), 
                      ncol=dim(independent.vars)[2])
  ind_var <- matrix(coefficients$a, 
                    nrow=length(coefficients$a), 
                    ncol=dim(independent.vars)[2])
  interact_var <- matrix(coefficients[,which_interact], 
                         nrow=length(coefficients[,which_interact]), 
                         ncol=dim(independent.vars)[2])
  phen_var <- matrix(coefficients$type, nrow=length(coefficients$type), 
                     ncol=dim(independent.vars)[2])
  
  # Compute the numerator of the prediction for each subject 
  # (sans covariate terms)
  ind_term <- independent.vars[coefficients$Analyte1,] * ind_var
  dep_term <- dependent.vars[coefficients$Analyte2,]
  pred_phenotype <- dep_term - (intercept + ind_term)

  # If there are covariates, include the covariate terms in the prediction
  # by subtracting them.
  covariates <- inputData$covar_matrix
  if(!is.null(covar)) {
    all_cov_terms<- lapply(covar, function(cov_name){
      coef_cov_name <- colnames(coefficients)[grepl(cov_name, colnames(coefficients))]
      cov_cov_name <- colnames(covariates)[grepl(cov_name, colnames(covariates))]
      this_covariate_term <- NULL
      
      # If the term is numeric, simply multiply.
      if(is.numeric(covariates[,cov_cov_name])){
        this_coefficient_mat <- matrix(coefficients[,coef_cov_name], 
                                       nrow=length(coefficients[,coef_cov_name]), 
                                       ncol=dim(covariates)[1])
        this_covariate_mat <- t(matrix(covariates[,cov_cov_name], 
                                       ncol=length(coefficients[,coef_cov_name]), 
                                       nrow=dim(covariates)[1]))
        this_covariate_term <- this_coefficient_mat * this_covariate_mat
      }
      else if(length(coef_cov_name) == 1){
        second_part_of_name <- strsplit(coef_cov_name, cov_name)[[1]][2]
        this_covariate_term <- multiplyCovariate(coefficients, covariates, coef_cov_name, 
                                                 cov_cov_name, second_part_of_name)
      }
      # Add together the multiples of the one-hot-encoded terms.
      else{
        these_covariate_terms <-lapply(1:length(coef_cov_name), function(i){
          second_part_of_name <- strsplit(coef_cov_name[i], cov_name)[[1]][2]
          cova <- covariates
          for(j in 1:dim(cova)[2]){
            cova[,j] <- as.character(cova[,j])
          }
          cova[multi.which(cova != second_part_of_name)] <- "Other"
          current_covariate_term <- multiplyCovariate(coefficients, cova, coef_cov_name[i], 
                                   cov_cov_name, second_part_of_name)
          return(current_covariate_term)
        })
        this_covariate_term <- Reduce('+', these_covariate_terms)
      }
      return(this_covariate_term)
    })
    
    # Include the covariates in the numerator prediction.
    final_covariate_val <- Reduce('+', all_cov_terms)
    pred_phenotype <- pred_phenotype - final_covariate_val
  }

  # Calculate the denominator and divide.
  div_term <- independent.vars[coefficients$Analyte1,] * interact_var
  div_term <- div_term + phen_var
  pred_phenotype <- pred_phenotype / div_term 

  # For discrete phenotypes only, round the value.
  if(!is.numeric(inputData$p)){
    pred_phenotype[multi.which(pred_phenotype >= 1)] <- 1
    pred_phenotype[multi.which(pred_phenotype <= 0)] <- 0
    pred_phenotype <- round(pred_phenotype,digits=0)
  }
  pred_phenotype = as.data.frame(pred_phenotype)
  
  # Add analyte information for each prediction.
  pred_phenotype$to <- coefficients$Analyte2
  pred_phenotype$from <- coefficients$Analyte1

  return(pred_phenotype)
}

#' Given each significant pairwise model and the input data, predict the phenotype
#' for each sample. Recall that IntLIM models take the following form:
#' m ~ beta0 + beta1(g) + beta2(phenotype) + beta3(g:phenotype) + beta4...n(covariates)
#' Therefore, to predict phenotype given the betas learned by IntLIM, we use the
#' following model:
#' p ~ (m - (beta0 + beta1(g) + beta4...n(covariates)) / (beta2 + beta3(g))
#' @param inputResults A list of IntLimResults objects. Each object must include
#'  model and processing results (output of ProcessResultsAllFolds()). All results must
#'  include learned covariate weights (i.e. must be run with save.covar.pvals = TRUE)
#' @param inputData MultiDataSet object (output of CreateCrossValFolds()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param stype The phenotype (outcome) to predict. This can be either a categorical
#' or numeric outcome.
#' @param covar The clinical covariates to include in the model. These should be the same
#' covariates that were included when running the IntLIM linear models.
#' @param testing Boolean indicating whether this is testing data. FALSE by default.
#' @param independentVarType The independent variable type ("gene" or "metabolite")
#' @param outcome The outcome type ("gene" or "metabolite")
#' @export
RunPairwisePredictionAllFolds <- function(inputResults, inputData, stype=NULL, covar=NULL, 
                                  testing = FALSE, independentVarType, outcome){
  # Get all predictions.
  all_preds <- lapply(1:length(inputResults[[1]]), function(i){
    predval <- lapply(1:length(inputResults), function(j){
      input <- inputData[[i]]$training
      if(testing == TRUE){
        input <- inputData[[i]]$testing
      }
      return(RunPairwisePrediction(inputResults = inputResults[[j]][[i]], 
                                 inputData = input,
                                 stype=stype, covar=covar,
                                 independent.var.type = independentVarType[[j]][[i]],
                                 outcome = outcome[[j]][[i]]))
    })
    return(do.call(rbind, predval))
  })
  
  # Assign names.
  names(all_preds) <- unlist(lapply(1:length(inputResults), function(i){
    return(paste("Fold", i, sep = "_"))
  }))
  
  # Return
  return(all_preds)
}

#' Multiply a covariate with its learned coefficients. This is straightforward for
#' numeric covariates, but requires some conversion for categorical covariates.
#' @param coefficients The coefficients of the selected covariates.
#' @param covariates The covariate values.
#' @param coef_cov_name The covariate name in the coefficient data
#' @param cov_cov_name The covariate name in the covariate data
#' @param level1 The level of the vactor which equals "1" in the one-hot encoding.
multiplyCovariate <- function(coefficients, covariates, coef_cov_name, cov_cov_name, level1){
  this_coefficient_mat <- matrix(coefficients[,coef_cov_name], 
                                 nrow=length(coefficients[,coef_cov_name]), 
                                 ncol=dim(covariates)[1])
  this_covariate_mat <- t(matrix(covariates[,cov_cov_name], 
                                 ncol=length(coefficients[,coef_cov_name]), 
                                 nrow=dim(covariates)[1]))
  original_mat <- this_covariate_mat
  this_covariate_df <- as.data.frame(t(this_covariate_mat))
  for(i in 1:dim(this_covariate_df)[2]){
    this_covariate_df[,i] <- as.numeric(as.factor(this_covariate_df[,i]))
  }
  this_covariate_mat <- as.matrix(t(this_covariate_df))
  this_covariate_mat[multi.which(original_mat != level1)] <- 0
  this_covariate_mat[multi.which(original_mat == level1)] <- 1
  this_covariate_term <- this_coefficient_mat * this_covariate_mat
  return(this_covariate_term)
}

#' Given a graph and a phenotype prediction for each significant pair, generate a 
#' new graph with the phenotype predictions included as the edge weights.
#' @param predictions A matrix of predictions. Each signficant pair of analytes 
#' results in a prediction for each subject.
#' @param coRegulationGraph An igraph object. This graph is the co-regulation graph
#' generated using IntLIM analysis of analyte pairs.
#' @export
ProjectPredictionsOntoGraph <- function(predictions, coRegulationGraph){
  
  # Convert graph into data frame.
  edges <- igraph::as_data_frame(coRegulationGraph, what = "edges")
  
  # Define continuous color based on prediction.
  pal <- grDevices::colorRampPalette(c("limegreen", "purple"))

  # Modify n copies of graph, where n is the number of subjects.
  new_graphs <- lapply(1:(dim(predictions)[2]-2), function(i){

    # Modify graph weight.
    subject_graph <- data.frame(edges)
    subject_graph$weight <- predictions[,i]
    
    # Modify color.
    bin_count <- 100
    intervals <- seq(range(predictions[,i])[1], range(predictions[,i])[2], 
                     by = (range(predictions[,i])[2] - range(predictions[,i])[1])
                     / (bin_count - 1))
    subject_color_scale <- findInterval(predictions[,i], intervals)
    subject_graph$color <- pal(bin_count + 1)[subject_color_scale]
    
    # Modify node properties.
    node_df <- igraph::as_data_frame(coRegulationGraph, what = "vertices")
    
    # Create graph.
    final_graph = igraph::graph_from_data_frame(subject_graph, vertices = node_df)
    return(final_graph)
  })
  names(new_graphs)<-colnames(predictions)[1:(length(colnames(predictions))-2)]
  return(new_graphs)
}

#' This is a wrapper for ProjectPredictionsOntoGraph for all folds of data.
#' @param predictions_list A list of lists of matrices of predictions. For each matrix,
#' each signficant pair of analytes results in a prediction for each subject. The
#' list may consist of: gene-gene model predictions, metabolite-metabolite model
#' predictions, metabolite-gene model predictions, and/or gene-metabolite model
#' predictions.
#' @param coRegulationGraphs A list of igraph objects. This graph is the co-regulation graph
#' generated using IntLIM analysis of analyte pairs.
#' @export
ProjectPredictionsOntoGraphAllFolds <- function(predictions_list, coRegulationGraphs){

  prediction_graphs <- lapply(1:length(coRegulationGraphs), function(i){
    ProjectPredictionsOntoGraph(predictions = predictions_list[[i]], 
                                coRegulationGraph = coRegulationGraphs[[i]])
  })
  return(prediction_graphs)
}