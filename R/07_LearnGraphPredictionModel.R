#' 
#' Format the input for graph-based learning. This input consists of:
#' 1. The Laplacian of a line graph built from the co-regulation graphs, where 
#' each node corresponds to a pair of analytes.
#' 2. A prediction value for each node of the line graph, for each sample X.
#' 3. The true prediction values Y for each sample X.
#' @include internalfunctions.R
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param predictionGraphs A list of igraph objects, each of which includes
#' predictions for each edge.
#' @param coregulationGraph An igraph object containing the coregulation graph.
#' @param stype.class The class of the outcome ("numeric" or "categorical")
#' @param edgeTypeList List containing one or more of the following to include
#' in the line graph:
#' - "shared.outcome.analyte"
#' - "shared.independent.analyte"
#' - "analyte.chain"
#' @export
formatInput <- function(predictionGraphs, coregulationGraph,
                                inputData, stype.class, edgeTypeList){
  
  # Extract edge-wise predictions.
  predictions_by_node <- lapply(names(predictionGraphs), function(sampName){
    df_predictions <- igraph::as_data_frame(predictionGraphs[[sampName]])
    node_names <- paste(make.names(df_predictions$to), make.names(df_predictions$from),
                        sep = "__")
    df_predictions_new <- data.frame(Node = node_names, Weight = df_predictions$weight)
    return(df_predictions_new)
  })
  names(predictions_by_node) <- names(predictionGraphs)
  predicted_weights_only <- lapply(predictions_by_node, function(pred){
    return(pred$Weight)
  })
  predictions_flattened <- t(data.frame(predicted_weights_only))
  colnames(predictions_flattened) <- predictions_by_node[[1]]$Node
  
  # Convert co-regulation graph into a line graph. Return the adjacency matrix.
  # If edges were not connected by nodes in the original graph, they may be
  # removed from the line graph. Remove these from the predictions_by_node df.
  A <- CreateLineGraph(predictionsByEdge = predictions_by_node[[1]],
                       graphWithPredictions = predictionGraphs[[1]],
                       edgeTypeList = edgeTypeList)
  predictions_flattened <- predictions_flattened[,colnames(A)]

  # Add self-loops.
  A_tilde <- as.matrix(A)
  diag(A_tilde) <- 1
  
  # Extract the diagonal (in degree) and raise to the negative half power.
  diags1 <- colSums(A_tilde)
  diags_neg_half <- 1 / sqrt(diags1)
  D_tilde_neg_half1 <- t(matrix(rep(diags_neg_half,length(diags_neg_half)),
                              nrow = length(diags_neg_half)))
  A_hat1 <- D_tilde_neg_half1 * A_tilde
  rm(D_tilde_neg_half1)
  
  # Extract the diagonal (out degree) and raise to the negative half power.
  diags2 <- rowSums(A_tilde)
  rm(A_tilde)
  diags_neg_half <- 1 / sqrt(diags2)
  D_tilde_neg_half2 <- t(matrix(rep(diags_neg_half,length(diags_neg_half)),
                              nrow = length(diags_neg_half)))

  # Obtain the final matrix. Note that we modify the matrix multiplication
  # problem to obtain an elementwise multiplication problem
  # because it speeds up computation.
  A_hat <- A_hat1 * D_tilde_neg_half2
  
  # Obtain the predictions.
  Y <- inputData$p
  if(stype.class == "factor"){
    Y <- as.numeric(Y)-1
  }
  # if(length(dim(input_data)) > 0){
  #   if(length(dim(input_data)) == 0 || dim(input_data)[1] == 1){
  #     Y <- unlist(input_data[stype])
  #   }else{
  #     Y <- input_data[names(predictions_by_node),stype]
  #   }
  # }
  names(Y) <- names(predictions_by_node)
  
  # Create a ModelInput object and return it.
  newModelInput <- methods::new("ModelInput", A.hat=A_hat, node.wise.prediction=t(predictions_flattened),
                       true.phenotypes=Y, outcome.type=stype.class, 
                       coregulation.graph=igraph::get.adjacency(coregulationGraph, sparse = FALSE), 
                       line.graph=as.matrix(A))
  return(newModelInput)
}

#' 
#' A wrapper for the formatInput class.
#' @param inputData List of MultiDataSet objects (output of CreateCrossValFolds()) 
#' with gene expression,
#' metabolite abundances, and associated meta-data
#' @param predictionGraphs A list of igraph objects, each of which includes
#' predictions for each edge.
#' @param coregulationGraph An igraph object containing the coregulation graph.
#' @param stype.class The class of the outcome ("numeric" or "categorical")
#' @param edgeTypeList List containing one or more of the following to include
#' in the line graph:
#' - "shared.outcome.analyte"
#' - "shared.independent.analyte"
#' - "analyte.chain"
#' @param testing A boolean indicating whether the testing data is to be used
#' (as opposed to the training data). Default is FALSE.
#' @export
formatInputAllFolds <- function(predictionGraphs, coregulationGraph,
                                inputData, stype.class, edgeTypeList,
                                testing = FALSE){
  # Apply to all folds.
  return(lapply(1:length(predictionGraphs), function(i){
    # Select whether input data is training or testing.
    input <- inputData[[i]]$training
    if(testing == TRUE){
      input <- inputData[[i]]$testing
    }
    return(formatInput(predictionGraphs=predictionGraphs[[i]], 
                       coregulationGraph=coregulationGraph[[i]], 
                       inputData = input, stype.class = stype.class, 
                       edgeTypeList = edgeTypeList))
  }))
}

#' #' Create the graph pooling filter, given the adjacency matrix of the input graph.
#' #' @param modelInputs An object of type "ModelInputs".
#' #' @param k The output dimensionality of the filter.
#' #' @param poolType One of "mean", "median", "max", or "min".
#' #' @export
#' CreatePoolingFilter <- function(modelInputs, k, poolType){
#'   
#'   # Perform hierarchical clustering.
#'   hier <- doHierarchicalClustering(modelInputs)
#'   clusters_as_sets <- initializeClusters(hier)
#'   hier <- initializeMergeDataFrame(hier)
#'  
#'   # Find the clusters.
#'   clusters <- findKClusters(modelInputs=modelInputs, clusters=clusters_as_sets,
#'                            hClustResults=hier, k=k, allClusters={}, allVariances = {})
#'   cluster_names <- sort(names(clusters))
#'   
#'   # Arrange cluster mappings in matrix.
#'   mappings <- matrix(0, ncol = k, nrow = length(igraph::V(g)))
#'   cluster_length <- matrix(0, length(unique(clusters)))
#'   for(i in 1:length(cluster_names)){
#'     which_in_cluster <- which(rownames(modelInputs@line.graph) %in% 
#'                                 unlist(clusters[[cluster_names[i]]]))
#'     mappings[which_in_cluster, i] <- 1
#'     cluster_length[i] <- length(which_in_cluster)
#'   }
#'   
#'   # Return pooling filter.
#'   newPoolingFilter <- methods::new("PoolingFilter", filter=mappings, filter.type=poolType,
#'                                    cluster.sizes=cluster_length, individual.filters=list())
#'   return(newPoolingFilter)
#' }

#' Create the graph pooling filter, given the adjacency matrix of the input graph.
#' @param modelInputs An object of type "ModelInputs".
#' @param k The output dimensionality of the filter.
#' @param poolType One of "mean", "median", "max", or "min".
#' @export
CreatePoolingFilterKMeans <- function(modelInputs, k, poolType){
  
  # Extract the normalized graph Laplacian.
  graph <- modelInputs@A.hat
  
  # Find the eigenvectors.
  eigenvecs <- eigen(graph, symmetric=TRUE)
  eigenvecs_k <- eigenvecs$vectors[,(ncol(eigenvecs$vectors)-k+1):ncol(eigenvecs$vectors)]
  
  # Cluster.
  kmeans_result <- stats::kmeans(eigenvecs_k, centers=k)
  
  # Arrange cluster mappings in matrix.
  mappings <- matrix(0, ncol = k, nrow = dim(graph)[1])
  cluster_length <- matrix(0, length(unique(kmeans_result$cluster)))
  for(i in 1:k){
    which_in_cluster <- which(kmeans_result$cluster == i)
    mappings[which_in_cluster, i] <- 1
    cluster_length[i] <- length(which_in_cluster)
  }

  
  # Return pooling filter.
  newPoolingFilter <- methods::new("PoolingFilter", filter=mappings, filter.type=poolType,
                         cluster.sizes=cluster_length, individual.filters=list())

  # Return.
  return(newPoolingFilter)
}

#' Create a line graph given the original, unweighted graph. In a line graph,
#' edges are nodes, and edges connected by a node are edges.
#' @param predictionsByEdge Prediction levels corresponding to each edge.
#' @param graphWithPredictions Original graph data frame.
#' @param edgeTypeList List containing one or more of the following:
#' - "shared.outcome.analyte"
#' - "shared.independent.analyte"
#' - "analyte.chain"
CreateLineGraph <- function(predictionsByEdge, graphWithPredictions, edgeTypeList){
  # Step 1: Define vertices.
  line_graph_vertices <- predictionsByEdge$Node
  list_to_add <- list()
  
  # Step 2: Identify edges that share "to" nodes.
  if("shared.outcome.analyte" %in% edgeTypeList){
    list_to_add[[length(list_to_add)+1]] <- FindEdgesSharingNodes(predictionsByEdge = predictionsByEdge,
                                          graphWithPredictions = graphWithPredictions,
                                          nodeType1 = "to", nodeType2 = "to")
  }

  # Step 3: Identify edges that share "from" nodes.
  if("shared.independent.analyte" %in% edgeTypeList){
    list_to_add[[length(list_to_add)+1]] <- FindEdgesSharingNodes(predictionsByEdge = predictionsByEdge,
                                                                graphWithPredictions = graphWithPredictions,
                                                                nodeType1 = "from", nodeType2 = "from")
  }

  # Step 4: For each "to" node, identify edges that share a "from" node.
  if("analyte.chain" %in% edgeTypeList){
    list_to_add[[length(list_to_add)+1]] <- FindEdgesSharingNodes(predictionsByEdge = predictionsByEdge,
                                                                graphWithPredictions = graphWithPredictions,
                                                                nodeType1 = "to", nodeType2 = "from")
  }

  # Step 6: Concatenate.
  shared_df <- list_to_add[[1]]
  if(length(list_to_add) > 1){
    shared_df <- do.call(rbind, list_to_add)
  }
  
  # Step 7: Convert to graph.
  line_graph <- igraph::graph_from_data_frame(shared_df)
  
  # Step 8: Convert to adjacency matrix.
  line_graph_A <- igraph::as_adj(line_graph)
  
  # Return.
  return(line_graph_A)
}

#' Find edges that share nodes and add them to a data frame.
#' @param predictionsByEdge Prediction levels corresponding to each edge.
#' @param graphWithPredictions Original graph data frame.
#' @param nodeType1 Either "to" or "from".
#' @param nodeType2 Either "to or "from".
FindEdgesSharingNodes <- function(predictionsByEdge, graphWithPredictions, nodeType1,
                                  nodeType2){
  # Convert predictions to data frame.
  graph_df <- igraph::as_data_frame(graphWithPredictions)
  
  # Find shared nodes.
  nodes <- unique(graph_df[,nodeType1])
  to_shared <- lapply(1:length(nodes), function(i){
    
    # Find all line graph nodes starting with or ending with the analyte in question.
    node <- nodes[i]
    combs <- NULL
    set_with_1 <- predictionsByEdge$Node[which(graph_df[,nodeType1] == node)]
    set_with_2 <- predictionsByEdge$Node[which(graph_df[,nodeType2] == node)]

    # If there are multiple line graph nodes including this analyte, return them.
    combs <- expand.grid(set_with_1, set_with_2)
    combs$Var1 <- as.character(combs$Var1)
    combs$Var2 <- as.character(combs$Var2)
    combs <- combs[which(combs$Var1 != combs$Var2),]
    colnames(combs) <- c("to", "from")
    return(combs)
  })
  
  # Concatenate all shared nodes.
  line_graph_df <- do.call(rbind, to_shared)

  return(line_graph_df)
}

#' Find edges that share nodes and add them to a data frame.
#' @param modelInputs An object of the ModelInput class.
#' @param poolingFilter A matrix that pools convolution results.
#' @param iterations Maximum number of iterations.
#' @param convergenceCutoff Cutoff for convergence.
#' @param outcomeType One of either "character" or "numeric"
#' @param learningRate Learning rate to use during training
#' @param activationType Activation function. May be "softmax", "sigmoid", 
#' or "tanh".
#' @param weightsAfterPooling Whether to include the weights after the pooling
#' @param optimizationType Type of optimization. May be "BGD", "momentum",
#' "adagrad", or "adam".
#' operation (as opposed to before). Must be TRUE or FALSE.
#' @export
InitializeGraphLearningModel <- function(modelInputs, poolingFilter, iterations,
                                         convergenceCutoff, outcomeType, learningRate,
                                         activationType, weightsAfterPooling,
                                         optimizationType){
  
  # Initialize data frame with maximum number of iterations.
  weights_count <- dim(modelInputs@A.hat)[1]
  wt_name <- rownames(modelInputs@node.wise.prediction)
  if(weightsAfterPooling == TRUE){
    weights_count <- dim(poolingFilter@filter)[2]
    wt_name <- rep(1:dim(poolingFilter@filter)[2])
  }
  tracking.frame <- as.data.frame(matrix(-1, nrow = iterations, 
                                         ncol = 2 + (2 * weights_count)))
  tracking.frame.cnames <- c("Iteration", "Error")
  tracking.frame.cnames <- c(tracking.frame.cnames, paste("Weight", wt_name, sep = "_"))
  tracking.frame.cnames <- c(tracking.frame.cnames, paste("Gradient", wt_name, sep = "_"))
  colnames(tracking.frame) <- tracking.frame.cnames
  tracking.frame$Error[1] <- .Machine$double.xmax
  tracking.frame$Iteration[1] <- 0

  # Initialize weights with uniform distribution.
  max_phen <- max(modelInputs@true.phenotypes)
  num_nodes <- dim(modelInputs@node.wise.prediction)[1]
  weights <- as.matrix(rep(1 / num_nodes, num_nodes))
  tracking.frame[1,3:(2+weights_count)] <- weights
  
  # Initialize and return results.
  newModelResults <- methods::new("ModelResults", model.input=modelInputs, pooling.filter=poolingFilter,
                        iteration.tracking=tracking.frame, max.iterations=iterations,
                        convergence.cutoff=convergenceCutoff, learning.rate=learningRate,
                        previous.weights=as.matrix(rep(0,length(weights))), 
                        current.weights=as.matrix(weights),
                        current.gradient=as.matrix(rep(-1,length(weights))),
                        previous.momentum=as.matrix(rep(0,length(weights))),
                        previous.update.vector=as.matrix(rep(0,length(weights))),
                        sum.square.gradients=as.matrix(rep(0,length(weights))),
                        current.iteration=0, activation.type=activationType,
                        weights.after.pooling=weightsAfterPooling,
                        optimization.type=optimizationType)
  return(newModelResults)
}

#' Run least squares optimization to optimize the weights
#' @param modelInput An object of the ModelInput class.
#' @param convolutions Number of convolutions.
#' @param maskPercentile The sum of error deltas will be computed over all
#' predictors. The lowest maskPercentile of predictors will be used to ensure
#' that reliable predictors are selected. By default, this value is set to 1,
#' which means that all predictors will be used.
#' @param nrmseCutoff This parameter allows you to use only predictors that
#' have NRMSE below a given cutoff. If NULL, a cutoff is not used. Default is NULL.
#' @param corrCutoff This parameter allows you to use only predictors that
#' have correlation above a given cutoff. If NULL, a cutoff is not used. Default is NULL.
#' @param ridgeLambda The lambda value to be used in ridge regression. The
#' higher the value of lambda, the more the sum of weights is minimized.
#' This encourages smaller weights. Default is 0.
#' @export
RunLeastSquaresOptimization <- function(modelInput, convolutions,
                                        maskPercentile = 1,
                                        nrmseCutoff = NULL,
                                        corrCutoff = NULL,
                                        ridgeLambda = 0){
  # Compute the value of Y.
  A.hat <- modelInput@A.hat
  X <- modelInput@node.wise.prediction
  Y.formula <- X
  
  # Convolve the input.
  if(convolutions > 0){
    
    # Modify the convolutional matrix if more than
    # one convolution is desired.
    if(convolutions > 1){
      for(c in 2:convolutions){
        A.hat.old <- A.hat
        A.hat <- A.hat.old %*% A.hat.old
      }
    }
    
    # Do convolution.
    Y.formula.list <- lapply(1:dim(X)[2], function(i){
      return(A.hat %*% X[,i])
    })
    Y.formula <- do.call(cbind, Y.formula.list)
  }
  
  # Transpose the Y values.
  Y.formula <- t(Y.formula)
  
  # Filter by NRMSE if applicable.
  which_less_nrmse <- seq(1:dim(Y.formula)[2])
  if(!is.null(nrmseCutoff)){
    
    nrmse <- unlist(lapply(1:dim(Y.formula)[2], function(i){
      return(ComputeNRMSE(modelInput@true.phenotypes, 
                          Y.formula[,i]))
    }))
    which_less <- which(nrmse <= nrmseCutoff)
    which_less_nrmse <- which_less
    Y.formula <- Y.formula[,which_less]
  }
  
  # Filter by correlation if applicable.
  which_less_corr <- which_less_nrmse
  if(!is.null(corrCutoff)){
    correlation <- unlist(lapply(1:dim(Y.formula)[2], function(i){
      return(stats::cor(modelInput@true.phenotypes, Y.formula[,i], method = "spearman"))
    }))
    which_more <- which(correlation >= corrCutoff)
    which_less_corr <- which_less_nrmse[which_more]
    Y.formula <- Y.formula[,which_more]
  }
  
  # Solve Least Squares. Theta = (A^TA)^-1 * A^Tb
  # Because computing the inverse is very slow for large matrices, we follow
  # the advice of John Cook ("Don't Invert That Matrix") and solve it as:
  # (A^TA)Theta = A^Tb. See https://www.r-bloggers.com/2015/07/dont-invert-that-matrix-why-and-how/
  #Theta <- solve(t(Y.formula) %*% Y.formula, t(Y.formula) %*% modelInput@true.phenotypes)
  
  # Use SVD to solve.
  decomp <- svd(Y.formula)
  Theta_nonzero <- rep(0, length(decomp$d))
  if(length(decomp$d) > 1){
    S_inv <- solve(diag(decomp$d^2) + ridgeLambda * diag(rep(1, length(decomp$d))))
    Theta_nonzero <- decomp$v %*% S_inv %*% diag(decomp$d) %*% t(decomp$u) %*% 
      as.matrix(modelInput@true.phenotypes)
  }else{
    S_inv <- decomp$v * 1 / ((decomp$d^2) + ridgeLambda)
    Theta_nonzero <- decomp$v * S_inv * decomp$d * t(decomp$u) %*%
      as.matrix(modelInput@true.phenotypes)
  }
  
  # Build final theta, including zeros for predictors not within the cutoff.
  Theta <- rep(0, dim(modelInput@node.wise.prediction)[1])
  Theta[which_less_corr] <- Theta_nonzero
  names(Theta) <- rownames(modelInput@node.wise.prediction)
  
  return(Theta)
}

#' Wrapper for RunLeastSquaresOptimization.
#' @param modelInput A list of objects of the ModelInput class.
#' @param convolutions Number of convolutions to perform.
#' @param maskPercentile The sum of error deltas will be computed over all
#' predictors. The lowest maskPercentile of predictors will be used to ensure
#' that reliable predictors are selected. By default, this value is set to 1,
#' which means that all predictors will be used.
#' @param ridgeLambda The lambda value to be used in ridge regression. The
#' higher the value of lambda, the more the sum of weights is minimized.
#' This encourages smaller weights. Default is 0.
#' @param nrmseCutoff This parameter allows you to use only predictors that
#' have NRMSE below a given cutoff. If NULL, a cutoff is not used. Default is NULL.
#' @param corrCutoff This parameter allows you to use only predictors that
#' have correlation above a given cutoff. If NULL, a cutoff is not used. Default is NULL.
#' @export
RunLeastSquaresOptimizationAllFolds <- function(modelInput, convolutions,
                                                maskPercentile = 1,
                                                nrmseCutoff = NULL,
                                                corrCutoff = NULL,
                                                ridgeLambda = 0){
  result <- lapply(1:length(modelInput), function(i){
    return(RunLeastSquaresOptimization(modelInput=modelInput[[i]],
                                       convolutions=convolutions,
                                       maskPercentile = maskPercentile,
                                       nrmseCutoff = nrmseCutoff,
                                       corrCutoff = corrCutoff,
                                       ridgeLambda = ridgeLambda))
  })
  names(result) <- paste("Fold", 1:length(modelInput), sep = "_")
  return(result)
}

#' Train the graph learning model, using the specifications in the ModelResults.
#' class and storing the results in the ModelResults class.
#' @param modelResults An object of the ModelResults class.
#' @param pooling Whether or not to perform pooling during training.
#' @param convolution Whether or not to perform convolution during training.
#' @param stochastic Whether or not training should be stochastic.
#' @param ridgeRegressionWeight The hyperparameter weight assigned
#' to the ridge regression parameter (often referred to as lambda in the
#' literature)
#' @param varianceWeight The hyperparameter weight assigned to the difference
#' in variances between Y and the predicted value of Y.
#' @param verbose Whether to print results as you run the model.
#' @export
TrainGraphLearningModel <- function(modelResults, pooling, convolution,
                                    stochastic, ridgeRegressionWeight,
                                    varianceWeight, verbose = TRUE){
  
  # Start the first iteration and calculate a dummy weight delta.
  modelResults@current.iteration <- 1
  modelResults@iteration.tracking$Iteration[modelResults@current.iteration+1]<-
    modelResults@current.iteration
  weight.delta <- sqrt(sum((modelResults@current.weights - modelResults@previous.weights)^2))
  
  # Repeat the training process for all iterations, until the maximum is reached
  # or until convergence.
  while(modelResults@current.iteration < (modelResults@max.iterations - 1)
        && (weight.delta > modelResults@convergence.cutoff)){
    # For stochastic training, permute the samples, then compute the gradient one
    # sample at a time.
    # For batch training, compute the gradient over all samples.
    if(stochastic == TRUE){
      # Permute samples.
      perm_samples <- sample(1:dim(modelResults@model.input@node.wise.prediction)[2],
                             dim(modelResults@model.input@node.wise.prediction)[2])
      for(i in perm_samples){
        # Do training iteration for each sample.
        newModelResults <- modelResults
        newModelResults@model.input@node.wise.prediction <- 
          as.matrix(modelResults@model.input@node.wise.prediction[,i])
        newModelResults@model.input@true.phenotypes <- 
          modelResults@model.input@true.phenotypes[i]
        newModelResults <- DoSingleTrainingIteration(newModelResults, modelResults@current.iteration,
                                                  pooling, convolution, ridgeRegressionWeight,
                                                  varianceWeight)
        # Update weights and gradient in the model results according to the
        # results of this sample.
        modelResults@current.weights <- newModelResults@current.weights
        modelResults@previous.weights <- newModelResults@previous.weights
        modelResults@current.gradient <- newModelResults@current.gradient
        modelResults@outcome.prediction <- newModelResults@outcome.prediction
        modelResults@previous.momentum <- newModelResults@previous.momentum
        modelResults@previous.update.vector <- newModelResults@previous.update.vector
        modelResults@iteration.tracking <- newModelResults@iteration.tracking
      }
      # Compute the prediction error over all samples.
      Y.pred <- DoPrediction(modelResults, modelResults@current.iteration,
                                                           pooling, convolution)
      if(modelResults@model.input@outcome.type == "categorical"){
        modelResults@iteration.tracking$Error[modelResults@current.iteration+1] <- 
          ComputeClassificationError(modelResults@model.input@true.phenotypes, Y.pred)
      }else{
        modelResults@iteration.tracking$Error[modelResults@current.iteration+1] <- 
          ComputeNRMSE(modelResults@model.input@true.phenotypes, Y.pred)
      }
      
    # This is the batch case.  
    }else{
      modelResults <- DoSingleTrainingIteration(modelResults, modelResults@current.iteration,
                                                pooling, convolution, ridgeRegressionWeight,
                                                varianceWeight)
    }
    
    # Update the iteration.
    modelResults@current.iteration <- modelResults@current.iteration + 1
    modelResults@iteration.tracking$Iteration[modelResults@current.iteration+1]<-
      modelResults@current.iteration
    
    # Print the weight delta and error.
    weight.delta <- sqrt(sum((modelResults@current.weights - modelResults@previous.weights)^2))
    if(modelResults@current.iteration %% 1 == 0 && verbose == TRUE){
      print(paste("iteration", modelResults@current.iteration, ": weight delta is", weight.delta,
                  "and error is", 
                  modelResults@iteration.tracking$Error[modelResults@current.iteration]))
    }
  }
  # If we exited before the maximum number of iterations, remove the rest of the
  # tracking data.
  if(modelResults@current.iteration < modelResults@max.iterations){
    modelResults@iteration.tracking <- modelResults@iteration.tracking[1:modelResults@current.iteration,]
  }
  return(modelResults)
}
  
#' Predict Y given current weights.
#' @param modelResults An object of the ModelResults class.
#' @param iteration The current iteration.
#' @param pooling Whether or not to pool the weights.
#' @param convolution Whether or not to perform convolution.
DoPrediction <- function(modelResults, iteration, pooling, convolution){
  # Propagate forward.
  A.hat <- modelResults@model.input@A.hat
  X <- modelResults@model.input@node.wise.prediction
  Theta.old <- matrix(rep(modelResults@current.weights, dim(X)[2]), ncol = dim(X)[2])
  Y.pred <- X
  # If convolution is to be performed, perform convolution.
  if(convolution == TRUE){
    Y.pred.list <- lapply(1:dim(X)[2], function(i){
      return(A.hat %*% X[,i])
    })
    Y.pred <- do.call(cbind, Y.pred.list)
  }
  # If pooling is to be performed, multiply by pool either after or before weights
  # are learned, respectively.
  if(pooling == TRUE){
    if(modelResults@weights.after.pooling == TRUE){
      S_all <- modelResults@pooling.filter@individual.filters
      if(modelResults@current.iteration == 1){
        S_all <- CreateFilter(poolingFilter = modelResults@pooling.filter, A.hat = A.hat, 
                              X = X)
        modelResults@pooling.filter@individual.filters <- S_all
      }
      Y.pred <- unlist(lapply(1:length(S_all), function(i){
        return(sum(t(Y.pred[,i]) %*% S_all[[i]] * Theta.old[,i]))
      }))
    }else{
      S_all <- AdjustFilter(poolingFilter = modelResults@pooling.filter, A.hat = A.hat, 
                            X = X, Theta = Theta.old)
      modelResults@pooling.filter@individual.filters <- S_all
      Y.pred <- unlist(lapply(1:length(S_all), function(i){
        return(sum(t(Y.pred[,i] * Theta.old[,i]) %*% S_all[[i]]))
      }))
    }
  }else{
    Y.pred <- unlist(lapply(1:(dim(Theta.old)[2]), function(i){
      return(sum(t(Y.pred[,i] * Theta.old[,i])))
    }))
  }
  
  # Use activation function if output is of a character type. Note that all character
  # types are converted into factors, and since only binary factors are accepted by
  # the package, the values will be 1 (for the alphanumerically lowest level) and 2
  # (for the alphanumerically highest level).
  if(modelResults@model.input@outcome.type == "categorical"){
    if(modelResults@activation.type == "softmax"){
      Y.pred <- round(SoftmaxWithCorrection(Y.pred))
    }else if(modelResults@activation.type == "tanh"){
      Y.pred <- round(TanhWithCorrection(Y.pred))
    }else{
      Y.pred <- round(SigmoidWithCorrection(Y.pred))
    }
  }else{
    Y.pred <- Y.pred / dim(Theta.old)[1]
  }
  
  # Return the prediction.
  return(Y.pred)
}

#' Train the graph learning model, using the specifications in the ModelResults
#' class and storing the results in the ModelResults class.T
#' @param modelResults An object of the ModelResults class.
#' @param iteration The current iteration.
#' @param pooling Whether or not to pool the weights.
#' @param convolution Whether or not to perform convolution.
#' @param ridgeRegressionWeight The hyperparameter weight assigned
#' to the ridge regression parameter (often referred to as lambda in the
#' literature)
#' @param varianceWeight The hyperparameter weight assigned to the difference
#' in variances between Y and the predicted value of Y.
DoSingleTrainingIteration <- function(modelResults, iteration, pooling, convolution,
                                      ridgeRegressionWeight, varianceWeight){
  # Predict Y.
  A.hat <- modelResults@model.input@A.hat
  X <- modelResults@model.input@node.wise.prediction
  Theta.old <- matrix(rep(modelResults@current.weights, dim(X)[2]), ncol = dim(X)[2])
  Y.pred <- DoPrediction(modelResults, iteration, pooling, convolution)
  modelResults@outcome.prediction <- Y.pred
  
  # Backpropagate and calculate the error.
  Theta.new <- Theta.old
  error <- modelResults@iteration.tracking$Error[iteration-1]
  modelResults <- BackpropagateSingleLayer(modelResults, iteration, convolution,
                                           pooling, ridgeRegressionWeight,
                                           varianceWeight)
  if(modelResults@model.input@outcome.type == "categorical"){
    modelResults@iteration.tracking$Error[iteration+1] <- 
      ComputeClassificationError(modelResults@model.input@true.phenotypes, Y.pred)
  }else{
    modelResults@iteration.tracking$Error[iteration+1] <- 
      ComputeNRMSE(modelResults@model.input@true.phenotypes, Y.pred)
  }
  
  # Modify the model results and return.
  return(modelResults)
}

#' Run a prediction on new data using the graph learning model, and compute
#' the absolute error delta values.
#' @param weights The list of all learned weights.
#' @param convolutions Number of convolutions to perform. This should be the
#' same number of convolutions used in learning the weights.
#' @param modelInput A list of objects of the ModelInput class.
#' @param minimum The minimum prediction value to allow. Default is NULL.
#' @param maximum The maximum prediction value to allow. Default is NULL.
#' @export
GetErrorDeltas <- function(weights, modelInput, convolutions, minimum = NULL, 
                           maximum = NULL){
  # Convolve the input.
  Y.formula <- modelInput@node.wise.prediction
  if(convolutions > 0){
    
    # Modify the convolutional matrix if more than
    # one convolution is desired.
    if(convolutions > 1){
      for(c in 2:convolutions){
        modelInput@A.hat <- modelInput@A.hat %*% modelInput@A.hat
      }
    }
    
    # Do convolution.
    Y.formula <- modelInput@A.hat %*% Y.formula
  }
  deltas <- unlist(lapply(1:length(modelInput@true.phenotypes),function(j){
    solution <- sum(Y.formula[,j] * weights)
    # Adjust to fit minimum and maximum.
    if(!is.null(minimum)){
      solution[which(solution < minimum)] <- minimum
    }
    if(!is.null(maximum)){
      solution[which(solution > maximum)] <- maximum
    }  
    phenotype <- modelInput@true.phenotypes[j]
    return(unlist(unname(abs(solution - phenotype))))
  }))
  names(deltas) <- names(modelInput@true.phenotypes)
  return(deltas)
}

#' Run a prediction on new data using the graph learning model, and compute
#' the absolute error delta values.
#' @param weights The list of all learned weights.
#' @param convolutions Number of convolutions to perform. This should be the
#' same number of convolutions used in learning the weights.
#' @param modelInput A list of objects of the ModelInput class.
#' @param minimum The minimum prediction value to allow. Default is NULL.
#' @param maximum The maximum prediction value to allow. Default is NULL.
#' @export
GetErrorDeltasAllFolds <- function(weights, modelInput, convolutions, minimum = NULL, 
                                   maximum = NULL){
  return(lapply(1:length(modelInput), function(i){
    # Get error deltas.
    return(GetErrorDeltas(weights=weights[[i]], modelInput=modelInput[[i]], 
                          convolutions=convolutions, minimum=minimum,
                          maximum=maximum))
  }))
}
#' Run a prediction on new data using the graph learning model.
#' @param modelResults An object of the ModelResults class.
#' @param pooling Whether or not to pool the weights.
#' @param convolution Whether or not to perform convolution.
#' @param testInput An object of the ModelInput class.
#' @export
PredictTesting <- function(modelResults, pooling, convolution, testInput){
  # Propagate forward.
  A.hat <- modelResults@model.input@A.hat
  X <- testInput@node.wise.prediction
  Theta.old <- matrix(rep(modelResults@current.weights, dim(X)[2]), ncol = dim(X)[2])
  Y.pred <- X
  # If convolution is to be performed, perform convolution.
  if(convolution == TRUE){
    Y.pred.list <- lapply(1:dim(X)[2], function(i){
      return(A.hat %*% X[,i])
    })
    Y.pred <- do.call(cbind, Y.pred.list)
  }
  # If pooling is to be performed, multiply by pool either after or before weights
  # are learned, respectively.
  if(pooling == TRUE){
    if(modelResults@weights.after.pooling == TRUE){
      S_all <- modelResults@pooling.filter@individual.filters
      if(modelResults@current.iteration == 1){
        S_all <- CreateFilter(poolingFilter = modelResults@pooling.filter, A.hat = A.hat, 
                              X = X)
        modelResults@pooling.filter@individual.filters <- S_all
      }
      Y.pred <- unlist(lapply(1:length(S_all), function(i){
        return(sum(t(Y.pred[,i]) %*% S_all[[i]] * Theta.old[,i]))
      }))
    }else{
      S_all <- AdjustFilter(poolingFilter = modelResults@pooling.filter, A.hat = A.hat, 
                            X = X, Theta = Theta.old)
      modelResults@pooling.filter@individual.filters <- S_all
      Y.pred <- unlist(lapply(1:length(S_all), function(i){
        return(sum(t(Y.pred[,i] * Theta.old[,i]) %*% S_all[[i]]))
      }))
    }
  }else{
    Y.pred <- unlist(lapply(1:(dim(Theta.old)[2]), function(i){
      return(sum(t(Y.pred[,i] * Theta.old[,i])))
    }))
  }
  
  
  # Use activation function if output is of a character type. Note that all character
  # types are converted into factors, and since only binary factors are accepted by
  # the package, the values will be 1 (for the alphanumerically lowest level) and 2
  # (for the alphanumerically highest level).
  if(modelResults@model.input@outcome.type == "categorical"){
    if(modelResults@activation.type == "softmax"){
      Y.pred <- round(SoftmaxWithCorrection(Y.pred))
    }else if(modelResults@activation.type == "tanh"){
      Y.pred <- round(TanhWithCorrection(Y.pred))
    }else{
      Y.pred <- round(SigmoidWithCorrection(Y.pred))
    }
  }else{
    Y.pred <- Y.pred / dim(Theta.old)[1]
  }
  
  # Calculate error.
  error <- 1
  if(modelResults@model.input@outcome.type == "categorical"){
    error <- ComputeClassificationError(testInput@true.phenotypes, Y.pred)
  }else{
    error <- ComputeNRMSE(testInput@true.phenotypes, Y.pred)
  }
  
  # Modify the model results and return.
  return(list("Y.pred" = Y.pred, "Error" = error))
}

#' Compute classification error.
#' @param true.Y The true phenotype of each sample.
#' @param pred.Y The predicted phenotype of each sample.
ComputeClassificationError <- function(true.Y, pred.Y){
  # Find false and true positives and negatives.
  FP <- length(intersect(which(true.Y == 1), which(pred.Y == 2)))
  TP <- length(intersect(which(true.Y == 2), which(pred.Y == 2)))
  FN <- length(intersect(which(true.Y == 2), which(pred.Y == 1)))
  TN <- length(intersect(which(true.Y == 1), which(pred.Y == 1)))

  # Compute error and return.
  error <- (FP + FN) / (FP + FN + TP + TN)
  return(error)
}

#' Compute the normalized root mean squared error.
#' @param true.Y The true phenotype of each sample.
#' @param pred.Y The predicted phenotype of each sample.
#' @export
ComputeNRMSE <- function(true.Y, pred.Y){
  RMSD <- sqrt(sum((true.Y - pred.Y)^2) / length(true.Y))
  NRMSE <- RMSD
  if(length(true.Y) > 1){
    NRMSE <- RMSD / (max(true.Y) - min(true.Y))
  }
  return(NRMSE)
}

#' Run a prediction on new data using the graph learning model, and compute
#' the NRMSE values.
#' @param weights The list of all learned weights.
#' @param convolutions Number of convolutions to perform. This should be the
#' same number of convolutions used in learning the weights.
#' @param modelInput A list of objects of the ModelInput class.
#' @param minimum The minimum prediction value to allow. Default is NULL.
#' @param maximum The maximum prediction value to allow. Default is NULL.
#' @export
GetNRMSE <- function(weights, modelInput, convolutions, minimum = NULL, 
                           maximum = NULL){
  # Convolve the input.
  Y.formula <- modelInput@node.wise.prediction
  if(convolutions > 0){
    
    # Modify the convolutional matrix if more than
    # one convolution is desired.
    if(convolutions > 1){
      for(c in 2:convolutions){
        modelInput@A.hat <- modelInput@A.hat %*% modelInput@A.hat
      }
    }
    
    # Do convolution.
    Y.formula <- modelInput@A.hat %*% Y.formula
  }
  solutions <- unlist(lapply(1:length(modelInput@true.phenotypes),function(j){
    solution <- sum(Y.formula[,j] * weights)
    return(unlist(solution))
  }))
  
  NRMSE <- ComputeNRMSE(modelInput@true.phenotypes, solutions)
  return(NRMSE)
}

#' Run a prediction on new data using the graph learning model, and compute
#' the correlations.
#' @param weights The list of all learned weights.
#' @param convolutions Number of convolutions to perform. This should be the
#' same number of convolutions used in learning the weights.
#' @param modelInput A list of objects of the ModelInput class.
#' @param minimum The minimum prediction value to allow. Default is NULL.
#' @param maximum The maximum prediction value to allow. Default is NULL.
#' @export
GetCorrelation <- function(weights, modelInput, convolutions, minimum = NULL, 
                     maximum = NULL){
  # Convolve the input.
  Y.formula <- modelInput@node.wise.prediction
  if(convolutions > 0){
    
    # Modify the convolutional matrix if more than
    # one convolution is desired.
    if(convolutions > 1){
      for(c in 2:convolutions){
        modelInput@A.hat <- modelInput@A.hat %*% modelInput@A.hat
      }
    }
    
    # Do convolution.
    Y.formula <- modelInput@A.hat %*% Y.formula
  }
  solutions <- unlist(lapply(1:length(modelInput@true.phenotypes),function(j){
    solution <- sum(Y.formula[,j] * weights)
    # Adjust to fit minimum and maximum.
    if(!is.null(minimum)){
      solution[which(solution < minimum)] <- minimum
    }
    if(!is.null(maximum)){
      solution[which(solution > maximum)] <- maximum
    }  
    return(unlist(solution))
  }))
  
  corr <- stats::cor(modelInput@true.phenotypes, solutions, method = "spearman")
  return(corr)
}

#' Wrapper for GetNRMSE.
#' @param weights The list of all learned weights.
#' @param convolutions Number of convolutions to perform. This should be the
#' same number of convolutions used in learning the weights.
#' @param modelInput A list of objects of the ModelInput class.
#' @param minimum The minimum prediction value to allow. Default is NULL.
#' @param maximum The maximum prediction value to allow. Default is NULL.
#' @export
GetNRMSEAllFolds <- function(weights, modelInput, convolutions, minimum = NULL, 
                                   maximum = NULL){
  return(lapply(1:length(modelInput), function(i){
    # Get error deltas.
    return(GetNRMSE(weights=weights[[i]], modelInput=modelInput[[i]], 
                          convolutions=convolutions, minimum=minimum,
                          maximum=maximum))
  }))
}