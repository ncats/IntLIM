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
#' @param stype The outcome
#' @param stype.class The class of the outcome ("numeric" or "categorical")
#' @export
formatTrainingInput <- function(predictionGraphs, inputData, stype, stype.class){
  
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
                       graphWithPredictions = predictionGraphs[[1]])
  predictions_flattened <- predictions_flattened[,colnames(A)]

  # Add self-loops.
  A_tilde <- as.matrix(A)
  diag(A_tilde) <- 1
  
  # Extract the diagonal and raise to the negative half power.
  diags <- colSums(A_tilde)
  diags_neg_half <- 1 / sqrt(diags)
  D_tilde_neg_half <- matrix(0, nrow = nrow(A_tilde), ncol = ncol(A_tilde))
  diag(D_tilde_neg_half) <- diags_neg_half
  
  # Obtain the final matrix.
  A_hat <- D_tilde_neg_half %*% A_tilde %*% D_tilde_neg_half
  
  # Obtain the predictions.
  input_data <- inputData@phenoData$expression$main@data
  Y <- input_data[names(predictions_by_node),stype]
  if(stype.class == "categorical"){
    Y <- as.numeric(as.factor(Y))
  }
  names(Y) <- names(predictions_by_node)
  
  # Create a ModelInput object and return it.
  newModelInput <- methods::new("ModelInput", A.hat=A_hat, node.wise.prediction=t(predictions_flattened),
                       true.phenotypes=Y, outcome.type=stype.class)
  return(newModelInput)
}

#' Create the graph pooling filter, given the adjacency matrix of the input graph.
#' @param modelInputs An object of type "ModelInputs".
#' @param k The output dimensionality of the filter.
#' @param poolType One of "mean", "median", "max", or "min".
#' @export
CreatePoolingFilter <- function(modelInputs, k, poolType){
  
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
CreateLineGraph <- function(predictionsByEdge, graphWithPredictions){
  # Step 1: Define vertices.
  line_graph_vertices <- predictionsByEdge$Node
  # Step 2: Identify edges that share "to" nodes.
  to_shared_df <- FindEdgesSharingNodes(predictionsByEdge = predictionsByEdge,
                                        graphWithPredictions = graphWithPredictions,
                                        nodeType1 = "to", nodeType2 = "to")
  # Step 3: Identify edges that share "from" nodes.
  from_shared_df <- FindEdgesSharingNodes(predictionsByEdge = predictionsByEdge,
                                        graphWithPredictions = graphWithPredictions,
                                        nodeType1 = "from", nodeType2 = "from")
  # Step 4: For each "to" node, identify edges that share a "from" node.
  to_from_shared_df <- FindEdgesSharingNodes(predictionsByEdge = predictionsByEdge,
                                          graphWithPredictions = graphWithPredictions,
                                          nodeType1 = "to", nodeType2 = "from")
  # Step 5: For each "from" node, identify edges that share a "to" node.
  from_to_shared_df <- FindEdgesSharingNodes(predictionsByEdge = predictionsByEdge,
                                             graphWithPredictions = graphWithPredictions,
                                             nodeType1 = "from", nodeType2 = "to")
  # Step 6: Concatenate.
  shared_df <- do.call(rbind, list(to_shared_df, from_shared_df, to_from_shared_df,
                                   from_to_shared_df))
  # Step 7: Convert to graph.
  line_graph <- igraph::graph_from_data_frame(shared_df, directed = FALSE)
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
  graph_df <- igraph::as_data_frame(graphWithPredictions)
  to_shared <- lapply(unique(graph_df[,nodeType1]), function(node){
    combs <- NULL
    set_starting_with <- predictionsByEdge$Node[which(graph_df[,nodeType1] == node)]
    if(length(set_starting_with) > 1){
      if(nodeType1 == nodeType2){
        combs <- as.data.frame(t(utils::combn(set_starting_with, 2)))
        colnames(combs) <- c("to", "from")
      }else{
        set_ending_with <- predictionsByEdge$Node[which(graph_df[,nodeType2] == node)]
        set_starting_with_ext <- unlist(lapply(set_starting_with, function(n){
          return(rep(n, length(set_ending_with)))
        }))
        combs <- data.frame(to = set_starting_with_ext,
                            from = rep(set_ending_with, length(set_starting_with)))
      }
    }
    return(combs)
  })
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
  weights <- as.matrix(stats::runif(weights_count))
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

#' Train the graph learning model, using the specifications in the ModelResults.
#' class and storing the results in the ModelResults class.
#' @param modelResults An object of the ModelResults class.
#' @export
TrainGraphLearningModel <- function(modelResults){
  
  # Start the first iteration and calculate a dummy weight delta.
  modelResults@current.iteration <- 1
  modelResults@iteration.tracking$Iteration[modelResults@current.iteration+1]<-
    modelResults@current.iteration
  weight.delta <- sqrt(sum((modelResults@current.weights - modelResults@previous.weights)^2))
  
  # Repeat the training process for all iterations, until the maximum is reached
  # or until convergence.
  while(modelResults@current.iteration < modelResults@max.iterations - 1
        && (weight.delta > modelResults@convergence.cutoff) || 
        (modelResults@current.iteration < 10)){
    modelResults <- DoSingleTrainingIteration(modelResults, modelResults@current.iteration)
    modelResults@current.iteration <- modelResults@current.iteration + 1
    modelResults@iteration.tracking$Iteration[modelResults@current.iteration+1]<-
      modelResults@current.iteration
    weight.delta <- sqrt(sum((modelResults@current.weights - modelResults@previous.weights)^2))
    if(modelResults@current.iteration %% 1 == 0){
      print(paste("iteration", modelResults@current.iteration, ": weight delta is", weight.delta,
                  "and error is", 
                  modelResults@iteration.tracking$Error[modelResults@current.iteration]))
    }
  }
  if(modelResults@current.iteration < modelResults@max.iterations){
    modelResults@iteration.tracking <- modelResults@iteration.tracking[1:modelResults@current.iteration,]
  }
  return(modelResults)
}
  
#' Train the graph learning model, using the specifications in the ModelResults
#' class and storing the results in the ModelResults class.
#' @param modelResults An object of the ModelResults class.
#' @param iteration The current iteration.
DoSingleTrainingIteration <- function(modelResults, iteration){
  # Propagate forward.
  A.hat <- modelResults@model.input@A.hat
  X <- modelResults@model.input@node.wise.prediction
  Theta.old <- matrix(rep(modelResults@current.weights, dim(X)[2]), ncol = dim(X)[2])
  if(modelResults@weights.after.pooling == TRUE){
    S_all <- modelResults@pooling.filter@individual.filters
    if(modelResults@current.iteration == 1){
      S_all <- CreateFilter(poolingFilter = modelResults@pooling.filter, A.hat = A.hat, 
                            X = X)
      modelResults@pooling.filter@individual.filters <- S_all
    }
    Y.pred <- unlist(lapply(1:length(S_all), function(i){
      return(sum(t(A.hat %*% X[,i]) %*% S_all[[i]] * Theta.old[,i]))
    }))
  }else{
    S_all <- AdjustFilter(poolingFilter = modelResults@pooling.filter, A.hat = A.hat, 
                          X = X, Theta = Theta.old)
    modelResults@pooling.filter@individual.filters <- S_all
    Y.pred <- unlist(lapply(1:length(S_all), function(i){
      return(sum(t(A.hat %*% X[,i] * Theta.old[,i]) %*% S_all[[i]]))
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
  }
  modelResults@outcome.prediction <- Y.pred
  
  # Backpropagate and calculate the error.
  Theta.new <- Theta.old
  error <- modelResults@iteration.tracking$Error[iteration-1]
  modelResults <- BackpropagateSingleLayer(modelResults, iteration)
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
ComputeNRMSE <- function(true.Y, pred.Y){
  RMSD <- sqrt(sum((true.Y - pred.Y)^2) / length(true.Y))
  NRMSE <- RMSD / mean(true.Y)
  return(NRMSE)
}