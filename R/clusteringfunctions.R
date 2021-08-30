#' Find k clusters, optimizing based on the variance of the node
#' values in each cluster (averaged over all samples). This is a recursive function.
#' @param modelInputs An object of the modelInputs class.
#' @param clusters A list of clusters, where each cluster is a set of nodes.
#' @param hClustResults The results of hierarchical clustering as performed using
#' the hclust() function.
#' @param k The number of clusters desired.
#' @param allClusters All clusters that have been explored.
#' @param allVariances All variances that have been computed.
findKClusters <- function(modelInputs, clusters, hClustResults, k, allClusters,
                          allVariances){
  # Initialize clusters to return.
  updated_clusters <- clusters
  
  # Recursive case - include the lowest variance aggregation.
  #while(length(updated_clusters)>length(clusters)-10){
  while(length(updated_clusters)>k){
    
    # Get lowest variance parent and update explored clusters and variances
    # accordingly.
    parentResults <- getLowestVarianceParent(clusters=updated_clusters, modelInputs=modelInputs,
                                             merges=hClustResults$merge, labels=hClustResults$labels,
                                             allClusters=allClusters, allVariances=allVariances)
    p <- parentResults$parent
    p_name <- parentResults$parent_name
    to_remove <- parentResults$to_remove
    hClustResults$merge <- parentResults$merges
    allClusters <- parentResults$allClusters
    allVariances <- parentResults$allVariances

    # Update the final clusters.
    l1 <- length(unlist(updated_clusters[[as.character(to_remove[1])]]))
    l2 <- length(unlist(updated_clusters[[as.character(to_remove[2])]]))
    updated_clusters <- updated_clusters[setdiff(names(updated_clusters), to_remove)]
    updated_clusters[p_name] <- p
    l3 <- length(unlist(updated_clusters[[as.character(p_name)]]))
    print(paste("Replacing clusters", to_remove[1], "(of size", paste0(l1, ") and"),
                to_remove[2], "(of size", paste0(l2, ") with cluster"), p_name,
                "(of size", paste0(l3, " and avg variance ", allVariances[p_name], "),"),
                length(updated_clusters), "remaining"))
  }
  
  # Return updated clusters.
  return(updated_clusters)
}

#' Perform hierarchical clustering of the unweighted graph.
#' @param modelInputs An object of type "ModelInputs".
#' @export
doHierarchicalClustering <- function(modelInputs){
  # Find graph distances between nodes.
  g <- igraph::graph_from_adjacency_matrix(modelInputs@line.graph)
  dist <- igraph::distances(g)
  dist[which(!is.finite(dist))] <- NA
  max_dist <- max(dist[which(!is.na(dist))])
  dist[which(is.na(dist))] <- 2 * max_dist
  
  # Perform hierarchical clustering.
  hier <- stats::hclust(stats::as.dist(dist), method = "average")
  return(hier)
}

#' Initialize a list of clusters for the graph, where each node is included
#' in its own cluster.
#' @param hierarchicalClustering A data frame of merges from the hierarchical clustering.
#' This data frame is output by the function hclust().
#' @export
initializeClusters <- function(hierarchicalClustering){
  # Initialize clusters to be each individual node.
  clusters_as_sets <- lapply(hierarchicalClustering$labels, function(node){
    return(unique(list(node)))
  })
  names(clusters_as_sets) <- paste("C", (-1) * 1:length(hierarchicalClustering$labels), sep = "_")
  
  # Return clusters.
  return(clusters_as_sets)
}

#' Initialize the data frame of merges to include individual cluster and parent
#' cluster information.
#' @param hierarchicalClustering A data frame of merges from the hierarchical clustering.
#' This data frame is output by the function hclust().
#' @export
initializeMergeDataFrame <- function(hierarchicalClustering){
  # Initialize clusters to be each individual node.
  clusters_as_sets <- lapply(hierarchicalClustering$labels, function(node){
    return(unique(list(node)))
  })
  names(clusters_as_sets) <- paste("C", (-1) * 1:length(hierarchicalClustering$labels), sep = "_")
  
  # Create a data frame that includes merge information. This data frame will be
  # modified later.
  merge_df <- as.data.frame(hierarchicalClustering$merge)
  colnames(merge_df) <- c("Merged", "Into")
  merge_df$Merged <- paste("C", merge_df$Merged, sep = "_")
  merge_df$Into <- paste("C", merge_df$Into, sep = "_")
  merge_df$Parent <- paste("C", rownames(merge_df), sep = "_")
  rownames(merge_df) <- merge_df$Parent
  hierarchicalClustering$merge <- merge_df
  
  # Return hierarchical clustering with modified merge DF.
  return(hierarchicalClustering)
}

#' Find the variance of each cluster in a hierarchical clustering.
#' @param modelInputs An object of type "ModelInputs".
#' @export
findHierarchicalClusteringVariance <- function(modelInputs){
  
  # Perform hierarchical clustering.
  hier <- doHierarchicalClustering(modelInputs)
  clusters_as_sets <- initializeClusters(hier)
  hier <- initializeMergeDataFrame(hier)
  temp_merge_df <- hier@merge
  
  # Initialize variance.
  allVariances <- rep(0, length(clusters_as_sets))
  names(allVariances) <- names(clusters_as_sets)
  
  # Find the variance at each merge.
  while(dim(temp_merge_df)[2] > 0){
    
    # Get lowest variance parent and update explored clusters and variances
    # accordingly.
    # parentResults <- getVarianceAllParents(clusters=clusters_as_sets, modelInputs=modelInputs,
    #                                          merges=hClustResults$merge, labels=hClustResults$labels,
    #                                          allClusters=allClusters, allVariances=allVariances)
    # temp_merge_df <- parentResults$merges
    # clusters_as_sets <- parentResults$allClusters
    # allVariances <- parentResults$allVariances
  }
}
  
# Find the lowest variance parent of any of the current clusters.
#' @param modelInputs An object of the modelInputs class.
#' @param clusters A list of clusters, where each cluster is a set of nodes.
#' @param merges The merge slot of hierarchical clustering results as performed using
#' the hclust() function.
#' @param labels The labels for the numeric nodes in the hierarchical clustering
#' algorithm.
#' @param allClusters A list of all clusters investigated, where each cluster is
#' a set of nodes.
#' @param allVariances A list of all variances investigated, where each variance 
#' refers to a single cluster.
getLowestVarianceParent <- function(clusters, modelInputs, merges, labels,
                                    allClusters, allVariances){

  # Find the location where each cluster is merged.
  merge_indices <- unlist(lapply(1:length(names(clusters)), function(i){
    name <- names(clusters)[i]
    return(which(merges$Merged == name))
  }))
  
  # Find the merged "parent" of the cluster being merged and the cluster
  # it is merged with.
  parents <- lapply(1:length(names(clusters)), function(i){
    
    # Get the names of the parent and child clusters.
    name <- names(clusters)[i]
    idx_of_merge <- merge_indices[i]
    parent_name <- rownames(merges)[idx_of_merge]
    being_merged <- merges$Merged[idx_of_merge]
    merged_into <- merges$Into[idx_of_merge]
    
    # Merge the parent and child clusters if the cluster being merged
    # with has already been visited.
    parent_cluster <- NA
    if(merged_into %in% names(clusters)){
      clust1 <- clusters[[being_merged]]
      clust2 <- clusters[[merged_into]]
      parent_cluster <- list(unlist(c(clust1, clust2)))
    }
    
    # Return the parent.
    return(parent_cluster)
  })
  
  # Assign parent names and filter parents that include not previously seen clusters.
  names(parents) <- rownames(merges)[merge_indices]
  parents <- parents[which(!is.na(parents))]
  
  # Update the clusters that have been explored.
  for(parent_name in names(parents)){
    allClusters[parent_name] <- parents[parent_name]
  }

  # Find lowest variance parent.
  variances <- lapply(names(parents), function(parent_name){
    return(computeVariance(clusterName=parent_name, cluster=parents[parent_name], 
                           allVariances=allVariances, modelInputs=modelInputs))
  })
  allVariances[names(parents)] <- variances 
  lowest_var_parent <- parents[which.min(variances)]
  lowest_parent_name <- names(parents)[which.min(variances)]
  
  # List the names of clusters to remove.
  to_rm <- c(merges[lowest_parent_name, "Merged"],
             merges[lowest_parent_name, "Into"])

  # Remove the parent merge from the list of merges.
  merges <- merges[setdiff(rownames(merges), lowest_parent_name),]

  # Return all relevant information.
  parentResults <- list(parent=lowest_var_parent, parent_name=lowest_parent_name, 
                        to_remove=to_rm, merges=merges, allClusters=allClusters, 
                     allVariances=allVariances)
}

# Compute the variance within a cluster.
#' @param modelInputs An object of the modelInputs class.
#' @param cluster A set of nodes.
#' @param clusterName The name assigned to a set of nodes.
#' @param allVariances A list of all variances investigated, where each variance 
#' refers to a single cluster.
computeVariance <- function(cluster, clusterName, modelInputs, allVariances){
  avg_variance <- -1
  # If the variance has already been computed, return it.
  if(clusterName %in% names(allVariances)){
    avg_variance <- allVariances[clusterName]
  }else{
    # Find the variance for each sample.
    variances <- unlist(lapply(colnames(modelInputs@node.wise.prediction), function(samp){
      vals <- unlist(modelInputs@node.wise.prediction[unlist(cluster), samp])
      mean_val <- rep(mean(vals), length(vals))
      sq_diff_vals <- (vals - mean_val)^2
      return(sum(sq_diff_vals)/length(vals))
    }))
    avg_variance <- stats::median(variances)
  }
  return(unname(unlist(avg_variance)))
}