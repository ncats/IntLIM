#' Run backpropagation for a single layer. In backpropagation, a gradient is
#' is computed by taking partial derivatives for each of the model weights. Given
#' a learning rate, the weights are adjusted according to the gradient.
#' @param modelResults An object of the ModelResults class.
#' @param iteration The current iteration.
#' @param convolution Whether or not to perform convolution
#' @param pooling Whether or not to perform pooling
#' @param ridgeRegressionWeight The hyperparameter weight assigned
#' to the ridge regression parameter (often referred to as lambda in the
#' literature)
#' @param varianceWeight The hyperparameter weight assigned to the difference
#' in variances between Y and the predicted value of Y.
BackpropagateSingleLayer <- function(modelResults, iteration, convolution, pooling,
                                     ridgeRegressionWeight, varianceWeight){
  # Calculate gradient.
  gradient <- computeGradientSingleLayer(modelResults, convolution, pooling,
                                         ridgeRegressionWeight, varianceWeight)
  
  # Update gradient.
  modelResults@current.gradient <- as.matrix(gradient)
  modelResults@iteration.tracking[iteration+1,
                                  which(grepl("Gradient", 
                                              colnames(modelResults@iteration.tracking)))] <- gradient
  
  # Set current weights, previous weights, and gradient.
  # Reference for the optimization algorithms: https://arxiv.org/abs/1609.04747
  modelResults@previous.weights <- modelResults@current.weights
  if(modelResults@optimization.type == "BGD"){
    # Batch Gradient Descent
    modelResults@current.weights <- modelResults@previous.weights - 
      modelResults@learning.rate * modelResults@current.gradient
  }else if(modelResults@optimization.type == "momentum"){
    # Momentum
    gamma <- 0.9
    update.vec <- gamma * modelResults@previous.update.vector + 
      modelResults@learning.rate * modelResults@current.gradient
    if(modelResults@current.iteration == 1){
      update.vec <- modelResults@learning.rate * modelResults@current.gradient
    }
    modelResults@current.weights <- modelResults@previous.weights - update.vec
    modelResults@previous.update.vector <- update.vec
  }else if(modelResults@optimization.type == "adagrad"){
    # Adagrad
    G <- modelResults@sum.square.gradients
    epsilon <- 0.00000001
    update.vec <- (modelResults@learning.rate / sqrt(G+epsilon)) * modelResults@current.gradient
    if(modelResults@current.iteration == 1){
      update.vec <- modelResults@learning.rate * modelResults@current.gradient
    }
    modelResults@current.weights <- modelResults@previous.weights - update.vec
    modelResults@sum.square.gradients <- modelResults@sum.square.gradients +
      (modelResults@current.gradient^2)
  }else if(modelResults@optimization.type == "adam"){
    # ADAM
    beta1 <- 0.9
    beta2 <- 0.999
    epsilon <- 0.00000001
    m <- beta1 * modelResults@previous.momentum + 
      (1-beta1) * modelResults@current.gradient
    v <- beta2 * modelResults@previous.update.vector +
      (1-beta2) * (modelResults@current.gradient^2)
    m.hat <- m / (1 - (beta1)^modelResults@current.iteration)
    v.hat <- v / (1 - (beta2)^modelResults@current.iteration)
    update.vec <- (modelResults@learning.rate * m.hat) / (sqrt(v.hat)+epsilon)
    modelResults@current.weights <- modelResults@previous.weights - update.vec
    modelResults@previous.momentum <- m
    modelResults@previous.update.vector <- v
  }else if(modelResults@optimization.type == "newton"){
    # Calculate Hessian and update weights.
    hessian <- computeHessianSingleLayer(modelResults, convolution, pooling)
    update <- MASS::ginv(hessian) %*% modelResults@current.gradient
    modelResults@current.weights <- modelResults@previous.weights - 
      (modelResults@learning.rate * update)
  }
  modelResults@iteration.tracking[iteration+1,
                                  which(grepl("Weight", 
                                              colnames(modelResults@iteration.tracking)))] <- modelResults@current.weights
  
  # Return.
  return(modelResults)
}

#' Compute the Hessian for a single layer neural network.
#' Note that this function is currently not implemented for categorical outcomes.
#' @param modelResults An object of the ModelResults class.
#' @param convolution Whether or not to perform convolution
#' @param pooling Whether or not to perform pooling
computeHessianSingleLayer <- function(modelResults, convolution, pooling){
  # Components for derivative.
  A.hat <- modelResults@model.input@A.hat
  X <- modelResults@model.input@node.wise.prediction
  Theta.old <- matrix(rep(modelResults@current.weights, dim(X)[2]), ncol = dim(X)[2])
  Y <- modelResults@model.input@true.phenotypes
  
  # Convolve X.
  conv <- lapply(1:dim(X)[2], function(i){
    ret_val <- X[,i]
    if(convolution == TRUE){
      ret_val <- A.hat %*% X[,i]
    }
    return(ret_val)
  })
  
  # Compute filter.
  S <- modelResults@pooling.filter@filter
  S_all <- modelResults@pooling.filter@individual.filters
  
  # Compute the convolution values.
  convolution.val <- lapply(1:length(conv), function(i){
    return(conv[[i]])
  })
  if(pooling == TRUE){
    convolution.val <- lapply(1:length(S_all), function(i){
      S.flat <- rowSums(S_all[[i]])
      ret_val <- conv[[i]] * S.flat
      if(modelResults@weights.after.pooling == TRUE){
        ret_val <- t(conv[[i]]) %*% S_all[[i]]
      }
      return(ret_val)
    })
  }
  
  # Create array, where each matrix is one sample. These will be summed
  # together to obtain the Hessian.
  hessian.single.samples <- lapply(1:dim(X)[2], function(i){
    mat.X <- matrix(rep(X[,i], dim(X)[1]), nrow = dim(X)[1], ncol = dim(X)[1])
    return(mat.X * t(mat.X))
  })
  
  # Sum these together to obtain the Hessian.
  hessian <- hessian.single.samples[[1]]
  if(length(hessian.single.samples)>1){
    for(i in 2:length(hessian.single.samples)){
      hessian <- hessian + hessian.single.samples[[i]]
    }
  }
  
  # Return Hessian.
  return(hessian)
}

#' Compute the gradient for a single layer neural network.
#' @param modelResults An object of the ModelResults class.
#' @param convolution Whether or not to perform convolution
#' @param pooling Whether or not to perform pooling
#' @param ridgeRegressionWeight The hyperparameter weight assigned
#' to the ridge regression parameter (often referred to as lambda in the
#' literature)
#' @param varianceWeight The hyperparameter weight assigned to the difference
#' in variances between Y and the predicted value of Y.
computeGradientSingleLayer <- function(modelResults, convolution, pooling,
                                       ridgeRegressionWeight, varianceWeight){
  # Components for derivative.
  A.hat <- modelResults@model.input@A.hat
  X <- modelResults@model.input@node.wise.prediction
  Theta.old <- matrix(rep(modelResults@current.weights, dim(X)[2]), ncol = dim(X)[2])
  Y <- modelResults@model.input@true.phenotypes

  # Convolve X.
  conv <- lapply(1:dim(X)[2], function(i){
    ret_val <- X[,i]
    if(convolution == TRUE){
      ret_val <- A.hat %*% X[,i]
    }
    return(ret_val)
  })

  # Compute filter.
  S <- modelResults@pooling.filter@filter
  S_all <- modelResults@pooling.filter@individual.filters
  Y.pred <- rep(0, dim(S)[2])
  if(pooling == TRUE){
    if(modelResults@weights.after.pooling == TRUE){
      Y.pred <- unlist(lapply(1:length(S_all), function(i){
        return(sum(t(conv[[i]]) %*% S_all[[i]] * Theta.old[,i]))
      }))
    }else{
      Y.pred <- unlist(lapply(1:length(S_all), function(i){
        return(sum(t(conv[[i]] * Theta.old[,i]) %*% S_all[[i]]))
      }))
    }
  }else{
    Y.pred <- unlist(lapply(1:length(conv), function(i){
      return(sum(t(conv[[i]] * Theta.old[,i])))
    }))
  }

  # Compute activation function.
  activation <- as.matrix(Y.pred)
  d.act.d.pred <- as.matrix(rep(1, length(Y.pred)))
  if(modelResults@model.input@outcome.type == "categorical"){
    if(modelResults@activation.type == "softmax"){
      activation <- SoftmaxWithCorrection(Y.pred)
      # Set maximum and minimum values to prevent infinite and infinitecimal values.
      exp.Y.pred <- exp(Y.pred)
      exp.Y.pred[which(exp.Y.pred > 10000000)] <- 100000000
      exp.Y.pred[which(exp.Y.pred < 0.00000001)] <- 0.000000001
      d.act.d.pred <- (exp.Y.pred * (exp.Y.pred - sum(exp.Y.pred)))/(sum(exp.Y.pred)^2)
    }else if(modelResults@activation.type == "tanh"){
      activation <- TanhWithCorrection(Y.pred)
      d.act.d.pred <- 2 * (1 - (tanh(2*(Y.pred-1.5)))^2)
    }else if(modelResults@activation.type == "sigmoid"){
      activation <- SigmoidWithCorrection(Y.pred)
      exp.Y.pred <- exp(1-Y.pred)
      exp.Y.pred[which(exp.Y.pred > 10000000)] <- 100000000
      exp.Y.pred[which(exp.Y.pred < 0.00000001)] <- 0.000000001
      d.act.d.pred <- -1 * exp.Y.pred/((1+exp.Y.pred)^2)
    }else{
      stop(paste("Invalid activation type", modelResults@activation.type))
    }
  }

  # Compute components of gradient based on activation function.
  diff <- Y - activation

  # Partial derivative of the loss (gradient). This is obtained via the chain
  # rule. 
  # dLoss/dTheta = -dLoss/dSumDiff * dSumDiff/dTheta 
  # = -Sum(dSquaredDiff/dDiff * dAct/dPred * dPred/dTheta) 
  d.pred.d.theta.list <- lapply(1:length(conv), function(i){
    return(conv[[i]])
  })
  if(pooling == TRUE){
    d.pred.d.theta.list <- lapply(1:length(S_all), function(i){
      S.flat <- rowSums(S_all[[i]])
      ret_val <- conv[[i]] * S.flat
      if(modelResults@weights.after.pooling == TRUE){
        ret_val <- t(conv[[i]]) %*% S_all[[i]]
      }
      return(ret_val)
    })
  }
  if(modelResults@weights.after.pooling == TRUE && pooling == TRUE){
    d.pred.d.theta <- do.call(rbind, d.pred.d.theta.list)
    const.terms <- t(matrix(rep(2 * diff * d.act.d.pred, dim(S)[2]), ncol = dim(S)[2]))
  }else{
    d.pred.d.theta <- do.call(cbind, d.pred.d.theta.list)
    const.terms <- t(matrix(rep(2 * diff * d.act.d.pred, dim(X)[1]), ncol = dim(X)[1]))
  }
  d.loss.d.theta.to.sum <- d.pred.d.theta * const.terms
  d.loss.d.theta <- -1 * rowSums(d.loss.d.theta.to.sum)
  gradient <- d.loss.d.theta

  # Add regularization term (ridge regression).
  ridge <- 2 * ridgeRegressionWeight * Theta.old[,1]
  
  # Add variance minimization term.
  # Note that we want to minimize (Var(Y)-Var(Y.pred))^2. The first derivative
  # is: 2(Var(Y)-Var(Y.pred))*(-d(Var(Y.pred*X))) = 
  # 2(Var(Y)-Var(Y.pred))*(-Mean(2(Y.pred-Mean(Y.pred))*(Y.pred.i-Mean(Y.pred.i))))
  # = 2(Var(Y)-Var(Y.pred))*(Mean(2(Mean(Y.pred.i)-Y.pred.i)*(X.i-Mean(X.i))))
  # Note that this term is only defined for batch training, not stochastic, as
  # we cannot compute a variance on a single training sample.
  variance.diff <- 2 * stats::var(Y) - stats::var(Y.pred)
  d.std.dev.Y.pred.list <- lapply(1:length(Y), function(i){
    Y.pred.stdev <- 2 * (mean(Y.pred) - Y.pred[i])
    Y.pred.stdev.vec <- rep(Y.pred.stdev, dim(d.pred.d.theta)[1])
    weightless.pred.stdev <- (d.pred.d.theta[,i] - rowMeans(d.pred.d.theta))/2
    return(Y.pred.stdev * weightless.pred.stdev)
  })
  d.std.dev.y.pred <- do.call(rbind, d.std.dev.Y.pred.list)
  mean.std.dev.Y.pred <- colMeans(d.std.dev.y.pred)
  variance <- varianceWeight * variance.diff * mean.std.dev.Y.pred

  # Return gradient.
  gradient <- gradient + ridge + variance
  return(gradient)
}

#' Apply the Softmax function to a numeric vector input. Note that we expect inputs
#' to be between 1 and 2, not 0 and 1. We therefore correct accordingly by subtracting
#' 1 from the input and adding 1 to the output.
#' @param input A vectorized numeric input.
SoftmaxWithCorrection <- function(input){
  numerator <- exp(input - 1)
  numerator[which(numerator > 10000000)] <- 100000000
  numerator[which(numerator < 0.00000001)] <- 0.000000001
  denominator <- sum(numerator)
  output <- (numerator/denominator)+1
  return(output)
}

#' Apply the Tanh function to a numeric vector input. Note that we expect inputs
#' to be between 1 and 2, not -1 and 1. We therefore correct accordingly.
#' @param input A vectorized numeric input.
TanhWithCorrection <- function(input){
  return(1.5 + tanh(2 * (input - 1.5)) / 2)
}

#' Apply the Sigmoid function to a numeric vector input. Note that we expect inputs
#' to be between 1 and 2, not -1 and 1. We therefore correct accordingly.
#' @param input A vectorized numeric input.
SigmoidWithCorrection <- function(input){
  exp.input <- exp(1-input)
  exp.input[which(exp.input > 10000000)] <- 100000000
  exp.input[which(exp.input < 0.00000001)] <- 0.000000001
  return(1 + (1 / (1 + exp.input)))
}

#' Adjust the filter such that it reflects the pooling operation performed.
#' @param poolingFilter A PoolingFilter object.
#' @param A.hat Convolutional operator.
#' @param X Input matrix from last training iteration.
#' @param Theta Weight matrix from last training iteration.
AdjustFilter <- function(poolingFilter, A.hat, X, Theta){
  S <- poolingFilter@filter
  S_all <- lapply(1:dim(X)[2], function(i){
    return(S)
  })

  # For mean filters, divide each column by the count of features and
  # replicate for each sample.
  if(poolingFilter@filter.type == "mean"){
    S_all <- lapply(1:dim(X)[2], function(i){
      return(S / colSums(S))
    })
  # For min, max, or median filters, select the feature meeting the criteria
  # for each sample and return.
  }else{
    S_all <- lapply(1:dim(X)[2], function(i){
      to_mult <- t(A.hat %*% X[,i] * Theta[,i])
      S.copy <- matrix(0, nrow = dim(S)[1], ncol = dim(S)[2])
      for(i in 1:dim(S)[2]){
        to_mult_filt <- to_mult[which(S[,i]!=0)]
        sel_samp <- 0
        if(poolingFilter@filter.type == "min"){
          sel_samp <- which.min(to_mult_filt)
        }else if(poolingFilter@filter.type == "max"){
          sel_samp <- which.max(to_mult_filt)
        }else if(poolingFilter@filter.type == "median"){
          if (length(to_mult_filt) %% 2 != 0) {
            sel_samp <- which(to_mult_filt == stats::median(to_mult_filt))
          } else if (length(to_mult_filt) %% 2 == 0) {
            a = sort(to_mult_filt)[c(length(to_mult_filt)/2, length(to_mult_filt)/2+1)]
            sel_samp <- c(which(to_mult_filt == a[1]), which(to_mult_filt == a[2]))
          }
          sel_samp <- sel_samp[1]
        }
        S.copy[sel_samp, i] <- 1
      }
      return(S.copy)
    })
  }
  return(S_all)
}

#' Create the initial filter to reflect the pooling operation. This function is
#' for scenarios where the weights are applied after pooling. Note that min 
#' and max pooling do not make sense in this scenario.
#' @param poolingFilter A PoolingFilter object.
#' @param A.hat Convolutional operator.
#' @param X Input matrix from last training iteration.
CreateFilter <- function(poolingFilter, A.hat, X){
  S <- poolingFilter@filter
  S_all <- lapply(1:dim(X)[2], function(i){
    return(S)
  })
  
  # For mean filters, divide each column by the count of features and
  # replicate for each sample.
  if(poolingFilter@filter.type == "mean"){
    S_all <- lapply(1:dim(X)[2], function(i){
      return(S / colSums(S))
    })
    # For min, max, or median filters, select the feature meeting the criteria
    # for each sample and return.
  }else if(poolingFilter@filter.type == "median"){
    S_all <- lapply(1:dim(X)[2], function(i){
      to_mult <- t(A.hat %*% X[,i])
      S.copy <- matrix(0, nrow = dim(S)[1], ncol = dim(S)[2])
      for(i in 1:dim(S)[2]){
        to_mult_filt <- to_mult[which(S[,i]!=0)]
        
        # Get median.
        if (length(to_mult_filt) %% 2 != 0) {
          sel_samp <- which(to_mult_filt == stats::median(to_mult_filt))
        } else if (length(to_mult_filt) %% 2 == 0) {
          a = sort(to_mult_filt)[c(length(to_mult_filt)/2, length(to_mult_filt)/2+1)]
          sel_samp <- c(which(to_mult_filt == a[1]), which(to_mult_filt == a[2]))
          sel_samp <- sel_samp[1]
        }
        S.copy[sel_samp, i] <- 1
      }
      return(S.copy)
    })
  }else{
    stop(paste("Pooling type", poolingFilter@filter.type, "is not defined when weights
         are applied after pooling."))
  }
  return(S_all)
}
