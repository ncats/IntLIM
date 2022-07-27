#' @name RemovePlusInCovars
#' @title RemovePlusInCovars
# If there is a '+' in the covariate names, this will cause a problem
# because the covariates will be collapsed using '+'. For example, if
# covariate names are 'age', 'sex', and 'stage3+4', the formula will
# become 'f ~ a + type + a:type + age + sex + stage3 + 4'. We fix
# this problem by changing all '+' to 'plus'.
#' @param covar vector of additional vectors to consider
#' @param sampleDataColnames vector of column names, which is a superset
#' of the covar vector.
#' @return list containing two elements: new covariates and new column names.
RemovePlusInCovars <- function(covar=c(), sampleDataColnames){
  # Find which covariates have plus signs.
  which_plus <- which(grepl("+", covar, fixed = TRUE) == TRUE)
  oldCovars <- covar
  
  # Replace each plus sign with "plus".
  covar <- unlist(lapply(1:length(covar), function(i){
    retval <- covar[i]
    if(i %in% which_plus){
      splitOnPlus <- strsplit(covar[i], "+", fixed = TRUE)
      retval <- paste(splitOnPlus[[1]], collapse = "plus")
    }
    return(retval)
  }))
  
  # Replace the plus signs in the sampleMetaData column names as well.
  sampleDataColnames <- unlist(lapply(1:length(sampleDataColnames),
                                                     function(i){
                                                       retval <- sampleDataColnames[i]
                                                       if(retval %in% oldCovars){
                                                         whichMatch <- which(oldCovars == retval)
                                                         retval <- covar[whichMatch]
                                                       }
                                                       return(retval)
                                                     }))
  # Return new values.
  return(list(covar = covar, sampleDataColnames = sampleDataColnames))
}

#' Function that runs linear models and returns interaction p-values.
#'
#' @include AllClasses.R
#'
#' @param incommon Named list (output of 
#' FilterData()) with analyte levels, 
#' and associated meta-data
#' @param outcome '1' or '2' must be set as outcome/independent variable
#' (default is '1')
#' @param independentVariable '1' or '2' must be set as outcome/independent variable
#' @param type vector of sample type (by default, it will be used in the interaction term).
#' Only 2 categories are currently supported.
#' @param covar vector of additional vectors to consider
#' @param continuous boolean to indicate whether the data is continuous or discrete
#' @param save.covar.pvals boolean to indicate whether or not to save the p-values of all covariates,
#' which can be analyzed later but will also lengthen computation time.
#' (rather than interaction terms).
#' @param keep.highest.pval boolean to indicate whether or not to remove the 
#' pair with the highest p-value across two duplicate models (e.g. m1~m2 and m2~m1)
#' @param suppressWarnings whether or not to suppress warnings.
RunLM <- function(incommon, outcome=1, independentVariable = 2, type="", covar=c(), 
                  continuous=FALSE, save.covar.pvals = FALSE, keep.highest.pval = FALSE,
                  suppressWarnings = FALSE) {
  
  # Initialize types 1 and 2 and warning list.
  type1 <- incommon@analyteType1
  type2 <- incommon@analyteType2
  
  # Initialize covariate matrix, type, and messages.
  covarMatrix <- as.matrix(incommon@sampleMetaData[,covar])
  stype <- type
  mymessage<-list()
  if(length(covar)>0){
    # If one of the covariate names contains a plus sign, change it because this
    # will cause problems when constructing the formula.
    adjNames <- RemovePlusInCovars(covar, colnames(incommon@sampleMetaData))
    covar <- adjNames$covar
    colnames(incommon@sampleMetaData) <- adjNames$sampleDataColnames
    
    # Convert covariates to matrix. Ensure that matrix is one-hot encoded.
    f <- paste('~ 0 + ', paste(covar, collapse = ' + '))
    dat <- incommon@sampleMetaData[,covar]
    if(length(covar) == 1){
      dat <- data.frame(V1 = incommon@sampleMetaData[,covar])
      colnames(dat) <- covar[1]
    }
    covarMatrix <- stats::model.matrix(stats::as.formula(f), data = dat)
    covar <- colnames(covarMatrix)
    
    # Since names will be changed now, we need to remove plus signs again. For example,
    # if we have a variable 'treatment' that has values 'med1', 'med2', and 'med1+med2',
    # the column names will now include 'treatmentmed1' and 'treatmentmed1+med2'.
    adjNames <- RemovePlusInCovars(covar, colnames(covarMatrix))
    covar <- adjNames$covar
    colnames(covarMatrix) <- adjNames$sampleDataColnames
  }

  # Find all standard deviations.
  type1sd <- as.numeric(apply(type1,1,function(x){stats::sd(as.numeric(x),na.rm=TRUE)}))
  type2sd <- as.numeric(apply(type2,1,function(x){stats::sd(as.numeric(x),na.rm=TRUE)}))
  covarsd <- as.numeric(apply(covarMatrix,2,function(x){
    return(stats::sd(x,na.rm=TRUE))}))
  if(methods::is(stype, "character")){
    stype <- as.numeric(as.factor(stype))
  }else if(methods::is(stype, "factor")){
    stype <- as.numeric(stype)
  }
  stypesd <- stats::sd(stype,na.rm=TRUE)

  # If the standard deviation of the phenotype is zero, then stop.
  if(stypesd == 0){
    stop("stype variable has a standard deviation of zero. Cannot run.")
  }
  # If the standard deviation of analyte type 1 is zero, then remove and add a warning.
  if(length(which(type1sd==0))>0) {
    toremove <- c(which(type1sd==0))
    namestoremove <- rownames(type1)[toremove]
    type1 <- type1[-toremove,]
    mymessage[[length(mymessage)+1]] <- paste("Removed",length(toremove),"analytes of",
                                              "type 1 that had",
                                              "a standard deviation of 0:",
                                              namestoremove)
  }
  # If the standard deviation of analyte type 2 is zero, then remove and add a warning.
  if(length(which(type2sd==0))>0) {
    toremove <- c(which(type2sd==0))
    namestoremove <- rownames(type2)[toremove]
    type2 <- type2[-toremove,]
    mymessage[[length(mymessage)+1]] <-paste("Removed",length(toremove),"analytes of",
                                             "type 2 that had",
                                             "a standard deviation of 0:",
                                             namestoremove)
  }
  # If the standard deviation of a covariate is zero, then remove and add a warning.
  if(length(which(covarsd==0))>0) {
    toremove <- c(which(covarsd==0))
    namestoremove <- colnames(covarMatrix)[toremove]
    namestokeep <- colnames(covarMatrix)[-toremove]
    rows <- rownames(covarMatrix)
    covarMatrix <- as.matrix(covarMatrix[,-toremove])
    colnames(covarMatrix) <- namestokeep
    rownames(covarMatrix) <- rows
    covar <- namestokeep
    mymessage[[length(mymessage)+1]] <-paste("Removed",length(toremove),
                                             "covariates that had",
                                             "a standard deviation of 0:",
                                             namestoremove)
  }
  mat.list <- getStatsAllLM(outcome = outcome, independentVariable = independentVariable,
                            type1 = type1, type2 = type2, type =
                              type, covar = covar, covarMatrix = covarMatrix,
                            continuous = continuous, save.covar.pvals = save.covar.pvals,
                            remove.tri = keep.highest.pval, suppressWarnings = suppressWarnings)
  if(length(mat.list[["warnings"]])>0){
    for(i in 1:length(mat.list[["warnings"]])){
      mymessage[[length(mymessage)+1]] <- mat.list[["warnings"]][[i]]
    }
  }
  mat.list <- mat.list[["list"]]

  # Add covariates.
  intLimResultsCovar <- ""
  if(!is.null(covar)){
    intLimResultsCovar <- covar
  }

  # Create object to return.
  myres <- methods::new('IntLimResults',
                        interaction.pvalues=mat.list$mat.pvals,
                        interaction.adj.pvalues = mat.list$mat.pvalsadj,
                        interaction.coefficients=mat.list$mat.coefficients,
                        model.rsquared = mat.list$mat.rsquared,
                        covariate.pvalues = mat.list$covariate.pvals,
                        covariate.coefficients = mat.list$covariate.coefficients,
                        warnings=mymessage, covar = intLimResultsCovar)

  return(myres)
}

#' Function that runs linear models for analyte vs. all analytes of the other type
#'
#' @include AllClasses.R
#'
#' @param form LM formulat (typically m~g+t+g:t)
#' @param clindata data frame with 1st column: expression of one analyte; 2nd column
#' sample type (e.g. cancer/non-cancer)
#' @param arraydata matrix of analyte values
#' @param analytename name of independent analyte in the model
#' @param suppressWarnings whether or not to suppress warnings
getstatsOneLM <- function(form, clindata, arraydata, analytename, suppressWarnings = FALSE) {
  #array data is one analyte type
  #clindata is the other analyte type
  call=match.call()
  YY <- t(arraydata)                      # the data matrix
  #mean of analytes across all samples
  EY <- apply(YY, 2, mean)                # its mean vector
  #sum of squares after centering
  SYY <- apply(YY, 2, function(y) {sum(y^2)}) - nrow(YY)*EY^2     # sum of squares after centering
  # The first "Y" column is added as a "dummy" Y-value. This is needed
  # for stats::model.matrix to generate a design matrix.
  clindata <- data.frame(y=YY[,1], clindata)
  dimnames(clindata)[[2]][1] <- 'Y'
  # design matrix
  X <- stats::model.matrix(form, clindata)
  N = dim(X)[1]
  p <- dim(X)[2]
  # Create a contrast matrix.
  XtX <- t(X) %*% X
  ixtx <- MASS::ginv(XtX)
  # Use the pseudoinverse if the inverse cannot be found.
  # Print out correlated covariates in this case.
  
  # Initialize warnings. We will later remove if needed.
  pinv_message <- paste("Using pseudoinverse for", analytename)
  cutoff = 0.9
  covariate_msg1 <- paste("The following covariates have correlation >", cutoff, ":")
  covariate_msg2 <- paste("The following covariates have correlation <", -1 * cutoff, ":")
  cormat <- stats::cor(X[,which(colnames(X) != "(Intercept)")])
  cormat[lower.tri(cormat, diag = TRUE)] <- 0
  if(length(which(cormat > cutoff)) > 0){
    which_greater <- multi.which(cormat > cutoff)
    for(i in 1:nrow(which_greater)){
      if(i == 1){
        covariate_msg1 <- paste(covariate_msg1, paste0("(",colnames(cormat)[which_greater[i,1]],
                                                       ", ", 
                                                       colnames(cormat)[which_greater[i,2]], 
                                                       ")"))
      }else{
        covariate_msg1 <- paste(covariate_msg1, paste0("(",colnames(cormat)[which_greater[i,1]],
                                                       ", ", 
                                                       colnames(cormat)[which_greater[i,2]], 
                                                       ")"), sep = ", ")
      }
    }
  }
  if(length(which(cormat < -1 * cutoff)) > 0){
    which_less <- multi.which(cormat < -1 * cutoff)
    for(i in 1:nrow(which_less)){
      if(i == 1){
        covariate_msg2 <- paste(covariate_msg2, paste0("(",colnames(cormat)[which_less[i,1]],
                                                       ", ", 
                                                       colnames(cormat)[which_less[i,2]], 
                                                       ")"))
      }else{
        covariate_msg2 <- paste(covariate_msg2, paste0("(",colnames(cormat)[which_less[i,1]],
                                                       ", ", 
                                                       colnames(cormat)[which_less[i,2]], 
                                                       ")"), sep = ", ")
      }
    }
  }
  warnings <- c(pinv_message, covariate_msg1, covariate_msg2)
  
  tryCatch({
    ixtx <- solve(XtX)
    warnings <- list()
  }, error=function(e){
    if(suppressWarnings == FALSE){
      warning(pinv_message)
    }
    if(length(which(cormat > cutoff)) > 0 && suppressWarnings == FALSE){
      warning(covariate_msg1)
    }
    if(length(which(cormat < -1 * cutoff)) > 0 && suppressWarnings == FALSE){
      warning(covariate_msg2)
    }
  })
  bhat <- NULL
  
  bhat <- ixtx %*% t(X) %*% YY            # Use the pseudo-inverse to estimate the parameters
  yhat <- X %*% bhat                      # Figure out what is predicted by the model
  # Now we partition the sum-of-square errors
  rdf <- ncol(X)-1                        # number of parameters in the model
  edf <- nrow(YY)-rdf-1                   # additional degrees of freedom
  errors <- YY - yhat                     # difference between observed and model predictions
  sse <- apply(errors^2, 2, sum)  # sum of squared errors over the samples
  mse <- sse/edf                  # mean squared error
  ssr <- SYY - sse                        # regression error
  msr <- ssr/rdf                  # mean regression error
  fval <- msr/mse                 # f-test for the overall regression
  pfval <- 1-stats::pf(fval, rdf, edf)           # f-test p-values
  
  stderror.coeff <- sapply(mse,function(x){sqrt(diag(ixtx)*x)})
  t.coeff <- bhat/stderror.coeff
  p.val.coeff <- 2*stats::pt(-abs(t.coeff),df = (N-p))
  y.dev <- lapply(1:(dim(YY)[2]), function(i){
    return(YY[,i]-EY[i])
  })
  var.y <- unlist(lapply(y.dev, function(i){
    return(sum(i^2))
  }))
  r.squared <- 1 - (sse / var.y)
  rownames(bhat) <- colnames(X)
  rownames(p.val.coeff) <- colnames(X)
  return(list("mlin" = list(coefficients=bhat,
              p.value.coeff = p.val.coeff, # interaction p-value
              r.squared.val = r.squared),# r-squared value
              "warnings" = warnings))
}

#' Function that runs Linear Models for all analytes
#' @include AllClasses.R
#' @param outcome '1' or '2' must be set as outcome/independent variable
#' @param independentVariable '1' or '2' must be set as outcome/independent variable
#' @param type1 Analyte type 1 dataset
#' @param type2 Analyte type 2 dataset
#' @param type vector of sample type (by default, it will be used in the interaction term).
#' Only 2 categories are currently supported.
#' @param covar vector of additional vectors to consider
#' @param covarMatrix covariate matrix
#' @param continuous indicate whether data is discrete (FALSE) or continuous (TRUE)
#' @param save.covar.pvals boolean to indicate whether or not to save the p-values of all covariates,
#' which can be analyzed later but will also lengthen computation time.
#' @param remove.tri boolean to indicate whether or not to remove the 1-1
#' or 2-2 pair with the highest p-value across two duplicate models (e.g. m1~m2 and m2~m1)
#' @param suppressWarnings whether or not to suppress warnings
#' @return list of matrices (interaction.pvalues, interaction.adj.pvalues, interaction.coefficients)
getStatsAllLM <- function(outcome, independentVariable, type1, type2, type, covar, covarMatrix, 
                          continuous, save.covar.pvals, remove.tri = FALSE,
                          suppressWarnings = FALSE) {
  outcomeArrayData <- NULL
  independentArrayData <- NULL
  num <- NULL
  
  # Get array data based on outcome.
  if (outcome==1) {
    outcomeArrayData <- data.frame(type1)
  } else if (outcome == 2){
    outcomeArrayData <- data.frame(type2)
  }
  
  # Get array data based on independent variable type.
  if (independentVariable==1) {
    independentArrayData <- data.frame(type1)
    num <- nrow(type1)
  } else if (independentVariable == 2){
    num <- nrow(type2)
    independentArrayData <- data.frame(type2)
  }

  # Set up formula and interaction term.
  form.add <- "Y ~ a + type + a:type"
  interactionTerm <- "a:type"
  
  # Set numprog.
  numprog <- round(num*0.1)
  
  # Add covariates to the formula.
  if (length(covar)>0) {
    len.covar <- length(covar)
    for (i in 1:len.covar) {
      form.add <- paste(form.add, '+', covar[i])
    }
  }
  
  # Initialize stats to collect.
  list.pvals <- list()
  list.coefficients <- list()
  list.rsquared <- list()
  list.covariate.pvals <- list()
  list.covariate.coefficients <- list()
  
  # Run each model.
  warnings <- list()
  for (i in 1:num) {
    
    # Set up clinical data.
    a <- as.numeric(independentArrayData[i, ])
    if (is.null(covar)) {
      clindata <- data.frame(a, type)
    } else {
      clindata <- data.frame(a, type, covarMatrix)
      colnames(clindata)[3:ncol(clindata)] <- covar
    }
    
    # Change type for continuous data (factor to numeric)
    if(continuous){
      clindata[2] <- lapply(clindata[2], as.character)
      clindata[2] <- lapply(clindata[2], as.numeric)
    }

    # Run all models for this independent analyte.
    mlin <- getstatsOneLM(stats::as.formula(form.add), clindata = clindata,
                          arraydata = outcomeArrayData, 
                          analytename = rownames(independentArrayData)[i],
                          suppressWarnings = suppressWarnings)
    warnings <- c(warnings, mlin[["warnings"]])
    mlin <- mlin[["mlin"]]
    term.pvals <- rownames(mlin$p.value.coeff)
    
    # Return the primary p-values and coefficients.
    index.interac <- grep(interactionTerm, term.pvals)
    term.coefficient <- rownames(mlin$coefficients)
    index.coefficient <-  grep(interactionTerm, term.coefficient)
    p.val.vector <- as.vector(mlin$p.value.coeff[index.interac,])
    coefficient.vector <- as.vector(mlin$coefficients[index.coefficient,])
    
    term.rsquared <- rownames(mlin$r.squared.val)
    r.squared.vector <- as.vector(mlin$r.squared.val)
    
    if (numprog != 0){
      if (i %% numprog == 0) {
        progX <- round(i/num*100)
        message(paste(progX,"% complete"))
      }
    }
    list.pvals[[i]] <-  p.val.vector
    list.coefficients[[i]] <- coefficient.vector
    list.rsquared[[i]] <- r.squared.vector
    
    if(save.covar.pvals == TRUE){
      # Save covariate p-values.
      covariate.pvals <- lapply(term.pvals, function(covariate){
        return(mlin$p.value.coeff[covariate,])
      })
      covariate.pvals.df <- do.call("cbind", covariate.pvals)
      rownames(covariate.pvals.df) <- paste(rownames(independentArrayData)[i], 
                                            rownames(covariate.pvals.df),
                                            sep="__")
      colnames(covariate.pvals.df) <- term.pvals
      list.covariate.pvals[[i]] <- covariate.pvals.df
      
      # Save covariate coefficients.
      covariate.coefficients <- lapply(term.coefficient, function(covariate){
        return(mlin$coefficients[covariate,])
      })
      covariate.coefficients.df <- do.call("cbind", covariate.coefficients)
      rownames(covariate.coefficients.df) <- paste(rownames(independentArrayData)[i],
                                                   rownames(covariate.coefficients.df), 
                                                   sep="__")
      colnames(covariate.coefficients.df) <- term.pvals
      list.covariate.coefficients[[i]] <- covariate.coefficients.df
    }
  }
  
  # Convert the stats into matrix form.
  mat.pvals <- do.call(rbind, list.pvals)
  mat.coefficients <- do.call(rbind, list.coefficients)
  mat.rsquared <- do.call(rbind, list.rsquared)
  covariate.pvals <- do.call(rbind, list.covariate.pvals)
  covariate.coefficients <- do.call(rbind, list.covariate.coefficients)
  
  # Adjust p-values.
  row.pvt <- dim(mat.pvals)[1]
  col.pvt <- dim(mat.pvals)[2]
  myps <- as.vector(mat.pvals)
  mypsadj <- stats::p.adjust(myps, method = 'fdr')
  mat.pvalsadj <- matrix(mypsadj, row.pvt, col.pvt)
  
  # Assign names to results.
  if (outcome==1) {
    colnames(mat.pvals) <- colnames(mat.pvalsadj) <- rownames(type1)
    colnames(mat.coefficients) <- rownames(type1)
    colnames(mat.rsquared) <- rownames(type1)
  } else if (outcome==2) {
    colnames(mat.pvals) <- colnames(mat.pvalsadj) <- rownames(type2)
    colnames(mat.coefficients) <- rownames(type2)
    colnames(mat.rsquared) <- rownames(type2)
  }
  if(independentVariable == 1){
    rownames(mat.pvals) <- rownames(mat.pvalsadj) <- rownames(type1)
    rownames(mat.coefficients) <- rownames(type1)
    rownames(mat.rsquared) <- rownames(type1)
  }else if (independentVariable == 2){
    rownames(mat.pvals) <- rownames(mat.pvalsadj) <- rownames(type2)
    rownames(mat.coefficients) <- rownames(type2)
    rownames(mat.rsquared) <- rownames(type2)
  }
  
  # Remove the triangular matrix if applicable.
  if(remove.tri == TRUE && outcome != independentVariable){
    warning("Cannot remove triangular matrix if the independent and outcome variables are",
            "different analyte types. Not removing.")
    remove.tri = FALSE
  }
  if(remove.tri == TRUE)
  {
    # Find locations where the upper triangular value is higher than the lower triangular.
    mat.pvals.t <- t(mat.pvals)
    mat.pvalsadj.t <- t(mat.pvalsadj)
    mat.coefficients.t <- t(mat.coefficients)
    mat.rsquared.t <- t(mat.rsquared)
    where_upper_triangular_higher <- which(mat.pvals > t(mat.pvals))
    
    # Replace those locations with values from the lower triangular.
    mat.pvals[where_upper_triangular_higher] <- mat.pvals.t[where_upper_triangular_higher]
    mat.pvalsadj[where_upper_triangular_higher] <- mat.pvalsadj.t[where_upper_triangular_higher]
    mat.coefficients[where_upper_triangular_higher] <- mat.coefficients.t[where_upper_triangular_higher]
    mat.rsquared[where_upper_triangular_higher] <- mat.rsquared.t[where_upper_triangular_higher]
    covariate.coefficients <- do.call(rbind, list.covariate.coefficients)
    
    # Remove the highest p-values and store the result in the upper triangular if required.
    # Otherwise, remove the diagonal if applicable.
    should_remove = unlist(lapply(rownames(covariate.pvals), function(name){
      ret_val = FALSE
      pieces = strsplit(name, "__")[[1]]
      a_type = colnames(covariate.pvals)[which(grepl("a:", colnames(covariate.pvals), 
                                                     fixed = TRUE) == TRUE)]
      pval_1 = covariate.pvals[name, a_type]
      pval_2 = covariate.pvals[paste(pieces[2], pieces[1], sep = "__"),a_type]
      if(pieces[1] == pieces[2]){
        ret_val = TRUE
      }
      else if(pval_2 < pval_1){
        ret_val = TRUE
      }
      else if(pval_2 == pval_1){
        pos_pval_1 = which(rownames(covariate.pvals)[1] == name)
        pos_pval_2 = which(rownames(covariate.pvals)[1] == paste(pieces[2], 
                                                                 pieces[1], 
                                                                 sep = "__"))
        if(pos_pval_1 > pos_pval_2){
          ret_val = TRUE
        }
      }
      return(ret_val)
    }))
    
    if(save.covar.pvals == TRUE){
      for(i in 1:length(colnames(covariate.pvals))){
        covariate.pvals[which(should_remove == TRUE),i] <- NA
        covariate.coefficients[which(should_remove == TRUE),i] <- NA
      }
    }
    mat.pvals[lower.tri(mat.pvals,diag=TRUE)] <- NA
    mat.pvalsadj[lower.tri(mat.pvalsadj,diag=TRUE)] <- NA
    mat.coefficients[lower.tri(mat.coefficients,diag=TRUE)] <- NA
    mat.rsquared[lower.tri(mat.rsquared,diag=TRUE)] <- NA
  }
  else if (outcome == independentVariable){
    # Only remove the diagonal.
    diag(mat.pvals) <- NA
    diag(mat.pvalsadj) <- NA
    diag(mat.coefficients) <- NA
    diag(mat.rsquared) <- NA
    if(save.covar.pvals == TRUE){
      pieces_list = lapply(rownames(covariate.pvals), function(name){
        split = strsplit(name, "__")[[1]]
        return(data.frame(from = split[1], to = split[2]))
      })
      pieces = do.call(rbind, pieces_list)
      for(i in 1:length(colnames(covariate.pvals))){
        covariate.pvals[which(pieces$from == pieces$to),i] <- NA
        covariate.coefficients[which(pieces$from == pieces$to),i] <- NA
      }
    }
  }
  
  # Add the matrices to a final list.
  list.mat <- list()
  list.final <- list()
  list.mat[["mat.pvals"]] <- as.matrix(mat.pvals)
  list.mat[["mat.pvalsadj"]] <- as.matrix(mat.pvalsadj)
  list.mat[["mat.coefficients"]] <- as.matrix(mat.coefficients)
  list.mat[["mat.rsquared"]] <- as.matrix(mat.rsquared)
  list.mat[["covariate.pvals"]] <- as.data.frame(covariate.pvals)
  list.mat[["covariate.coefficients"]] <- as.data.frame(covariate.coefficients)
  list.final[["list"]] <- list.mat
  list.final[["warnings"]] <- warnings
  return(list.final)
}

#' Function that gets numeric cutoffs from percentile
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient (default bottom 10 percent (high negative coefficients) and top 10 percent (high positive coefficients))
#' @param tofilter dataframe for percentile filtering
#' @return vector with numeric cutoffs
getQuantileForInteractionCoefficient<-function(tofilter, interactionCoeffPercentile){
  
  if(interactionCoeffPercentile > 1 || interactionCoeffPercentile < 0) {
    stop("interactionCoeffPercentile parameter must be between 0 and 1")
  }
  
  #get top and bottom cutoffs (need highest positive and highest negative coeffs)
  tofilter_abs = abs(tofilter)
  abs_cutoff = as.numeric(stats::quantile(tofilter_abs, interactionCoeffPercentile, na.rm = TRUE))
  
  return(c(0 - abs_cutoff, abs_cutoff))
  
}

#' A which for multidimensional arrays.
#' Mark van der Loo 16.09.2011
#' 
#' @name multi.which
#' @param A Boolean function defined over a matrix
#' @return vector with numeric cutoffs
multi.which <- function(A){
  if ( is.vector(A) ) return(which(A))
  d <- dim(A)
  T.mat <- which(A) - 1
  nd <- length(d)
  t( sapply(T.mat, function(t){
    I <- integer(nd)
    I[1] <- t %% d[1]
    sapply(2:nd, function(j){
      I[j] <<- (t %/% prod(d[1:(j-1)])) %% d[j]
    })
    I
  }) + 1 )
 }