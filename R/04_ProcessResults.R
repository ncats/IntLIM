#' Retrieve significant pairs, based on adjusted p-values, interaction coefficient
#' percentile, and R^2 values.
#'
#' @include internalfunctions.R
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData IntLimData object (output of ReadData()) with analyte levels
#'  and associated meta-data
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param coeffPercentile percentile cutoff for absolute value of coefficient
#' @param rsquaredCutoff cutoff for lowest r-squared value
#' @param coefficient Coefficient to filter by. Default is interaction coefficient.
#' Other options are "stype" (phenotype) and "analyte" (independent analyte).
#' @return IntLimResults object with model results.
#' @export
ProcessResults <- function(inputResults,
				inputData,
				pvalcutoff=0.05,
				coeffPercentile=0,
				rsquaredCutoff = 0,
				coefficient = "interaction"){
  
  if(coefficient == "interaction" || 
     (coefficient == "stype" && length(inputResults@covariate.pvalues) > 0) ||
     (coefficient == "analyte" && length(inputResults@covariate.pvalues) > 0)){
    # Check input types.
    if(!methods::is(inputResults, "IntLimResults")){
      stop("Results must be an IntLIMResults object")
    }
    if(!methods::is(inputData, "IntLimData")){
      stop("Data must be an IntLIMData object")
    }
    
    # Check parameter ranges.
    if(pvalcutoff < 0 || pvalcutoff > 1){
      stop("P-value must be between 0 and 1")
    }
    if(rsquaredCutoff < 0 || rsquaredCutoff > 1){
      stop("R-squared value must be between 0 and 1")
    }
    if(coeffPercentile < 0 || coeffPercentile > 1){
      stop("Interaction coefficient percentile must be between 0 and 1")
    }
    
    # Store all results in shorter variables.
    mydat <-inputResults@interaction.adj.pvalues
    mydat.interac <- inputResults@interaction.coefficients
    mydat.rsq <- inputResults@model.rsquared 
    mydat.covar.coeff <- inputResults@covariate.coefficients
    
    # Set input and output type.
    indvarname <- "IndependentVariable"
    outcomename <- "Outcome"
    indtype <- NULL
    outcometype <- NULL
    if(inputResults@independent.var.type == 1 && inputResults@outcome == 2){
      indtype <- inputData@analyteType1
      outcometype <- inputData@analyteType2
    }else if(inputResults@independent.var.type == 2 && inputResults@outcome == 1){
      indtype <- inputData@analyteType2
      outcometype <- inputData@analyteType1
    }else if(inputResults@independent.var.type == 2 && inputResults@outcome == 2){
      indtype <- inputData@analyteType2
      outcometype <- inputData@analyteType2
    }else if(inputResults@independent.var.type == 1 && inputResults@outcome == 1){
      indtype <- inputData@analyteType1
      outcometype <- inputData@analyteType1
    }else{
      stop(paste("Independent variable and outcome must both",
                 "be either 1 or 2."))
    }
    
    # Throw error if type is not present in original data.
    if(length(indtype) == 0){
      stop("Independent data type is not present in original data")
    }
    if(length(outcometype) == 0){
      stop("Outcome type is not present in original data")
    }
    
    p <- inputData@sampleMetaData[,inputResults@stype]
    
    # Call continuous function if applicable.
    if((length(unique(p)) !=2) && (inputResults@continuous == 1)){
      filtResults <- ProcessResultsContinuous(inputResults,
                                              coeffPercentile,
                                              pvalcutoff, rsquaredCutoff,
                                              coefficient)
    } else if (length(unique(p)) !=2){
      stop(paste(
        "IntLim requires two categories only for correlation analysis. Make sure the column",
        inputResults@stype, "only has two unique values or is continuous"))
    }else{
      gp1 <- which(p == unique(p)[1])
      gp2 <- which(p == unique(p)[2])
      indgp1 <- t(indtype[rownames(mydat),gp1])
      outgp1 <- t(outcometype[colnames(mydat),gp1])
      indgp2 <- t(indtype[rownames(mydat),gp2])
      outgp2 <- t(outcometype[colnames(mydat),gp2])
      
      # Create melted matrices to filter by inputs.
      whichTypeCol <- which(grepl("^type", colnames(inputResults@covariate.adj.pvalues)) == TRUE)
      finmydat <- reshape2::melt(mydat)
      if(coefficient == "analyte"){
        finmydat <- data.frame(Var1 = finmydat$Var1,
                               Var2 = finmydat$Var2,
                               value = inputResults@covariate.adj.pvalues$a)
      }else if(coefficient == "stype"){
        finmydat <- data.frame(Var1 = finmydat$Var1,
                                  Var2 = finmydat$Var2,
                                  value = inputResults@covariate.adj.pvalues[,whichTypeCol])
      }
      colnames(finmydat) = c(indvarname, outcomename, "PAdjVal")
      finmydat[,indvarname] = as.character(finmydat[,indvarname])
      finmydat[,outcomename] = as.character(finmydat[,outcomename])
      finmydat.coef <- reshape2::melt(mydat.interac)
      if(coefficient == "analyte"){
        finmydat.coef <- data.frame(Var1 = finmydat.coef$Var1,
                               Var2 = finmydat.coef$Var2,
                               value = inputResults@covariate.coefficients$a)
      }else if(coefficient == "stype"){
        finmydat.coef <- data.frame(Var1 = finmydat.coef$Var1,
                                    Var2 = finmydat.coef$Var2,
                                    value = inputResults@covariate.coefficients[,whichTypeCol])
      }
      colnames(finmydat.coef) = c(indvarname, outcomename, "Coef")
      finmydat.coef[,indvarname] = as.character(finmydat.coef[,indvarname])
      finmydat.coef[,outcomename] = as.character(finmydat.coef[,outcomename])
      finmydat.rsq <- reshape2::melt(mydat.rsq)
      colnames(finmydat.rsq) = c(indvarname, outcomename, "rsquared")
      finmydat.rsq[,indvarname] = as.character(finmydat.rsq[,indvarname])
      finmydat.rsq[,outcomename] = as.character(finmydat.rsq[,outcomename])
      
      # Filter by coefficient.
      if(coeffPercentile > 0){
        first_half <- getQuantileForCoefficient(finmydat.coef$Coef,
                                                coeffPercentile)[1]
        second_half <- getQuantileForCoefficient(finmydat.coef$Coef,
                                                 coeffPercentile)[2]
        keepers_first <- which(finmydat.coef$Coef > second_half)
        keepers_second <- which(finmydat.coef$Coef < first_half)
        keepers <- c(keepers_first, keepers_second)
        finmydat <- finmydat[keepers,]
        finmydat.coef <- finmydat.coef[keepers,]
        finmydat.rsq <- finmydat.rsq[keepers,]
      }
      # Filter by p-value cutoff.
      if(pvalcutoff != 1) { #(no filtering)
        keepers2 <- which(finmydat$PAdjVal <= pvalcutoff)
        finmydat <- finmydat[keepers2,]
        finmydat.coef <- finmydat.coef[keepers2,]
        finmydat.rsq <- finmydat.rsq[keepers2,]
      }
      # Filter by r-squared.
      keepers3 <- which(finmydat.rsq$rsquared >= rsquaredCutoff)
      filtResults <- data.frame(as.character(finmydat[keepers3,indvarname]),
                                as.character(finmydat[keepers3,outcomename]))
      colnames(filtResults) <- c(indvarname, outcomename)
      
      # Melt the adjusted p-values and interaction terms.
      cornames = NULL
      adjp <- reshape2::melt(inputResults@interaction.adj.pvalues)
      p <-  reshape2::melt(inputResults@interaction.pvalues)
      interact <- reshape2::melt(inputResults@interaction.coefficients)
      cornames <- paste(as.character(filtResults[,indvarname]),
                        as.character(filtResults[,outcomename]), sep = "__")
      rownames(p) <- paste(as.character(p[,1]),as.character(p[,2]), sep = "__")
      rownames(adjp) <- paste(as.character(adjp[,1]),as.character(adjp[,2]), sep = "__")
      rownames(interact) <- paste(as.character(interact[,1]),as.character(interact[,2]), sep = "__")
      outp <- p[cornames,]
      outpadj <- adjp[cornames,]
      outinteract <- interact[cornames,]
      
      # Combine all results.
      filtResults <- cbind(filtResults,
                           finmydat.coef$Coef[keepers3],
                           finmydat.rsq$rsquared[keepers3],
                           outp$value,
                           outpadj$value)
      
      # Filter the covariate coefficients, if appropriate, and add them to the frame.
      if(length(mydat.covar.coeff) > 0){
        original_names = setdiff(rownames(outp),
                                 rownames(mydat.covar.coeff)[which(!is.na(mydat.covar.coeff[,"(Intercept)"]))])
        flipped_names = unlist(lapply(original_names, function(name){
          analytes = strsplit(name, "__")[[1]]
          new_name = paste(analytes[2], analytes[1], sep = "__")
        }))
        mydat.covar.coeff[original_names,] = mydat.covar.coeff[flipped_names,]
        mydat.covar.coeff <- mydat.covar.coeff[which(!is.na(mydat.covar.coeff[,1])),]
        mydat.covar.coeff <- mydat.covar.coeff[rownames(outp),]
        filtResults <- cbind(filtResults, mydat.covar.coeff)
      }
      coefTypeName <- "interaction_coeff"
      if(coefficient == "analyte"){
        coefTypeName <- "analyte_coeff"
      }else if(coefficient == "stype"){
        coefTypeName <- "phenotype_coeff"
      }
      colnames(filtResults)[3:6]<-c(coefTypeName, "rsquared", "Pval","FDRadjPval")
    }
    
    # Set first two column names.
    colnames(filtResults)[1:2] <- c("Analyte1", "Analyte2")
    
    # If outcome isn't "type", change it.
    which_int <- which(grepl("a:", colnames(filtResults), fixed = TRUE) == TRUE)
    colnames(filtResults)[which_int] <- "a:type"
    which_type <- which(lapply(colnames(filtResults), function(name){
      retval <- FALSE
      if(substr(name, 1, 4) == "type"){
        retval <- TRUE
      }
      return(retval)
    }) == TRUE)
    colnames(filtResults)[which_type] <- "type"
    
    # Set the row names.
    rownames(filtResults) <- paste(as.character(filtResults$Analyte1),
                                   as.character(filtResults$Analyte2), sep = "__")
    
    # Print and return the results.
    message(paste(nrow(filtResults), 'pairs found given cutoffs'))
    return(filtResults)
  }else{
    stop(paste("Valid coefficients for filtering are 'interaction', 'stype', and 'analyte'.",
          "If 'stype' or 'analyte' are chosen, the model results must include covariate p-values",
          "and coefficients. To obtain these, set save.covar.pvals = TRUE in the RunIntLim() function."))
  }
  
}

#' Retrieve significant pairs, based on adjusted p-values, interaction coefficient
#' percentile, and R^2 values for continuous models
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param pvalCutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param coeffPercentile percentile cutoff for absolute value of coefficient
#' @param rsquaredCutoff cutoff for lowest r-squared value
#' @param coefficient Coefficient to filter by. Default is interaction coefficient.
#' Other options are "stype" (phenotype) and "analyte" (independent analyte).
#' Other options are "stype" and "analyte".
#' @return A data frame with the following columns for each pair of analytes:
#' "Analyte1", "Analyte2", "interaction_coeff", "Pval", "FDRadjPval", and "rsquared".
#' Optionally, coefficients for each covariate may also be included.
#' @export
ProcessResultsContinuous<- function(inputResults,
                         coeffPercentile=0.1,
                         pvalCutoff=0.05,
                         rsquaredCutoff=0.0,
                         coefficient = "interaction"){

  if(!methods::is(inputResults, "IntLimResults")) {
    stop("input data is not a IntLim class")
  }

  # Melt matrices to generate pairwise data frames.
  format_coeff = reshape2::melt(inputResults@interaction.coefficients)
  format_pval = reshape2::melt(inputResults@interaction.pvalues)
  format_adjp = reshape2::melt(inputResults@interaction.adj.pvalues)
  format_rsquared = reshape2::melt(inputResults@model.rsquared)
  whichTypeCol <- which(grepl("^type", colnames(inputResults@covariate.coefficients)) == TRUE)
  if(coefficient == "analyte"){
    format_coeff <- data.frame(Var1 = format_coeff$Var1,
                               Var2 = format_coeff$Var2,
                               value = inputResults@covariate.coefficients$a)
    format_adjp <- data.frame(Var1 = format_coeff$Var1,
                               Var2 = format_coeff$Var2,
                               value = inputResults@covariate.adj.pvalues$a)
  }else if(coefficient == "stype"){
    format_coeff <- data.frame(Var1 = format_coeff$Var1,
                               Var2 = format_coeff$Var2,
                               value = inputResults@covariate.coefficients[,whichTypeCol])
    format_adjp <- data.frame(Var1 = format_coeff$Var1,
                              Var2 = format_coeff$Var2,
                              value = inputResults@covariate.adj.pvalues[,whichTypeCol])
  }

  # Set rownames for melted frames.
  rownames(format_coeff) <- paste(as.character(format_coeff[,1])
                                  ,as.character(format_coeff[,2]),
                                  sep = "__")
  rownames(format_pval) <- paste(as.character(format_pval[,1]),
                                 as.character(format_pval[,2]), 
                                 sep = "__")
  rownames(format_adjp) <- paste(as.character(format_adjp[,1]),
                                 as.character(format_adjp[,2]), 
                                 sep = "__")
  rownames(format_rsquared) <- paste(as.character(format_rsquared[,1]),
                                     as.character(format_rsquared[,2]), 
                                     sep = "__")
  tofilter = cbind(format_coeff, format_pval$value, 
                   format_adjp$value, format_rsquared$value)
  coefTypeName <- "interaction_coeff"
  if(coefficient == "analyte"){
    coefTypeName <- "analyte_coeff"
  }else if(coefficient == "stype"){
    coefTypeName <- "phenotype_coeff"
  }
  colnames(tofilter)<-c("Analyte1", "Analyte2", coefTypeName, "Pval","FDRadjPval", "rsquared")

  #get top and bottom cutoffs (need highest positive and highest negative coeffs)
  first_half = getQuantileForCoefficient(tofilter[,coefTypeName], 
                                                    coeffPercentile)[1]
  second_half = getQuantileForCoefficient(tofilter[,coefTypeName], 
                                                     coeffPercentile)[2]
  #sort
  tofilter_sortedbycoeff <- tofilter[order(tofilter[,coefTypeName]),]

  #filter by coefficient
  filtered_by_coeff = tofilter_sortedbycoeff[tofilter_sortedbycoeff[,coefTypeName]>=
                                               second_half | 
                                               tofilter_sortedbycoeff[,coefTypeName]
                                             <= first_half,]

  #filter by pvalue
  filtered_by_pval = filtered_by_coeff[filtered_by_coeff$FDRadjPval <= pvalCutoff,]
  
  #filter by r-squared value
  filtered_by_rsquared = filtered_by_pval[filtered_by_pval$rsquared >= rsquaredCutoff,]

  #remove NA's.
  filtered_no_na = filtered_by_rsquared[which(!is.na(filtered_by_rsquared$Pval)),]

  # Add the covariate coefficients, if appropriate.
  filt.results <- filtered_no_na
  mydat.covar.coeff <- inputResults@covariate.coefficients
  if(length(mydat.covar.coeff) > 0){
    original_names = setdiff(rownames(filtered_no_na), 
                             rownames(mydat.covar.coeff)[which(!is.na(mydat.covar.coeff[,"(Intercept)"]))])
    flipped_names = unlist(lapply(original_names, function(name){
      analytes = strsplit(name, "__")[[1]]
      new_name = paste(analytes[2], analytes[1], sep = "__")
    }))
    mydat.covar.coeff[original_names,] = mydat.covar.coeff[flipped_names,]
    mydat.covar.coeff <- mydat.covar.coeff[which(!is.na(mydat.covar.coeff[,1])),]
    mydat.covar.coeff <- mydat.covar.coeff[rownames(filtered_no_na),]
    filt.results <- cbind(filtered_no_na, mydat.covar.coeff)
  }

  #place in objec to return
  return(filt.results)
}