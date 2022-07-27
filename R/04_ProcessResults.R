#' Retrieve significant pairs, based on adjusted p-values.
#' For each pair that is statistically significant, calculate the
#' correlation within group1 (e.g. cancer) and the correlation within group2 (e.g.
#' non-cancer).  Users can then remove pairs with a difference in correlations between
#' groups 1 and 2 less than a user-defined threshold.
#'
#' @include internalfunctions.R
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with analyte levels
#'  and associated meta-data
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient 
#' (default bottom 10 percent (high negative coefficients) and top 10 percent 
#' (high positive coefficients))
#' @param rsquaredCutoff cutoff for lowest r-squared value
#' @return IntResults object with model results (now includes correlations)
#' @export
ProcessResults <- function(inputResults,
				inputData,
				pvalcutoff=0.05,
				interactionCoeffPercentile=0,
				rsquaredCutoff = 0){
  
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
  if(interactionCoeffPercentile < 0 || interactionCoeffPercentile > 1){
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
  	                           interactionCoeffPercentile,
  	                           pvalcutoff, rsquaredCutoff)
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
      finmydat <- reshape2::melt(mydat)
      colnames(finmydat) = c(indvarname, outcomename, "PAdjVal")
      finmydat[,indvarname] = as.character(finmydat[,indvarname])
      finmydat[,outcomename] = as.character(finmydat[,outcomename])
      finmydat.interac <- reshape2::melt(mydat.interac)
      colnames(finmydat.interac) = c(indvarname, outcomename, "InteracCoef")
      finmydat.interac[,indvarname] = as.character(finmydat.interac[,indvarname])
      finmydat.interac[,outcomename] = as.character(finmydat.interac[,outcomename])
      finmydat.rsq <- reshape2::melt(mydat.rsq)
      colnames(finmydat.rsq) = c(indvarname, outcomename, "rsquared")
      finmydat.rsq[,indvarname] = as.character(finmydat.rsq[,indvarname])
      finmydat.rsq[,outcomename] = as.character(finmydat.rsq[,outcomename])

      # Filter by interaction coefficient.
      if(interactionCoeffPercentile > 0){
        first_half <- getQuantileForInteractionCoefficient(mydat.interac,
                                                             interactionCoeffPercentile)[1]
        second_half <- getQuantileForInteractionCoefficient(mydat.interac,
                                                            interactionCoeffPercentile)[2]
        keepers_first <- which(finmydat.interac$InteracCoef > second_half)
        keepers_second <- which(finmydat.interac$InteracCoef < first_half)
        keepers <- c(keepers_first, keepers_second)
        finmydat <- finmydat[keepers,]
        finmydat.interac <- finmydat.interac[keepers,]
        finmydat.rsq <- finmydat.rsq[keepers,]
      }
      # Filter by p-value cutoff.
      if(pvalcutoff != 1) { #(no filtering)
    		keepers2 <- which(finmydat$PAdjVal <= pvalcutoff)
    		finmydat <- finmydat[keepers2,]
    		finmydat.interac <- finmydat.interac[keepers2,]
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
    	                                  finmydat.interac$InteracCoef[keepers3],
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
    	colnames(filtResults)[3:6]<-c("interaction_coeff", "rsquared", "Pval","FDRadjPval")
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
}

#' Retrieve significant pairs (aka filter out nonsignificant pairs) 
#' based on value of analyte:type interaction coefficient from linear model
#'
#' @param inputResults IntLimResults object with model results: output of RunIntLim
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient
#' default bottom 10 percent (high negative coefficients) and top 10 percent (high positive coefficients)
#' @param pvalCutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param rsquaredCutoff cutoff of R-squared value for filtering (default 0, no filtering)
#' @return A data frame with the following columns for each pair of analytes:
#' "Analyte1", "Analyte2", "interaction_coeff", "Pval", "FDRadjPval", and "rsquared".
#' Optionally, coefficients for each covariate may also be included.
#' @export
ProcessResultsContinuous<- function(inputResults,
                         interactionCoeffPercentile=0.1,
                         pvalCutoff=0.05,
                         rsquaredCutoff=0.0){

  if(!methods::is(inputResults, "IntLimResults")) {
    stop("input data is not a IntLim class")
  }

  # Melt matrices to generate pairwise data frames.
  format_coeff = reshape2::melt(inputResults@interaction.coefficients)
  format_pval = reshape2::melt(inputResults@interaction.pvalues)
  format_adjp = reshape2::melt(inputResults@interaction.adj.pvalues)
  format_rsquared = reshape2::melt(inputResults@model.rsquared)

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
  colnames(tofilter) = c("Analyte1", "Analyte2", "interaction_coeff", "Pval","FDRadjPval", "rsquared")

  #get top and bottom cutoffs (need highest positive and highest negative coeffs)
  first_half = getQuantileForInteractionCoefficient(tofilter$interaction_coeff, 
                                                    interactionCoeffPercentile)[1]
  second_half = getQuantileForInteractionCoefficient(tofilter$interaction_coeff, 
                                                     interactionCoeffPercentile)[2]
  #sort
  tofilter_sortedbycoeff <- tofilter[order(tofilter$interaction_coeff),]

  #filter by coefficient
  filtered_by_coeff = tofilter_sortedbycoeff[tofilter_sortedbycoeff$interaction_coeff>=
                                               second_half | 
                                               tofilter_sortedbycoeff$interaction_coeff 
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