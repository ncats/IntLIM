#' Retrieve significant gene-metabolite pairs, based on adjusted p-values.
#' For each gene-metabolite pair that is statistically significant, calculate the
#' correlation within group1 (e.g. cancer) and the correlation within group2 (e.g.
#' non-cancer).  Users can then remove pairs with a difference in correlations between
#' groups 1 and 2 less than a user-defined threshold.
#'
#' @include internalfunctions.R
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param diffcorr cutoff of differences in correlations for filtering (default 0.5)
#' @param corrtype spearman or pearson or other parameters allowed by cor() function 
#' (default spearman)
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient 
#' (default bottom 10 percent (high negative coefficients) and top 10 percent 
#' (high positive coefficients))
#' @param treecuts user-selected number of clusters (of gene-metabolite pairs) 
#' to cut the tree into
#' @param rsquaredCutoff cutoff for lowest r-squared value
#' @return IntResults object with model results (now includes correlations)
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,mydata,treecuts=2)
#' }
#' @export
ProcessResults <- function(inputResults,
				inputData,
				pvalcutoff=0.05,
				diffcorr=0.5,
				corrtype="spearman",
				interactionCoeffPercentile=0.5,
				rsquaredCutoff = 0.0,
				treecuts = 0){
  independent.var.type = "gene"
	if(inputResults@outcome == "metabolite") {
		mydat <-inputResults@interaction.adj.pvalues
		mydat.interac <- inputResults@interaction.coefficients
		mydat.rsq <- inputResults@model.rsquared 
	}
		#mydat <- reshape2::melt(inputResults@interaction.adj.pvalues)}
	else if (inputResults@outcome == "gene") {
		mydat <-t(inputResults@interaction.adj.pvalues)
		mydat.interac <- t(inputResults@interaction.coefficients)
		mydat.rsq <- t(inputResults@model.rsquared)
		independent.var.type = "metabolite"
	}

	incommon <- getCommon(inputData,inputResults@stype)
	p <- incommon$p


	if(length(unique(p)) !=2 & is.null(diffcorr) & is.null(corrtype)) {
	  inputResults = ProcessResultsContinuous(inputResults,
  	                           interactionCoeffPercentile,
  	                           pvalcutoff, rsquaredCutoff, 
  	                           outcome=inputResults@outcome, 
  	                           independent.var.type = independent.var.type)
	} else if (length(unique(p)) !=2){
  	  stop(paste("IntLim requires two categories only for correlation analysis.  
  	             Make sure the column",inputResults@stype,"only has two unique 
  	             values or set diffcorr and corrtype to null to switch to interaction 
  	             coefficient analysis"))
	}  else{
    	gene <- incommon$gene
    	metab <- incommon$metab
    	gp1 <- which(p == unique(p)[1])
    	cor1.m <- stats::cor(t(gene[rownames(mydat),gp1]),t(metab[colnames(mydat),gp1]),
    	                     method=corrtype)
    	gp2 <- which(p == unique(p)[2])
      cor2.m <- stats::cor(t(gene[rownames(mydat),gp2]),t(metab[colnames(mydat),gp2]),
                           method=corrtype)
      
      # Create melted matrices to filter by inputs.
      fincor1 <- reshape2::melt(cor1.m)
      colnames(fincor1) = c("Gene", "Metabolite", "Cor1")
      fincor1$Gene = as.character(fincor1$Gene)
      fincor1$Metabolite = as.character(fincor1$Metabolite)
      fincor2 <- reshape2::melt(cor2.m)
      colnames(fincor2) = c("Gene", "Metabolite", "Cor2")
      fincor2$Gene = as.character(fincor2$Gene)
      fincor2$Metabolite = as.character(fincor2$Metabolite)
      finmydat <- reshape2::melt(mydat)
      colnames(finmydat) = c("Gene", "Metabolite", "PAdjVal")
      finmydat$Gene = as.character(finmydat$Gene)
      finmydat$Metabolite = as.character(finmydat$Metabolite)
      finmydat.interac <- reshape2::melt(mydat.interac)
      colnames(finmydat.interac) = c("Gene", "Metabolite", "InteracCoef")
      finmydat.interac$Gene = as.character(finmydat.interac$Gene)
      finmydat.interac$Metabolite = as.character(finmydat.interac$Metabolite)
      finmydat.rsq <- reshape2::melt(mydat.rsq)
      colnames(finmydat.rsq) = c("Gene", "Metabolite", "Rsq")
      finmydat.rsq$Gene = as.character(finmydat.rsq$Gene)
      finmydat.rsq$Metabolite = as.character(finmydat.rsq$Metabolite)
      
      if(!is.null(interactionCoeffPercentile)){
        if(interactionCoeffPercentile > 0){
          first_half <- getQuantileForInteractionCoefficient(mydat.interac, 
                                                             interactionCoeffPercentile)[1]
          second_half <- getQuantileForInteractionCoefficient(mydat.interac, 
                                                              interactionCoeffPercentile)[2]
          keepers_first <- which(finmydat.interac$InteracCoef > second_half)
          keepers_second <- which(finmydat.interac$InteracCoef < first_half)
          keepers <- c(keepers_first, keepers_second)
          fincor1 = fincor1[keepers,]
          fincor2 = fincor2[keepers,]
          finmydat <- finmydat[keepers,]
          finmydat.interac <- finmydat.interac[keepers,]
          finmydat.rsq <- finmydat.rsq[keepers,]
        }
        
      }
      if(pvalcutoff == 1) { #(no filtering)
    		genenames <- as.character(fincor1$Gene)
    		metabnames <- as.character(fincor2$Metabolite)
    	} else {
    		keepers2 <- which(finmydat$PAdjVal <= pvalcutoff)
    		fincor1 <- fincor1[keepers2,]
    		fincor2 <- fincor2[keepers2,]
    		finmydat <- finmydat[keepers2,]
    		finmydat.interac <- finmydat.interac[keepers2,]
    		finmydat.rsq <- finmydat.rsq[keepers2,]
    		genenames <- as.character(fincor1$Gene)
    		metabnames <- as.character(fincor2$Metabolite)
    	}
      if(diffcorr > 0){
        mydiffcor <- abs(fincor1$Cor1-fincor2$Cor2)
        keepers3 <- which(mydiffcor >= diffcorr)
        fincor1 <- fincor1[keepers3,]
        fincor2 <- fincor2[keepers3,]
        finmydat <- finmydat[keepers3,]
        finmydat.interac <- finmydat.interac[keepers3,]
        finmydat.rsq <- finmydat.rsq[keepers3,]
        genenames <- as.character(fincor1$Gene)
        metabnames <- as.character(fincor2$Metabolite)
      }
      keepers4 <- which(finmydat.rsq$Rsq >= rsquaredCutoff)
    	inputResults@filt.results <- data.frame(metab=metabnames[keepers4],
    		gene=genenames[keepers4])
    	inputResults@filt.results <- cbind(inputResults@filt.results,
    	                                   fincor1$Cor1[keepers4],
    	                                   fincor2$Cor2[keepers4])
    	colnames(inputResults@filt.results)[3:4]=paste0(setdiff(as.character(unlist(unique(p))),""),"_cor")
    	diff.corr <- inputResults@filt.results[,4] - inputResults@filt.results[,3]

    	inputResults@filt.results <- cbind(inputResults@filt.results, diff.corr)
    	if(inputResults@outcome == "metabolite") {
                    adjp <- reshape2::melt(inputResults@interaction.adj.pvalues)
    		p <-  reshape2::melt(inputResults@interaction.pvalues)
    	} else if (inputResults@outcome == "gene") {
                    adjp <- reshape2::melt(t(inputResults@interaction.adj.pvalues))
    		p <- reshape2::melt(t(inputResults@interaction.pvalues))
    	}

    	cornames <- paste(as.character(inputResults@filt.results[,"metab"]),
    	                  as.character(inputResults@filt.results[,"gene"]))
    	rownames(p) <- paste(as.character(p[,2]),as.character(p[,1]))
    	rownames(adjp) <- paste(as.character(adjp[,2]),as.character(adjp[,1]))
    	outp <- p[cornames,]
    	outpadj <- adjp[cornames,]
    	outinteract <- (reshape2::melt(t(inputResults@interaction.coefficients)))[cornames,]

    	inputResults@filt.results = cbind(inputResults@filt.results,
    	                                  finmydat.interac$InteracCoef[keepers4],
    	                                  finmydat.rsq$Rsq[keepers4],
    	                                  outp$value, 
    	                                  outpadj$value)
    	colnames(inputResults@filt.results)[6:9]=c("interaction_coeff", "Rsq", "Pval","FDRadjPval")
	}

	if (treecuts > 0){

	hc.rows<- stats::hclust(stats::dist(inputResults@filt.results[,c(3,4)]))
	cluster <- stats::cutree(hc.rows, k = treecuts)

	inputResults@filt.results = cbind(inputResults@filt.results, cluster)


	}

print(paste(nrow(inputResults@filt.results), 'gene-metabolite pairs found given cutoffs'))
return(inputResults)
}

#' Retrieve significant metabolite-metabolite pairs, based on adjusted p-values.
#' For each metabolite-metabolite pair that is statistically significant, calculate the
#' correlation within group1 (e.g. cancer) and the correlation within group2 (e.g.
#' non-cancer).  Users can then remove pairs with a difference in correlations between
#' groups 1 and 2 less than a user-defined threshold.
#'
#' @include internalfunctions.R
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression
#' ,metabolite abundances, and associated meta-data
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param diffcorr cutoff of differences in correlations for filtering (default 0.5)
#' @param corrtype spearman or pearson or other parameters allowed by cor() function 
#' (default spearman)
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient 
#' (default bottom 10 percent (high negative coefficients) and top 10 percent 
#' (high positive coefficients))
#' @param treecuts user-selected number of clusters (of gene-metabolite pairs) 
#' to cut the tree into
#' @param rsquaredCutoff cutoff for lowest r-squared value
#' @return IntResults object with model results (now includes correlations)
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,mydata,treecuts=2)
#' }
#' @export
ProcessResultsMetabolitePairs <- function(inputResults,
                           inputData,
                           pvalcutoff=0.05,
                           diffcorr=0.5,
                           corrtype="spearman",
                           interactionCoeffPercentile=0.5,
                           rsquaredCutoff=0.0,
                           treecuts = 0){
  
  mydat <-inputResults@interaction.adj.pvalues
  mydat.interac <- inputResults@interaction.coefficients
  mydat.rsq <- inputResults@model.rsquared
  
  mytypes <- names(Biobase::assayData(inputData))
  incommon <- NULL
  if(any(mytypes == "expression")){
    incommon <- getCommon(inputData,inputResults@stype)
  }else{
    incommon <- formatSingleOmicInput(inputData,inputResults@stype, type = "metabolite")
  }
  p <- incommon$p
  metab <- incommon$metab 
  
  if(length(unique(p)) !=2 & is.null(diffcorr) & is.null(corrtype)) {
    inputResults = ProcessResultsContinuous(inputResults,
                                            interactionCoeffPercentile,
                                            pvalcutoff,
                                            rsquaredCutoff,
                                            outcome="metabolite", 
                                            independent.var.type = "metabolite")
  } else if (length(unique(p)) !=2){
    stop(paste("IntLim requires two categories only for correlation analysis. 
               Make sure the column",inputResults@stype,"only has two unique 
               values or set diffcorr and corrtype to null to switch to interaction 
               coefficient analysis"))
  }  else{
    gp1 <- which(p == unique(p)[1])
    cor1.m <- cor(t(metab[rownames(mydat),gp1]), method = corrtype)
    cor1.m[lower.tri(cor1.m,diag=TRUE)] <- NA
    gp2 <- which(p == unique(p)[2])
    cor2.m <- cor(t(metab[rownames(mydat),gp2]), method = corrtype)
    cor2.m[lower.tri(cor2.m,diag=TRUE)] <- NA
    finmydat <- reshape2::melt(mydat)
    colnames(finmydat) = c("Metab1", "Metab2", "PAdjVal")
    finmydat$Metab1 = as.character(finmydat$Metab1)
    finmydat$Metab2 = as.character(finmydat$Metab2)
    fincor1 <- reshape2::melt(cor1.m)
    colnames(fincor1) = c("Metab1", "Metab2", "Cor1")
    fincor1$Metab1 = as.character(fincor1$Metab1)
    fincor1$Metab2 = as.character(fincor1$Metab2)
    fincor2 <- reshape2::melt(cor2.m)
    colnames(fincor2) = c("Metab1", "Metab2", "Cor2")
    fincor2$Metab1 = as.character(fincor2$Metab1)
    fincor2$Metab2 = as.character(fincor2$Metab2)
    finmydat.interac <- reshape2::melt(mydat.interac)
    colnames(finmydat.interac) = c("Metab1", "Metab2", "InteracCoef")
    finmydat.interac$Metab1 = as.character(finmydat.interac$Metab1)
    finmydat.interac$Metab2 = as.character(finmydat.interac$Metab2)
    finmydat.rsq <- reshape2::melt(mydat.rsq)
    colnames(finmydat.rsq) = c("Metab1", "Metab2", "Rsq")
    finmydat.rsq$Metab1 = as.character(finmydat.rsq$Metab1)
    finmydat.rsq$Metab2 = as.character(finmydat.rsq$Metab2)
    
    if(!is.null(interactionCoeffPercentile)){
      if(interactionCoeffPercentile > 0){
        first_half <- getQuantileForInteractionCoefficient(mydat.interac, 
                                                           interactionCoeffPercentile)[1]
        second_half <- getQuantileForInteractionCoefficient(mydat.interac, 
                                                            interactionCoeffPercentile)[2]
        keepers_first <- which(finmydat.interac$InteracCoef > second_half)
        keepers_second <- which(finmydat.interac$InteracCoef < first_half)
        keepers <- c(keepers_first, keepers_second)
        fincor1 = fincor1[keepers,]
        fincor2 = fincor2[keepers,]
        finmydat <- finmydat[keepers,]
        finmydat.interac <- finmydat.interac[keepers,]
        finmydat.rsq <- finmydat.rsq[keepers,]
        metab1names <- as.character(fincor1$Metab1)
        metab2names <- as.character(fincor2$Metab2)
      }
    }
    if(pvalcutoff == 1) { #(no filtering)
      metab1names <- as.character(fincor1$Metab1)
      metab2names <- as.character(fincor2$Metab2)
    } else {
      keepers2 <- which(finmydat$PAdjVal <= pvalcutoff)
      fincor1 <- fincor1[keepers2,]
      fincor2 <- fincor2[keepers2,]
      finmydat.interac <- finmydat.interac[keepers2,]
      finmydat.rsq <- finmydat.rsq[keepers2,]
      finmydat <- finmydat[keepers2,]
      metab1names <- as.character(fincor1$Metab1)
      metab2names <- as.character(fincor2$Metab2)
    }
    if(diffcorr > 0){
      mydiffcor = abs(fincor1$Cor1-fincor2$Cor2)
      keepers3 <- which(mydiffcor >= diffcorr)
      fincor1 <- fincor1[keepers3,]
      fincor2 <- fincor2[keepers3,]
      finmydat.interac <- finmydat.interac[keepers3,]
      finmydat.rsq <- finmydat.rsq[keepers3,]
      finmydat <- finmydat[keepers3,]
      metab1names <- as.character(fincor1$Metab1)
      metab2names <- as.character(fincor2$Metab2)
    }
    myrsquared <- finmydat.rsq$Rsq
    keepers4 <- which(myrsquared >= rsquaredCutoff)
    inputResults@filt.results <- data.frame(metab1=metab1names[keepers4],
                                            metab2=metab2names[keepers4])
    inputResults@filt.results <- cbind(inputResults@filt.results,
                                       fincor1$Cor1[keepers4],
                                       fincor2$Cor2[keepers4])

    colnames(inputResults@filt.results)[3:4]=paste0(setdiff(as.character(unlist(unique(p))),""),"_cor")
    diff.corr <- inputResults@filt.results[,4] - inputResults@filt.results[,3]
    inputResults@filt.results <- cbind(inputResults@filt.results, diff.corr)
    adjp <- reshape2::melt(inputResults@interaction.adj.pvalues)
    p <-  reshape2::melt(inputResults@interaction.pvalues)
    cornames <- paste(as.character(inputResults@filt.results[,"metab2"]),
                      as.character(inputResults@filt.results[,"metab1"]))
    rownames(p) <- paste(as.character(p[,2]),as.character(p[,1]))
    rownames(adjp) <- paste(as.character(adjp[,2]),as.character(adjp[,1]))
    outp <- p[cornames,]
    outpadj <- adjp[cornames,]
    
    inputResults@filt.results = cbind(inputResults@filt.results,
                                      finmydat.interac$InteracCoef[keepers4],
                                      finmydat.rsq$Rsq[keepers4],
                                      outp$value, outpadj$value)
    colnames(inputResults@filt.results)[6:9]=c("interaction_coeff", "Rsq", "Pval","FDRadjPval")
  }
  
  if (treecuts > 0){
    
    hc.rows<- stats::hclust(stats::dist(inputResults@filt.results[,c(3,4)]))
    cluster <- stats::cutree(hc.rows, k = treecuts)
    
    inputResults@filt.results = cbind(inputResults@filt.results, cluster)
    
    
  }
  print(paste(nrow(inputResults@filt.results), 'metabolite-metabolite pairs found given cutoffs'))
  return(inputResults)
}

#' Retrieve significant gene-gene pairs, based on adjusted p-values.
#' For each gene-gene pair that is statistically significant, calculate the
#' correlation within group1 (e.g. cancer) and the correlation within group2 (e.g.
#' non-cancer).  Users can then remove pairs with a difference in correlations between
#' groups 1 and 2 less than a user-defined threshold.
#'
#' @include internalfunctions.R
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,metabolite abundances, and associated meta-data
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param diffcorr cutoff of differences in correlations for filtering (default 0.5)
#' @param corrtype spearman or pearson or other parameters allowed by cor() function (default spearman)
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient (default bottom 10 percent (high negative coefficients) and top 10 percent (high positive coefficients))
#' @param treecuts user-selected number of clusters (of gene-gene pairs) to cut the tree into
#' @param rsquaredCutoff cutoff for lowest r-squared value
#' @return IntResults object with model results (now includes correlations)
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,mydata,treecuts=2)
#' }
#' @export
ProcessResultsGenePairs <- function(inputResults,
                                    inputData,
                                    pvalcutoff=0.05,
                                    diffcorr=0.5,
                                    corrtype="spearman",
                                    interactionCoeffPercentile=0.5,
                                    rsquaredCutoff = 0.0,
                                    treecuts = 0){
  
  mydat <-inputResults@interaction.adj.pvalues
  mydat.interac <- inputResults@interaction.coefficients
  mydat.rsq <- inputResults@model.rsquared
  
  mytypes <- names(Biobase::assayData(inputData))
  incommon <- NULL
  if(any(mytypes == "metabolite")){
    incommon <- getCommon(inputData,inputResults@stype)
  }else{
    incommon <- formatSingleOmicInput(inputData,inputResults@stype, type = "expression")
  }
  p <- incommon$p
  gene <- incommon$gene 
  
  if(length(unique(p)) !=2 & is.null(diffcorr) & is.null(corrtype)) {
    inputResults = ProcessResultsContinuous(inputResults,
                                            interactionCoeffPercentile,
                                            pvalcutoff,
                                            rsquaredCutoff,
                                            outcome = "gene",
                                            independent.var.type = "gene")
  } else if (length(unique(p)) !=2){
    stop(paste("IntLim requires two categories only for correlation analysis. 
               Make sure the column",inputResults@stype,"only has two unique 
               values or set diffcorr and corrtype to null to switch to interaction 
               coefficient analysis"))
  }  else{
    gp1 <- which(p == unique(p)[1])
    cor1.g <- cor(t(gene[rownames(mydat),gp1]), method = corrtype)
    cor1.g[lower.tri(cor1.g,diag=TRUE)] <- NA
    gp2 <- which(p == unique(p)[2])
    cor2.g <- cor(t(gene[rownames(mydat),gp2]), method = corrtype)
    cor2.g[lower.tri(cor2.g,diag=TRUE)] <- NA
    finmydat <- reshape2::melt(mydat)
    colnames(finmydat) = c("Gene1", "Gene2", "PAdjVal")
    finmydat$Gene1 = as.character(finmydat$Gene1)
    finmydat$Gene2 = as.character(finmydat$Gene2)
    fincor1 <- reshape2::melt(cor1.g)
    colnames(fincor1) = c("Gene1", "Gene2", "Cor1")
    fincor1$Gene1 = as.character(fincor1$Gene1)
    fincor1$Gene2 = as.character(fincor1$Gene2)
    fincor2 <- reshape2::melt(cor2.g)
    colnames(fincor2) = c("Gene1", "Gene2", "Cor2")
    fincor2$Gene1 = as.character(fincor2$Gene1)
    fincor2$Gene2 = as.character(fincor2$Gene2)
    finmydat.interac <- reshape2::melt(mydat.interac)
    colnames(finmydat.interac) = c("Gene1", "Gene2", "InteracCoef")
    finmydat.interac$Gene1 = as.character(finmydat.interac$Gene1)
    finmydat.interac$Gene2 = as.character(finmydat.interac$Gene2)
    finmydat.rsq <- reshape2::melt(mydat.rsq)
    colnames(finmydat.rsq) = c("Gene1", "Gene2", "Rsq")
    finmydat.rsq$Gene1 = as.character(finmydat.rsq$Gene1)
    finmydat.rsq$Gene2 = as.character(finmydat.rsq$Gene2)
    
    if(!is.null(interactionCoeffPercentile)){
      if(interactionCoeffPercentile > 0){
        first_half <- getQuantileForInteractionCoefficient(mydat.interac, 
                                                           interactionCoeffPercentile)[1]
        second_half <- getQuantileForInteractionCoefficient(mydat.interac, 
                                                            interactionCoeffPercentile)[2]
        keepers_first <- which(finmydat.interac$InteracCoef > second_half)
        keepers_second <- which(finmydat.interac$InteracCoef < first_half)
        keepers <- c(keepers_first, keepers_second)
        fincor1 = fincor1[keepers,]
        fincor2 = fincor2[keepers,]
        finmydat <- finmydat[keepers,]
        finmydat.interac <- finmydat.interac[keepers,]
        finmydat.rsq <- finmydat.rsq[keepers,]
        gene1names <- as.character(fincor1$Gene1)
        gene2names <- as.character(fincor2$Gene2)
      }
    }
    if(pvalcutoff == 1) { #(no filtering)
      gene1names <- as.character(fincor1$Gene1)
      gene2names <- as.character(fincor2$Gene2)
    } else {
      keepers2 <- which(finmydat$PAdjVal <= pvalcutoff)
      fincor1 <- fincor1[keepers2,]
      fincor2 <- fincor2[keepers2,]
      finmydat.interac <- finmydat.interac[keepers2,]
      finmydat.rsq <- finmydat.rsq[keepers2,]
      finmydat <- finmydat[keepers2,]
      gene1names <- as.character(fincor1$Gene1)
      gene2names <- as.character(fincor2$Gene2)
    }
    if(diffcorr > 0){
      mydiffcor = abs(fincor1$Cor1-fincor2$Cor2)
      keepers3 <- which(mydiffcor >= diffcorr)
      fincor1 <- fincor1[keepers3,]
      fincor2 <- fincor2[keepers3,]
      finmydat.interac <- finmydat.interac[keepers3,]
      finmydat.rsq <- finmydat.rsq[keepers3,]
      finmydat <- finmydat[keepers3,]
      gene1names <- as.character(fincor1$Gene1)
      gene2names <- as.character(fincor2$Gene2)
    }
    myrsquared <- finmydat.rsq$Rsq
    keepers4 <- which(myrsquared >= rsquaredCutoff)
    inputResults@filt.results <- data.frame(gene1=gene1names[keepers4],
                                            gene2=gene2names[keepers4])
    inputResults@filt.results <- cbind(inputResults@filt.results,
                                       fincor1$Cor1[keepers4],
                                       fincor2$Cor2[keepers4])
    
    colnames(inputResults@filt.results)[3:4]=paste0(setdiff(as.character(unlist(unique(p))),""),"_cor")
    diff.corr <- inputResults@filt.results[,4] - inputResults@filt.results[,3]
    inputResults@filt.results <- cbind(inputResults@filt.results, diff.corr)
    adjp <- reshape2::melt(inputResults@interaction.adj.pvalues)
    p <-  reshape2::melt(inputResults@interaction.pvalues)
    cornames <- paste(as.character(inputResults@filt.results[,"gene2"]),
                      as.character(inputResults@filt.results[,"gene1"]))
    rownames(p) <- paste(as.character(p[,2]),as.character(p[,1]))
    rownames(adjp) <- paste(as.character(adjp[,2]),as.character(adjp[,1]))
    outp <- p[cornames,]
    outpadj <- adjp[cornames,]
    
    inputResults@filt.results = cbind(inputResults@filt.results,
                                      finmydat.interac$InteracCoef[keepers4],
                                      finmydat.rsq$Rsq[keepers4],
                                      outp$value, outpadj$value)
    colnames(inputResults@filt.results)[6:9]=c("interaction_coeff", "Rsq", "Pval","FDRadjPval")
  }
  
  if (treecuts > 0){
    
    hc.rows<- stats::hclust(stats::dist(inputResults@filt.results[,c(3,4)]))
    cluster <- stats::cutree(hc.rows, k = treecuts)
    
    inputResults@filt.results = cbind(inputResults@filt.results, cluster)
    
    
  }
  print(paste(nrow(inputResults@filt.results), 'gene-gene pairs found given cutoffs'))
  return(inputResults)
}

#' Retrieve significant gene-metabolite / metabolite-metabolite pairs (aka filter out nonsignificant pairs) based on value of gene:type interaction coefficient from linear model
#'
#' @param inputResults IntLimResults object with model results: output of RunIntLim
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient
#' default bottom 10 percent (high negative coefficients) and top 10 percent (high positive coefficients)
#' @param pvalCutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param rsquaredCutoff cutoff of R-squared value for filtering (default 0, no filtering)
#' @param outcome 'metabolite' or 'gene' must be set as outcome/independent variable
#' (default is 'metabolite')
#' @param independent.var.type 'metabolite' or 'gene' must be set as independent variable
#' (default is 'metabolite')
#' @export
ProcessResultsContinuous<- function(inputResults,
                         interactionCoeffPercentile=0.1,
                         pvalCutoff=0.05,
                         rsquaredCutoff=0.0,
                         outcome = "metabolite",
                         independent.var.type = "gene"){

  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }

  #merge and properly name all data to return
  gene_metabolite_format_coeff = reshape2::melt(inputResults@interaction.coefficients)
  gene_metabolite_format_pval = reshape2::melt(inputResults@interaction.pvalues)
  gene_metabolite_format_adjp = reshape2::melt(inputResults@interaction.adj.pvalues)
  gene_metabolite_format_rsquared = reshape2::melt(inputResults@model.rsquared)
  tofilter = cbind(gene_metabolite_format_coeff, gene_metabolite_format_pval$value, 
                   gene_metabolite_format_adjp$value, gene_metabolite_format_rsquared$value)
  if((outcome == "metabolite" && independent.var.type == "gene") || (outcome == "gene" &&
                                                               independent.var.type == "metabolite")){
    colnames(tofilter) = c("gene", "metab", "interaction_coeff", "Pval","FDRadjPval", "rsquared")
  }
  else if(outcome == "metabolite" && independent.var.type == "metabolite"){
    colnames(tofilter) = c("metab1", "metab2", "interaction_coeff", "Pval","FDRadjPval", "rsquared")
  }
  else if(outcome == "gene" && independent.var.type == "gene"){
    colnames(tofilter) = c("gene1", "gene2", "interaction_coeff", "Pval","FDRadjPval", "rsquared")
  }
  else{
    stop("Error! outcome and independent.var.type must each be one of the following: gene, metabolite.")
  }

  #get top and bottom cutoffs (need highest positive and highest negative coeffs)
  first_half = getQuantileForInteractionCoefficient(tofilter$interaction_coeff, 
                                                    interactionCoeffPercentile)[1]
  second_half = getQuantileForInteractionCoefficient(tofilter$interaction_coeff, 
                                                     interactionCoeffPercentile)[2]

  #sort
  tofilter_sortedbycoeff <- tofilter[order(tofilter$interaction_coeff),]

  #filter by coefficient
  filtered_by_coeff = tofilter_sortedbycoeff[tofilter_sortedbycoeff$interaction_coeff>
                                               second_half | 
                                               tofilter_sortedbycoeff$interaction_coeff 
                                             < first_half,]

  #filter by pvalue
  filtered_by_pval = filtered_by_coeff[filtered_by_coeff$FDRadjPval < pvalCutoff,]
  
  #filter by r-squared value
  filtered_by_rsquared = filtered_by_pval[filtered_by_pval$rsquared >= rsquaredCutoff,]
  
  #remove NA's.
  filtered_no_na = filtered_by_rsquared[which(!is.na(filtered_by_rsquared$Pval)),]

  #place in objec to return
  inputResults@filt.results = filtered_no_na
  return(inputResults)

}

#' Create results table, which includes significant gene:metabolite pairs, associated p-values,
#' and correlations in each category evaluated.
#'
#' @param inputResults IntLimResults object with model results (output of ProcessResults())
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,mydata)
#' mytable <- CreateResultsTable(myres)
#' }
#' @export
#   CreateResultsTable <- function(inputResults) {
#        a<-inputResults@corr
#        a$cordiff<-round(abs(a[,3]-a[,4]),3)
#        a[,3]<-round(a[,3],2)
#        a[,4]<-round(a[,4],2)
#        p <- padj <- c()
#        if(inputResults@outcome=="metabolite") {
#                for (i in 1:nrow(a)) {
#                        g <- which(rownames(inputResults@interaction.pvalues) == a$gene[i])
#                        m <- which(colnames(inputResults@interaction.pvalues) == a$metab[i])
#                        if(length(g)==0 || length(m)==0) {p<-c(p,NA);padj<-c(padj,NA)} else {
#                                p <- c(p,inputResults@interaction.pvalues[g,m])
#				padj <- c(padj,inputResults@interaction.adj.pvalues[g,m])
#                              padj <- c(padj,inputResults@interaction.adj.pvalues[a$gene[i],a$metab[i]])
#                        }
#                }
#        } else if (inputResults@outcome=="gene") {
#               for (i in 1:nrow(a)) {
#                        g <- which(rownames(inputResults@interaction.pvalues) == a$gene[i])
#                        m <- which(colnames(inputResults@interaction.pvalues) == a$metab[i])
#                        if(length(g)==0 || length(m)==0) {p<-c(p,NA)} else {
 #                               p <- c(p,inputResults@interaction.pvalues[a$metab[i],a$gene[i]])
#				p <- c(p,inputResults@interaction.pvalues[g,m])
#                                padj <- c(padj,inputResults@interaction.adj.pvalues[m,g])
#                        }
#                }
#        }
#        else {stop("Outcome should be either 'metabolite' or 'gene'")}
#        a$pval <- p
#        a$adjpval <- padj
#        table<-a[order(a$adjpval,decreasing = TRUE),]
#	rownames(table) <- NULL
#        return(table)
#   }
#
#
