#' Get some stats after reading in data
#'
#' @import magrittr
#' @import highcharter
#'
#' @include MultiDataSet_extendedfunctions.R
#'
#' @param inputData IntLimObject output of ReadData()
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @return a highcharter object
#'
#' @examples
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' PlotDistributions(mydata)
#' @export
PlotDistributions <- function(inputData,viewer=T, palette="Set1"){
  . <- c()
  if (length(palette) == 2) {
    cols <- c(palette)
  }
  else if (length(palette) == 1) {
    cols <- RColorBrewer::brewer.pal(3, palette)[1:2]
  }
  else {
    stop("palette must either be an RColorBrewer palette or a vector of hex colors of size 2")
  }
  mytypes <- names(Biobase::assayData(inputData))
  g <- NULL
  m <- NULL
  p <- NULL
  boxplotOptions <- list(
    fillColor = '#ffffff',
    lineWidth = 2,
    medianColor = '#000000',
    medianWidth = 2,
    stemColor = '#000000',
    stemDashStyle = 'dot',
    stemWidth = 1,
    whiskerColor = '#000000',
    whiskerLength = '20%',
    whiskerWidth = 3)
  if(any(mytypes=="expression")){
    mygene <- as.data.frame(Biobase::assayDataElement(inputData[["expression"]],'exprs'))
    toplot <- suppressMessages(reshape2::melt(mygene))
    df <- dplyr::tibble(value = toplot$value, by = toplot$variable) %>% dplyr::group_by_at("by") %>%
      dplyr::do(data = grDevices::boxplot.stats(.$value))
    bxps <- purrr::map(df$data, "stats")
    outs <- purrr::map2_df(seq(nrow(df)), df$data, function(x, y) {
      if (length(y$out) > 0)
        d <- dplyr::tibble(x = x - 1, y = y$out)
      else d <- dplyr::tibble()
      d
    })
    outs <- data.frame(outs, 'z' = colnames(mygene)[outs$x + 1])
    z <- outs$z
    # To try to get the gene names of outliers, would have to go back and get the gene names from original data frame and put htem in outs$color
    
    g <- highcharter::highchart(width = 750, height = 750 ) %>%
      highcharter::hc_title(text = "Gene Expression",
                            style = list(color = '#2E1717',
                                         fontWeight = 'bold', fontSize = "20px")) %>%
      highcharter::hc_plotOptions(
        boxplot = boxplotOptions
      ) %>%
      hc_add_series(data = bxps,name = "Gene Expression", type="boxplot",color=cols[1],showInLegend=FALSE) %>%
      highcharter::hc_add_series(data=list_parse(outs),name = "Gene Expression",
                                 type="scatter",color=cols[1],showInLegend=FALSE,
                                 tooltip = list(headerFormat = "", pointFormat = "{point.z} <br/> {point.y}",
                                                showInLegend = FALSE)) %>%
      highcharter::hc_yAxis(title = list(text = "log(expression)",
                                         style = list(fontSize = "13px")),
                            labels = list(format = "{value}")) %>%
      highcharter::hc_xAxis(labels="", categories = colnames(mygene)) %>%
      highcharter::hc_tooltip(valueDecimals = 2) %>%
      highcharter::hc_exporting(enabled = TRUE)
  }
  if(any(mytypes=="metabolite")){
    mymetab <- Biobase::assayDataElement(inputData[["metabolite"]],'metabData')
    toplot <- suppressMessages(reshape2::melt(mymetab))
    df <- dplyr::data_frame(value = toplot$value, by = toplot$variable) %>%
      dplyr::group_by_at("by") %>%
      dplyr::do(data = grDevices::boxplot.stats(.$value))
    bxps <- purrr::map(df$data, "stats")
    outs <- purrr::map2_df(seq(nrow(df)), df$data, function(x, y) {
      if (length(y$out) > 0)
        d <- dplyr::data_frame(x = x - 1, y = y$out)
      else d <- dplyr::data_frame()
      d
    })
    outs <- data.frame(outs, 'z' = colnames(mymetab)[outs$x + 1])
    z <- outs$z
    
    m <- highcharter::highchart(width = 750, height = 750 ) %>%
      highcharter::hc_title(text = "Metabolite Levels",
                            style = list(color = '#2E1717',
                                         fontWeight = 'bold', fontSize = "20px")) %>%
      highcharter::hc_plotOptions(
        boxplot = boxplotOptions
      ) %>%
      highcharter::hc_add_series(data = bxps,name = "Metabolite Levels",
                                 type="boxplot",color=cols[2],showInLegend=FALSE) %>%
      highcharter::hc_add_series(data=list_parse(outs),name = "Metabolite Levels",
                                 type="scatter",color=cols[2],showInLegend=FALSE,tooltip = list(headerFormat = "", pointFormat = "{point.z} <br/> {point.y}",
                                                                                                showInLegend = FALSE)) %>%
      
      highcharter::hc_yAxis(title = list(text = "log(abundances)",
                                         style = list(fontSize = "13px")),
                            labels = list(format = "{value}")) %>%
      highcharter::hc_xAxis(labels="", categories = colnames(mymetab)) %>%
      highcharter::hc_tooltip(valueDecimals = 2) %>%
      highcharter::hc_exporting(enabled = TRUE)
  }
  if(!is.null(g) & !is.null(m)){
    if (viewer == TRUE) {
      p <-
        htmltools::browsable(highcharter::hw_grid(g, m, ncol = 2, rowheight = 550))
    }
    else {
      p <- highcharter::hw_grid(g, m)
    }
  } else if(!is.null(g)){
    if (viewer == TRUE) {
      p <-
        htmltools::browsable(highcharter::hw_grid(g, ncol = 1, rowheight = 550))
    }
    else {
      p <- highcharter::hw_grid(g)
    }
  } else if(!is.null(m)){
    if (viewer == TRUE) {
      p <-
        htmltools::browsable(highcharter::hw_grid(m, ncol = 1, rowheight = 550))
    }
    else {
      p <- highcharter::hw_grid(m)
    }
  }
  
  return(p)
}

#' PCA plots of data for QC
#'
#' @import magrittr
#' @import highcharter
#'
#' @include MultiDataSet_extendedfunctions.R
#'
#' @param inputData IntLimObject output of ReadData()
#' @param stype category to color-code by (can be more than two categories)
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @param common whether or not samples that are in common between the metabolite and gene expression datasets should be plotted (T/F); default is TRUE
#' @return a highcharter object
#'
#' @examples
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' PlotPCA(mydata,stype = "PBO_vs_Leukemia")
#' @export
PlotPCA <- function(inputData,viewer=T,stype=NULL,common=T,
        palette = "Set1") {

  mytypes <- names(Biobase::assayData(inputData))
  biobase_pdata <- NULL
  if(any(mytypes == "metabolite")){
    biobase_pdata <- Biobase::pData(inputData[["metabolite"]])
  }else{
    biobase_pdata <- Biobase::pData(inputData[["expression"]])
  }

  if(is.null(stype) || is.numeric(biobase_pdata[,stype]) == TRUE) {
		warning("The resulting PCA plot is not color-coded because you did not provide 
		        a categorical variable in 'stype'")
		mytype <- NULL
  } else if (length(intersect(colnames(biobase_pdata),stype))!=1) {
		stop(paste0("You provided ",stype, "as your stype variable but it does not exist in your data"))
  } else {
  	mytype <- as.character(biobase_pdata[,stype])
    numcateg <- length(unique(mytype))
    if(length(palette) >= 2) {
      cols <- palette
    } else {
      if(numcateg == 1) {
         if(length(palette)==1) {cols <- RColorBrewer::brewer.pal(3, palette)[1]
         } else {stop("palette should be an RColorBrewer palette or a vector of colors")}
      } else if (numcateg == 2) {
          if(length(palette)==1) {cols <- RColorBrewer::brewer.pal(3, palette)[1:2]
          } else {stop("palette should be an RColorBrewer palette or a vector of colors")}
      } else if (numcateg > 2) {
          if(length(palette)==1) {cols <- RColorBrewer::brewer.pal(numcateg, palette)
          } else {stop("palette should be an RColorBrewer palette or a vector of colors")}
      } else {stop("There are no values in your 'stype' column")}
               }
    }
  p <- NULL
  pg <- NULL
  pm <- NULL

	if(common==T) {
	  if(!any(mytypes=="expression") || !any(mytypes == "metabolite")){
	    stop("A dataset not containing both expression and metabolite data cannot run
	         with 'common' set to TRUE. Set 'common' to FALSE.")
	  } else {
	    if(is.null(stype) || is.numeric(biobase_pdata[,stype]) == TRUE) {
  			incommon <- getCommon(inputData)
  			mygene <- incommon$gene
  			gpca <- stats::prcomp(t(mygene),center=T,scale=F)
  			mymetab <- incommon$metab
  			mpca <- stats::prcomp(t(mymetab),center=T,scale=F)
  			gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),color=rep("blue",nrow(gpca$x)))
  			mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),color=rep("blue",nrow(mpca$x)))
        gds <- list_parse(gtoplot)
        pg <- highcharter::highchart(width = 350, height = 350 )
        pg <- pg %>% highcharter::hc_add_series(data=gds,type="scatter",
                                  tooltip = list(headerFormat="",
                                  pointFormat=paste("{point.label}","{point.z}")),
                                  showInLegend=FALSE)
        mds <- list_parse(mtoplot)
        pm <- highcharter::highchart(width = 350, height = 350)
        pm <- pm %>% highcharter::hc_add_series(data=mds,type="scatter",
                                  tooltip = list(headerFormat="",
                                  pointFormat=paste("{point.label}","{point.z}")),
                                  showInLegend=FALSE)
	    } else {
  			incommon <- getCommon(inputData,stype)
  			mygene <- incommon$gene
  			mymetab <- incommon$metab
  			alltype <- incommon$p
  			uniqtypes <- unique(alltype)
  			mycols <- as.character(alltype)
  			for (i in 1:numcateg) {
  				mycols[which(alltype==uniqtypes[i])] <- cols[i]
  			}
  			gpca <- stats::prcomp(t(mygene),center=T,scale=F)
  			mpca <- stats::prcomp(t(mymetab),center=T,scale=F)
  			gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),label=alltype,color=mycols)
  			mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),label=alltype,color=mycols)
  			mds <- list_parse(mtoplot)
        gds <- list_parse(gtoplot)
        pg <- highcharter::highchart(width = 350, height = 350)
  			pm <- highcharter::highchart(width = 350, height = 350)
        for (i in 1:length(uniqtypes)) {
          mytype <- unique(alltype)[i]
          gds <- list_parse(gtoplot[which(gtoplot$label==mytype),])
          pg <- pg %>% highcharter::hc_add_series(data=gds,type="scatter",
                                          name=mytype,color=cols[which(alltype==mytype)[1]],tooltip = list(headerFormat="",
                                          pointFormat=paste("{point.label}","{point.z}")),
                                          showInLegend=TRUE)
          mds <- list_parse(mtoplot[which(mtoplot$label==mytype),])
          pm <- pm %>% highcharter::hc_add_series(data=mds,type="scatter",
                                          name=mytype,color=cols[which(alltype==mytype)[1]],tooltip = list(headerFormat="",
                                          pointFormat=paste("{point.label}","{point.z}")),
                                          showInLegend=TRUE)
        }
	    }
		}
	} else { # common == F
	  
	  if(!is.null(stype) || is.numeric(biobase_pdata[,stype]) == TRUE) {
      # Compute PC's.
      gtypes <- NULL
      mtypes <- NULL
      uniqtypes <- NULL
      if(any(mytypes == "expression")){
        mygene <- as.data.frame(Biobase::assayDataElement(inputData[["expression"]],'exprs'))
        gpca <- stats::prcomp(t(mygene),center=T,scale=F)
        gtypes <- as.character(Biobase::pData(inputData[["expression"]])[,stype])
        uniqtypes <- unique(gtypes)
      }
      if(any(mytypes == "metabolite")){
        mymetab <- Biobase::assayDataElement(inputData[["metabolite"]],'metabData')
        mpca <- stats::prcomp(t(mymetab),center=T,scale=F)
        mtypes <- as.character(Biobase::pData(inputData[["metabolite"]])[,stype])
        uniqtypes <- unique(mtypes)
      }
      
      # Take the union of types.
      if(any(mytypes=="expression") && any(mytypes == "metabolite")){
	      uniqtypes <- unique(c(mtypes,gtypes))
      }
	      
      # Set up plots.
	    if(any(mytypes == "expression")){
	      gcols <- as.character(gtypes)
	      for (i in 1:numcateg) {
	        gcols[which(gtypes==uniqtypes[i])] <- cols[i]
	      }
	      # Deal with missing values or ""
	      if(length(which(gtypes==""))>0) {
	        gcols[which(gtypes=="")]="grey"
	        gtypes[which(gtypes=="")]="NA"
	      }
	      if(length(which(is.na(gtypes)))>0) {
	        gcols[which(is.na(gtypes))]="grey"
	      }
	      gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),label=gtypes,color=gcols)
	      gds <- list_parse(gtoplot)
	      pg <- highcharter::highchart(width = 350, height = 350 )
	      for (i in 1:length(unique(gtypes))) {
	        mytype <- unique(gtypes)[i]
	        gds <- list_parse(gtoplot[which(gtoplot$label==mytype),])
	        pg <- pg %>% highcharter::hc_add_series(data=gds,type="scatter",
	                                                name=mytype,color=cols[which(gtypes==mytype)[1]],tooltip = list(headerFormat="",
	                                                                                                                pointFormat=paste("{point.label}","{point.z}")),
	                                                showInLegend=TRUE)
	      }
	    }
      if(any(mytypes == "metabolite")){
        mcols <- as.character(mtypes)
        for (i in 1:numcateg) {
          mcols[which(mtypes==uniqtypes[i])] <- cols[i]
        }
        if (length(which(mtypes==""))>0) {
          mcols[which(mtypes=="")]="grey"
          mtypes[which(mtypes=="")]="NA"
        }
        if(length(which(is.na(mtypes)))>0) {
          mcols[which(is.na(mtypes))]="grey"
        }
        mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),label=mtypes,color=mcols)
        mds <- list_parse(mtoplot)
        pm <- highcharter::highchart(width = 350, height = 350 )
        for (i in 1:length(unique(mtypes))) {
          mytype <- unique(mtypes)[i]
          mds <- list_parse(mtoplot[which(mtoplot$label==mytype),])
          pm <- pm %>% highcharter::hc_add_series(data=mds,type="scatter",
                                                  name=mytype,color=cols[which(mtypes==mytype)[1]],tooltip = list(headerFormat="",
                                                                                                                  pointFormat=paste("{point.label}","{point.z}")),
                                                  showInLegend=TRUE)
        }
      }
	    
	  } else { #stype is null
	    if(any(mytypes == "expression")){
	      gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),label="",color=rep("blue",nrow(gpca$x)))
	      gds <- list_parse(gtoplot)
	      pg <- highcharter::highchart(width = 350, height = 350 )
	      pg <- pg %>% highcharter::hc_add_series(data=gds,type="scatter",
	                                              tooltip = list(headerFormat="",
	                                                             pointFormat=paste("{point.label}","{point.z}")),
	                                              showInLegend=FALSE)
	    }
	    if(any(mytypes == "metabolite")){
	      mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),label="",color=rep("blue",nrow(mpca$x)))
	      mds <- list_parse(mtoplot)
	      pm <- highcharter::highchart(width = 350, height = 350)
	      pm <- pm %>% highcharter::hc_add_series(data=mds,type="scatter",
	                                              tooltip = list(headerFormat="",
	                                                             pointFormat=paste("{point.label}","{point.z}")),
	                                              showInLegend=FALSE)
	    }
	  }
	}
  # end common == F

  # Set up plots.
  if(any(mytypes == "metabolite")){
    mpercvar=round((mpca$sdev)^2 / sum(mpca$sdev^2)*100,2)
    pm <- pm %>% highcharter::hc_title(text="PCA of metabolites") %>%
      highcharter::hc_xAxis(title=list(text=paste0("PC1:",round(mpercvar[1],1),"%"))) %>%
      highcharter::hc_yAxis(title=list(text=paste0("PC2:",round(mpercvar[2],2),"%"))) %>%
      hc_chart(zoomType = "xy")
  }
  if(any(mytypes == "expression")){
    gpercvar=round((gpca$sdev)^2 / sum(gpca$sdev^2)*100,2)
    pg <- pg %>% highcharter::hc_title(text="PCA of genes") %>%
      highcharter::hc_xAxis(title=list(text=paste0("PC1:",round(gpercvar[1],1),"%"))) %>%
      highcharter::hc_yAxis(title=list(text=paste0("PC2:",round(gpercvar[2],2),"%"))) %>%
      hc_chart(zoomType = "xy")
  }
  p <- NULL
  if(any(mytypes=="expression") && any(mytypes == "metabolite")){
    if (viewer == TRUE) {
      p <-htmltools::browsable(highcharter::hw_grid(pg, pm, ncol = 2, rowheight = 550))
    } else {
      p <- highcharter::hw_grid(pg, pm)
    }
  } else if(any(mytypes=="expression")){
    if (viewer == TRUE) {
      p <-htmltools::browsable(highcharter::hw_grid(pg, ncol = 1, rowheight = 550))
    } else {
      p <- highcharter::hw_grid(pg)
    }
  } else if(any(mytypes=="metabolite")){
    if (viewer == TRUE) {
      p <-htmltools::browsable(highcharter::hw_grid(pm, ncol = 1, rowheight = 550))
    } else {
      p <- highcharter::hw_grid(pm)
    }
  }
  
  return(p)

}



#' Visualize the distribution of unadjusted p-values from linear models
#'
#' @include IntLimResults_extendedfunctions.R
#'
#' @param IntLimResults output of RunIntLim()
#' @param breaks the number of breaks to use in histogram (see hist() documentation for more details)
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' DistPvalues(myres)
#' }
#' @importFrom graphics boxplot par
#' @export
DistPvalues<- function(IntLimResults,breaks=100) {

    hist(IntLimResults@interaction.pvalues,breaks=breaks,
	main="Histogram of Interaction P-values")
}

#' Visualize the distribution of unadjusted p-values for all covariates
#' from linear models using a bar chart.
#'
#' @include IntLimResults_extendedfunctions.R
#'
#' @param IntLimResults output of RunIntLim()
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' DistPvalues(myres)
#' }
#' @export
PValueBoxPlots<- function(IntLimResults) {
  if(length(IntLimResults@covariate.pvalues) == 0){
    print("Error! You must set save.covar.pvals to TRUE when running IntLIM to run PValueBoxPlots")
  }else{
    par(mar=c(8, 4.1, 4.1, 2.1))
    boxplot(IntLimResults@covariate.pvalues, las = 3, ylim = c(0,1), ylab = "P-Value")
  }
}

#' Visualize the distribution of unadjusted p-values from linear models
#'
#' @include IntLimResults_extendedfunctions.R
#'
#' @param IntLimResults output of RunIntLim()
#' @param breaks the number of breaks to use in histogram (see hist() documentation for more details)
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' DistPvalues(myres)
#' }
#' @export
DistRSquared<- function(IntLimResults,breaks=100) {
  
  hist(IntLimResults@model.rsquared,breaks=breaks,
       main="Histogram of Interaction R-Squared Values")
}


#' Plot correlation heatmap
#'
#' @import magrittr
#'
#' @param inputResults IntLimResults object (output of ProcessResults())
#' @param inputData IntLimObject output of ReadData() or FilterData()
#' @param top_pairs cutoff of the top pairs, sorted by adjusted p-values, to be plotted (plotting more than 1200 can take some time) (default: 1200)
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @param treecuts number of clusters (of gene-metabolite pairs) to cut the tree into for color-coding
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @return a highcharter object
#'@param static allows user to decide whether heatmap is interactive or static
#'@param html.file allows user to specify file path to output heatmap onto (used for non-static heatmaply objects)
#'@param pdf.file allows user to specify file path to output heatmap onto (used for static heatmap.2 objects)
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,mydata)
#' CorrHeatmap(myres)
#' }
#' @export
CorrHeatmap <- function(inputResults,inputData, viewer=T,top_pairs=1200,treecuts=2, palette = NULL, static = FALSE,
                        html.file=NULL, pdf.file=NULL) {
type <- cor <- c()

	if(nrow(inputResults@filt.results)==0) {
		stop("Make sure you run ProcessResults before making the heatmap")
	}
  mytypes <- names(Biobase::assayData(inputData))
  if(any(mytypes == "expression") && any(mytypes == "metabolite")){
    incommon <- getCommon(inputData,inputResults@stype)
  }else if(any(mytypes == "metabolite")){
    incommon <- formatSingleOmicInput(inputData,inputResults@stype, type = "metabolite")
  }else if(any(mytypes == "expression")){
    incommon <- formatSingleOmicInput(inputData,inputResults@stype, type = "expression")
  }
  p <- incommon$p
  if(length(unique(p)) !=2){
    stop("CorrHeatmap requires 2 discrete phenotypes. Do not run with continuous phenotypes.")
  }
  else{
    allres <- inputResults@filt.results
    if(nrow(allres)>top_pairs) {
      allp <- inputResults@filt.results[,"FDRadjPval"]
      allres <- allres[order(allp,decreasing=F)[1:top_pairs],]
    }
    
    toplot <- data.frame(name=paste(allres[,1],allres[,2],sep=" vs "),
                         allres[,3:4])
    suppressMessages(
      meltedtoplot <- tidyr::gather(
        toplot,
        type,cor,colnames(toplot)[2],colnames(toplot)[3]))
    
    #all possible values of X (type) and Y (name)
    theXAxis <- as.character(meltedtoplot[, "type"])
    theYAxis <- as.character(meltedtoplot[, "name"])
    
    #unique values of X and Y
    theUniqueY <- as.character(unique(theYAxis))
    theUniqueX <- as.character(unique(theXAxis))
    
    # Substitute words with position on the meatrix
    for (i in 1:length(theUniqueY)){
      num <- which(theYAxis == theUniqueY[i])
      theYAxis[num] <- i
    }
    for (i in 1:length(theUniqueX)) {
      num <- which(theXAxis == theUniqueX[i])
      theXAxis[num] <- i
    }
    # New package heatmaply here
    type <- unique(meltedtoplot[,'type'])
    num <- nrow(meltedtoplot[meltedtoplot[,'type'] == type[1],])
    heat_data <- matrix(data = 0, nrow =num,ncol = 2)
    row.names(heat_data) <- meltedtoplot[1:num,1]
    colnames(heat_data) <- gsub("_cor","",c(type[1],type[2]))
    heat_data[,1] <- meltedtoplot[1:num,3]
    
    heat_data[,2] <- meltedtoplot[-1:-num,3]
    if (is.null(palette)){
      palette=grDevices::colorRampPalette(c("#D01C8B", "#F1B6DA", "#F7F7F7", "#B8E186", "#4DAC26")) (255)[255:1]
    }
    
    if(static == FALSE){
      hm <- heatmaply::heatmaply(heat_data,main = "Correlation heatmap",
                                 k_row = treecuts,#k_col = 2,
                                 margins = c(80,5),
                                 dendrogram = "row",
                                 y_axis_font_size ="1px",
                                 colors = palette,
                                 key.title = 'Correlation \n differences',
                                 file=html.file)
      hm
      
      if(!is.null(pdf.file)){
        
        hmr <- heatmaply::heatmapr(heat_data,main = "Correlation heatmap",
                                   k_row = treecuts,#k_col = 2,
                                   margins = c(80,5),
                                   dendrogram = "row",
                                   y_axis_font_size ="1px",
                                   colors = palette,
                                   key.title = 'Correlation \n differences' )
        
        row_dend = hmr$rows
        grDevices::pdf(file=pdf.file, width=12, height=6.3)
        gplots::heatmap.2(heat_data,main = "Correlation \n heatmap",
                          dendrogram = "row",
                          col = palette,
                          density.info = 'none',
                          key.title = 'Correlation \n differences',
                          labRow = rep('',nrow(heat_data)),
                          cexCol = 0.05 + 0.25/log10(ncol(heat_data)),
                          trace = 'none', Rowv = row_dend)
        grDevices::dev.off()
      }
      return(hm)
    }else{
      
      hmr <- heatmaply::heatmapr(heat_data,main = "Correlation heatmap",
                                 k_row = treecuts,#k_col = 2,
                                 margins = c(80,5),
                                 dendrogram = "row",
                                 y_axis_font_size ="1px",
                                 colors = palette,
                                 key.title = 'Correlation \n differences' )
      
      row_dend = hmr$rows
      gplots::heatmap.2(heat_data,main = "Correlation \n heatmap",
                        dendrogram = "row",
                        col = palette,
                        density.info = 'none',
                        key.title = 'Correlation \n differences',
                        labRow = rep('',nrow(heat_data)),
                        cexCol = 0.05 + 0.25/log10(ncol(heat_data)),
                        trace = 'none', Rowv = row_dend)
      
      if(!is.null(pdf.file)){
        grDevices::pdf(file=pdf.file, width=12, height=6.3)
        gplots::heatmap.2(heat_data,main = "Correlation \n heatmap",
                          dendrogram = "row",
                          col = palette,
                          density.info = 'none',
                          key.title = 'Correlation \n differences',
                          labRow = rep('',nrow(heat_data)),
                          cexCol = 0.05 + 0.25/log10(ncol(heat_data)),
                          trace = 'none', Rowv = row_dend)
        grDevices::dev.off()
      }
      
      
    }
    
    
    if(!is.null(html.file) & static==TRUE){
      hm.html.out <- heatmaply::heatmaply(heat_data,main = "Correlation heatmap",
                                          k_row = treecuts,#k_col = 2,
                                          margins = c(80,5),
                                          dendrogram = "row",
                                          y_axis_font_size ="1px",
                                          colors = palette,
                                          key.title = 'Correlation \n differences',
                                          file=html.file)
    }
  }
		

}

#' scatter plot of gene-metabolite pairs (based on user selection)
#'
#' @import magrittr
#' @import highcharter
#'
#' @param inputData IntLimObject output of ReadData() or FilterData()
#' @param stype category to color-code by
##' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
##' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param geneName string of select geneName
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @param metabName string of select metabName
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' PlotGMPair(mydata,stype="PBO_vs_Leukemia","DLG4","(p-Hydroxyphenyl)lactic acid")
#'
#' }
#' @export
PlotGMPair<- function(inputData,stype=NULL,geneName,metabName,palette = "Set1",
	viewer=T) {

      if(is.null(stype)) {
	stop("Users must define stype which defines the categories to be compared (e.g. tumor vs non-tumor).  This could be the same parameter that was used to run RunIntLim()")
	}
      if (length(palette) == 2) {
        cols <- c(palette)
      } else if (length(palette) == 1) {
        cols <- RColorBrewer::brewer.pal(3, palette)[1:2]
      } else {
        stop("palette must either be an RColorBrewer palette or a vector of hex colors of size 2")
      }

   if (class(inputData) != "MultiDataSet") {
        stop("input data is not a MultiDataSet class")
}

    incommon <- getCommon(inputData,stype)

	if(is.null(stype)) {
                stop("A category to colorcode by (e.g. stype) must be provided")
        } else if (length(intersect(colnames(Biobase::pData(inputData[["metabolite"]])),stype))!=1) {
                stop(paste0("You provided ",stype, "as your stype variable but it does not exist in your data"))
	} else {
                mytypes <- incommon$p
        }

    gene<-incommon$gene
    if(length(which(rownames(gene)==geneName))>0) {
	    sGene<-gene[geneName,]
    } else {
	stop(paste0("The gene ",geneName," was not found in your data"))
    }

    metab<-incommon$metab
    if(length(which(rownames(metab)==metabName))>0) {
    	sMetab<-as.numeric(metab[metabName,])
    } else {
	stop(paste0("The metabolite ",metabName," was not found in your data"))
    }

    if(length(unique(mytypes))!=2) {
	stop(paste0("The group selected, '",stype,"', should only contain two different categories"))
    }

    mycols <- as.character(mytypes)
    mycols[which(mytypes==unique(mytypes)[1])] <- cols[1]
    mycols[which(mytypes==unique(mytypes)[2])] <- cols[2]

    data<-data.frame(x=sGene,y=sMetab,z=colnames(gene),label=mytypes,color=mycols)
    
    # Get points to draw the lines for each phenotype by hand

    uniqtypes=as.character(unique(mytypes))

    # Starting with phenotype 1, get min and max x values constrained to the values of y
    # The reason we do this, is because the lines do not necessary need to go out to the max or min of x, particularly
    # when slopes are really steep (abline does this automatically but not highcharter)
    getLinePoints <- function(data,mytypes, uniqtypes, currenttype) {
    	y=data$y[which(data$label==uniqtypes[currenttype])]; x=data$x[which(data$label==uniqtypes[currenttype])]
	min <- min(data$x[which(mytypes==uniqtypes[currenttype])])
    	max <- max(data$x[which(mytypes==uniqtypes[currenttype])])

    	m1<-stats::glm(y ~ x)
    	line1<-data.frame(x=c(max,min),
	y=c(stats::predict(m1,data.frame(x=c(max,min)))))
	return(data.frame(x=c(max,min), y=c(stats::predict(m1,data.frame(x=c(max,min))))))
    }

    line1 <- getLinePoints(data,mytypes,uniqtypes,currenttype=1)
    line2 <- getLinePoints(data,mytypes, uniqtypes, currenttype=2)

    ds <- highcharter::list_parse(data)

        hc <- highcharter::highchart(width = 350, height = 350 ) %>%
                highcharter::hc_title(text=paste(geneName,' vs. ', metabName, sep = '')) %>%
                highcharter::hc_xAxis(title=list(text=geneName)) %>%
                highcharter::hc_yAxis(title=list(text=metabName)) %>%
                hc_chart(zoomType = "xy") %>%
                highcharter::hc_add_series(data=ds,type="scatter",#col=cols[1],
                        tooltip = list(headerFormat="",
                          pointFormat=paste("{point.label}","{point.z}")),
                        showInLegend=FALSE)

    hc <- hc %>%
        highcharter::hc_add_series(name = uniqtypes[1],
		data=line1,type='line',#name=sprintf("regression line %s",type1),
		color = cols[1],enableMouseTracking=FALSE,marker=FALSE) %>%
        highcharter::hc_add_series(name = uniqtypes[2],
		data=line2,type='line',#name=sprintf("regression line %s",type2),
		color = cols[2],enableMouseTracking=FALSE,marker=FALSE)

    hc
}

#' scatter plot of metabolite-gene pairs (based on user selection)
#'
#' @import magrittr
#'
#' @param inputData IntLimObject output of ReadData() or FilterData()
#' @param stype category to color-code by
##' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
##' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param geneName string of select geneName
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @param metabName string of select metabName
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' PlotGMPair(mydata,stype="PBO_vs_Leukemia","DLG4","(p-Hydroxyphenyl)lactic acid")
#'
#' }
#' @export
PlotMGPair<- function(inputData,stype=NULL,metabName,geneName,palette = "Set1",
                      viewer=T) {
  
  if(is.null(stype)) {
    stop("Users must define stype which defines the categories to be compared (e.g. tumor vs non-tumor).  This could be the same parameter that was used to run RunIntLim()")
  }
  if (length(palette) == 2) {
    cols <- c(palette)
  } else if (length(palette) == 1) {
    cols <- RColorBrewer::brewer.pal(3, palette)[1:2]
  } else {
    stop("palette must either be an RColorBrewer palette or a vector of hex colors of size 2")
  }
  
  if (class(inputData) != "MultiDataSet") {
    stop("input data is not a MultiDataSet class")
  }
  
  incommon <- getCommon(inputData,stype)
  
  if(is.null(stype)) {
    stop("A category to colorcode by (e.g. stype) must be provided")
  } else if (length(intersect(colnames(Biobase::pData(inputData[["metabolite"]])),stype))!=1) {
    stop(paste0("You provided ",stype, "as your stype variable but it does not exist in your data"))
  } else {
    mytypes <- incommon$p
  }
  
  gene<-incommon$gene
  if(length(which(rownames(gene)==geneName))>0) {
    sGene<-gene[geneName,]
  } else {
    stop(paste0("The gene ",geneName," was not found in your data"))
  }
  
  metab<-incommon$metab
  if(length(which(rownames(metab)==metabName))>0) {
    sMetab<-as.numeric(metab[metabName,])
  } else {
    stop(paste0("The metabolite ",metabName," was not found in your data"))
  }
  
  if(length(unique(mytypes))!=2) {
    stop(paste0("The group selected, '",stype,"', should only contain two different categories"))
  }
  
  mycols <- as.character(mytypes)
  mycols[which(mytypes==unique(mytypes)[1])] <- cols[1]
  mycols[which(mytypes==unique(mytypes)[2])] <- cols[2]
  
  data<-data.frame(x=sMetab,y=sGene,z=colnames(gene),label=mytypes,color=mycols)
  
  # Get points to draw the lines for each phenotype by hand
  
  uniqtypes=as.character(unique(mytypes))
  
  # Starting with phenotype 1, get min and max x values constrained to the values of y
  # The reason we do this, is because the lines do not necessary need to go out to the max or min of x, particularly
  # when slopes are really steep (abline does this automatically but not highcharter)
  getLinePoints <- function(data,mytypes, uniqtypes, currenttype) {
    y=data$y[which(data$label==uniqtypes[currenttype])]; x=data$x[which(data$label==uniqtypes[currenttype])]
    min <- min(data$x[which(mytypes==uniqtypes[currenttype])])
    max <- max(data$x[which(mytypes==uniqtypes[currenttype])])
    
    m1<-stats::glm(y ~ x)
    line1<-data.frame(x=c(max,min),
                      y=c(stats::predict(m1,data.frame(x=c(max,min)))))
    return(data.frame(x=c(max,min), y=c(stats::predict(m1,data.frame(x=c(max,min))))))
  }
  
  line1 <- getLinePoints(data,mytypes,uniqtypes,currenttype=1)
  line2 <- getLinePoints(data,mytypes, uniqtypes, currenttype=2)
  
  ds <- highcharter::list_parse(data)
  #cols=c("blue","pink")
  
  hc <- highcharter::highchart(width = 350, height = 350 ) %>%
    highcharter::hc_title(text=paste(metabName,' vs. ', geneName, sep = '')) %>%
    highcharter::hc_xAxis(title=list(text=metabName)) %>%
    highcharter::hc_yAxis(title=list(text=geneName)) %>%
    hc_chart(zoomType = "xy") %>%
    highcharter::hc_add_series(data=ds,type="scatter",#col=cols[1],
                               tooltip = list(headerFormat="",
                                              pointFormat=paste("{point.label}","{point.z}")),
                               showInLegend=FALSE)
  
  hc <- hc %>%
    highcharter::hc_add_series(name = uniqtypes[1],
                               data=line1,type='line',#name=sprintf("regression line %s",type1),
                               color = cols[1],enableMouseTracking=FALSE,marker=FALSE) %>%
    highcharter::hc_add_series(name = uniqtypes[2],
                               data=line2,type='line',#name=sprintf("regression line %s",type2),
                               color = cols[2],enableMouseTracking=FALSE,marker=FALSE)
  
  hc
}

#' scatter plot of metabolite-metabolite pairs (based on user selection)
#'
#' @import magrittr
#' @import highcharter
#'
#' @param inputData IntLimObject output of ReadData() or FilterData()
#' @param stype category to color-code by
##' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
##' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param metab1Name string of select metab1Name
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @param metab2Name string of select metab2Name
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' PlotMMPair(mydata,stype="PBO_vs_Leukemia","arginine","(p-Hydroxyphenyl)lactic acid")
#'
#' }
#' @export
PlotMMPair<- function(inputData,stype=NULL,metab1Name,metab2Name,palette = "Set1",
                      viewer=T) {
  
  if(is.null(stype)) {
    stop("Users must define stype which defines the categories to be compared (e.g. tumor vs non-tumor).  This could be the same parameter that was used to run RunIntLim()")
  }
  if (length(palette) == 2) {
    cols <- c(palette)
  } else if (length(palette) == 1) {
    cols <- RColorBrewer::brewer.pal(3, palette)[1:2]
  } else {
    stop("palette must either be an RColorBrewer palette or a vector of hex colors of size 2")
  }
  
  if (class(inputData) != "MultiDataSet") {
    stop("input data is not a MultiDataSet class")
  }
  
  mytypes <- names(Biobase::assayData(inputData))
  if(any(mytypes == "expression")){
    incommon <- getCommon(inputData,stype)
  }else{
    incommon <- formatSingleOmicInput(inputData,stype, type = "metabolite")
  }
  
  if(is.null(stype)) {
    stop("A category to colorcode by (e.g. stype) must be provided")
  } else if (length(intersect(colnames(Biobase::pData(inputData[["metabolite"]])),stype))!=1) {
    stop(paste0("You provided ",stype, "as your stype variable but it does not exist in your data"))
  } else {
    mytypes <- incommon$p
  }
  
  metab<-incommon$metab
  if(length(which(rownames(metab)==metab1Name))>0) {
    sMetab1<-as.numeric(metab[metab1Name,])
  } else {
    stop(paste0("The metabolite ",metab1Name," was not found in your data"))
  }
  
  metab<-incommon$metab
  if(length(which(rownames(metab)==metab2Name))>0) {
    sMetab2<-as.numeric(metab[metab2Name,])
  } else {
    stop(paste0("The metabolite ",metab2Name," was not found in your data"))
  }
  
  if(length(unique(mytypes))!=2) {
    stop(paste0("The group selected, '",stype,"', should only contain two different categories"))
  }
  
  mycols <- as.character(mytypes)
  mycols[which(mytypes==unique(mytypes)[1])] <- cols[1]
  mycols[which(mytypes==unique(mytypes)[2])] <- cols[2]
  
  data<-data.frame(x=sMetab1,y=sMetab2,z=colnames(metab),label=mytypes,color=mycols)

  # Get points to draw the lines for each phenotype by hand
  
  uniqtypes=as.character(unique(mytypes))
  
  # Starting with phenotype 1, get min and max x values constrained to the values of y
  # The reason we do this, is because the lines do not necessary need to go out to the max or min of x, particularly
  # when slopes are really steep (abline does this automatically but not highcharter)
  getLinePoints <- function(data,mytypes, uniqtypes, currenttype) {
    y=data$y[which(data$label==uniqtypes[currenttype])]; x=data$x[which(data$label==uniqtypes[currenttype])]
    min <- min(data$x[which(mytypes==uniqtypes[currenttype])])
    max <- max(data$x[which(mytypes==uniqtypes[currenttype])])
    
    m1<-stats::glm(y ~ x)
    line1<-data.frame(x=c(max,min),
                      y=c(stats::predict(m1,data.frame(x=c(max,min)))))
    return(data.frame(x=c(max,min), y=c(stats::predict(m1,data.frame(x=c(max,min))))))
  }
  
  line1 <- getLinePoints(data,mytypes,uniqtypes,currenttype=1)
  line2 <- getLinePoints(data,mytypes, uniqtypes, currenttype=2)
  
  ds <- highcharter::list_parse(data)

  #cols=c("blue","pink")
  
  hc <- highcharter::highchart(width = 350, height = 350 ) %>%
    highcharter::hc_title(text=paste(metab1Name,' vs. ', metab2Name, sep = '')) %>%
    highcharter::hc_xAxis(title=list(text=metab1Name)) %>%
    highcharter::hc_yAxis(title=list(text=metab2Name)) %>%
    hc_chart(zoomType = "xy") %>%
    highcharter::hc_add_series(data=ds,type="scatter",#col=cols[1],
                               tooltip = list(headerFormat="",
                                              pointFormat=paste("{point.label}","{point.z}")),
                               showInLegend=FALSE)
  
  
  hc <- hc %>%
    highcharter::hc_add_series(name = uniqtypes[1],
                               data=line1,type='line',#name=sprintf("regression line %s",type1),
                               color = cols[1],enableMouseTracking=FALSE,marker=FALSE) %>%
    highcharter::hc_add_series(name = uniqtypes[2],
                               data=line2,type='line',#name=sprintf("regression line %s",type2),
                               color = cols[2],enableMouseTracking=FALSE,marker=FALSE)
  
  hc
}

#' scatter plot of metabolite-metabolite pairs (based on user selection)
#'
#' @import magrittr
#'
#' @param inputData IntLimObject output of ReadData() or FilterData()
#' @param stype category to color-code by
##' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
##' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param gene1Name string of select gene1Name
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @param gene2Name string of select gene2Name
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' PlotMMPair(mydata,stype="PBO_vs_Leukemia","arginine","(p-Hydroxyphenyl)lactic acid")
#'
#' }
#' @export
PlotGGPair<- function(inputData,stype=NULL,gene1Name,gene2Name,palette = "Set1",
                      viewer=T) {
  
  if(is.null(stype)) {
    stop("Users must define stype which defines the categories to be compared (e.g. tumor vs non-tumor).  This could be the same parameter that was used to run RunIntLim()")
  }
  if (length(palette) == 2) {
    cols <- c(palette)
  } else if (length(palette) == 1) {
    cols <- RColorBrewer::brewer.pal(3, palette)[1:2]
  } else {
    stop("palette must either be an RColorBrewer palette or a vector of hex colors of size 2")
  }
  
  if (class(inputData) != "MultiDataSet") {
    stop("input data is not a MultiDataSet class")
  }
  biobase_pdata <- NULL
  mytypes <- names(Biobase::assayData(inputData))
  if(any(mytypes == "metabolite")){
    incommon <- getCommon(inputData,stype)
    biobase_pdata <- Biobase::pData(inputData[["metabolite"]])
  }else{
    incommon <- formatSingleOmicInput(inputData,stype, type = "expression")
    biobase_pdata <- Biobase::pData(inputData[["expression"]])
  }
  
  if(is.null(stype)) {
    stop("A category to colorcode by (e.g. stype) must be provided")
  } else if (length(intersect(colnames(biobase_pdata),stype))!=1) {
    stop(paste0("You provided ",stype, "as your stype variable but it does not exist in your data"))
  } else {
    mytypes <- incommon$p
  }
  
  gene<-incommon$gene
  if(length(which(rownames(gene)==gene1Name))>0) {
    sGene1<-as.numeric(gene[gene1Name,])
  } else {
    stop(paste0("The gene ",gene1Name," was not found in your data"))
  }
  
  gene<-incommon$gene
  if(length(which(rownames(gene)==gene2Name))>0) {
    sGene2<-as.numeric(gene[gene2Name,])
  } else {
    stop(paste0("The gene ",gene2Name," was not found in your data"))
  }
  
  if(length(unique(mytypes))!=2) {
    stop(paste0("The group selected, '",stype,"', should only contain two different categories"))
  }
  
  mycols <- as.character(mytypes)
  mycols[which(mytypes==unique(mytypes)[1])] <- cols[1]
  mycols[which(mytypes==unique(mytypes)[2])] <- cols[2]
  
  data<-data.frame(x=sGene1,y=sGene2,z=colnames(gene),label=mytypes,color=mycols)
  
  # Get points to draw the lines for each phenotype by hand
  
  uniqtypes=as.character(unique(mytypes))
  
  # Starting with phenotype 1, get min and max x values constrained to the values of y
  # The reason we do this, is because the lines do not necessary need to go out to the max or min of x, particularly
  # when slopes are really steep (abline does this automatically but not highcharter)
  getLinePoints <- function(data,mytypes, uniqtypes, currenttype) {
    y=data$y[which(data$label==uniqtypes[currenttype])]; x=data$x[which(data$label==uniqtypes[currenttype])]
    min <- min(data$x[which(mytypes==uniqtypes[currenttype])])
    max <- max(data$x[which(mytypes==uniqtypes[currenttype])])
    
    m1<-stats::glm(y ~ x)
    line1<-data.frame(x=c(max,min),
                      y=c(stats::predict(m1,data.frame(x=c(max,min)))))
    return(data.frame(x=c(max,min), y=c(stats::predict(m1,data.frame(x=c(max,min))))))
  }
  
  line1 <- getLinePoints(data,mytypes,uniqtypes,currenttype=1)
  line2 <- getLinePoints(data,mytypes, uniqtypes, currenttype=2)
  
  ds <- highcharter::list_parse(data)
  
  hc <- highcharter::highchart(width = 350, height = 350 ) %>%
    highcharter::hc_title(text=paste(gene1Name,' vs. ', gene2Name, sep = '')) %>%
    highcharter::hc_xAxis(title=list(text=gene1Name)) %>%
    highcharter::hc_yAxis(title=list(text=gene2Name)) %>%
    hc_chart(zoomType = "xy") %>%
    highcharter::hc_add_series(data=ds,type="scatter",#col=cols[1],
                               tooltip = list(headerFormat="",
                                              pointFormat=paste("{point.label}","{point.z}")),
                               showInLegend=FALSE)
  
  
  hc <- hc %>%
    highcharter::hc_add_series(name = uniqtypes[1],
                               data=line1,type='line',#name=sprintf("regression line %s",type1),
                               color = cols[1],enableMouseTracking=FALSE,marker=FALSE) %>%
    highcharter::hc_add_series(name = uniqtypes[2],
                               data=line2,type='line',#name=sprintf("regression line %s",type2),
                               color = cols[2],enableMouseTracking=FALSE,marker=FALSE)
  
  hc
}


#' 'volcano' plot (difference in correlations vs p-values)
#' of all gene-metabolite pairs
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' @param nrpoints number of points to be plotted in lowest density areas (see 'smoothScatter' documentation for more detail)
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param diffcorr cutoff of differences in correlations for filtering (default 0.5)
#' @return a smoothScatter plot
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' pvalCorrVolcano(inputResults=myres,inputData=mydata)
#' }
#' @export
pvalCorrVolcano <- function(inputResults, inputData,nrpoints=10000,diffcorr=0.5,pvalcutoff=0.05){
   if (class(inputData) != "MultiDataSet") {
        stop("input data is not a MultiDataSet class")
	}
    if(class(inputResults) != "IntLimResults") {
	stop("input data is not a IntLim class")
    }
  incommon <- getCommon(inputData,inputResults@stype)
  p <- incommon$p
  
  if (length(unique(p)) !=2){
    stop(paste("pvalCorrVolcano is invalid for continuous outcomes and outcomes
               with more than two categories."))
  }
    volc.results <- IntLIM::ProcessResults(inputResults,  inputData, diffcorr = 0, pvalcutoff = 1)
    volc.table <- volc.results@filt.results
    Corrdiff <- volc.table[,4] - volc.table[,3]
    pval <- -log10(volc.table$FDRadjPval)
    graphics::smoothScatter(x = Corrdiff, pval, xlab = 'Difference in Correlation between Phenotypes',
		ylab = '-log10(FDR-adjusted p-value)', nrpoints=nrpoints,
                main = 'Volcano Plot')
    graphics::abline(h=-log10(pvalcutoff),lty=2,col="blue")
    graphics::abline(v=c(diffcorr,-diffcorr),lty=2,col="blue")
}

#' 'volcano' plot (difference in correlations vs p-values)
#' of all metabolite-metabolite pairs
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with metabolite abundance,
#' @param nrpoints number of points to be plotted in lowest density areas (see 'smoothScatter' documentation for more detail)
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param diffcorr cutoff of differences in correlations for filtering (default 0.5)
#' @return a smoothScatter plot
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' pvalCorrVolcano(inputResults=myres,inputData=mydata)
#' }
#' @export
pvalCorrVolcanoMetabolitePairs <- function(inputResults, inputData,nrpoints=10000,diffcorr=0.5,pvalcutoff=0.05){
  if (class(inputData) != "MultiDataSet") {
    stop("input data is not a MultiDataSet class")
  }
  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }
  mytypes <- names(Biobase::assayData(inputData))
  incommon <- NULL
  if(any(mytypes == "expression")){
    incommon <- getCommon(inputData,inputResults@stype)
  }else{
    incommon <- formatSingleOmicInput(inputData,inputResults@stype, type = "metabolite")
  }
  
  p <- incommon$p
  
  if (length(unique(p)) !=2){
    stop(paste("pvalCorrVolcano is invalid for continuous outcomes and outcomes
               with more than two categories."))
  }
  volc.results <- IntLIM::ProcessResults(inputResults,  inputData, diffcorr = 0, pvalcutoff = 1)
  volc.table <- volc.results@filt.results
  Corrdiff <- volc.table[,4] - volc.table[,3]
  pval <- -log10(volc.table$FDRadjPval)
  graphics::smoothScatter(x = Corrdiff, pval, xlab = 'Difference in Correlation between Phenotypes',
                          ylab = '-log10(FDR-adjusted p-value)', nrpoints=nrpoints,
                          main = 'Volcano Plot')
  graphics::abline(h=-log10(pvalcutoff),lty=2,col="blue")
  graphics::abline(v=c(diffcorr,-diffcorr),lty=2,col="blue")
}

#' 'volcano' plot (difference in correlations vs p-values)
#' of all gene-gene pairs
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with gene abundance,
#' @param nrpoints number of points to be plotted in lowest density areas (see 'smoothScatter' documentation for more detail)
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param diffcorr cutoff of differences in correlations for filtering (default 0.5)
#' @return a smoothScatter plot
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' pvalCorrVolcano(inputResults=myres,inputData=mydata)
#' }
#' @export
pvalCorrVolcanoGenePairs <- function(inputResults, inputData,nrpoints=10000,
                                     diffcorr=0.5,pvalcutoff=0.05){
  if (class(inputData) != "MultiDataSet") {
    stop("input data is not a MultiDataSet class")
  }
  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }
  mytypes <- names(Biobase::assayData(inputData))
  incommon <- NULL
  if(any(mytypes == "metabolite")){
    incommon <- getCommon(inputData,inputResults@stype)
  }else{
    incommon <- formatSingleOmicInput(inputData,inputResults@stype, type = "expression")
  }
  p <- incommon$p
  
  if (length(unique(p)) !=2){
    stop(paste("pvalCorrVolcano is invalid for continuous outcomes and outcomes
               with more than two categories."))
  }
  volc.results <- IntLIM::ProcessResults(inputResults,  inputData,
                                                        diffcorr = 0, pvalcutoff = 1)
  volc.table <- volc.results@filt.results
  Corrdiff <- volc.table[,4] - volc.table[,3]
  pval <- -log10(volc.table$FDRadjPval)
  graphics::smoothScatter(x = Corrdiff, pval, xlab = 'Difference in Correlation between Phenotypes',
                          ylab = '-log10(FDR-adjusted p-value)', nrpoints=nrpoints,
                          main = 'Volcano Plot')
  graphics::abline(h=-log10(pvalcutoff),lty=2,col="blue")
  graphics::abline(v=c(diffcorr,-diffcorr),lty=2,col="blue")
}

#' Graphs a scatterplot of gene-metabolite pairs vs. the interaction coefficient
#' for the gene-metabolite pair
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient 
#' (default bottom 10 percent (high negative coefficients) and top 10 percent 
#' (high positive coefficients))
#' @param percentageToPlot percentage of points to plot (the points will be 
#' randomly selected) -- plotting all points will likely overwhelm plotting function.
#' @param independent.var.type type of analyte used as the independent variable 
#' ("gene" or "metabolite")
#' @param outcome type of analyte used as the outcome/dependent variable ("gene"
#' or "metabolite")
#' @return a scatterplot
#'
#' @export
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' InteractionCoefficientGraph(inputResults=myres)
#' }
InteractionCoefficientGraph<-function(inputResults,
                                      interactionCoeffPercentile=0.10,
                                      percentageToPlot = 0.01, 
                                      independent.var.type = "gene",
                                      outcome = "metabolite"){


    if(class(inputResults) != "IntLimResults") {
      stop("input data is not a IntLim class")
    }

    #merge and properly name all data to return
    gene_metabolite_format_coeff = reshape2::melt(inputResults@interaction.coefficients)
    gene_metabolite_format_pval = reshape2::melt(inputResults@interaction.pvalues)
    gene_metabolite_format_adjp = reshape2::melt(inputResults@interaction.adj.pvalues)
    tofilter = cbind(gene_metabolite_format_coeff, gene_metabolite_format_pval$value, 
                     gene_metabolite_format_adjp$value)
    colnames(tofilter) = c("gene", "metab", "interaction_coeff", "Pval","FDRadjPval")


    #get top and bottom cutoffs (need highest positive and highest negative coeffs)
    first_half = getQuantileForInteractionCoefficient(tofilter$interaction_coeff, 
                                                      interactionCoeffPercentile)[1]
    second_half = getQuantileForInteractionCoefficient(tofilter$interaction_coeff, 
                                                       interactionCoeffPercentile)[2]


    toplot = data.frame(tofilter$interaction_coeff)
    colnames(toplot) = c("interaction_coeff")
    toplot$adjpval = tofilter$FDRadjPval
    toplot_sort = toplot[order(toplot$interaction_coeff),]
    colnames(toplot_sort) = c("interaction_coeff", "adjpval")
    toplot_sort$color = "black"
    toplot_sort$color[(toplot_sort$interaction_coeff > second_half | toplot_sort$
                         interaction_coeff <first_half)]="red"
    randomize = function(x) sample(1:nrow(toplot_sort),x,replace=F)
    random_rows_to_keep = sort(randomize(nrow(toplot_sort)*percentageToPlot))
    toplot_sort = toplot_sort[random_rows_to_keep,]
    if((independent.var.type == "gene" && outcome == "metabolite") || 
       (independent.var.type == "metabolite" && outcome == "gene"))
    {
      plot(1:length(toplot_sort$interaction_coeff),toplot_sort$interaction_coeff, 
           col=toplot_sort$color, xlab = "Gene Metabolite Pairs", ylab = 
             "Interaction Coefficient", pch=16)
    }else if(independent.var.type == "metabolite" && outcome == "metabolite"){
      toplot_sort = toplot_sort[which(!is.na(toplot_sort$interaction_coeff)),]
      plot(1:length(toplot_sort$interaction_coeff),toplot_sort$interaction_coeff, 
           col=toplot_sort$color, xlab = "Metabolite Pairs", ylab = 
             "Interaction Coefficient", pch=16)
    }
    else if(independent.var.type == "gene" && outcome == "gene"){
      toplot_sort = toplot_sort[which(!is.na(toplot_sort$interaction_coeff)),]
      plot(1:length(toplot_sort$interaction_coeff),toplot_sort$interaction_coeff, 
           col=toplot_sort$color, xlab = "Gene Pairs", ylab = 
             "Interaction Coefficient", pch=16)
    }
    else{
      stop("Error! outcome and independent.var.type must each be one of the following: gene, metabolite.")
    }
    
}


#' Creates a dataframe of the marginal effect of phenotype
#'
#' @import margins
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' @param metaboliteOfInterest metabolite in gene-metabolite pair
#' @param geneOfInterest gene in gene-metabolite pair
#' @param continuous whether or not the outcome is continuous (TRUE or FALSE)
#' @return dataframe for further analysis
#' @export
MarginalEffectsGraphDataframe<-function(inputResults, inputData, geneOfInterest, 
                                        metaboliteOfInterest, continuous){

  if (class(inputData) != "MultiDataSet") {
    stop("input data is not a MultiDataSet class")
  }
  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }

  #get covariates
  covariates = as.character(inputResults@covar$covariate)
  covariates_class = as.character(inputResults@covar$class.var)

  #get dataframes
  if(length(covariates) == 0){
    incommon <- getCommon(inputData,inputResults@stype,covar=covariates)
  }else{
    incommon <- getCommon(inputData,inputResults@stype,covar=covariates,class.covar=covariates_class)
  }
  pheno <- incommon$p
  gene <- incommon$gene
  metab <- incommon$metab

  #get one gene an metabolite
  gene_data = gene[geneOfInterest,]
  metab_data = metab[metaboliteOfInterest,]

  #Add gene, phenotype and metabolite data for glm
  forglm  = data.frame(row.names = 1:length(gene_data))
  forglm$g = gene_data
  forglm$type = pheno
  forglm$Y = as.numeric(metab_data)


  if (!is.null(covariates)) {

    #Add all covariates to dataframe for glm()
    i=3
    for(each in covariates){
      names = colnames(forglm)
      i = i+1
      forglm[,i] = incommon$covar_matrix[,each]
      colnames(forglm) = c(names, each)
    }
  }
  return(forglm)
}

#' Creates a dataframe of the marginal effect of phenotype
#'
#' @import margins
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' @param metaboliteOfInterest1 outcome metabolite in metabolite pair
#' @param metaboliteOfInterest2 independent metabolite in metabolite pair
#' @return dataframe for further analysis
#' @export
MarginalEffectsGraphDataframeMetabolitePairs<-function(inputResults, inputData, 
                                                       metaboliteOfInterest1, 
                                                       metaboliteOfInterest2){
  
  if (class(inputData) != "MultiDataSet") {
    stop("input data is not a MultiDataSet class")
  }
  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }
  
  #get covariates
  covariates = as.character(inputResults@covar$covariate)
  covariates_class = as.character(inputResults@covar$class.var)
  
  #get dataframes
  mytypes <- names(Biobase::assayData(inputData))
  incommon <- NULL
  if(length(covariates) == 0){
    if(any(mytypes == "expression")){
      incommon <- getCommon(inputData,inputResults@stype,covar=covariates)
    } else {
      incommon <- formatSingleOmicInput(inputData,inputResults@stype,covar=covariates, 
                                        type = "metabolite")
    }
    
  }else{
    if(any(mytypes == "expression")){
      incommon <- getCommon(inputData,inputResults@stype,covar=covariates,class.covar=covariates_class)
    } else {
      incommon <- formatSingleOmicInput(inputData,inputResults@stype,covar=covariates,
                                        class.covar=covariates_class, type = "metabolite")
    }
  }
  pheno <- incommon$p
  metab <- incommon$metab
  
  #get two metabolites
  metab1_data = metab[metaboliteOfInterest1,]
  metab2_data = metab[metaboliteOfInterest2,]
  
  #Add phenotype and metabolite data for glm
  forglm  = data.frame(row.names = 1:length(metab1_data))
  forglm$g = as.numeric(metab1_data)
  forglm$type = pheno
  forglm$Y = as.numeric(metab2_data)
  
  
  if (!is.null(covariates)) {
    
    #Add all covariates to dataframe for glm()
    i=3
    for(each in covariates){
      names = colnames(forglm)
      i = i+1
      forglm[,i] = incommon$covar_matrix[,each]
      colnames(forglm) = c(names, each)
    }
  }
  return(forglm)
}

#' Creates a dataframe of the marginal effect of phenotype
#'
#' @import margins
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' @param geneOfInterest1 outcome gene in gene pair
#' @param geneOfInterest2 independent gene in gene pair
#' @return dataframe for further analysis
#' @export
MarginalEffectsGraphDataframeGenePairs<-function(inputResults, inputData, 
                                                 geneOfInterest1, 
                                                 geneOfInterest2){
  
  if (class(inputData) != "MultiDataSet") {
    stop("input data is not a MultiDataSet class")
  }
  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }
  
  #get covariates
  covariates = as.character(inputResults@covar$covariate)
  covariates_class = as.character(inputResults@covar$class.var)
  
  #get dataframes
  mytypes <- names(Biobase::assayData(inputData))
  if(length(covariates) == 0){
    if(any(mytypes == "metabolite")){
      incommon <- getCommon(inputData,inputResults@stype,covar=covariates)
    } else {
      incommon <- formatSingleOmicInput(inputData,inputResults@stype,covar=covariates,
                                        type = "expression")
    }
  }else{
    if(any(mytypes == "metabolite")){
      incommon <- getCommon(inputData,inputResults@stype,covar=covariates)
    } else {
      incommon <- formatSingleOmicInput(inputData,inputResults@stype,covar=covariates,
                                        class.covar=covariates_class, type = "expression")
    }
  }
  pheno <- incommon$p
  gene <- incommon$gene
  
  #get two metabolites
  gene1_data = gene[geneOfInterest1,]
  gene2_data = gene[geneOfInterest2,]
  
  #Add phenotype and metabolite data for glm
  forglm  = data.frame(row.names = 1:length(gene1_data))
  forglm$g = as.numeric(gene1_data)
  forglm$type = pheno
  forglm$Y = as.numeric(gene2_data)
  
  
  if (!is.null(covariates)) {
    
    #Add all covariates to dataframe for glm()
    i=3
    for(each in covariates){
      names = colnames(forglm)
      i = i+1
      forglm[,i] = incommon$covar_matrix[,each]
      colnames(forglm) = c(names, each)
    }
  }
  return(forglm)
}

#' Creates a dataframe of the marginal effect of phenotype
#' @param dataframe from MarginalEffectsGraphDataframe
#' @param title for graph
#' @return values used for graphing
#' @export
MarginalEffectsGraph<-function(dataframe, title){
  form = "Y ~ g + type + g:type"
  if (ncol(dataframe) > 3) {

    covariates = colnames(dataframe)[4:ncol(dataframe)]
    #Add all covariates to formula for glm()
    for(i in 1:length(covariates)){
      form <- paste(form, '+', covariates[i])
    }
  }
  model = stats::glm(formula = form, data=dataframe)
  tryCatch({
    cplot(model, "type", data = dataframe, what = "prediction", main = title)
  }, error = function(cond){
    print("Could not plot the data. Check to see whether your outcome is continuous.
          Only categorical outcomes are valid for this function.")
  })
  return(model)

}


#' histogram of gene-metabolite pairs
#' depending upon metabolite or gene
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim() and ProcessResults())
#' @param type 'metabolite' or 'gene'.  'metabolite' set as default
#' @param breaks Number of breaks selected for histogram
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(inputResults=myres,inputData=mydata)
#' HistogramGMPairs(inputResults=myres)
#' }
#' @export
HistogramGMPairs <- function(inputResults, type = 'metabolite', breaks = 50){

  x <- inputResults@filt.results
  
  if(is.null(x)){
      stop('Please run ProcessResults() before inputting into HistogramGMPairs')
  }
  if (type == 'metabolite'){
    metab.pairs <- data.frame(table(x$Metabolite))
    metab.pairs.number <- as.vector(metab.pairs$Freq)
    hist(metab.pairs.number, breaks = breaks, main = "Number of gene-metabolite 
         pairs based on metabolite", xlab = 'Gene-metabolite pairs based on metabolite')
  }else if (type == 'gene'){
    gene.pairs <- data.frame(table(x$Gene))
    gene.pairs.number <- as.vector(gene.pairs$Freq)
    str(x)
    str(x$Gene)
    hist(gene.pairs.number, main = "Number of gene-metabolite pairs based on gene", 
         breaks = breaks, xlab = 'Gene-metabolite pairs based on gene')
  }else{
      stop("Only two valid types:  gene or metabolite.  Invalid type entered")
  }
}

#' histogram of metabolite-metabolite pairs
#' depending upon metabolite
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim() and ProcessResults())
#' @param type 'outcome' or 'independent'.  'outcome' set as default
#' @param breaks Number of breaks selected for histogram
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(inputResults=myres,inputData=mydata)
#' HistogramGMPairs(inputResults=myres)
#' }
#' @export
HistogramMMPairs <- function(inputResults, type = 'outcome', breaks = 50){
  
  x <- inputResults@filt.results
  if(is.null(x)){
    stop('Please run ProcessResults() before inputting into HistogramMMPairs')
  }
  if (type == 'outcome'){
    metab1.pairs <- data.frame(table(x$Metab1))
    metab1.pairs.number <- as.vector(metab1.pairs$Freq)
    hist(metab1.pairs.number, breaks = breaks, 
         main = "Number of metabolite pairs based on outcome metabolite", 
         xlab = 'Metabolite pairs based on outcome metabolite')
  }else if (type == 'independent'){
    
    metab2.pairs <- data.frame(table(x$Metab2))
    metab2.pairs.number <- as.vector(metab2.pairs$Freq)
    hist(metab2.pairs.number, 
         main = "Number of metabolite pairs based on independent metabolite", 
         breaks = breaks, xlab = 'Metabolite pairs based on independent metabolite')
  }else{
    stop("Only two valid types:  outcome or independent.  Invalid type entered")
  }
}

#' histogram of gene-gene pairs
#' depending upon gene
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim() 
#' and ProcessResults())
#' @param type 'outcome' or 'independent'.  'outcome' set as default
#' @param breaks Number of breaks selected for histogram
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(inputResults=myres,inputData=mydata)
#' HistogramGMPairs(inputResults=myres)
#' }
#' @export
HistogramGGPairs <- function(inputResults, type = 'outcome', breaks = 50){
  
  x <- inputResults@filt.results
  if(is.null(x)){
    stop('Please run ProcessResults() before inputting into HistogramGGPairs')
  }
  if (type == 'outcome'){
    gene1.pairs <- data.frame(table(x$Gene1))
    gene1.pairs.number <- as.vector(gene1.pairs$Freq)
    hist(gene1.pairs.number, breaks = breaks, 
         main = "Number of gene pairs based on outcome gene", 
         xlab = 'Gene pairs based on outcome gene')
  }else if (type == 'independent'){
    
    gene2.pairs <- data.frame(table(x$Gene2))
    gene2.pairs.number <- as.vector(gene2.pairs$Freq)
    hist(gene2.pairs.number, 
         main = "Number of gene pairs based on independent gene", 
         breaks = breaks, xlab = 'Gene pairs based on independent gene')
  }else{
    stop("Only two valid types:  outcome or independent.  Invalid type entered")
  }
}

#' For each sample, save a plot of phenotype predictions from each edge in the
#' IntLIM graph. The edge weight is the phenotype prediction
#' using the pair of analytes connected to the edge. Weights correspond to color
#' in the graph. A color bar is shown for reference, and the true phenotype is
#' listed at the top of the plot.
#' @param graphWithPredictions An igraph object. This graph is the co-regulation graph
#' generated using IntLIM analysis of analyte pairs. Weights correspond to phenotype
#' predictions.
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param stype The phenotype of interest. This should correspond to a column in the
#' input data.
#' @param dirName The name of the directory where the output images will be saved.
#' @param continuous A boolean indicating whether the phenotype is continuous (TRUE)
#' or discrete (FALSE).
#' @export
SaveGraphPredictionPlots <- function(graphWithPredictions, inputData, stype, dirName,
                                     continuous = TRUE){
  # Create folder.
  dir.create(dirName)
  
  # Save each graph.
  for(name in names(graphWithPredictions)){
    
    # Connect to the file.
    grDevices::png(paste0(paste(dirName, make.names(name), sep = "\\"), ".png"))
    
    # Extract graph for the subject of interest.
    g <- graphWithPredictions[[name]]
    
    # Set up the layout and margins.
    graphics::layout(t(1:2),widths=c(5.5,1.5))
    par(mar=c(0,0,2,3))
    
    # Plot the graph.
    plot(g, layout = igraph::layout.fruchterman.reingold, vertex.label = NA)
    
    # Add the true phenotype and subject ID.
    inputDat <- inputData@phenoData$expression$main@data
    true_phen <- inputDat[which(rownames(inputDat) == name), stype]
    graphics::text(x = -1, y = 1.2, paste0(name, " (true phenotype is ", true_phen, ")"), pos = 4)
    
    # Add information about conversion to factors for discrete data.
    if(continuous == FALSE){
      inputDataPhen <- as.factor(inputDat[,stype])
      graphics::text(x = -1, y = 1.1, paste(levels(inputDataPhen)[1], " <= 0"), pos = 4)
      graphics::text(x = -1, y = 1.0, paste(levels(inputDataPhen)[2], " >= 1"), pos = 4)
    }
    
    # Add the color bar.
    color <- igraph::edge_attr(g, name = "color")[order(igraph::edge_attr(g, name = "weight"))]
    labs <- igraph::edge_attr(g, name = "weight")[order(igraph::edge_attr(g, name = "weight"))]
    lab_quants <- seq(min(labs), max(labs), by = (max(labs)-min(labs))/5)
    graphics::image(y=1:100,z=t(1:100), col=color, axes=FALSE, main="Prediction", cex.main=.8)
    graphics::axis(4,cex.axis=0.8, at = seq(0, 100, by = 20), labels = format(as.list(lab_quants), 
                                                                    digits=0, 
                                                                    scientific=FALSE), 
         las = 1)
    
    # Close the file connection.
    grDevices::dev.off()
  }
}

#' Plot the graph with edges colored by weight in the final outcome.
#' @param graph The co-regulation graph.
#' @param results A modelResults object.
#' @export
PlotGraphWeights <- function(graph, results){
  
  # Set up the layout and margins.
  graphics::layout(t(1:2),widths=c(5.5,1.5))
  par(mar=c(0,0,2,3))
  
  # Match the weights to graph edges.
  g <- igraph::as_data_frame(graph)
  g$to <- make.names(g$to)
  g$from <- make.names(g$from)
  weights <- results@current.weights
  if(results@weights.after.pooling == TRUE){
    S <- results@pooling.filter@filter
    weights <- t(matrix(rep(weights, dim(S)[1]), ncol = dim(S)[1]))
    sum_S <- colSums(S)
    S.weighted <- S * weights / sum_S
    S.flat <- rowSums(S.weighted)
    weights <- S.flat
  }
  names(weights) <- rownames(results@model.input@node.wise.prediction)
  weights_by_edge_name <- lapply(1:dim(g)[1], function(edge){
    forwards <- paste(g$from[edge], g$to[edge], sep = "__")
    backwards <- paste(g$to[edge], g$from[edge], sep = "__")
    which_weight <- union(which(names(weights) == forwards), 
                          which(names(weights) == backwards))
    the_weight <- NA
    if(length(which_weight) > 0){
      the_weight <- weights[which_weight]
    }
    return(the_weight)
  })
  
  # Add the weights to the data frame.
  g$weight <- unlist(weights_by_edge_name)
  g <- g[which(!is.na(g$weight)),]
  
  # Map weights to colors.
  pal <- grDevices::colorRampPalette(c("blue", "red"))(100)
  range_weight <- range(g$weight)
  color_scale <- pal[findInterval(g$weight, seq(range_weight[1], range_weight[2], 
                                            length.out = length(pal)+1), all.inside = TRUE)]
  g$color <- color_scale

  # Plot the graph.
  new_graph <- igraph::graph_from_data_frame(g, directed = FALSE)
  plot(new_graph, layout = igraph::layout.fruchterman.reingold, vertex.label = NA,
       vertex.size = 3)
  
  # Add the color bar.
  color <- igraph::edge_attr(new_graph, name = "color")[order(igraph::edge_attr(new_graph, 
                                                                                name = "weight"))]
  labs <- igraph::edge_attr(new_graph, name = "weight")[order(igraph::edge_attr(new_graph, 
                                                                                name = "weight"))]
  lab_quants <- seq(min(labs), max(labs), by = (max(labs)-min(labs))/5)
  graphics::image(y=1:100,z=t(1:100), col=color, axes=FALSE, main="Weight", cex.main=.8)
  graphics::axis(4,cex.axis=0.8, at = seq(0, 100, by = 20), labels = format(as.list(lab_quants), 
                                                                            digits=0, 
                                                                            scientific=FALSE), 
                 las = 1)
}

#' Plot the graph as a heatmap with edges colored by weight in the final outcome.
#' @param graph The co-regulation graph.
#' @param results A modelResults object.
#' @export
PlotGraphWeightsHeatmap <- function(graph, results){
  # Match the weights to graph edges.
  g <- igraph::as_data_frame(graph)
  g$to <- make.names(g$to)
  g$from <- make.names(g$from)
  weights <- results@current.weights
  if(results@weights.after.pooling == TRUE){
    S <- results@pooling.filter@filter
    weights <- t(matrix(rep(weights, dim(S)[1]), ncol = dim(S)[1]))
    sum_S <- colSums(S)
    S.weighted <- S * weights / sum_S
    S.flat <- rowSums(S.weighted)
    weights <- S.flat
  }
  names(weights) <- rownames(results@model.input@node.wise.prediction)
  weights_by_edge_name <- lapply(1:dim(g)[1], function(edge){
    forwards <- paste(g$from[edge], g$to[edge], sep = "__")
    backwards <- paste(g$to[edge], g$from[edge], sep = "__")
    which_weight <- union(which(names(weights) == forwards), 
                          which(names(weights) == backwards))
    the_weight <- NA
    if(length(which_weight) > 0){
      the_weight <- weights[which_weight]
    }
    return(the_weight)
  })
  
  # Add the weights to the data frame.
  g$weight <- unlist(weights_by_edge_name)
  g <- g[which(!is.na(g$weight)),]
  
  # Map weights to colors.
  pal <- grDevices::colorRampPalette(c("blue", "red"))(100)
  range_weight <- range(g$weight)
  color_scale <- pal[findInterval(g$weight, seq(range_weight[1], range_weight[2], 
                                                length.out = length(pal)+1), all.inside = TRUE)]
  g$color <- color_scale
  
  # Create an adjacency matrix.
  new_graph <- igraph::graph_from_data_frame(g, directed = FALSE)
  heat <- igraph::as_adjacency_matrix(new_graph, sparse = FALSE, attr = "weight")
  
  # Plot the heatmap.
  heatmap(heat)
}
