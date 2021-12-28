#' Get some stats after reading in data
#'
#' @param inputData IntLimObject output of ReadData()
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @return a highcharter object
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
  if(length(inputData@analyteType1)>0){
    type1Data <- inputData@analyteType1
    toplot <- suppressMessages(reshape2::melt(type1Data))
    df <- dplyr::tibble(value = toplot$value, by = toplot$Var2)
    df <- dplyr::group_by_at(df, "by")
    df <- dplyr::do(df, data = grDevices::boxplot.stats(.$value))
    bxps <- purrr::map(df$data, "stats")
    outs <- purrr::map2_df(seq(nrow(df)), df$data, function(x, y) {
      if (length(y$out) > 0)
        d <- dplyr::tibble(x = x - 1, y = y$out)
      else d <- dplyr::tibble()
      d
    })
    outs <- data.frame(outs, 'z' = colnames(type1Data)[outs$x + 1])
    z <- outs$z
    # To try to get the analyte names of outliers, would have to go back and get 
    # the analyte names from original data frame and put htem in outs$color
    
    g <- highcharter::highchart(width = 750, height = 750 )
    g <- highcharter::hc_title(g, text = "Analyte Type 1 Levels",
                            style = list(color = '#2E1717',
                                         fontWeight = 'bold', fontSize = "20px"))
    g <- highcharter::hc_plotOptions(g, boxplot = boxplotOptions)
    g <- highcharter::hc_add_series(g, data = bxps,name = "Analyte Type 1 Levels", type="boxplot",
                    color=cols[1],showInLegend=FALSE)
    g <- highcharter::hc_add_series(g, data=highcharter::list_parse(outs),name = "Analyte Type 1 Levels",
                                 type="scatter",color=cols[1],showInLegend=FALSE,
                                 tooltip = list(headerFormat = "",
                                                pointFormat = "{point.z} <br/> {point.y}",
                                                showInLegend = FALSE))
    g <- highcharter::hc_yAxis(g, title = list(text = "Levels",
                                         style = list(fontSize = "13px")),
                            labels = list(format = "{value}"))
    g <- highcharter::hc_xAxis(g, labels="", categories = colnames(type1Data))
    g <- highcharter::hc_tooltip(g, valueDecimals = 2)
    g <- highcharter::hc_exporting(g, enabled = TRUE)
  }
  if(length(inputData@analyteType2)>0){
    type2Data <- inputData@analyteType2
    toplot <- suppressMessages(reshape2::melt(t(type2Data)))
    df <- dplyr::data_frame(value = toplot$value, by = toplot$Var1)
    df <- dplyr::group_by_at(df, "by")
    df <- dplyr::do(df, data = grDevices::boxplot.stats(.$value))
    bxps <- purrr::map(df$data, "stats")
    outs <- purrr::map2_df(seq(nrow(df)), df$data, function(x, y) {
      if (length(y$out) > 0)
        d <- dplyr::data_frame(x = x - 1, y = y$out)
      else d <- dplyr::data_frame()
      d
    })
    outs <- data.frame(outs, 'z' = colnames(type2Data)[outs$x + 1])
    z <- outs$z
    
    m <- highcharter::highchart(width = 750, height = 750 )
    m <- highcharter::hc_title(m, text = "Analyte Type 2 Levels",
                            style = list(color = '#2E1717',
                                         fontWeight = 'bold', fontSize = "20px"))
    m <- highcharter::hc_plotOptions(m, boxplot = boxplotOptions)
    m <- highcharter::hc_add_series(m, data = bxps,name = "Analyte Type 2 Levels",
                                 type="boxplot",color=cols[2],showInLegend=FALSE)
    m <- highcharter::hc_add_series(m, data=highcharter::list_parse(outs),name = "Analyte Type 2 Levels",
                                 type="scatter",color=cols[2],showInLegend=FALSE,
                                 tooltip = list(headerFormat = "", 
                                                pointFormat = "{point.z} <br/> {point.y}",
                                                showInLegend = FALSE))
    m <- highcharter::hc_yAxis(m, title = list(text = "Levels",
                                         style = list(fontSize = "13px")),
                            labels = list(format = "{value}"))
    m <- highcharter::hc_xAxis(m, labels="", categories = colnames(type2Data))
    m <- highcharter::hc_tooltip(m, valueDecimals = 2)
    m <- highcharter::hc_exporting(m, enabled = TRUE)
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
#' @param inputData IntLimObject output of ReadData()
#' @param stype category to color-code by (can be more than two categories)
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @return a highcharter object
#' @export
PlotPCA <- function(inputData,viewer=T,stype="",palette = "Set1") {
  
  if(is.numeric(inputData@sampleMetaData[,stype]) == TRUE) {
    warning("The resulting PCA plot is not color-coded because you did not provide 
		        a categorical variable in 'stype'")
    mytype <- NULL
  } 
  else {
    mytype <- as.character(inputData@sampleMetaData[,stype])
    numcateg <- length(unique(mytype))
    if(length(palette) >= 2) {
      cols <- palette
    } 
    else {
      if(numcateg == 1) {
        if(length(palette)==1) {
          cols <- RColorBrewer::brewer.pal(3, palette)[1]
        }
        else {
          stop("palette should be an RColorBrewer palette or a vector of colors")
        }
      } 
      else if (numcateg == 2) {
        if(length(palette)==1) {
          cols <- RColorBrewer::brewer.pal(3, palette)[1:2]
        } 
        else {
          stop("palette should be an RColorBrewer palette or a vector of colors")
        }
      } 
      else if (numcateg > 2) {
        if(length(palette)==1) {
          cols <- RColorBrewer::brewer.pal(numcateg, palette)
        } 
        else {
          stop("palette should be an RColorBrewer palette or a vector of colors")
        }
      } 
      else {
        stop("There are no values in your 'stype' column")
      }
    }
  }
  p <- NULL
  pg <- NULL
  pm <- NULL
  
  if(length(inputData@analyteType1)>0 && length(inputData@analyteType2)>0){
    if(is.numeric(inputData@sampleMetaData[,stype]) == TRUE) {
      mpca <- stats::prcomp(t(inputData@analyteType1),center=T,scale=F)
      gpca <- stats::prcomp(t(inputData@analyteType2),center=T,scale=F)
      gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),color=rep("blue",nrow(gpca$x)))
      mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),color=rep("blue",nrow(mpca$x)))
      gds <- highcharter::list_parse(gtoplot)
      pg <- highcharter::highchart(width = 350, height = 350 )
      pg <- highcharter::hc_add_series(pg, data=gds,type="scatter",
                                              tooltip = list(headerFormat="",
                                                             pointFormat=paste("{point.label}","{point.z}")),
                                              showInLegend=FALSE)
      mds <- highcharter::list_parse(mtoplot)
      pm <- highcharter::highchart(width = 350, height = 350)
      pm <- highcharter::hc_add_series(pm, data=mds,type="scatter",
                                              tooltip = list(headerFormat="",
                                                             pointFormat=paste("{point.label}","{point.z}")),
                                              showInLegend=FALSE)
    } else {
      type1 <- inputData@analyteType1
      type2 <- inputData@analyteType2
      alltype <- inputData@sampleMetaData[,stype]
      uniqtypes <- unique(alltype)
      mycols <- as.character(alltype)
      for (i in 1:numcateg) {
        mycols[which(alltype==uniqtypes[i])] <- cols[i]
      }
      gpca <- stats::prcomp(t(type1),center=T,scale=F)
      mpca <- stats::prcomp(t(type2),center=T,scale=F)
      gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),label=alltype,color=mycols)
      mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),label=alltype,color=mycols)
      mds <- highcharter::list_parse(mtoplot)
      gds <- highcharter::list_parse(gtoplot)
      pg <- highcharter::highchart(width = 350, height = 350)
      pm <- highcharter::highchart(width = 350, height = 350)
      for (i in 1:length(uniqtypes)) {
        mytype <- unique(alltype)[i]
        gds <- highcharter::list_parse(gtoplot[which(gtoplot$label==mytype),])
        pg <- highcharter::hc_add_series(pg, data=gds,type="scatter",
                                                name=mytype,
                                         color=cols[which(alltype==mytype)[1]],
                                         tooltip = list(headerFormat="",
                                                        showInLegend=TRUE))
        mds <- highcharter::list_parse(mtoplot[which(mtoplot$label==mytype),])
        pm <- highcharter::hc_add_series(pm, data=mds,type="scatter",
                                                name=mytype,
                                                color=cols[which(alltype==mytype)[1]],
                                         tooltip = list(headerFormat="",
                                                        pointFormat=paste("{point.label}","{point.z}")),
                                                showInLegend=TRUE)
      }
    }
  }
  
  # Set up plots.
  if(length(inputData@analyteType1)>0){
    mpercvar=round((mpca$sdev)^2 / sum(mpca$sdev^2)*100,2)
    pm <- highcharter::hc_title(pm, text="PCA of analyte type 2")
    pm <- highcharter::hc_xAxis(pm, title=list(text=paste0("PC1:",round(mpercvar[1],1),"%")))
    pm <- highcharter::hc_yAxis(pm, title=list(text=paste0("PC2:",round(mpercvar[2],2),"%")))
    pm <- highcharter::hc_chart(pm, zoomType = "xy")
  }
  if(length(inputData@analyteType2)>0){
    gpercvar=round((gpca$sdev)^2 / sum(gpca$sdev^2)*100,2)
    pg <- highcharter::hc_title(pg, text="PCA of analyte type 1")
    pg <- highcharter::hc_xAxis(pg, title=list(text=paste0("PC1:",round(gpercvar[1],1),"%")))
    pg <- highcharter::hc_yAxis(pg, title=list(text=paste0("PC2:",round(gpercvar[2],2),"%")))
    pg <- highcharter::hc_chart(pg, zoomType = "xy")
  }
  p <- NULL
  if(length(inputData@analyteType1)>0 && length(inputData@analyteType2)>0){
    if (viewer == TRUE) {
      p <-htmltools::browsable(highcharter::hw_grid(pg, pm, ncol = 2, rowheight = 550))
    } else {
      p <- highcharter::hw_grid(pg, pm)
    }
  } else if(length(inputData@analyteType1)>0){
    if (viewer == TRUE) {
      p <-htmltools::browsable(highcharter::hw_grid(pg, ncol = 1, rowheight = 550))
    } else {
      p <- highcharter::hw_grid(pg)
    }
  } else if(length(inputData@analyteType2)>0){
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
#' @param adjusted Whether or not to plot adjusted p-values. If TRUE (default),
#' adjusted p-values are plotted. If FALSE, unadjusted p-values are plotted.
#' @export
DistPvalues<- function(IntLimResults,breaks=100,adjusted = TRUE) {

  if(adjusted == FALSE){
    hist(IntLimResults@interaction.pvalues,breaks=breaks,
         main="Histogram of Interaction P-values")
  }else{
    hist(IntLimResults@interaction.adj.pvalues,breaks=breaks,
         main="Histogram of Adjusted Interaction P-values")
  }

}

#' Visualize the distribution of unadjusted p-values for all covariates
#' from linear models using a bar chart.
#'
#' @include IntLimResults_extendedfunctions.R
#'
#' @param IntLimResults output of RunIntLim()
#' @export
PValueBoxPlots<- function(IntLimResults) {
  if(length(IntLimResults@covariate.pvalues) == 0){
    print("Error! You must set save.covar.pvals to TRUE when running IntLIM to run PValueBoxPlots")
  }else{
    graphics::par(mar=c(8, 4.1, 4.1, 2.1))
    graphics::boxplot(IntLimResults@covariate.pvalues, las = 3, ylim = c(0,1), ylab = "P-Value")
  }
}

#' Visualize the distribution of unadjusted p-values from linear models
#'
#' @include IntLimResults_extendedfunctions.R
#'
#' @param IntLimResults output of RunIntLim()
#' @param breaks the number of breaks to use in histogram (see hist() documentation for more details)
#' @export
DistRSquared<- function(IntLimResults,breaks=100) {
  
  hist(IntLimResults@model.rsquared,breaks=breaks,
       main="Histogram of Interaction R-Squared Values")
}

#' Returns the clusters found using CorrHeatmap.
#'
#' @param inputResults Data frame (output of ProcessResults())
#' @param treecuts number of clusters (of pairs) to cut the tree 
#' into for color-coding
#' @return A data frame including the independent and outcome analytes in each
#' for each pair and the cluster to which that pair belongs.
#' @export
GetCorrClusters <- function(inputResults,treecuts=2) {
  type <- cor <- c()
  clusterdata <- NULL
  
  if(is.null(inputResults)){
    stop('Please run ProcessResults() before inputting into HistogramPairs')
  }
  
  # Stop if not two discrete phenotypes.
  if(colnames(inputResults)[3] == "interaction_coeff"){
    stop("GetCorrClusters requires 2 discrete phenotypes. Do not run with continuous phenotypes.")
  }
  else{
    allres <- inputResults
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
    
    # hclust is used under the hood of heatmaply.
    dend <- stats::hclust(stats::dist(heat_data, method = "euclidean"))
    clusters <- stats::cutree(dend, k = treecuts)
    
    # Split clusters.
    analyte1 <- unlist(lapply(names(clusters), function(c){
      return(strsplit(c, " vs ")[[1]][1])
    }))
    analyte2 <- unlist(lapply(names(clusters), function(c){
      return(strsplit(c, " vs ")[[1]][2])
    }))
    
    # Formulate cluster data.
    clusterdata <- data.frame(IndependentAnalyte = analyte1, OutcomeAnalyte = analyte2,
                           Cluster = clusters)
    
    # Add average correlation for each phenotype.
    for(phen in colnames(heat_data)){
      clusterdata[,paste0(phen, "Corr")] <- heat_data[,phen]
    }
  }
  return(clusterdata)
}
#' Plot correlation heatmap
#'
#' @param inputResults Data frame (output of ProcessResults())
#' @param top_pairs cutoff of the top pairs, sorted by adjusted p-values, to be 
#' plotted (plotting more than 1200 can take some time) (default: 1200)
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @param treecuts number of clusters (of pairs) to cut the tree into for color-coding
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors. "Set1" is the default.
#' @return a highcharter object
#'@param static allows user to decide whether heatmap is interactive or static
#'@param html.file allows user to specify file path to output heatmap onto (used for non-static heatmaply objects)
#'@param pdf.file allows user to specify file path to output heatmap onto (used for static heatmap.2 objects)
#' @export
CorrHeatmap <- function(inputResults,viewer=T,top_pairs=1200,treecuts=2, 
                        palette = "Set1", static = FALSE,
                        html.file=NA, pdf.file=NA) {
  
  # Stop if not two discrete phenotypes.
  if(colnames(inputResults)[3] == "interaction_coeff"){
    stop("CorrHeatmap requires 2 discrete phenotypes. Do not run with continuous phenotypes.")
  }
  type <- cor <- c()

	if(nrow(inputResults)==0) {
		stop("Make sure you run ProcessResults before making the heatmap")
	}
  else{
    allres <- inputResults
    if(nrow(allres)>top_pairs) {
      allp <- inputResults[,"FDRadjPval"]
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
    
    if(static == FALSE){
      if(is.na(pdf.file) && is.na(html.file)){
        hm <- heatmaply::heatmaply(heat_data,main = "Correlation heatmap",
                                   k_row = treecuts,#k_col = 2,
                                   margins = c(80,5),
                                   dendrogram = "row",
                                   y_axis_font_size ="1px",
                                   colors = palette,
                                   key.title = 'Correlation \n differences')
        hm
      }
      else if(!is.na(pdf.file)){
        
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
      if(is.na(pdf.file) && is.na(html.file)){
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
      }
      else{  
        if(!is.na(pdf.file)){
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
        if(!is.na(html.file) & static==TRUE){
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
  }
}

#' scatter plot of pairs (based on user selection)
#'
#' @param inputData IntLimObject output of ReadData() or FilterData()
#' @param stype category to color-code by
##' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
##' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @param outcomeAnalyteOfInterest outcome analyte in pair
#' @param independentAnalyteOfInterest independent analyte in pair
#' @param outcome '1' or '2' must be set as outcome/independent variable
#' @param independentVariable '1' or '2' must be set as outcome/independent variable
#' @export
PlotPair<- function(inputData,stype,outcome,independentVariable, independentAnalyteOfInterest, 
                    outcomeAnalyteOfInterest, palette = "Set1",	viewer=T) {
  # Convert names.
  name_outcomeAnalyteOfInterest <- make.names(outcomeAnalyteOfInterest)
  name_independentAnalyteOfInterest <- make.names(independentAnalyteOfInterest)
  stype <- make.names(stype)
  
  if(is.null(stype)) {
    stop("Users must define stype which defines the categories to be compared 
         (e.g. tumor vs non-tumor).  This could be the same parameter that was 
         used to run RunIntLim()")
  }
  if (length(palette) == 2) {
    cols <- c(palette)
  } else if (length(palette) == 1) {
    cols <- RColorBrewer::brewer.pal(3, palette)[1:2]
  } else {
    stop("palette must either be an RColorBrewer palette or a vector of hex colors of size 2")
  }
  
  # Extract the outcome and independent data.
  outcomeData <- NULL
  independentData <- NULL
  sOutcome <- NULL
  sIndependent <- NULL
  if(outcome == 1){
    outcomeData <- inputData@analyteType1
  }else if(outcome == 2){
    outcomeData <- inputData@analyteType2
  }
  if(independentVariable == 1){
    independentData <- inputData@analyteType1
  }else if(independentVariable == 2){
    independentData <- inputData@analyteType2
  }
  
  # Check that analytes of interest are found in data.
  if(length(which(rownames(outcomeData)==name_outcomeAnalyteOfInterest))>0) {
    sOutcome<-as.numeric(outcomeData[name_outcomeAnalyteOfInterest,])
  } else {
    stop(paste0("The analyte ",outcomeAnalyteOfInterest," was not found in your data"))
  }
  if(length(which(rownames(independentData)==name_independentAnalyteOfInterest))>0) {
    sIndependent<-as.numeric(independentData[name_independentAnalyteOfInterest,])
  } else {
    stop(paste0("The analyte ",independentAnalyteOfInterest," was not found in your data"))
  }
  
  # Check that data only contains two categories.
  if(length(unique(inputData@sampleMetaData[,stype]))!=2) {
    stop(paste0("The group selected, '",stype,"', should only contain two different categories"))
  }
  
  # Set up data.
  mycols <- as.character(inputData@sampleMetaData[,stype])
  mycols[which(inputData@sampleMetaData[,stype]==
                 unique(inputData@sampleMetaData[,stype])[1])] <- cols[1]
  mycols[which(inputData@sampleMetaData[,stype]==
                 unique(inputData@sampleMetaData[,stype])[2])] <- cols[2]
  
  data<-data.frame(x=sIndependent,y=sOutcome,z=colnames(independentData),
                   label=inputData@sampleMetaData[,stype],color=mycols)

  # Get points to draw the lines for each phenotype by hand
  
  uniqtypes=as.character(unique(inputData@sampleMetaData[,stype]))
  
  # Starting with phenotype 1, get min and max x values constrained to the values of y
  # The reason we do this, is because the lines do not necessary need to go out to the 
  # max or min of x, particularly
  # when slopes are really steep (abline does this automatically but not highcharter)
  mytypes <- inputData@sampleMetaData[,stype]
  getLinePoints <- function(data,mytypes, uniqtypes, currenttype) {
    y=data$y[which(data$label==uniqtypes[currenttype])]; 
    x=data$x[which(data$label==uniqtypes[currenttype])]
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
  
  hc <- highcharter::highchart(width = 350, height = 350 )
  hc <- highcharter::hc_title(hc, text=paste(independentAnalyteOfInterest,' vs. ', 
                                       outcomeAnalyteOfInterest, sep = ''))
  hc <- highcharter::hc_xAxis(hc, title=list(text=independentAnalyteOfInterest))
  hc <- highcharter::hc_yAxis(hc, title=list(text=outcomeAnalyteOfInterest))
  hc <- highcharter::hc_chart(hc, zoomType = "xy")
  hc <- highcharter::hc_add_series(hc, data=ds,type="scatter",#col=cols[1],
                               tooltip = list(headerFormat="",
                                              pointFormat=paste("{point.label}","{point.z}")),
                               showInLegend=FALSE)
  
  hc <- highcharter::hc_add_series(hc, name = uniqtypes[1],
                               data=line1,type='line',#name=sprintf("regression line %s",type1),
                               color = cols[1],enableMouseTracking=FALSE,marker=FALSE)
  hc <- highcharter::hc_add_series(hc, name = uniqtypes[2],
                               data=line2,type='line',#name=sprintf("regression line %s",type2),
                               color = cols[2],enableMouseTracking=FALSE,marker=FALSE)
  
  hc
}


#' 'volcano' plot (difference in correlations vs p-values)
#' of all pairs
#'
#' @param inputResults Data frame with model results (output of ProcessResults())
#' @param inputData Named list (output of 
#' FilterData()) with analyte levels 
#' and associated meta-data
#' @param nrpoints number of points to be plotted in lowest density areas (see 'smoothScatter' documentation for more detail)
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param diffcorr cutoff of differences in correlations for filtering (default 0.5)
#' @return a smoothScatter plot
#' @export
pvalCorrVolcano <- function(inputResults, inputData,nrpoints=10000,diffcorr=0.5,pvalcutoff=0.05){
    if(class(inputResults) != "IntLimResults") {
	stop("input data is not a IntLim class")
    }
  p <- inputData@sampleMetaData[,inputResults@stype]
  
  if (length(unique(p)) !=2){
    stop(paste("pvalCorrVolcano is invalid for continuous outcomes and outcomes
               with more than two categories."))
  }
    volc.results <- IntLIM::ProcessResults(inputResults,  inputData, diffcorr = 0, pvalcutoff = 1)
    volc.table <- volc.results
    Corrdiff <- volc.table[,4] - volc.table[,3]
    pval <- -log10(volc.table$FDRadjPval)
    graphics::smoothScatter(x = Corrdiff, pval, xlab = 'Difference in Correlation between Phenotypes',
		ylab = '-log10(FDR-adjusted p-value)', nrpoints=nrpoints,
                main = 'Volcano Plot')
    graphics::abline(h=-log10(pvalcutoff),lty=2,col="blue")
    graphics::abline(v=c(diffcorr,-diffcorr),lty=2,col="blue")
}

#' Makes an UpSet plot showing the filtered pairs of analytes found in each fold.
#' This plot should only be made for cross-validation data.
#' @param inputResults List of outputs of ProcessResultsAllFolds(), each of which
#' is a list of IntLIMResults.
#' @return an UpSet plot
PlotFoldOverlapUpSet<-function(inputResults){
  sig_list <- lapply(1:length(inputResults), function(i){
    return(paste(inputResults[[i]][,1], inputResults[[i]][,2], sep = "_"))
  })
  sig_mat <- ComplexHeatmap::list_to_matrix(sig_list)
  colnames(sig_mat) <- lapply(1:length(inputResults), function(i){
    return(paste("Fold", i, sep = "_"))
  })
  comb_mat <- ComplexHeatmap::make_comb_mat(sig_mat)
  ComplexHeatmap::UpSet(comb_mat)
}

#' Graphs a scatterplot of pairs vs. the interaction coefficient
#' for the pair
#' @param inputResults Data frame with model results (output of ProcessResults())
#' @param interactionCoeffPercentile percentile cutoff for interaction coefficient 
#' (default bottom 10 percent (high negative coefficients) and top 10 percent 
#' (high positive coefficients))
#' @param percentageToPlot percentage of points to plot (the points will be 
#' randomly selected) -- plotting all points will likely overwhelm plotting function.
#' @param independent.var.type type of analyte used as the independent variable 
#' ("1" or "2")
#' @param outcome type of analyte used as the outcome/dependent variable ("1"
#' or "2")
#' @return a scatterplot
#'
#' @export
InteractionCoefficientGraph<-function(inputResults,
                                      interactionCoeffPercentile=0.10,
                                      percentageToPlot = 0.01, 
                                      independent.var.type = 1,
                                      outcome = 2){


    if(class(inputResults) != "IntLimResults") {
      stop("input data is not a IntLim class")
    }

    #merge and properly name all data to return
    format_coeff = reshape2::melt(inputResults@interaction.coefficients)
    format_pval = reshape2::melt(inputResults@interaction.pvalues)
    format_adjp = reshape2::melt(inputResults@interaction.adj.pvalues)
    tofilter = cbind(format_coeff, format_pval$value, 
                     format_adjp$value)
    colnames(tofilter) = c("analyte1", "analyte2", "interaction_coeff", "Pval","FDRadjPval")


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
    if(independent.var.type == outcome){
      toplot_sort = toplot_sort[which(!is.na(toplot_sort$interaction_coeff)),]
    }
    if((independent.var.type != 1 && independent.var.type != 2) ||
       (outcome != 1 && outcome != 2)){
      stop("Error! outcome and independent.var.type must each be one of the following: 1, 2")
    }else{
      plot(1:length(toplot_sort$interaction_coeff),toplot_sort$interaction_coeff, 
           col=toplot_sort$color, xlab = "AnalytePairs", ylab = 
             "Interaction Coefficient", pch=16)
    }
}


#' Creates a dataframe of the marginal effect of phenotype
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData Named list (output of 
#' FilterData()) with analyte levels 
#' and associated meta-data
#' @param outcomeAnalyteOfInterest outcome analyte in pair
#' @param independentAnalyteOfInterest independent analyte in pair
#' @param continuous whether or not the outcome is continuous (TRUE or FALSE)
#' @param outcome '1' or '2' must be set as outcome/independent variable
#' @param independentVariable '1' or '2' must be set as outcome/independent variable
#' @return dataframe for further analysis
#' @export
MarginalEffectsGraphDataframe<-function(inputResults, inputData, independentAnalyteOfInterest, 
                                        outcomeAnalyteOfInterest, continuous, outcome,
                                        independentVariable){
  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }
  
  #get covariates
  covariates = as.character(inputResults@covar$covariate)
  covariates_class = as.character(inputResults@covar$class.var)
  
  #get dataframes
  pheno <- inputData@sampleMetaData[,inputResults@stype]
  outcomeAnalytes <- NULL
  independentAnalytes <- NULL
  if(outcome == 1){
    outcomeAnalytes <- inputData@analyteType1
  }else if(outcome == 2){
    outcomeAnalytes <- inputData@analyteType2
  }
  if(independentVariable == 1){
    independentAnalytes <- inputData@analyteType1
  }else if(independentVariable == 2){
    independentAnalytes <- inputData@analyteType2
  }
  
  #get one pair
  outcomeData = as.numeric(outcomeAnalytes[make.names(outcomeAnalyteOfInterest),])
  independentData = as.numeric(independentAnalytes[make.names(independentAnalyteOfInterest),])
  
  #Add analyte and phenotype data for glm
  forglm  = data.frame(row.names = 1:length(outcomeData))
  forglm$g = outcomeData
  forglm$type = pheno
  forglm$Y = as.numeric(independentData)
  
  
  if (!is.null(covariates)) {
    
    #Add all covariates to dataframe for glm()
    i=3
    for(each in covariates){
      names = colnames(forglm)
      i = i+1
      forglm[,i] = inputData@sampleMetaData[,each]
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
    margins::cplot(model, "type", data = dataframe, what = "prediction", main = title)
  }, error = function(cond){
    print("Could not plot the data. Check to see whether your outcome is continuous.
          Only categorical outcomes are valid for this function.")
  })
  return(model)

}


#' histogram of analyte pairs
#' depending upon independent or outcome analyte
#'
#' @param inputResults Data frame with model results (output of ProcessResults())
#' @param type 'independent' or 'outcome'.  'outcome' set as default
#' @param breaks Number of breaks selected for histogram
#' @export
HistogramPairs <- function(inputResults, type = 'outcome', breaks = 50){

  x <- inputResults
  
  if(is.null(x)){
      stop('Please run ProcessResults() before inputting into HistogramPairs')
  }
  if (type == 'outcome'){
    pairs <- data.frame(table(x$Analyte1))
    pairs.number <- as.vector(pairs$Freq)
    hist(pairs.number, breaks = breaks, main = "Number of analyte 
         pairs based on outcome analyte", xlab = 'Analyte pairs based on outcome analyte')
  }else if (type == 'independent'){
    pairs <- data.frame(table(x$Analyte2))
    pairs.number <- as.vector(pairs$Freq)
    hist(pairs.number, main = "Number of analyte pairs based on independent
         variable analyte", 
         breaks = breaks, xlab = 'Analyte pairs based on independent variable analyte')
  }else{
      stop("Only two valid types:  outcome or independent.  Invalid type entered")
  }
}