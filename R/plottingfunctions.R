#' Get some stats after reading in data
#'
#' @import magrittr
#' @import highcharter
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
  if("gene" %in% names(inputData)){
    mygene <- inputData$gene
    toplot <- suppressMessages(reshape2::melt(mygene))
    df <- dplyr::tibble(value = toplot$value, by = toplot$Var2) %>% dplyr::group_by_at("by") %>%
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
  if("metab" %in% names(inputData)){
    mymetab <- inputData$metab
    toplot <- suppressMessages(reshape2::melt(t(mymetab)))
    df <- dplyr::data_frame(value = toplot$value, by = toplot$Var1) %>%
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
#' @param inputData IntLimObject output of ReadData()
#' @param stype category to color-code by (can be more than two categories)
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @return a highcharter object
#' @export
PlotPCA <- function(inputData,viewer=T,stype=NULL,palette = "Set1") {

  if(is.numeric(inputData$p) == TRUE) {
		warning("The resulting PCA plot is not color-coded because you did not provide 
		        a categorical variable in 'stype'")
		mytype <- NULL
  } else {
  	mytype <- as.character(inputData$p)
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

  if(!("gene" %in% names(inputData) && "metab" %in% names(inputData))){
    stop("A dataset not containing both expression and metabolite data cannot run
	         with 'common' set to TRUE. Set 'common' to FALSE.")
  } else {
    if(is.numeric(inputData$p) == TRUE) {
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
      mymetab <- inputData$metab
      mygene <- inputData$gene
      alltype <- inputData$p
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

  # Set up plots.
  if("metab" %in% names(inputData)){
    mpercvar=round((mpca$sdev)^2 / sum(mpca$sdev^2)*100,2)
    pm <- pm %>% highcharter::hc_title(text="PCA of metabolites") %>%
      highcharter::hc_xAxis(title=list(text=paste0("PC1:",round(mpercvar[1],1),"%"))) %>%
      highcharter::hc_yAxis(title=list(text=paste0("PC2:",round(mpercvar[2],2),"%"))) %>%
      hc_chart(zoomType = "xy")
  }
  if("gene" %in% names(inputData)){
    gpercvar=round((gpca$sdev)^2 / sum(gpca$sdev^2)*100,2)
    pg <- pg %>% highcharter::hc_title(text="PCA of genes") %>%
      highcharter::hc_xAxis(title=list(text=paste0("PC1:",round(gpercvar[1],1),"%"))) %>%
      highcharter::hc_yAxis(title=list(text=paste0("PC2:",round(gpercvar[2],2),"%"))) %>%
      hc_chart(zoomType = "xy")
  }
  p <- NULL
  if("gene" %in% names(inputData) && "metab" %in% names(inputData)){
    if (viewer == TRUE) {
      p <-htmltools::browsable(highcharter::hw_grid(pg, pm, ncol = 2, rowheight = 550))
    } else {
      p <- highcharter::hw_grid(pg, pm)
    }
  } else if("gene" %in% names(inputData)){
    if (viewer == TRUE) {
      p <-htmltools::browsable(highcharter::hw_grid(pg, ncol = 1, rowheight = 550))
    } else {
      p <- highcharter::hw_grid(pg)
    }
  } else if("metab" %in% names(inputData)){
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
#' @importFrom graphics boxplot par
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
#' @export
DistRSquared<- function(IntLimResults,breaks=100) {
  
  hist(IntLimResults@model.rsquared,breaks=breaks,
       main="Histogram of Interaction R-Squared Values")
}

#' Returns the clusters found using CorrHeatmap.
#'
#' @param inputResults Data frame (output of ProcessResults())
#' @param inputData Named list (output of 
#' FilterData()) with gene expression, metabolite abundances, 
#' and associated meta-data
#' @param top_pairs cutoff of the top pairs, sorted by adjusted p-values, to be plotted (plotting more than 1200 can take some time) (default: 1200)
#' @param treecuts number of clusters (of gene-metabolite pairs) to cut the tree into for color-coding
#' @return a highcharter object
#' @export
GetCorrClusters <- function(inputResults,inputData,top_pairs=1200,treecuts=2) {
  type <- cor <- c()
  clusters <- NULL
  
  if(nrow(inputResults)==0) {
    stop("Make sure you run ProcessResults before making the heatmap")
  }
  p <- inputData$p
  if(length(unique(p)) !=2){
    stop("GetCorrClusters requires 2 discrete phenotypes. Do not run with continuous phenotypes.")
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
    
    # hclust is used under the hood of heatmaply.
    dend <- stats::hclust(stats::dist(heat_data, method = "euclidean"))
    clusters <- stats::cutree(dend, k = treecuts)
  }
  return(clusters)
}
#' Plot correlation heatmap
#'
#' @import magrittr
#'
#' @param inputResults Data frame (output of ProcessResults())
#' @param inputData Named list (output of 
#' FilterData()) with gene expression, metabolite abundances, 
#' and associated meta-data
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
#' @export
CorrHeatmap <- function(inputResults,inputData,viewer=T,top_pairs=1200,treecuts=2, 
                        palette = NULL, static = FALSE,
                        html.file=NULL, pdf.file=NULL) {
  type <- cor <- c()

	if(nrow(inputResults)==0) {
		stop("Make sure you run ProcessResults before making the heatmap")
	}
  p <- inputData$p
  if(length(unique(p)) !=2){
    stop("CorrHeatmap requires 2 discrete phenotypes. Do not run with continuous phenotypes.")
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
#' @export
PlotGMPair<- function(inputData,stype,geneName,metabName,palette = "Set1",
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

  gene<-inputData$gene
  if(length(which(rownames(gene)==geneName))>0) {
    sGene<-gene[geneName,]
  } else {
    stop(paste0("The gene ",geneName," was not found in your data"))
  }

  metab<-inputData$metab
  if(length(which(rownames(metab)==metabName))>0) {
  	sMetab<-as.numeric(metab[metabName,])
  } else {
    stop(paste0("The metabolite ",metabName," was not found in your data"))
  }

  if(length(unique(inputData$p))!=2) {
    stop(paste0("The group selected, '",stype,"', should only contain two different categories"))
  }

  mycols <- as.character(inputData$p)
  mycols[which(inputData$p==unique(inputData$p)[1])] <- cols[1]
  mycols[which(inputData$p==unique(inputData$p)[2])] <- cols[2]

  data<-data.frame(x=sGene,y=sMetab,z=colnames(gene),label=inputData$p,color=mycols)
  
  # Get points to draw the lines for each phenotype by hand

  uniqtypes=as.character(unique(inputData$p))

  # Starting with phenotype 1, get min and max x values constrained to the values of y
  # The reason we do this, is because the lines do not necessary need to go out to the max or min of x, particularly
  # when slopes are really steep (abline does this automatically but not highcharter)
  mytypes <- inputData$p
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
#' @export
PlotMGPair<- function(inputData,stype,metabName,geneName,palette = "Set1",
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
  
  gene<-inputData$gene
  if(length(which(rownames(gene)==geneName))>0) {
    sGene<-gene[geneName,]
  } else {
    stop(paste0("The gene ",geneName," was not found in your data"))
  }
  
  metab<-inputData$metab
  if(length(which(rownames(metab)==metabName))>0) {
    sMetab<-as.numeric(metab[metabName,])
  } else {
    stop(paste0("The metabolite ",metabName," was not found in your data"))
  }
  
  if(length(unique(inputData$p))!=2) {
    stop(paste0("The group selected, '",stype,"', should only contain two different categories"))
  }
  
  mycols <- as.character(inputData$p)
  mycols[which(inputData$p==unique(inputData$p)[1])] <- cols[1]
  mycols[which(inputData$p==unique(inputData$p)[2])] <- cols[2]
  
  data<-data.frame(x=sMetab,y=sGene,z=colnames(gene),label=inputData$p,color=mycols)
  
  # Get points to draw the lines for each phenotype by hand
  
  uniqtypes=as.character(unique(inputData$p))
  
  # Starting with phenotype 1, get min and max x values constrained to the values of y
  # The reason we do this, is because the lines do not necessary need to go out to the max or min of x, particularly
  # when slopes are really steep (abline does this automatically but not highcharter)
  mytypes <- inputData$p
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
#' @export
PlotMMPair<- function(inputData,stype,metab1Name,metab2Name,palette = "Set1",
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

  metab<-inputData$metab
  if(length(which(rownames(metab)==metab1Name))>0) {
    sMetab1<-as.numeric(metab[metab1Name,])
  } else {
    stop(paste0("The metabolite ",metab1Name," was not found in your data"))
  }
  
  metab<-inputData$metab
  if(length(which(rownames(metab)==metab2Name))>0) {
    sMetab2<-as.numeric(metab[metab2Name,])
  } else {
    stop(paste0("The metabolite ",metab2Name," was not found in your data"))
  }
  
  if(length(unique(inputData$p))!=2) {
    stop(paste0("The group selected, '",stype,"', should only contain two different categories"))
  }
  
  mycols <- as.character(inputData$p)
  mycols[which(inputData$p==unique(inputData$p)[1])] <- cols[1]
  mycols[which(inputData$p==unique(inputData$p)[2])] <- cols[2]
  
  data<-data.frame(x=sMetab1,y=sMetab2,z=colnames(metab),label=inputData$p,color=mycols)

  # Get points to draw the lines for each phenotype by hand
  
  uniqtypes=as.character(unique(inputData$p))
  
  # Starting with phenotype 1, get min and max x values constrained to the values of y
  # The reason we do this, is because the lines do not necessary need to go out to the max or min of x, particularly
  # when slopes are really steep (abline does this automatically but not highcharter)
  mytypes <- inputData$p
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
  
  gene<-inputData$gene
  if(length(which(rownames(gene)==gene1Name))>0) {
    sGene1<-as.numeric(gene[gene1Name,])
  } else {
    stop(paste0("The gene ",gene1Name," was not found in your data"))
  }
  
  gene<-inputData$gene
  if(length(which(rownames(gene)==gene2Name))>0) {
    sGene2<-as.numeric(gene[gene2Name,])
  } else {
    stop(paste0("The gene ",gene2Name," was not found in your data"))
  }
  
  if(length(unique(inputData$p))!=2) {
    stop(paste0("The group selected, '",stype,"', should only contain two different categories"))
  }
  
  mycols <- as.character(inputData$p)
  mycols[which(inputData$p==unique(inputData$p)[1])] <- cols[1]
  mycols[which(inputData$p==unique(inputData$p)[2])] <- cols[2]
  
  data<-data.frame(x=sGene1,y=sGene2,z=colnames(gene),label=inputData$p,color=mycols)
  
  # Get points to draw the lines for each phenotype by hand
  
  uniqtypes=as.character(unique(inputData$p))
  
  # Starting with phenotype 1, get min and max x values constrained to the values of y
  # The reason we do this, is because the lines do not necessary need to go out to the max or min of x, particularly
  # when slopes are really steep (abline does this automatically but not highcharter)
  mytypes <- inputData$p
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
#' @param inputResults Data frame with model results (output of ProcessResults())
#' @param inputData Named list (output of 
#' FilterData()) with gene expression, metabolite abundances, 
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
  p <- inputData$p
  
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
  names(sig_list) <- names(inputResults)
  comb_mat <- ComplexHeatmap::make_comb_mat(ComplexHeatmap::list_to_matrix(sig_list))
  ComplexHeatmap::UpSet(comb_mat)
}

#' Graphs a scatterplot of gene-metabolite pairs vs. the interaction coefficient
#' for the gene-metabolite pair
#' @param inputResults Data frame with model results (output of ProcessResults())
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
#' @param inputData Named list (output of 
#' FilterData()) with gene expression, metabolite abundances, 
#' and associated meta-data
#' @param metaboliteOfInterest metabolite in gene-metabolite pair
#' @param geneOfInterest gene in gene-metabolite pair
#' @param continuous whether or not the outcome is continuous (TRUE or FALSE)
#' @return dataframe for further analysis
#' @export
MarginalEffectsGraphDataframe<-function(inputResults, inputData, geneOfInterest, 
                                        metaboliteOfInterest, continuous){
  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }

  #get covariates
  covariates = as.character(inputResults@covar$covariate)
  covariates_class = as.character(inputResults@covar$class.var)

  #get dataframes
  pheno <- inputData$p
  gene <- inputData$gene
  metab <- inputData$metab

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
      forglm[,i] = inputData$covar_matrix[,each]
      colnames(forglm) = c(names, each)
    }
  }
  return(forglm)
}

#' Creates a dataframe of the marginal effect of phenotype
#'
#' @import margins
#'
#' @param inputResults IntLIMResults with model results (output of RunIntLim())
#' @param inputData Data frame object (output of FilterData()) with gene expression,
#' @param metaboliteOfInterest1 outcome metabolite in metabolite pair
#' @param metaboliteOfInterest2 independent metabolite in metabolite pair
#' @return dataframe for further analysis
#' @export
MarginalEffectsGraphDataframeMetabolitePairs<-function(inputResults, inputData, 
                                                       metaboliteOfInterest1, 
                                                       metaboliteOfInterest2){
  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }
  
  #get covariates
  covariates = as.character(inputResults@covar$covariate)
  covariates_class = as.character(inputResults@covar$class.var)
  
  #get dataframes
  pheno <- inputData$p
  metab <- inputData$metab
  
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
      forglm[,i] = inputData$covar_matrix[,each]
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
#' @param inputData Named list (output of 
#' FilterData()) with gene expression, metabolite abundances, 
#' and associated meta-data
#' @param geneOfInterest1 outcome gene in gene pair
#' @param geneOfInterest2 independent gene in gene pair
#' @return dataframe for further analysis
#' @export
MarginalEffectsGraphDataframeGenePairs<-function(inputResults, inputData, 
                                                 geneOfInterest1, 
                                                 geneOfInterest2){
  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }
  
  #get covariates
  covariates = as.character(inputResults@covar$covariate)
  covariates_class = as.character(inputResults@covar$class.var)
  
  #get dataframes
  pheno <- inputData$p
  gene <- inputData$gene
  
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
      forglm[,i] = inputData$covar_matrix[,each]
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