#' Get some stats after reading in data
#'
#' @param inputData IntLimObject output of ReadData()
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param viewer whether the plot should be displayed in the RStudio viewer (TRUE) or
#' in Shiny/Knittr (FALSE)
#' @return a highcharter object
#' @export
PlotDistributions <- function(inputData,viewer=TRUE, palette="Set1"){
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
    # Compute boxplot statistics for analyte type 1.
    type1Data <- inputData@analyteType1
    toplot <- suppressMessages(reshape2::melt(type1Data))
    df <- data.frame(value = toplot$value, by = toplot$Var2)
    stats <- lapply(sort(unique(df$by)), function(grp){
      return(grDevices::boxplot.stats(df$value[which(df$by == grp)]))
    })
    bxps <- lapply(stats, function(stat){
      return(stat$stats)
    })
    
    # Construct output.
    outsList <- lapply(seq(length(unique(df$by))), function(x) {
      y <- stats[[x]]
      d <- data.frame()
      if (length(y$out) > 0)
        d <- data.frame(x = x - 1, y = y$out)
      return(d)
    })
    outs <- do.call(rbind, outsList)
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
    g <- highcharter::hc_xAxis(g, categories = colnames(type1Data))
    g <- highcharter::hc_tooltip(g, valueDecimals = 2)
    g <- highcharter::hc_exporting(g, enabled = TRUE)
  }
  if(length(inputData@analyteType2)>0){
    # Compute boxplot statistics for analyte type 2.
    type2Data <- inputData@analyteType2
    toplot <- suppressMessages(reshape2::melt(type2Data))
    df <- data.frame(value = toplot$value, by = toplot$Var2)
    stats <- lapply(sort(unique(df$by)), function(grp){
      return(grDevices::boxplot.stats(df$value[which(df$by == grp)]))
    })
    bxps <- lapply(stats, function(stat){
      return(stat$stats)
    })
    
    # Construct output.
    outsList <- lapply(seq(length(unique(df$by))), function(x) {
      y <- stats[[x]]
      d <- data.frame()
      if (length(y$out) > 0)
        d <- data.frame(x = x - 1, y = y$out)
      return(d)
    })
    outs <- do.call(rbind, outsList)
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
    m <- highcharter::hc_xAxis(m, categories = colnames(type2Data))
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
#' @param viewer whether the plot should be displayed in the RStudio viewer (TRUE) or
#' in Shiny/Knittr (FALSE)
#' @return a highcharter object
#' @export
PlotPCA <- function(inputData,viewer=TRUE,stype="",palette = "Set1") {
  
  if(is.numeric(inputData@sampleMetaData[,stype]) == TRUE) {
    mytype <- inputData@sampleMetaData[,stype]
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
    if(is.numeric(mytype) == TRUE) {
      mpca <- stats::prcomp(t(inputData@analyteType1),center=TRUE,scale=FALSE)
      gpca <- stats::prcomp(t(inputData@analyteType2),center=TRUE,scale=FALSE)
      
      # Set colors.
      bin_count <- 100
      # Make sure the spacing is even. We need to do this using seq.
      intervals <- seq(min(mytype), max(mytype),
                       by = (max(mytype) - min(mytype)) / (bin_count - 1))
      subject_color_scale <- findInterval(mytype, intervals)
      pal <- grDevices::colorRampPalette(c("#89CFF0", "#002366"))(bin_count+1)
      mycols <-pal[subject_color_scale]
      
      gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),color=mycols)
      mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),color=mycols)
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
      pm <- highcharter::hc_colorAxis(pm, min=min(mytype), max=max(mytype), 
                                      minColor = "#89CFF0", maxColor = "#002366")
    } else {
      type1 <- inputData@analyteType1
      type2 <- inputData@analyteType2
      alltype <- inputData@sampleMetaData[,stype]
      uniqtypes <- unique(alltype)
      mycols <- as.character(alltype)
      for (i in 1:uniqtypes) {
        mycols[which(alltype==uniqtypes[i])] <- cols[i]
      }
      gpca <- stats::prcomp(t(type1),center=TRUE,scale=FALSE)
      mpca <- stats::prcomp(t(type2),center=TRUE,scale=FALSE)
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
  }else if(length(inputData@analyteType1)>0){
    mpca <- stats::prcomp(t(inputData@analyteType1),center=TRUE,scale=FALSE)
    
    if(is.numeric(mytype) == TRUE) {
      # Set colors.
      bin_count <- 100
      # Make sure the spacing is even. We need to do this using seq.
      intervals <- seq(min(mytype), max(mytype),
                       by = (max(mytype) - min(mytype)) / (bin_count - 1))
      subject_color_scale <- findInterval(mytype, intervals)
      pal <- grDevices::colorRampPalette(c("#89CFF0", "#002366"))(bin_count+1)
      mycols <-pal[subject_color_scale]
      
      mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),color=mycols)
      mds <- highcharter::list_parse(mtoplot)
      pm <- highcharter::highchart(width = 350, height = 350)
      pm <- highcharter::hc_add_series(pm, data=mds,type="scatter",
                                       tooltip = list(headerFormat="",
                                                      pointFormat=paste("{point.label}","{point.z}")),
                                       showInLegend=FALSE)
      pm <- highcharter::hc_colorAxis(pm, min=min(mytype), max=max(mytype), 
                                      minColor = "#89CFF0", maxColor = "#002366")
    } else {
      type1 <- inputData@analyteType1
      alltype <- inputData@sampleMetaData[,stype]
      uniqtypes <- unique(alltype)
      mycols <- as.character(alltype)
      for (i in 1:numcateg) {
        mycols[which(alltype==uniqtypes[i])] <- cols[i]
      }
      mpca <- stats::prcomp(t(type1),center=TRUE,scale=FALSE)
      mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),label=alltype,color=mycols)
      mds <- highcharter::list_parse(mtoplot)
      pm <- highcharter::highchart(width = 350, height = 350)
      for (i in 1:length(uniqtypes)) {
        mytype <- unique(alltype)[i]
        mds <- highcharter::list_parse(mtoplot[which(mtoplot$label==mytype),])
        pm <- highcharter::hc_add_series(pm, data=mds,type="scatter",
                                         name=mytype,
                                         color=cols[which(alltype==mytype)[1]],
                                         tooltip = list(headerFormat="",
                                                        pointFormat=paste("{point.label}","{point.z}")),
                                         showInLegend=TRUE)
      }
    }
  }else if(length(inputData@analyteType2)>0){
    if(is.numeric(mytype) == TRUE) {
      # Set colors.
      bin_count <- 100
      # Make sure the spacing is even. We need to do this using seq.
      intervals <- seq(min(mytype), max(mytype),
                       by = (max(mytype) - min(mytype)) / (bin_count - 1))
      subject_color_scale <- findInterval(mytype, intervals)
      pal <- grDevices::colorRampPalette(c("#89CFF0", "#002366"))(bin_count+1)
      mycols <-pal[subject_color_scale]
      
      gpca <- stats::prcomp(t(inputData@analyteType2),center=TRUE,scale=FALSE)
      gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),color=mycols)
      gds <- highcharter::list_parse(gtoplot)
      pg <- highcharter::highchart(width = 350, height = 350 )
      pg <- highcharter::hc_add_series(pg, data=gds,type="scatter",
                                       tooltip = list(headerFormat="",
                                                      pointFormat=paste("{point.label}","{point.z}")),
                                       showInLegend=FALSE)
    } else {
      type2 <- inputData@analyteType2
      alltype <- inputData@sampleMetaData[,stype]
      uniqtypes <- unique(alltype)
      mycols <- as.character(alltype)
      for (i in 1:numcateg) {
        mycols[which(alltype==uniqtypes[i])] <- cols[i]
      }
      gpca <- stats::prcomp(t(type2),center=TRUE,scale=FALSE)
      gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),label=alltype,color=mycols)
      gds <- highcharter::list_parse(gtoplot)
      pg <- highcharter::highchart(width = 350, height = 350)
      for (i in 1:length(uniqtypes)) {
        mytype <- unique(alltype)[i]
        gds <- highcharter::list_parse(gtoplot[which(gtoplot$label==mytype),])
        pg <- highcharter::hc_add_series(pg, data=gds,type="scatter",
                                         name=mytype,
                                         color=cols[which(alltype==mytype)[1]],
                                         tooltip = list(headerFormat="",
                                                        showInLegend=TRUE))
      }
    }
  }
  if(length(inputData@analyteType1)>0){
    mpercvar=round((mpca$sdev)^2 / sum(mpca$sdev^2)*100,2)
    pm <- highcharter::hc_title(pm, text="PCA of analyte type 1")
    pm <- highcharter::hc_xAxis(pm, title=list(text=paste0("PC1:",round(mpercvar[1],1),"%")))
    pm <- highcharter::hc_yAxis(pm, title=list(text=paste0("PC2:",round(mpercvar[2],2),"%")))
    pm <- highcharter::hc_chart(pm, zoomType = "xy")
  }
  if(length(inputData@analyteType2)>0){
    gpca <- stats::prcomp(t(inputData@analyteType2),center=TRUE,scale=FALSE)
    gpercvar=round((gpca$sdev)^2 / sum(gpca$sdev^2)*100,2)
    pg <- highcharter::hc_title(pg, text="PCA of analyte type 2")
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
      p <-htmltools::browsable(highcharter::hw_grid(pm, ncol = 1, rowheight = 550))
    } else {
      p <- highcharter::hw_grid(pm)
    }
  } else if(length(inputData@analyteType2)>0){
    if (viewer == TRUE) {
      p <-htmltools::browsable(highcharter::hw_grid(pg, ncol = 1, rowheight = 550))
    } else {
      p <- highcharter::hw_grid(pg)
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
#' @return No return value, called for side effects
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
#' @return No return value, called for side effects
#' @export
PValueBoxPlots<- function(IntLimResults) {
  if(length(IntLimResults@covariate.pvalues) == 0){
    stop("Error! You must set save.covar.pvals to TRUE when running IntLIM to run PValueBoxPlots")
  }else{
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
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
#' @return No return value, called for side effects
#' @export
DistRSquared<- function(IntLimResults,breaks=100) {
  
  hist(IntLimResults@model.rsquared,breaks=breaks,
       main="Histogram of Interaction R-Squared Values")
}

#' A helper function for the PlotPair functions (i.e. the highcharter one and
#' the flat, base-R one).
#' @param inputData IntLimObject output of ReadData() or FilterData()
#' @param inputResults Data frame with model results (output of ProcessResults())
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param outcomeAnalyteOfInterest outcome analyte in pair
#' @param independentAnalyteOfInterest independent analyte in pair
#' @param outcome '1' or '2' must be set as outcome/independent variable
#' @param independentVariable '1' or '2' must be set as outcome/independent variable
#' @param stype Phenotype or outcome variable
BuildDataAndLines <- function(inputData,inputResults,outcome,independentVariable, 
                              independentAnalyteOfInterest, 
                              outcomeAnalyteOfInterest, palette = "Set1",	stype){
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
  
  # Return everything needed.
  return(list(data=data,uniqtypes=uniqtypes,line1=line1,line2=line2,cols=cols))
}
  
#' scatter plot of pairs (based on user selection)
#'
#' @param inputData IntLimObject output of ReadData() or FilterData()
#' @param inputResults Data frame with model results (output of ProcessResults())
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param viewer whether the plot should be displayed in the RStudio viewer (TRUE) or
#' in Shiny/Knittr (FALSE)
#' @param outcomeAnalyteOfInterest outcome analyte in pair
#' @param independentAnalyteOfInterest independent analyte in pair
#' @param outcome '1' or '2' must be set as outcome/independent variable
#' @param independentVariable '1' or '2' must be set as outcome/independent variable
#' @return No return value, called for side effects
#' @export
PlotPair<- function(inputData,inputResults,outcome,independentVariable, independentAnalyteOfInterest, 
                    outcomeAnalyteOfInterest, palette = "Set1",	viewer=TRUE) {
  
  # Set type.
  stype <- inputResults@stype
  
  # Check whether continuous or discrete.
  unique_stypes <- unique(inputData@sampleMetaData[,stype])
  
  # For continuous data, plot the marginal effects graph.
  if(length(unique_stypes) > 2){
    MarginalEffectsGraph(
      dataframe = MarginalEffectsGraphDataframe(inputResults = inputResults,
                                                inputData = inputData,
                                                outcomeAnalyteOfInterest = outcomeAnalyteOfInterest,
                                                independentAnalyteOfInterest = independentAnalyteOfInterest,
                                                outcome = outcome,
                                                independentVariable = independentVariable), 
      title = paste("Marginal Effects -", independentAnalyteOfInterest,
                    "and", outcomeAnalyteOfInterest), xlab = independentAnalyteOfInterest,
      ylab = outcomeAnalyteOfInterest)
  }else{
    
    # Get data.
    plotdata <- BuildDataAndLines(inputData,inputResults,outcome,independentVariable, 
                                  independentAnalyteOfInterest, 
                                  outcomeAnalyteOfInterest, palette,stype)
    data <- plotdata$data
    uniqtypes <- plotdata$uniqtypes
    line1 <- plotdata$line1
    line2 <- plotdata$line2
    cols <- plotdata$cols
    
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
}

#' scatter plot of pairs (based on user selection). This version does not use
#' highcharter and instead plots a base R plot.
#'
#' @param inputData IntLimObject output of ReadData() or FilterData()
#' @param inputResults Data frame with model results (output of ProcessResults())
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param outcomeAnalyteOfInterest outcome analyte in pair
#' @param independentAnalyteOfInterest independent analyte in pair
#' @param outcome '1' or '2' must be set as outcome/independent variable
#' @param independentVariable '1' or '2' must be set as outcome/independent variable
#' @return No return value, called for side effects
#' @export
PlotPairFlat<- function(inputData,inputResults,outcome,independentVariable, independentAnalyteOfInterest, 
                    outcomeAnalyteOfInterest, palette = "Set1") {
    
  # Set type.
  stype <- inputResults@stype
  
  # Check whether continuous or discrete.
  unique_stypes <- unique(inputData@sampleMetaData[,stype])
  
  # For continuous data, plot the marginal effects graph.
  if(length(unique_stypes) > 2){
    MarginalEffectsGraph(
      dataframe = MarginalEffectsGraphDataframe(inputResults = inputResults,
                                                inputData = inputData,
                                                outcomeAnalyteOfInterest = outcomeAnalyteOfInterest,
                                                independentAnalyteOfInterest = independentAnalyteOfInterest,
                                                outcome = outcome,
                                                independentVariable = independentVariable), 
      title = paste("Marginal Effects -", independentAnalyteOfInterest,
                    "and", outcomeAnalyteOfInterest), xlab = independentAnalyteOfInterest,
      ylab = outcomeAnalyteOfInterest)
  }else{  
    
    # Get data.
    plotdata <- BuildDataAndLines(inputData,inputResults,outcome,independentVariable, 
                                  independentAnalyteOfInterest, 
                                  outcomeAnalyteOfInterest, palette,stype)
    data <- plotdata$data
    uniqtypes <- plotdata$uniqtypes
    line1 <- plotdata$line1
    line2 <- plotdata$line2
    cols <- plotdata$cols
    
    # Plot
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    graphics::par(mar = c(5, 4, 4, 8), xpd = TRUE, pty="s")
    plot(data$x, data$y, col = data$color, xlab = independentAnalyteOfInterest,
         ylab = outcomeAnalyteOfInterest, main = paste(independentAnalyteOfInterest,
                                                       "vs.", outcomeAnalyteOfInterest),
         pch = 16)
    graphics::lines(line1, col = cols[1])
    graphics::lines(line2, col = cols[2])
    coord <- graphics::par("usr")
    graphics::legend(x = coord[2] * 1.05, y = coord[4], legend=c(uniqtypes[1],uniqtypes[2]), 
           col=c(cols[1],cols[2]), title="stype",lty=1,bg="transparent")
  }
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
#' @param coefPercentileCutoff cutoff of interaction coefficient percentile.
#' @return a smoothScatter plot
#' @export
pvalCoefVolcano <- function(inputResults, inputData,nrpoints=10000,pvalcutoff=0.05,
                            coefPercentileCutoff=0.9){
    if(!methods::is(inputResults, "IntLimResults")) {
	    stop("input data is not a IntLim class")
    }
  
    # Get the formatted results of processing, including all results (p-val <= 1)
    volc.table <- IntLIM::ProcessResults(inputResults,  inputData, pvalcutoff = 1)
    interaction_coeff <- volc.table$interaction_coeff
    pval <- -log10(volc.table$Pval)
    
    # Get the p-value cutoff using the FDR-adjusted cutoff.
    if(length(which(volc.table$FDRadjPval <= pvalcutoff)) == 0){
      stop(paste("No p-values meet the provided FDR-adjusted cutoff of", 
                 pvalcutoff, "- please choose a higher p-value cutoff."))
    }
    pvals_below_cutoff <- volc.table[which(volc.table$FDRadjPval <= pvalcutoff),]
    highest_pval_below_cutoff <- pvals_below_cutoff[which.max(pvals_below_cutoff$FDRadjPval), "Pval"]
    
    # Create the scatter plot.
    graphics::smoothScatter(x = interaction_coeff, pval, xlab = "Interaction Coefficient",
		ylab = '-log10(p-value)', nrpoints=nrpoints,
                main = 'Volcano Plot')
    
    # Plot cutoff lines.
    graphics::abline(h=-log10(highest_pval_below_cutoff),lty=2,col="blue")
    lower_line = getQuantileForInteractionCoefficient(interaction_coeff, 
                                                      coefPercentileCutoff)[1]
    upper_line = getQuantileForInteractionCoefficient(interaction_coeff, 
                                                       coefPercentileCutoff)[2]
    graphics::abline(v=c(lower_line,upper_line),lty=2,col="blue")
}

#' Makes an UpSet plot showing the filtered pairs of analytes found in each fold.
#' This plot should only be made for cross-validation data.
#' @param inputResults List of outputs of ProcessResultsAllFolds(), each of which
#' is a list of IntLIMResults.
#' @return an UpSet plot
#' @export
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


    if(!methods::is(inputResults, "IntLimResults")) {
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
    randomize = function(x) sample(1:nrow(toplot_sort),x,replace=FALSE)
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
MarginalEffectsGraphDataframe<-function(inputResults, inputData, independentAnalyteOfInterest, 
                                        outcomeAnalyteOfInterest, continuous, outcome,
                                        independentVariable){
  if(!methods::is(inputResults, "IntLimResults")) {
    stop("input data is not a IntLim class")
  }
  
  #get covariates
  covariates = as.character(inputResults@covar)

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
  forglm$Y = outcomeData
  forglm$type = pheno
  forglm$g = as.numeric(independentData)
  
  
  if (covariates != "") {
    
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
#' @param ylab outcome analyte in pair
#' @param xlab independent analyte in pair
#' @return values used for graphing
MarginalEffectsGraph<-function(dataframe, title, ylab, xlab){
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
    margins::cplot(model, x = "type", data = dataframe, what = "prediction", 
                   main = title)
  }, error = function(cond){
    stop(cond)
  })
  return(model)

}


#' histogram of analyte pairs
#' depending upon independent or outcome analyte
#'
#' @param inputResults Data frame with model results (output of ProcessResults())
#' @param type 'independent' or 'outcome'.  'outcome' set as default
#' @param breaks Number of breaks selected for histogram
#' @return No return value, called for side effects
#' @export
HistogramPairs <- function(inputResults, type = 'outcome', breaks = 50){

  x <- inputResults
  
  if(is.null(x)){
      stop('Please run ProcessResults() before inputting into HistogramPairs')
  }
  if (type == 'outcome'){
    pairs <- data.frame(table(as.character(x$Analyte1)))
    pairs.number <- as.vector(pairs$Freq)
    hist(pairs.number, breaks = breaks, main = "Number of analyte 
         pairs based on outcome analyte", xlab = 'Analyte pairs based on outcome analyte')
  }else if (type == 'independent'){
    pairs <- data.frame(table(as.character(x$Analyte2)))
    pairs.number <- as.vector(pairs$Freq)
    hist(pairs.number, main = "Number of analyte pairs based on independent
         variable analyte", 
         breaks = breaks, xlab = 'Analyte pairs based on independent variable analyte')
  }else{
      stop("Only two valid types:  outcome or independent.  Invalid type entered")
  }
}

#' Return the number of significant analytes and the number 
#' of permutations in which each analyte is significant.
#' If plot = TRUE, show a box plot of number of significant analytes over permutations, 
#' overlaid with the number of significant analytes in the original data.
#'
#' @param inputResults Data frame with model results (output of ProcessResults())
#' @param permResults An object of type PermutationResults (output of PermuteIntLIM())
#' @param plot Whether or not to show the boxplot. Default is TRUE.
#' @return A data frame that includes, for each permutation, the number of significant
#' pairs and the number of unique analytes of each analyte type within those pairs
#' @export
PermutationCountSummary <- function(inputResults, permResults, plot){
  
  # Prevent "visible binding for global variable" notes.
  Count <- Type <- NULL
  
  # Compute summary.
  pairCountDistrib <- permResults$numSigPairs$Num_Significant_Pairs
  analyte1 <- lapply(1:length(permResults$listOfSigPairs), function(i){
    return(unlist(lapply(permResults$listOfSigPairs[[i]], function(string){
      return(strsplit(string, split="__V__")[[1]][1])
    })))
  })
  analyte2 <- lapply(1:length(permResults$listOfSigPairs), function(i){
    return(unlist(lapply(permResults$listOfSigPairs[[i]], function(string){
      return(strsplit(string, split="__V__")[[1]][2])
    })))
  })
  analyte1CountDistrib <- unlist(lapply(1:length(analyte1), function(i){
    return(length(unique(analyte1[[i]])))
  }))
  analyte2CountDistrib <- unlist(lapply(1:length(analyte2), function(i){
    return(length(unique(analyte2[[i]])))
  }))
  countDistribs <- data.frame(Pairs = pairCountDistrib, 
                              Independent = analyte1CountDistrib,
                              Outcome = analyte2CountDistrib)

  # Set up data for input.
  significantCounts <- data.frame(Count=c(countDistribs$Pairs, 
                                          countDistribs$Independent,
                                          countDistribs$Outcome),
                                  Type=c(rep("Pair", nrow(permResults$numSigPairs)),
                                         rep("Independent.Variable", nrow(permResults$numSigPairs)),
                                         rep("Outcome", nrow(permResults$numSigPairs))))
  
  # Compute the same for the original data.
  PairCount <- nrow(inputResults)
  inputIndependentCount <- length(unique(inputResults$Analyte1))
  inputOutcomeCount <- length(unique(inputResults$Analyte2))

  # Make plot.
  cols <- c("Original.Data"="red")
  fills <- c("Permuted.Data"="black")
  if(plot == TRUE){
    plt <- ggplot2::ggplot(significantCounts, ggplot2::aes(x=Type, y=Count)) + 
      ggplot2::geom_violin(ggplot2::aes(fill = "Permuted.Data")) + 
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title.x=ggplot2::element_blank(), 
                     panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black")) + 
      ggplot2::annotation_logticks(sides = "l", scaled = TRUE) + 
      ggplot2::scale_y_log10() + 
      ggplot2::geom_boxplot(width=0.05, fill = "white") + 
      ggplot2::geom_segment(ggplot2::aes(x = 0.6, xend = 1.4, y = inputIndependentCount, 
                                         yend = inputIndependentCount, color = "Original.Data")) + 
      ggplot2::geom_segment(ggplot2::aes(x = 1.6, xend = 2.4, y = inputOutcomeCount, 
                                         yend = inputOutcomeCount, color = "Original.Data")) +
      ggplot2::geom_segment(ggplot2::aes(x = 2.6, xend = 3.4, y = PairCount, 
                                         yend = PairCount, color = "Original.Data")) +
      ggplot2::scale_colour_manual(name = "", values=cols) +
      ggplot2::scale_fill_manual(name = "", values=fills)
    print(plt)
  }
  
  # Return summary.
  return(countDistribs)
}

#' Return the number of significant analytes / pairs per permutation and the number 
#' of permutations in which each analyte is significant.
#' If plot = TRUE, show a box plot of number of significant analytes over permutations, 
#' overlaid with the number of significant analytes in the original data.
#'
#' @param inputResults Data frame with model results (output of ProcessResults())
#' @param permResults An object of type PermutationResults (output of PermuteIntLIM())
#' @param plot Whether or not to show the boxplot. Default is TRUE.
#' @return A data frame that includes each significant pair from the unpermuted 
#' data and the number of times that pair was significant in the permuted data.
#' @export
PermutationPairSummary <- function(inputResults, permResults, plot){
  
  # Prevent "visible binding for global variable" notes.
  Pair <- Perm.Count <- NULL
  
  # Compute pair significance counts.
  allSignificantPairs <- do.call(c, permResults$listOfSigPairs)
  summaryCount <- table(allSignificantPairs)
  
  myres.sig.pairs <- paste(inputResults$Analyte1, inputResults$Analyte2, sep = "__V__")
  original.pairs.count <- unlist(lapply(myres.sig.pairs, function(pair){
    freq <- 0
    if(pair %in% allSignificantPairs){
      freq <- summaryCount[which(names(summaryCount) == pair)]
    }
    return(freq)
  }))
  original.pairs.df <- data.frame(Pair=myres.sig.pairs, 
                                  Perm.Count=original.pairs.count)
  original.pairs.df <- original.pairs.df[order(-original.pairs.df$Perm.Count),]
  original.pairs.df$Pair <- order(-original.pairs.df$Perm.Count)
  plt <- ggplot2::ggplot(data=original.pairs.df, ggplot2::aes(x = Pair, y = Perm.Count))+
    ggplot2::geom_bar(stat = "identity", width = 1) + 
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::labs(x="Significant Pair (Original Data)",
                  y="Number of Permutations Where Significant") +
    ggplot2::ylim(0,nrow(permResults$numSigPairs))
  print(plt)

  
  # Return summary.
  original.pairs.df$Pair <- myres.sig.pairs[order(-original.pairs.df$Perm.Count)]
  return(original.pairs.df)
}