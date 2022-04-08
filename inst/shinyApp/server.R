# Set trace to FALSE.
options(shiny.trace=F)

# Called when the Shiny app is loaded
shinyServer(function(input, output, session) {  
  
  # Stop the app when the exit button is clicked.
  stop<-observe({
    if(input$sidebar == "stop") {
      stopApp(returnValue = invisible())
    }
  })
  
  # STEP 1
  
  # Function to upload the files.
  fixUploadedFilesNames <- function(x) {
    if (is.null(x)) {
      return()
    }
    oldNames = x$datapath
    newNames = file.path(dirname(x$datapath),x$name)
    file.rename(from = oldNames, to = newNames)
    x$datapath <- newNames
    x
  }
  
  # Display the output "filename" widget if all files are present, 
  # or print an error to the wiget if not.
  output$filename <- renderPrint({
    myFile <- fixUploadedFilesNames(input$file1) #fixing uploaded filenames.
    if (is.null(myFile)) {
      cat("Select all files by holding down Shift key and clicking 'Browse' above")
    } else if (length(myFile$name) !=6) {
      cat("Load all required files before continuing.")
    } else {
      paste(myFile$name,"- uploaded")
    }
  })

  # Display the 'id' text box widget for analyte type 1. 
  output$idChooseType1 <- renderUI({
    # If the file name has not been entered yet for Step 1, skip.
    if (is.null(input$file1$datapath)) {
    } else {
      # Read the "input.csv" file.
      myFile <- fixUploadedFilesNames(input$file1)
      indexfile = which(myFile$name == 'input.csv')
      if (indexfile == which(myFile$name == 'input.csv')) {
        file <- read.csv(myFile$datapath[indexfile])
        rownames(file) <- file$type
        
        # Ignore metadata.
        if (file["analyteType1MetaData", "filenames"] == "") {
          return()
        }
        
        # Create a text box where the user can enter their ID.
        textInput("analtyeType1id", "Analyte Type 1 ID", "")
      } else {
        tags$b(paste("Please rename input file as input.csv"))
      }
    }
  })
  
  # Display the 'id' text box widget for analyte type 2. 
  output$idChooseType2 <- renderUI({
    # If the file name has not been entered yet for Step 1, skip.
    if (is.null(input$file1)) {
    } else {
      # Read the "input.csv" file.
      myFile <- fixUploadedFilesNames(input$file1)
      indexfile = which(myFile$name == 'input.csv')
      if (indexfile == which(myFile$name == 'input.csv')) {
        file <- read.csv(myFile$datapath[indexfile])
        rownames(file) <- file$type
        
        # Ignore metadata.
        if (file["analyteType2MetaData", "filenames"] == "") {
          return()
        }
        
        # Create a text box where the user can enter their ID.
        textInput("analyteType2id", "Analyte Type 2 ID", "")
      } else {
        tags$b(paste("Please rename input file as input.csv"))
      }
    }
  })
  
  # Function to reading the data using the ReadData function.
  multiData <- eventReactive(input$run, {
    myFile <- fixUploadedFilesNames(input$file1)
    print(length(myFile))
    indexfile = which(myFile$name == 'input.csv')
    IntLIM::ReadData(req(myFile$datapath[[indexfile]]),
                     input$analyteType1id, input$analyteType2id)
  })
  
  # Display the 'stats' widget using the ShowStats function.
  output$stats<-renderDataTable({
        table<- as.data.frame(t(IntLIM::ShowStats(multiData())))
        colnames(table)<-"value"
        cbind(names=rownames(table),table)
  },options = list(dom = 't'))

  # Display the 'plot' widget using the PlotDistributions function.
  output$plot<-renderUI(
        IntLIM::PlotDistributions(multiData())
  )
    
  # STEP 2
  
  # Function to filter the data using the FilterData function.
  FmultiData<-eventReactive(input$run2,{
    if(input$run2==0){
      FmultiData<-multiData()
    }
    if(input$run2!=0){
      FmultiData<-IntLIM::FilterData(multiData(),
                                     analyteType1perc=input$analyteType1perc,
                                     analyteType2perc=input$analyteType2perc,
                                     analyteMiss=input$analyteMiss,
                                     cov.cutoff=input$cov.cutoff)
    }
    FmultiData
  },ignoreNULL=FALSE)
    
  # Display the 'Ostats' widget using the ShowStats function.
  output$Ostats<-renderDataTable({
    if(input$run2==0) return()
    table<- as.data.frame(t(IntLIM::ShowStats(multiData())))
    colnames(table)<-"value"
    cbind(names=rownames(table),table)
  },options = list(dom = 't'))
    
  # Display the 'Oplot' widget using the PlotDistributions function.
  output$Oplot<-renderUI({
    if(input$run2==0) return()
    IntLIM::PlotDistributions(multiData())
  })
  
  # Display the 'Fstats' widget using the ShowStats function on the
  # filtered data.
  output$Fstats<-renderDataTable({
    if(input$run2==0) return()
    table<- as.data.frame(t(IntLIM::ShowStats(FmultiData())))
    colnames(table)<-"value"
    cbind(names=rownames(table),table)
  },options = list(dom = 't'))
    
  # Display the 'Fplot' widget using the PlotDistributions function
  # on the filtered data.
  output$Fplot<-renderUI({
    if(input$run2==0) return()
    IntLIM::PlotDistributions(FmultiData())
  })
    
  # Save the file.
  output$downloadFdata <- downloadHandler(
    filename = "Filtered_data.zip",
    content = function(con) {
      IntLIM::OutputData(FmultiData(),con)
    }
  )
    
  # Step 3
  
  # Choose the outcome type using the 'choosestype' widget.
  output$choosestype <- renderUI({
    choice<-reactive({
      return(colnames(FmultiData()@sampleMetaData))
    })
    selectInput("stype", "Sample Type:", c(Choose='',choice()),selected = NULL)
  })
    
  # Function to obtain results using RunIntLim.  
  continuous<-reactive(input$continuous)
  myres <- eventReactive(input$run3,{
    shinyjs::html("text", "")
    IntLIM::RunIntLim(FmultiData(),stype=input$stype,outcome=1, independent.var.type=2,
                      continuous=continuous())
  })
  rsquared<-reactive(input$rsquared1)
  pvalcutoff<-reactive(input$pvalcutoff1)
  interactionCoeff<-reactive(input$interactionCoeff)
    
  # Display the 'volcanoPlot' widget using pvalCoefVolcano.
  output$volcanoPlot<-renderPlot({
    IntLIM::pvalCoefVolcano(myres(),FmultiData(),input$nrpoints,pvalcutoff(),
                            interactionCoeff())},
    height=500
  )
  
  # Display the 'Pdist' widget using DistPvalues.
  output$Pdist<-renderPlot({
    if(!is.null(myres())){
      IntLIM::DistPvalues(myres(),breaks = input$breaks, adjusted = FALSE)
    }
  })
    
  # Display the 'Ptext' widget, which provides information about the meaning
  # of the p-values.
  output$Ptext<-renderPrint(
      (paste("Distribution of unadjusted p-values (a peak close to zero",
      "suggests that there are significant analyte pairs that are found)."))
  )
    
  # Step 4
    
  # Process the results using the cutoffs.
  rsquared2<-reactive(input$rsquared2)
  pvalcutoff2<-reactive(input$pvalcutoff2)
  # We make the interaction coefficient reactive so that the interaction plot
  # does not show up until the user clicks "Run".
  interactionCoeff2 <- eventReactive(input$run4, {
    input$interactionCoeff2
  })
  myres2 <- eventReactive(input$run4,{
    IntLIM::ProcessResults(myres(),FmultiData(),pvalcutoff=pvalcutoff2(),
                           rsquared=rsquared2(),
                           interactionCoeffPercentile=interactionCoeff2())
  })
        
  # Display the 'betagraph' widget using the InteractionCoefficientGraph function.
  output$betagraph<-renderPlot({
    # Grabbing FmultiData() data leads to plotting being delayed until button is pressed.
    fmdat <- FmultiData()
    IntLIM::InteractionCoefficientGraph(inputResults = myres(),
                                        interactionCoeffPercentile=interactionCoeff2(),
                                        percentageToPlot = 0.01, 
                                        independent.var.type = 1,
                                        outcome = 2)
  })
  
  # Download the betagraph as a PDF.
  output$downloadplot1 <- downloadHandler(
    filename = "plot.pdf",
    content = function(file) {
      pdf(file, width=12, height=6.3)
      IntLIM::InteractionCoefficientGraph(myres(),interactionCoeff())
      dev.off()
    }
  )
    
  # Function to generate beta graph.
  beta_graph <- function() {
    IntLIM::InteractionCoefficientGraph(myres(),interactionCoeff())
  }

  # Download the results as a CSV file using the OutputResults function.
  output$downloadData <- downloadHandler(
    filename = "results.csv",
    content = function(con) {
      IntLIM::OutputResults(req(myres2()),con)
    }
  )
    
  # Step 5

  # List table of pairs in 'table' widget.
  output$table<-DT::renderDT(
    myres2()
  )
  
  # Function to select rows from possible pairs.
  scatterrows<-eventReactive(input$run5,{
      input$table_rows_selected
  })
  
  output$scatterplot<-renderPlot({
    # Obtain pairs.
    a<-as.matrix(scatterrows())
    pair1<-as.matrix(myres2()[a[1,],])
    independentAnalyteOfInterest1<-pair1[,"Analyte1"]
    outcomeAnalyteOfInterest1<-pair1[,"Analyte2"]
    
    # Plot.
    IntLIM::PlotPairFlat(FmultiData(),
                     myres(),
                     independentAnalyteOfInterest=independentAnalyteOfInterest1,
                     outcomeAnalyteOfInterest=outcomeAnalyteOfInterest1,
                     outcome = 1,
                     independentVariable = 2)
  })
  
  # Information and status boxes used throughout app.
    
  # Step 1 box
  output$statusbox1 <- renderInfoBox({
      if (is.null(input$file1)) {
          infoBox(
              "Status",
              "File Not Loaded Yet!",
              icon = icon("import", lib = "glyphicon"),
              color = "aqua",
              fill = TRUE
          )
      }
      else if (!is.null(input$file1)&&input$run==0) {
          infoBox(
              "Status",
              "Step 1 is Not Complete Yet!",

              "Press Run button",

              icon = icon("warning-sign", lib = "glyphicon"),
              color = "aqua",
              fill = TRUE
          )
      }
      else if (!input$run==0) {
          infoBox(
              "Status",
              HTML(paste("Data is loaded.",
                         "You can proceed to Step 2 (optional) or Step 3.",
                         sep = "<br/>")),
              icon = icon("thumbs-up", lib = "glyphicon"),
              color = "green", fill = TRUE
          )
      }
  })
    
  # Step 2 box
  output$statusbox2 <- renderInfoBox({
      if (input$analyteType1perc==0&&input$analyteType2perc==0) {
          infoBox(
              "Status",
              "Please provide input",
              icon = icon("flag", lib = "glyphicon"),
              color = "aqua",
              fill = TRUE
          )
      }
      else if (input$run2==0) {
          infoBox(
              "Status",
              "Press Run button",
              icon = icon("flag", lib = "glyphicon"),
              color = "aqua",
              fill = TRUE
          )
      }
      else if (!input$run2==0) {
          infoBox(
              "Status",
              HTML(paste("Data filtering is complete.",
                         "You can proceed to Step 3",
                         sep = "<br/>")),
              icon = icon("thumbs-up", lib = "glyphicon"),
              color = "green", fill = TRUE)
      }
  })
    
  # Step 3 box
  output$statusbox3 <- renderInfoBox({
      if (input$stype == "") {
          infoBox(
              "Status",
              "Please select your sample type",
              icon = icon("flag", lib = "glyphicon"),
              color = "aqua",
              fill = TRUE
          )
      }
      else if (input$run3==0) {
          infoBox(
              "Status",
              "Step 3 is Not Complete Yet!",
              "Press Run button",
              "This function can take several minutes, please be patient",
              icon = icon("warning-sign", lib = "glyphicon"),
              color = "aqua",
              fill = TRUE
          )
      }
      else if (!is.null(myres())) {
          infoBox(
              "Status",
              HTML(paste("IntLIM models are calculated.",
                         "You can proceed to Step 4",
                         sep = "<br/>")),
              icon = icon("thumbs-up", lib = "glyphicon"),
              color = "green", fill = TRUE
          )
      }
  })
    
  # Step 4 box
  output$statusbox4 <- renderInfoBox({
    if (input$run4==0) {
      infoBox(
        "Status",
        "Step 4 is Not Complete Yet!",
        "Press Run button",
        "This function will perform the data filtering for Step 5",
        icon = icon("warning-sign", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE
      )
    }
    else if (!is.null(myres2())) {
          infoBox(
              "Status",
              HTML(paste("Results processed.",
                         "You can proceed to Step 5",
                         sep = "<br/>")),
              icon = icon("thumbs-up", lib = "glyphicon"),
              color = "green", fill = TRUE
          )
      } 
  })
  
  # Step 5 box
  output$statusbox5 <- renderInfoBox({
      if (is.null(input$table_rows_selected)) {
          infoBox(
              "Status",
              "Please choose a pair of analytes by clicking the table",
              icon = icon("flag", lib = "glyphicon"),
              color = "aqua",
              fill = TRUE
          )
      }
      else if (input$run5==0) {
          infoBox(
              "Status",
              "Press Run",
              icon = icon("flag", lib = "glyphicon"),
              color = "aqua",
              fill = TRUE
          )
      }
      else if (!is.null(scatterrows())) {
          infoBox(
              "Status",
              HTML(paste("Scatter plot running complete!",
                         
                         sep = "<br/>")),
              icon = icon("thumbs-up", lib = "glyphicon"),
              color = "green", fill = TRUE
          )
      }
  })
})

