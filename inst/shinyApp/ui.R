headerbar <- shinydashboard::dashboardHeader(
    title = "IntLIM",
    titleWidth = 270,
    shinydashboard::dropdownMenu(
        type = "notifications",
        shinydashboard::notificationItem(
            text = "Plots might take some time to display",
            icon("truck"),
            status = "warning"
        )
    )
)

sidebar <- shinydashboard::dashboardSidebar(
    width = 270,
    shinydashboard::sidebarMenu(
        id="sidebar",
        shinydashboard::menuItem("About",
                 tabName = "about",
                 icon = icon("info")),
        shinydashboard::menuItem(
            "Load Data",
            tabName = "loaddata",
            icon = icon("folder-open"),
            badgeLabel = "step 1"
        ),
        shinydashboard::menuItem(
            "Filter Data (optional)",
            tabName = "Filterdata",
            icon = icon("bullseye"),
            badgeLabel = "step 2"
        ),
        shinydashboard::menuItem(
            "Run Linear Models",
            tabName = "RunLM",
            icon = icon("bullseye"),
            badgeLabel = "step 3"
        ),
        shinydashboard::menuItem(
            "Process result",
            tabName = "processresult",
            icon = icon("bullseye"),
            badgeLabel = "step 4"
        ),
        shinydashboard::menuItem(
            "Scatter plot",
            tabName = "scatterplot",
            icon = icon("bullseye"),
            badgeLabel = "step 5"
        ),
        shinydashboard::menuItem(
            "Exit app",
            tabName = "stop",
            icon = icon("sign-out"),
            badgeLabel = "exit"
        )
    )
)

body <- shinydashboard::dashboardBody(

    shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "about",
                             fluidRow(
                                             shinydashboard::box(
                                                width = 12,
                                                solidHeader = TRUE,
                                                includeMarkdown("README.md")
                                              )
                              )
        ),
        shinydashboard::tabItem(tabName = "loaddata",
                fluidRow(
                    shinydashboard::box(
                        title = strong("Load Data"),
                        width = 8,
                        solidHeader = TRUE,
                        h5("This step takes all the relevant CSV files as input, including the following (See About for more details):"),
                  			tags$ul(
                  				tags$li("input.csv (required): contains the names of all files input (See About)"),
                  				tags$li("analyteType1Data (required): rows are analytes of the first type, columns are samples; the first row is assumed to have sample ids and these ids should be unique; the first column is assumed to have feature ids and those should be unique."),
                  				tags$li("sampleMetaData (required): rows are samples, features are columns"),
                  				tags$li("analyteType2Data (required): rows are analytes of the second type, columns are samples; the first row is assumed to have sample ids and these ids should be unique; the first column is assumed to have feature ids and those should be unique."),
                  				tags$li("analyteType1MetaData (optional): rows are analytes of the first type, features are columns"),
                      				tags$li("analyteType2MetaData (optional): rows are analytes of the second type, features are columns")
                  			)
                    ),
                    shinydashboard::box(
                        width = 4,
                        infoBoxOutput("statusbox1", width = NULL)
                    )
                ),#end of info flow
                fluidRow(
                    shinydashboard::box(
                      fileInput("file1", "Select all CSV Files (input, analyteType1Data, analyteType2Data, and sampleMetaData are required)", multiple = TRUE, accept = ".csv"),
                        verbatimTextOutput('filename'),
                        uiOutput('idChooseType1'),
                        uiOutput('idChooseType2'),
                        hr(),
                        actionButton("run", "Run"),
                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                         tags$div("Loading...",id="loadmessage"))
                    ),
                    shinydashboard::box(
                        tags$b("Summary Statistics"),
                        pre(dataTableOutput('stats')),
                        tags$style(type="text/css", '#stats tfoot {display:none;}')
                    )
                ),#end of select file and stats flow
                fluidRow(
                    shinydashboard::box(
                        width = 12,
                        tags$b("Distribution of Input Data"),
                        pre(htmlOutput("plot"))
                    )
                )#end of plot flow
        ), # end tab loaddata
        shinydashboard::tabItem(tabName = "Filterdata",
                fluidRow(
                    shinydashboard::box(
                        title = strong("Filter Data (optional)") ,
                        width = 8,
                        solidHeader = TRUE,
                        h5("This step allows you to filter the data by a user-defined percentile cutoff."),
                        hr(),
                        numericInput("analyteType1perc", "percentile cutoff (0-1) for filtering analyte type 1 (e.g. remove analytes with mean values < cutoff):", 0, min = 0, max = 1),
                        numericInput("analyteType2perc", "percentile cutoff (0-1) for filtering analyte type 2 (e.g. remove analytes with mean values < cutoff):", 0, min = 0, max = 1),
                        numericInput("analyteMiss", "missing value percent cutoff (0-1) for filtering analytes (e.g. analytes with > % cutoff missing values will be removed)", 0,min=0,max=1),
                        numericInput("cov.cutoff", "coefficient of variation percentile cutoff for filtering analytes (e.g. analytes with cov > cutoff)", 0,min=0,max=1),
                        actionButton("run2", "Run")
                    ),
                    shinydashboard::box(
                        width=4,
                        infoBoxOutput("statusbox2", width = NULL),
                        downloadButton('downloadFdata', 'Download')

                    )
                ), # end filter option flow
                fluidRow(
                    shinydashboard::box(
                        tags$b("The statistic summary of origin data"),
                        dataTableOutput('Ostats'),
                        tags$style(type="text/css", '#Ostats tfoot {display:none;}')
                    ),
                    shinydashboard::box(
                        tags$b("The statistic summary of filtered data"),
                        dataTableOutput('Fstats'),
                        tags$style(type="text/css", '#Fstats tfoot {display:none;}')
                    )
                ),#end stats comparison flow
                fluidRow(
                    shinydashboard::box(
                        tags$b("The distribution of the origin data."),
                        uiOutput('Oplot')
                    ),
                    shinydashboard::box(
                        tags$b("Verify the distribution of the filtered data."),
                        uiOutput('Fplot')
                    )
                )#end plot comparison flow
        ),
        shinydashboard::tabItem(tabName = "RunLM",
                fluidRow(
                    shinydashboard::box(
                        title = strong("Run IntLIM") ,
                        width = 6,
                        height =800,
                        solidHeader = TRUE,
                        h5("This step performs the linear models for all combinations of analyte pairs and then plots distribution of p-values."),
                  			h5("The linear model performed is 'a_i ~ a_j + p + a_j:p' where "),
                  			tags$ul(
                  				tags$li("'a_i' is the outcome analyte level (may be of types 1 or 2)"),
                  				tags$li("a_j is the independent analyte level (may be of types 1 or 2"),
                  				tags$li("p is the phenotype (e.g. tumor vs non-tumor)"),
                  				tags$li("a_j:p is the interaction between phenotype and independent analyte level")
                  			),
                  			h5("A statistically significant p-value of the the interaction term a_j:p indicates that the analyte pair relationship is phenotype-specific. Please see manuscript for more details."),
                        hr(),
                        uiOutput('choosestype'),
                        numericInput("nrpoints", "number of points to be plotted in lowest density areas:", 10000, min = 0, max = 30000),
                        numericInput("pvalcutoff1","cutoff of FDR-adjusted p-value for filtering(0 - 1) :", 0.05, min = 0, max = 1),
                        numericInput("rsquared1", "cutoff of R^2 value for filtering (0-1):", 0.5, min = 0, max = 1),
                  		numericInput("interactionCoeff", "cutoff of interaction coefficient percentile (0-1):", 0.9, min = 0, max = 1),
                  		checkboxInput("continuous", "stype is continuous", value = FALSE, width = NULL),
                  		actionButton("run3", "Run")
                    ),
                    shinydashboard::box(
                        width = 6,
                        height =750,
                        numericInput("breaks","Breaks of histogram",100,min=10,max = 500),
                        plotOutput("Pdist"),
                        hr(),
                        textOutput("Ptext")
                    )
                ),
                fluidRow(
                    shinydashboard::box(
                        width = 15,
                        height = 550,
                        tags$head(tags$style(type="text/css", "
                        loadmessage {
                                         position: fixed;
                                         top: 0px;
                                         left: 0px;
                                         width: 100%;
                                         padding: 5px 0px 5px 0px;
                                         text-align: center;
                                         font-weight: bold;
                                         font-size: 100%;
                                         color: #000000;
                                         background-color: #CCFF66;
                                         z-index: 105;
                                         }
                                         ")),
                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                         tags$div("Loading...(It might take several minutes depending on the size of the dataset,please be patient!)",id="loadmessage")),
    
                        plotOutput("volcanoPlot")
                    ),
                    shinydashboard::box(
                      width = 4,
                      infoBoxOutput("statusbox3", width = NULL)
                    )
                )
        ),
        shinydashboard::tabItem(tabName = "processresult",
                fluidRow(
                    shinydashboard::box(
                        title = strong("Process the result") ,
                        width = 8,
                        height = 200,
                        solidHeader = TRUE,
                        h5("Process the results and filter pairs of analytes based on
                           adjusted p-values, R^2 values, and interaction coefficient cutoffs between the two groups being compared."),
                        h5("Then plot beta graph of significant gene-metabolite pairs by filling out parameters below and clicking 'Run'.")
                    ),
                    shinydashboard::box(
                        width = 4,
                        height = 200,
                        infoBoxOutput("statusbox4", width = NULL),
                        downloadButton('downloadData', 'Download')
                    )
                ),#end of info floww
                fluidRow(
                    shinydashboard::box(
                        width = 4,
                        uiOutput("numericChoice1"),
                        uiOutput("numericChoice2"),
                        numericInput("pvalcutoff2","cutoff of FDR-adjusted p-value for filtering(0 - 1) :", 0.05, min = 0, max = 1),
                        numericInput("rsquared2", "cutoff of R^2 value for filtering (0-1):", 0.5, min = 0, max = 1),
                        numericInput("interactionCoeff2", "cutoff of interaction coefficient percentile (0-1):", 0.9, min = 0, max = 1),
                        actionButton("run4", "Run")
                    ),
                    shinydashboard::box(
                        width = 8,
                        tags$head(tags$style(type="text/css", "
                          loadmessage {
                                           position: fixed;
                                           top: 0px;
                                           left: 0px;
                                           width: 100%;
                                           padding: 5px 0px 5px 0px;
                                           text-align: center;
                                           font-weight: bold;
                                           font-size: 100%;
                                           color: #000000;
                                           background-color: #CCFF66;
                                           z-index: 105;
                                           }
                                           ")),
                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                         tags$div("Loading...",id="loadmessage")),
                        downloadButton('downloadplot1',label='Download graph as pdf'),
                        plotOutput("betagraph")
                    )
                )
        ),
        shinydashboard::tabItem(tabName = "scatterplot",
                                fluidRow(
                                    shinydashboard::box(
                                        title = strong("Scatter plot") ,
                                        width = 8,
                                        solidHeader = TRUE,
                                        h5("This step presents the table of analyte pairs that are significant."),
                                        h5("You can plot the scatter plot of preferred analyte pairs by clicking table")
                                    ),
                                    shinydashboard::box(
                                        width = 4,
                                        infoBoxOutput("statusbox5", width = NULL)
                                    )
                                ),# end of info flow
                                fluidRow(
                                    shinydashboard::box(
                                        width = NULL,
                                        tags$b("Significant pairs"),
                                        DT::DTOutput('table'),
                                        hr(),
                                        actionButton("run5", "Run"),
                                        tags$head(tags$style(type="text/css", "
                                                             loadmessage {
                                                             position: fixed;
                                                             top: 0px;
                                                             left: 0px;
                                                             width: 100%;
                                                             padding: 5px 0px 5px 0px;
                                                             text-align: center;
                                                             font-weight: bold;
                                                             font-size: 100%;
                                                             color: #37649b;
                                                             background-color: #CCFF66;
                                                             z-index: 105;
                                                             }
                                                             ")),
                                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                         tags$div("Loading...",id="loadmessage"))
                                    )
                                ),#end of table flow
                                fluidRow(
                                    shinydashboard::box(
                                        width = 12,
                                        pre(plotOutput("scatterplot")),
                                    )
                                )#end of scatterplot flow
        )
    )
)
shinyUI(fluidPage(
    tags$head(includeHTML(("www/google_analytics.html"))),
    shinydashboard::dashboardPage(
        headerbar,
        sidebar,
        body
    )
)) # end shinUI
