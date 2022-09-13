#' run shiny app
#' @param port set port
#' @return No return value, called for side effects
#' @export
runIntLIMApp <- function(port="127.0.0.1") {
    appDir <- system.file("shinyApp", package = "IntLIM")
    if (appDir == "") {
        stop(" The ShinyApp directory was not found.
             Try re-installing `IntLIM`.",
             call. = FALSE)
    }

    shiny::runApp(appDir, display.mode = "normal", host = getOption("shiny.host", port))
}
