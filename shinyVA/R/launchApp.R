#' launches the shinyVA app
#'
#' @export launchApp
#'
#' @return shiny application object
#'
#' @example \dontrun {shinyVA::launchApp()}
#'
#'


# wrapper for shiny::shinyApp()
launchApp <- function() {
  shiny::shinyApp(ui = shinyAppUI, server = shinyAppServer)
}
