#' Launch openVA shiny app in internet browser
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(shinyVA)
#' launchApp()
#' }
launchApp <- function() {
  shiny::runApp(appDir = system.file('app', package = 'openVAapp'))
}
