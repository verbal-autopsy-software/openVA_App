#' Launch openVA shiny app in internet browser
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(openVAapp)
#' launchApp()
#' }
launchApp <- function() {
  origMaxRS <- getOption("shiny.maxRequestSize")
  origWidth <- getOption("width")
  old <- options(shiny.maxRequestSize = origMaxRS, width = origWidth)
  options(shiny.maxRequestSize = 100*1024^2, width = 100)
  shiny::runApp(appDir = system.file('app', package = 'openVAapp'))
  on.exit(options(old), add = TRUE)
}
