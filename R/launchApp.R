#' Launch App function
#'
#' @return
#' @export
#'
#' @examples
launchApp <- function(){
  appDir <- system.file("app", package = "biomeshiny.package")
  if (appDir == "") {
    stop("Could not find app Try re-installing `biomeshiny.package`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal", port = 4800, launch.browser = FALSE)
}
