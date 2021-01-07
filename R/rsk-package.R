#' rsk: A tool to work with pressure and temperature datasets from RBR
#' SQLite3 databases
#'
#'
#' You can learn about the rsk package in the vignettes:
#' \code{browseVignettes(package = "rsk")}
#'
#'
#'
#' @useDynLib rsk
#' @docType package
#' @aliases rsk-package
#'
#' @importFrom R6 R6Class
#'
#' @importFrom RSQLite SQLite
#' @importFrom RSQLite dbConnect
#' @importFrom RSQLite dbDisconnect
#' @importFrom RSQLite dbGetQuery
#'
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table setDT
#' @importFrom data.table melt
#'
#' @importFrom plotly plot_ly
#'
#' @importFrom Rcpp evalCpp
"_PACKAGE"
