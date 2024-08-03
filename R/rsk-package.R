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
#' @importFrom RSQLite dbSendQuery
#' @importFrom RSQLite dbListFields
#' @importFrom RSQLite dbListTables
#' @importFrom RSQLite dbWriteTable
#'
#' @importFrom DBI dbWriteTable
#' @importFrom DBI dbConnect
#' @importFrom DBI dbDisconnect
#'
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table as.data.table
#' @importFrom data.table setDT
#' @importFrom data.table melt
#' @importFrom data.table fwrite
#' @importFrom data.table set
#' @importFrom data.table setkey
#' @importFrom data.table setkeyv
#' @importFrom data.table setnames
#'
#' @importFrom plotly plot_ly
#' @importFrom collapse collapse
#'
#' @importFrom duckdb duckdb
#'
#' @importFrom tools file_ext
#'
#'
#' @importFrom Rcpp evalCpp
"_PACKAGE"
