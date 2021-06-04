#===============================================================================
#' @title to_duckdb
#'
#' @description convert .rsk to .duckdb database
#'
#' @author Jonathan Kennel \email{jkennel@uoguelph.ca}
#'
#'
#' @param input the path of the .rsk input file
#' @param output the path of the .duckdb output file
#'
#' @return NULL - copies an .rsk file to .duckdb
#'
#' @export
#'
#===============================================================================
to_duckdb <- function(input, output = NULL) {


  # generate the output file name
  if(is.null(output)) {
    output <- gsub('.rsk', '.duckdb', x)
  }


  # open databases
  db_in   <- RSQLite::dbConnect(RSQLite::SQLite(), input)
  db_out  <- duckdb::dbConnect(duckdb::duckdb(), dbdir = output)
  nm_tbls <- dbListTables(db_in)


  # do not copy the downloads and downsample tables
  nm_tbls <- nm_tbls[!grepl('downloads', nm_tbls)]
  nm_tbls <- nm_tbls[!grepl('downsample', nm_tbls)]


  # copy the tables from the .rsk to the .duckdb
  for(i in seq_along(nm_tbls)) {
    tbl_dat <- setDT(dbGetQuery(db_in, paste0("SELECT * FROM ", nm_tbls[i])))

    if(nm_tbls[i] == 'data') {
      tbl_dat[, tstamp := as.POSIXct(as.numeric(tstamp * 1e-3),
                                     origin = '1970-01-01', tz = 'UTC')]
    }

    DBI::dbWriteTable(db_out, nm_tbls[i], tbl_dat, overwrite = TRUE)
  }


  dbDisconnect(db_in, shutdown=TRUE)
  dbDisconnect(db_out, shutdown=TRUE)


  return(NULL)

}



# library(duckdb)
# library(RSQLite)
# library(data.table)
#
# x1 <- '/home/jonathankennel/Storage/data/rbr/rd45a 081871_20191118_1213.rsk'
# x2 <- '/home/jonathankennel/Storage/data/rbr/rd45a 081871_20191118_1213.duckdb'
# system.time(to_duck(x1,x2))
#
# system.time({
#   db <- dbConnect(duckdb(), x2)
#   a <- setDT(dbGetQuery(db, "SELECT * FROM data"))
#   dbDisconnect(db, shutdown=TRUE)
# })
# x1 <- '/home/jonathankennel/Storage/data/rbr/rd130 077623_20191119_1500.rsk'
# x2 <- '/home/jonathankennel/Storage/data/rbr/rd130 077623_20191119_1500.duckdb'
# system.time(to_duck(x1,x2))
#
# system.time({
#   db <- dbConnect(duckdb(), x2)
#   a <- setDT(dbGetQuery(db, "SELECT * FROM data"))
#   dbDisconnect(db, shutdown=TRUE)
# })
