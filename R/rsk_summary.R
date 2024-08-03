#===============================================================================
#' @title rsk_summary
#'
#' @author Jonathan Kennel \email{jkennel@uoguelph.ca}
#'
#' @description Quick summary of an .rsk file
#'
#'
#' @param file_name the *.rsk file to read
#'
#'
#' @return table with some information about the rsk file
#'
#' @export
#'
#===============================================================================
rsk_summary <- function(file_name) {

  db   <- dbConnect(SQLite(), file_name)

  if (.check_for_table(db, 'channels')) {
    channels <- collapse::qDT(dbGetQuery(db, "SELECT * FROM channels"))
  } else {
    channels <- NULL
  }

  if (.check_for_table(db, 'instruments')) {
    instruments  <- collapse::qDT(dbGetQuery(db, "SELECT * FROM instruments"))
  } else {
    instruments <- NULL
  }

  if (.check_for_table(db, 'continuous')) {
    continuous <- collapse::qDT(dbGetQuery(db, "SELECT * FROM continuous"))
  } else if (.check_for_table(db, 'schedules')) {
    continuous <- collapse::qDT(dbGetQuery(db, "SELECT * FROM schedules"))
  } else {
    continuous <- NULL
  }

  if (.check_for_table(db, 'data')) {

    start_id <- paste0("SELECT tstamp/1000 FROM data ORDER BY tstamp LIMIT 1")
    end_id   <- paste0("SELECT tstamp/1000 FROM data ORDER BY tstamp DESC LIMIT 1")
    start_id <- as.POSIXct(as.numeric(RSQLite::dbGetQuery(db, start_id)$tstamp),
                           origin = '1970-01-01',
                           tz = 'UTC')
    end_id   <- as.POSIXct(as.numeric(RSQLite::dbGetQuery(db, end_id)$tstamp),
                           origin = '1970-01-01',
                           tz = 'UTC')
    # count    <- paste0("SELECT COUNT(tstamp) FROM data")
    # count    <- RSQLite::dbGetQuery(db, count)[[1]]
  } else {
    start_id <- NULL
    end_id   <- NULL
  }

  dbDisconnect(db)


  return(
    data.table::data.table(file_path = file_name,
                         file_name = basename(file_name),
                         file_size = file.size(file_name),
                         serial = instruments[["serialID"]],
                         model = instruments[["model"]],
                         start_id,
                         end_id,
                         # n_records = count,
                         n_channels = nrow(channels),
                         measurement_interval = continuous[["samplingPeriod"]] / 1000
                         )

  )
}


.check_for_table <- function(db, tbl_name) {
  tbls <- RSQLite::dbListTables(db)
  tbl_name %in% tbls

}


# fn <- list.files('/media/jonathankennel/Seagate Expansion Drive/rbr_g360', full.names = TRUE)
# tmp <- list()
# for (i in seq_along(fn)) {
#   print(fn[i])
#   tmp[[i]] <- rsk_summary(fn[i])
# }
# tmp <- rbindlist(tmp, fill = TRUE)
# setkey(tmp, file_size)
