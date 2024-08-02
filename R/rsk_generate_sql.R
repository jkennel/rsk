#===============================================================================
#' @title rsk_generate_sql
#'
#' @author Jonathan Kennel \email{jkennel@uoguelph.ca}
#'
#' @description Generate SQL string
#'
#'
#' @param db the database connection
#'
#'
#' @return SQL string for a query
#'
#' @export
#'
#===============================================================================
rsk_generate_sql <- function(db, start = NULL, end = NULL, by = NULL) {


  nm_tbl <- RSQLite::dbListTables(db)
  if (!'data' %in% nm_tbl) {
    return('')
  }

  # get column names
  col_names <- DBI::dbListFields(db, 'data')

  # remove datasetID column if present
  col_names <- col_names[!grepl('datasetID', col_names)]

  #nm <- c(paste0('`', col_names, '`'))
  sql_base <- paste0('SELECT ', paste(col_names, collapse = ', '),
                     ' FROM data')

  # no subset
  if (is.null(start) & is.null(end) & is.null(by)) {
    return(sql_base)
  }

  # subset dates
  if (!is.null(start) & !is.null(end) & is.null(by)) {
    return(paste0(sql_base, ' WHERE tstamp BETWEEN ', as.numeric(start) * 1000, ' AND ', as.numeric(end) * 1000))
  }

  # get start and end
  # this speeds up short time period subsetting
  if (is.null(start)) {
    start_id <- paste0("SELECT tstamp FROM data ORDER BY tstamp LIMIT 1")
  } else {
    start_id <- paste0("SELECT tstamp FROM data WHERE tstamp >= ", as.numeric(start) * 1000, " ORDER BY tstamp LIMIT 1")
  }
  if(is.null(end)) {
    end_id <- paste0("SELECT tstamp FROM data ORDER BY tstamp DESC LIMIT 1")
  } else {
    end_id <- paste0("SELECT tstamp FROM data WHERE tstamp <= ", as.numeric(end) * 1000, " ORDER BY tstamp DESC LIMIT 1")
  }

  start_id <- as.numeric(RSQLite::dbGetQuery(db, start_id)$tstamp)
  end_id   <- as.numeric(RSQLite::dbGetQuery(db, end_id)$tstamp)

  n <- (end_id - start_id) %/% (by * 1000)

  # do first sub
  if (!is.null(start) | !is.null(end)){
    sql_str <- paste0('SELECT * FROM (',
                      sql_base,
                      ' WHERE tstamp BETWEEN ', start_id, ' AND ', end_id,
                      ') WHERE tstamp % ', by * 1000, ' = 0')

  } else {
    sql_str <- paste0(sql_base, paste0(' WHERE tstamp % ', by * 1000, ' = 0'))
  }

  return(sql_str)

}



#===============================================================================
#' @title rsk_generate_sql_times
#'
#' @author Jonathan Kennel \email{jkennel@uoguelph.ca}
#'
#' @description Generate SQL string to filter times
#'
#'
#' @return SQL string for a query
#'
#' @export
#'
#===============================================================================
rsk_generate_sql_times <- function(start = NULL, end = NULL, by = NULL, times = NULL) {

  # read data
  if (!is.null(start) & !is.null(end) & !is.null(by)) {

    if((end-start) / by > 1e5) {
      return("")
    }

    times <- seq.int(as.numeric(start), as.numeric(end), by = by)
  }

  # set times
  if(!is.null(times)) {
    if (length(times) > 1e5) {
      return("")
    }
    return(paste0(' WHERE tstamp IN (', paste(times*1000, collapse = ', '), ')'))
  }


  return("")
}



