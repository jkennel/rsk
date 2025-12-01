#' adjust_over_time_period
#'
#' @param x data.table with values column to adjust
#' @param start of period to use for adjustment
#' @param end of period to use for adjustment
#' @param id optional id to use for adjustment (baro, top transducer, etc.)
#'
#' @return data.table with adjusted values
#' @export
#'
adjust_over_time_period <- function(x,
                                    start,
                                    end,
                                    baro_file_name = NULL) {

  adj <- x[between(datetime, start, end),
               .(mean_val = mean(value, na.rm = TRUE)), by = file_name]

  if (is.null(baro_file_name)) {
    mean_val_all <- mean(adj$mean_val, na.rm = TRUE)
    adj[, sh := mean_val - mean_val_all]
  } else {
    adj[, sh := mean_val - mean_val[file_name == baro_file_name]]

  }

  x[adj, value_adj := value - sh, on = "file_name"]

}


#' add_baro_column
#' add a column of barometric pressure to the data.table
#'
#' @param x data.table with values column to adjust
#' @param baro_file_name file name for the baro data
#' @param remove_baro remove the baro data from the data.frame
#'
#' @return data.table with adjusted values
#' @export
#'
remove_baro <- function(x,
                        baro_file_name = NULL,
                        remove_baro = TRUE) {

  ba <- x[file_name == baro_file_name, .(datetime, baro = value)]
  if (remove_baro) {
    x  <- x[file_name != baro_file_name]
  }

  x[ba, on = "datetime"]


}
