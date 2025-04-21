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
adjust_over_time_period <- function(x, start, end, id = NULL) {

  air_adj <- x[between(datetime, start_air, end_air), .(mean_val = mean(value)), by = file_name]

  if (is.null(id)) {
    mean_val_all <- mean(air_adj$mean_val)
    air_adj[, sh := mean_val - mean_val_all]
  } else {
    air_adj[, sh := mean_val - mean_val[file_name == id]]

  }

  x[air_adj, value_adj := value - sh, on = "file_name"]

}
