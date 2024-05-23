#' rename_data
#'
#' Renames the RBR column heading names
#'
#' @param x data.table with rbr names to rename
#'
#' @return table with new simplified names
#' @export
#'
#' @examples
rename_data <- function(x) {
  data(rbr_channels)
  nms <- names(x)
  nms <- intersect(nms, rbr_channels$channel_name)
  nms_dt <- data.table(channel_name = nms)

  nms_new <- rbr_channels[nms_dt, on = "channel_name"][["field_name"]]

  setnames(x, nms, nms_new)

  x

}
