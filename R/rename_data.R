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

  type_desc <- data.table::fread(
    system.file("include/rbr_channels.csv", package = "rsk"))

  nms <- names(x)
  nms <- intersect(nms, type_desc$channel_name)
  nms_dt <- data.table(channel_name = nms)

  nms_new <- type_desc[nms_dt, on = "channel_name"][["field_name"]]
  # nms_new <- type_desc[channel_name %in% nms]$field_name
  setnames(x, nms, nms_new)

  x

}
