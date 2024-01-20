#' @title
#'
#' @description This data.frame contains wavegroups for different data time
#' spans.  The wavegroups should be subset prior to use and the 'time' column
#' provides guidelines based on your input time span.
#'
#' @format A \code{data.table} The columns are:
#' \describe{
#'  \item{\code{channel_name}}{the rbr channel name}
#'  \item{\code{field_name}}{long name of channel with type}
#'  \item{\code{parameter}}{long name of channel}
#'  \item{\code{units}}{the units of the measuremnt}
#'  \item{\code{type}}{whether it is a primary or compensated value}
#' }
#'
#' @keywords internal
#' @examples
#' utils::data(rbr_channels)
"rbr_channels"

# rbr_channels <- data.table::fread(
#   system.file("include/rbr_channels.csv", package = "rsk"))
# usethis::use_data(rbr_channels,
#                   internal = FALSE,
#                   overwrite = TRUE)
