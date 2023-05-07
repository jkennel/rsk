file_path <- "../../../Downloads/ruskin_instrument.log"
#' read_instrument_log
#'
#' @param file_path the path of the log file "ruskin_instrument.log"
#'
#' @return table containing the transducer settings
#' @export
#'
read_instrument_log <- function(file_path) {

  find_nearest_text <- function(txt, to_merge, to_find) {

    found_ids <- grep(to_find, txt, fixed = TRUE)
    ids <- data.table(found_ids, found_ids)

    setnames(ids, c(to_find, "found_ids"))
    setkeyv(ids, to_find)
    line_txt <- txt[unlist(ids[to_merge, roll = -Inf][["found_ids"]], use.names = FALSE)]
    trimws(sapply(line_txt, function(x) strsplit(x, ":")[[1]][2]))

  }

  txt <- readLines(file_path)
  headings <- grep("^^^", txt, fixed = TRUE)

  # found_instrument <- grep("FOUND INSTRUMENT", txt, fixed = TRUE)
  # end_of_logger <- grep("END OF LOGGER", txt, fixed = TRUE)
  # starting_instrument <- grep("STARTING INSTRUMENT", txt, fixed = TRUE)

  dat <- data.table(row_id = grep("STARTING INSTRUMENT", txt, fixed = TRUE))
  setkey(dat, row_id)
  dat_ids <- copy(dat)

  dat[, ruskin_version := find_nearest_text(txt, dat_ids, "RUSKIN VERSION")]
  dat[, model := find_nearest_text(txt, dat_ids, "Model")]
  dat[, serial := find_nearest_text(txt, dat_ids, "Serial")]
  dat[, version := find_nearest_text(txt, dat_ids, "Version")]
  dat[, clock := find_nearest_text(txt, dat_ids, "Clock")]
  dat[, host := find_nearest_text(txt, dat_ids, "Host")]
  dat[, state := find_nearest_text(txt, dat_ids, "State")]
  dat[, start := find_nearest_text(txt, dat_ids, "Start")]
  dat[, end := find_nearest_text(txt, dat_ids, "End")]
  dat[, gate := find_nearest_text(txt, dat_ids, "Gate")]
  dat[, mode := find_nearest_text(txt, dat_ids, "Mode")]
  dat[, period := find_nearest_text(txt, dat_ids, "Period")]
  dat[, internal := find_nearest_text(txt, dat_ids, "Internal")]
  dat[, external := find_nearest_text(txt, dat_ids, "External")]

  dat[, clock := as.POSIXct(clock, format = "%Y%m%d%H%M%S")]
  dat[, host := as.POSIXct(host, format = "%Y%m%d%H%M%S")]
  dat[, start := as.POSIXct(start, format = "%Y%m%d%H%M%S")]
  dat[, end := as.POSIXct(end, format = "%Y%m%d%H%M%S")]
  dat[, period := as.numeric(period)]
  dat[, battery_voltage := as.numeric(sub("V", "", internal))]

}






