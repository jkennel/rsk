#===============================================================================
#' @title read_rsk
#'
#' @author Jonathan Kennel \email{jkennel@uoguelph.ca}
#'
#' @description read an *.rsk file or collection of *.rsk files
#'
#'
#' @inheritParams read_rbr
#' @param file_name the *rsk file to read
#' @param return_data_table return a data.table of results?
#' @param include_params vector of parameter names contained in *rsk file to
#'   include in returned data.table.
#' @param simplify do you want the final temperature and pressure results only.
#'   only works when return_data_table is TRUE
#'
#'
#' @return table with some information about the rsk file
#'
#' @export
#'
#===============================================================================
read_rsk <- function(file_name,
                     return_data_table = FALSE,
                     include_params = NULL,
                     simplify = TRUE,
                     ...) {

  transducer_data <- list()

  for (i in seq_along(file_name)) {
    transducer_data[[i]] <- Rsk$new(file_name[i], ...)
  }

  if (!return_data_table) {
    return(transducer_data)
  }


  if (is.null(include_params)) {
    return(rbindlist(lapply(transducer_data, function(x) {
      data.table::melt(rename_data(x$data), id.vars = "datetime")
    })))
  }


  result <- list()
  for (i in seq_along(transducer_data)) {

    data_renamed <- rename_data(transducer_data[[i]][["data"]])

    # remove non-final columns
    if (simplify) {
      nms <- names(data_renamed)
      if ("pressure_compensated" %in% nms) {
        data_renamed[, pressure := NULL]
        setnames(data_renamed, "pressure_compensated", "pressure", skip_absent = TRUE)
      }
      if ("temperature" %in% nms & "temperature_onboard" %in% nms) {
        data_renamed[, temperature_onboard := NULL]
      }
      setnames(data_renamed, "temperature_onboard", "temperature", skip_absent = TRUE)
    }

    data <- data.table::melt(data_renamed, id.vars = "datetime")

    for (j in seq_along(include_params)) {
      data.table::set(data, j = include_params[j], value = transducer_data[[i]][[include_params[j]]])
    }

    result[[i]] <- data
  }

  return(data.table::rbindlist(result))

}

#===============================================================================
#' @title get_rsk_data
#'
#' @author Jonathan Kennel \email{jkennel@uoguelph.ca}
#'
#' @description read an *.rsk file or collection of *.rsk files
#'
#'
#' @inheritParams read_rbr
#' @param rsk the rsk collection
#'
#'
#' @return data.table of rsk results
#'
#' @export
#'
#===============================================================================
get_rsk_data <- function(rsk) {

  out <- list()

  for (i in seq_along(rsk)) {

    serial <- rsk[[i]]$serial
    file_name <- rsk[[i]]$file_name
    d <- rename_data(rsk[[i]]$data)
    d <- melt(d, id.vars = "datetime")
    d[, serial := serial]
    d[, file_name := file_name]

    out[[i]] <- d
  }

  rbindlist(out)
}
