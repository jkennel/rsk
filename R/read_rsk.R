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
                     ...) {

  transducer_data <- list()

  for (i in seq_along(file_name)) {
    transducer_data[[i]] <- Rsk$new(file_name[i], ...)
  }

  if (is.null(return_data_table)) {
    return(transducer_data)
  }


  if (is.null(include_params)) {
    rbindlist(lapply(transducer_data, function(x) {
      data.table::melt(rename_data(x$data), id.vars = "datetime")
    }))
  }


  result <- list()
  for (i in seq_along(transducer_data)) {
    data <- data.table::melt(rename_data(transducer_data[[i]][["data"]]),
                            id.vars = "datetime")
    for (j in seq_along(include_params)) {
      set(data, j = include_params[i], value = transducer_data[[i]][[include_params[i]]])
    }

    result[[i]] <- data
  }

  return(rbindlist(result))

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
