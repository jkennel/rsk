#===============================================================================
#' @title read_rsk
#'
#' @author Jonathan Kennel \email{jkennel@uoguelph.ca}
#'
#' @description read an *.rsk file or collection of *.rsk files
#'
#'
#' @param file_name the *rsk file to read
#' @param return_data_table return a data.table of results?
#' @param simplify_names simplify the final names?
#' @param include_params vector of parameter names contained in *rsk file to
#'   include in returned data.table.
#' @param ... arguments to pass to Rsk$new()
#'
#'
#' @return table with some information about the rsk file
#'
#' @export
#'
#===============================================================================
read_rsk <- function(file_name,
                     return_data_table = FALSE,
                     simplify_names = TRUE,
                     include_params = NULL,
                     ...) {

  transducer_data <- list()

  for (i in seq_along(file_name)) {
    transducer_data[[i]] <- Rsk$new(file_name[i], ...)
    if (simplify_names) {
      transducer_data[[i]]$simplify_names()
    }
  }


  if (!return_data_table) {
    return(transducer_data)
  }


  if (is.null(include_params)) {
    return(rbindlist(lapply(transducer_data, function(x) {
      data.table::melt(x$data, id.vars = "datetime")
    })))
  }


  result <- list()
  for (i in seq_along(transducer_data)) {


    data <- data.table::melt(transducer_data[[i]][["data"]], id.vars = "datetime")

    for (j in seq_along(include_params)) {
      data.table::set(data, j = include_params[j], value = transducer_data[[i]][[include_params[j]]])
    }

    result[[i]] <- data
  }

  return(data.table::rbindlist(result))

}
