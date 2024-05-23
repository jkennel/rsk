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
#'
#'
#' @return table with some information about the rsk file
#'
#' @export
#'
#===============================================================================
read_rsk <- function(file_name, ...) {

  transducer_data <- list()
  for (i in seq_along(file_name)) {
    transducer_data[[i]] <- Rsk::new(file_name[i], ...)
  }

  return(transducer_data)

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
