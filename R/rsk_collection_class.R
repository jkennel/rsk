#' @export
Rsk_collection <- R6Class("Rsk_collection",
                          public = list(

                            # database name
                            rsk_list     = NULL,
                            calibration_type = NULL,

                            air_calibration_start = NULL,
                            air_calibration_end = NULL,

                            blended_calibration_start = NULL,
                            blended_calibration_end = NULL,

                            shifts = NULL,

                          #' @description
                          #' Class to hold data from an .rsk or .duckdb file
                          #'
                          #' @param rsk_list
                          #'
                          #' @return
                          #' @export
                          #'
                          #' @examples
                          initialize = function(rsk_list
                          ) {
                            self$rsk_list <- rsk_list

                          },

                          calculate_shifts = function(
                          ) {

                            self$shifts <- data.table(file_name = self$get_file_names())

                            if (!is.null(self$air_calibration_start) & !is.null(self$air_calibration_end))
                            {

                            }

                            if (!is.null(self$blended_calibration_start) & !is.null(self$blended_calibration_end))
                            {

                            }

                          },

                          apply_shifts = function(
                          ) {

                          },

                          get_serial_numbers = function(
                          ) {
                            vapply(self$rsk_list, FUN = '[[', 'serial', FUN.VALUE = NA_real_)
                          },

                          get_file_names = function(
                          ) {
                            vapply(self$rsk_list, FUN = '[[', 'file_name', FUN.VALUE = NA_character_)
                          },

                          to_database = function(
                          ) {
                          })


)
