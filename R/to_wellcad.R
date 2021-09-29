#===============================================================================
#' @title to_wellcad
#'
#' @description convert .rsk to .duckdb database
#'
#' @author Jonathan Kennel \email{jkennel@uoguelph.ca}
#'
#'
#' @param input the path of the .rsk input file
#' @param output the path of the .wcd output file
#'
#' @return NULL - export an .rsk file to .wcd (WellCAD)
#'
#' @export
#'
#===============================================================================
# stub to be filled in
to_wellcad <- function(input, output = NULL, ...) {

  # generate the output file name
  if(is.null(output)) {
    output <- gsub('.rsk', '.csv', x)
  }

  return(NULL)


}
