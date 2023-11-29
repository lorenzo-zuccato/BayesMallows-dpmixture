#' @export
leap_and_shift <- function(rho, leap_size = NULL){
  if(is.null(leap_size)){
    leap_size <- max(1L, floor(length(rho) / 5))
  }

  return(leap_and_shiftR(rho, leap_size))
}
