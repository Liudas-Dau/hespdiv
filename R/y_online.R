#' Interpolate point on polygon
#' @param x a numeric vector sorted in increasing order (x coordinates)
#' @param y a numeric vector associated with x (y coordinates)
#' @param x3 a numeric value, in the range of x. y coordinate will be obtained
#' for this x coordiante by linear interpolation.
#' @author Liudas Daumantas
#' @noRd

.y.online <- function(x,y,x3) {
  ind<-.binar_search(x,x3)
  y[ind[1]]+(x3-x[ind[1]])*(y[ind[2]]-y[ind[1]])/(x[ind[2]]-x[ind[1]])
}
