#' Interpolate point on polygon
#' @param lowest_up_intervals data frame of filtered polygon segments are the
#' nearest to the y axis.
#' @param polygon data frame of a polygon
#' @author Liudas Daumantas
#' @import sp
#' @NoRd

.y.online <- function(x,y,x3) {
  ind<-.binar_search(x,x3)
  y[ind[1]]+(x3-x[ind[1]])*(y[ind[2]]-y[ind[1]])/(x[ind[2]]-x[ind[1]])
}
