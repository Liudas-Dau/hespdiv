#' Check whether points fall inside a polygon
#'
#' This function is identical replica of point.in.poly function
#' from 'sp' package
#' @param point.x numeric vector of point x coordinates
#' @param point.y numeric vector of point y coordinates
#' @param pol.x numeric vector of polygon x coordinates
#' @param pol.y numeric vector of polygon y coordinates
#' @param mode.checked default FALSE, used internally to save time when all the other argument are known to be of storage mode double
#' @return integer array; values are: 0: point is strictly exterior to pol; 1: point is strictly interior to pol; 2: point lies on the relative interior of an edge of pol; 3: point is a vertex of pol.
#' @references
#' Uses the C function InPoly(). InPoly is Copyright (c) 1998 by Joseph O'Rourke. It may be freely redistributed in its entirety provided that this copyright notice is not removed
#' Pebesma E, Bivand R (2005). “Classes and methods for spatial data in R.” _R News_, *5*(2), 9-13. <https://CRAN.R-project.org/doc/Rnews/>.
#' Bivand R, Pebesma E, Gomez-Rubio V (2013). _Applied spatial data analysis with R, Second edition_. Springer, NY. <https://asdar-book.org/>.
#' @examples
#'#Creating data.frame of a polygon
#'# open polygon:
#'.point.in.polygon(1:10,1:10,c(3,5,5,3),c(3,3,5,5))
#'# closed polygon:
#'.point.in.polygon(1:10,rep(4,10),c(3,5,5,3,3),c(3,3,5,5,3))
#' @useDynLib hespdiv, R_point_in_polygon_sp
#' @noRd
.point.in.polygon = function(point.x, point.y, pol.x, pol.y,
    mode.checked=FALSE) {
    if (mode.checked) res <- .Call(R_point_in_polygon_sp, point.x,
        point.y, pol.x, pol.y)
    else res <- .Call(R_point_in_polygon_sp, as.numeric(point.x),
        as.numeric(point.y), as.numeric(pol.x), as.numeric(pol.y))
    res
}
