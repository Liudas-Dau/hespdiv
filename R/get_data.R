#' Get data points that lie inside a polygon
#'
#' This function extracts points from a provided dataset that lie inside a given polygon.
#' @param polygon A data frame of 2 columns ("x","y") that contain coordinates of polygon vertices.
#' @param xy.dat A data frame, containing "x" and "y" columns. Other columns may also be present.
#' @return A filtered data frame.
#' @note This function excludes points that are strictly outside the given
#' polygon. Points lying on the edges or vertices of the
#' polygon (if the polygon is not closed) will be included in the
#' filtered data frame.
#' @author Liudas Daumantas
#' @examples #Creating data.frame of a polygon
#' poly<- data.frame(c(3.38,3.30,1.70,0.78,-0.06,-2.30,-2.94,-3.97,-1.61,-0.39,0.68,1.28,1.60,3.38),
#' c(-0.12,-0.31,-2.73,-3.22,-3.29,-2.19,-1.62,0.94,3.10,3.00,2.91,2.49,2.20,-0.12))
#'
#' #Creating a data set of points
#' xy.dat<-data.frame(X=runif(250,-4,4),Y=runif(250,-4,4))
#' plot(poly,type='l',xlab="X",ylab="Y")
#' points(xy.dat)
#'
#' #Extracting points that lie inside a polygon
#' points(get_data(poly,xy.dat),pch=19,col=2)
#' @export
get_data<-function(polygon, xy.dat){
  data.frame(xy.dat[.point.in.polygon(pol.x = polygon[,1],
                                   pol.y = polygon[,2],
                                   point.x = xy.dat$x,
                                   point.y = xy.dat$y)!=0,])
}

