#' Generate points on a perimeter
#'
#' This function generates regularly spaced points on a perimeter of a polygon. Their number is specified by arguments n.pts or dst.pts.
#' @param polygon A data frame of 2 columns (x,y) that contain coordinates of polygon vertices. Both, closed and open polygons, are accepted.
#' @param n.pts A number of regularly spaced points to be generated on a perimeter of a polygon. If n.pts is not specified, then it is calculated according to the argument dst.pts.
#' @param dst.pts A distance along a polygon perimeter between adjecent points to be generated on a perimeter of a polygon. If dst.pts is not specified, then it is calculated according to the argument n.pts.
#' @return A list of 2 elements:
#' \describe{
#'   \item{\code{per_pts}}{A data frame of 4 columns, providing the information about the generated points on a perimeter of a polygon. This data frame is used as an input in \code{\link{pair_pts}} function.}
#'   \itemize{
#'   \item \code{x} - X coordinates of generated points.
#'   \item \code{y} - Y coordinates of generated points.
#'   \item \code{ID} - An ID that reflects the relative location of a point along a perimeter of a polygon in relation to other generated points and polygon vertices.
#'   \item \code{segment.no} = A vector indicating the ID of a polygon segment on which a generated point is located. It helps to indentify points located on the same
#' polygon segment.
#'   }
#'   \item{\code{full.poly}}{ A data frame that contains coordinates of the provided polygon vertices and generated points. \code{coords[,"ID"]} can be used to
#' extract rows of generated points.This data frame is used as an input in \code{\link{curvial.split}} function.}
#' }
#' @note If both, n.pts and dst.pts, are specified, then points are generated according to n.pts.
#' @author Liudas Daumantas
#' @examples #Creating data frame of a polygon
#' poly<- data.frame(X=c(3.38,3.30,1.70,0.78,-0.06,-2.30,-2.94,-3.97,-1.61,
#' -0.39,0.68,1.28,1.60,3.38),
#'                   Y=c(-0.12,-0.31,-2.73,-3.22,-3.29,-2.19,-1.62,0.94,3.10,
#'                   3.00,2.91,2.49,2.20,-0.12))
#' plot(poly,type='o')
#' #Generating 10 points on a polygon perimeter
#' a<-perimeter_pts(poly,n.pts = 10)
#' #location of points
#' points(a[[1]][,-3],col=2,pch=19)
#' #ID of points
#' text(x=a[[1]][,1],y=a[[1]][,2],a[[1]][,3],-0.3)
#' @export


.perimeter_pts<-function (polygon,n.pts=NULL,dst.pts=NULL){
  if (all(polygon[1,]!=polygon[nrow(polygon),])){
    polygon<-close_poly(polygon)
  }
  x<-polygon[,1]
  y<-polygon[,2]
  perimetras<-0
  for( i in 1:c(length(x)-1)){
    perimetras<-perimetras+sqrt((x[i]-x[i+1])^2+(y[i]-y[i+1])^2)
  }
  if (!is.null(n.pts)){
    dst.pts<-perimetras/n.pts
  } else {
    if (!is.null(dst.pts)) {
    n.pts<-round(perimetras/dst.pts,0)
  } else {
    stop(print("Specify either n.pts or dst.pts"))
    }
  }
  passed<-0
  generated_x_pt<-x[1]
  generated_y_pt<-y[1]
  ID<-1
  i=2
  segment.no<-i
  x.poly<-x[1]
  y.poly<-y[1]
  lastx<-generated_x_pt
  lasty<-generated_y_pt
  while (length(generated_x_pt)<n.pts){
    m<-sqrt((lastx-x[i])^2+(lasty-y[i])^2)
    if (m+passed>=dst.pts){
      kof<-(dst.pts-passed)/m
      if (lastx > x[i]) {
        X<-lastx-(lastx-x[i])*kof
      } else {
        X<-lastx+(x[i]-lastx)*kof
      }
      if (lasty > y[i]){
        Y<-lasty-(lasty-y[i])*kof
      } else {
        Y<-lasty+(y[i]-lasty)*kof
      }
      generated_x_pt<-c(generated_x_pt,X)
      generated_y_pt<-c(generated_y_pt,Y)
      x.poly<-c(x.poly,X)
      y.poly<-c(y.poly,Y)
      ID<-c(ID,length(x.poly))
      lasty<-Y
      lastx<-X
      segment.no<-c(segment.no,i)
      passed<-0
    } else {
      passed<-passed+m
      x.poly<-c(x.poly,x[i])
      y.poly<-c(y.poly,y[i])
      lasty<-y[i]
      lastx<-x[i]
      i<-i+1
    }
  }

  coords<-data.frame(x=generated_x_pt,y=generated_y_pt,ID,segment.no)
  full.poly<-rbind(data.frame(x.poly,y.poly),data.frame(x.poly=x[i:length(x)],y.poly=y[i:length(x)]))
  return(list(per_pts=coords,full.poly))
}
