#' Pairing connectable polygon points
#'
#' This function pairs points that are located on a perimeter of a polygon and can be joined by a straight line that is strictly contained within the polygon.
#' @param per_pts a data frame of 4 columns, first item of list output of \code{\link{perimeter_pts}} function.
#' @param polygon A data frame containing coordinates of a polygon. The polygon can be either closed or open.
#' @return A data frame with 7 columns. Each row represents a connectable pair of points located on a perimeter of a polygon.
#' \itemize{
#'   \item \code{x1} - x coordinate of a point that has lower x coordinate.
#'   \item \code{y1} - y coordinate of a point that has lower x coordinate.
#'   \item \code{x2} - x coordinate of a point that has higher x coordinate.
#'   \item \code{y2} - y coordinate of a point that has higher x coordinate.
#'   \item \code{b} - slope of a straight line connecting the paired points.
#'   \item \code{id1} - relative location of a point that has lower x coordinate. This ID is derived from \code{per_pts} third column (\code{ID}).
#'   \item \code{id2} - relative location of a point that has higher x coordinate. This ID is derived from \code{per_pts} third column (\code{ID}).
#' }
#' @note If both points have the same x coordinate, a point with lower ID will be treated as the one with lower x coordinate.
#' @author Liudas Daumantas
#' @examples #Creating data.frame of a polygon
#' poly<-data.frame(X=c(3.38,3.30,1.70,0.78,-0.06,-2.30,-2.94,-3.97,-1.61,
#'                     -0.39,0.68,1.28,1.60,3.38),
#'                 Y=c(-0.12,-0.31,-2.73,-3.22,0,-2.19,-1.62,0.94,3.10,
#'                 3.00,2.91,2.49,2.20,-0.12))
#' plot(poly,type='o')
#' #Generating 10 points on a polygon perimeter
#' a<-perimeter_pts(poly,n.pts = 15)
#' #location of points
#' points(a[[1]][,-3],col=2,pch=19)
#' #ID of points
#' text(x=a[[1]][,1],y=a[[1]][,2],a[[1]][,3],-0.3)
#' #pairing points
#' b<-pair_pts(a[[1]],poly)
#' #all arrows point to the right and all lines are within the polygon.
#' arrows(x0=b[,1],y0=b[,2],x1=b[,3],y1=b[,4],length = 0.15,col=2)
#' @export

.pair_pts<-function(per_pts,polygon){
  pairs_pts<-data.frame(numeric())
  for (a in 1:(nrow(per_pts)-1)){
    k<-1
    while(per_pts[a+k,4]==per_pts[a,4] &(a+k)<nrow(per_pts)){
      k=k+1
    }
    if(c(a+k)<=c(nrow(per_pts))&per_pts[a+k,4]!=per_pts[a,4]){
      for (o in c(a+k):c(nrow(per_pts))){
        x<-per_pts[c(a,o),1]
        y<-per_pts[c(a,o),2]
        id.min<-which.min(x)
        b<-(y[-id.min]-y[id.min])/(x[-id.min]-x[id.min])
        line.x<-seq(x[id.min],x[-id.min],length.out = 100)
        line.y<-y[id.min]+b*(line.x-x[id.min])
        if(all(point.in.polygon(line.x,line.y,polygon[,1],polygon[,2])[-c(1,100)]==1)){
          pairs_pts<-rbind(pairs_pts,data.frame(x1=x[id.min],y1=y[id.min],x2=x[-id.min],y2=y[-id.min],
                                        b=b,id1=per_pts[(c(a,o)[id.min]),3],id2=per_pts[(c(a,o))[-id.min],3]))
        }
      }} else{
        break
      }}
  pairs_pts
}
