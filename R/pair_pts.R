#' Pairing connectable polygon points
#'
#' This function pairs points that are located on a perimeter of a polygon and can be joined by a straight line that is strictly contained within the polygon.
#' @param per_pts a data frame of 4 columns, first item of list output of \code{perimeter_pts} function.
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
#' a <- .perimeter_pts(poly,n.pts = 15)
#' #location of points
#' points(a[[1]][,-3],col=2,pch=19)
#' #ID of points
#' text(x=a[[1]][,1],y=a[[1]][,2],a[[1]][,3],-0.3)
#' #pairing points
#' b <- .pair_pts(a[[1]],poly)
#' #all arrows point to the right and all lines are within the polygon.
#' arrows(x0=b[,1],y0=b[,2],x1=b[,3],y1=b[,4],length = 0.15,col=2)
#' @noRd
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
        bad <- FALSE
        b1 <- (y[1]-y[2])/(x[1]-x[2])
        a1 <- y[1] - b1*x[1]
        mid.x <- mean(x)
        mid.y <- a1 + mid.x * b1
        if (.point.in.polygon(point.x = mid.x + 7e-14, point.y = mid.y+7e-14,
                                 pol.x = polygon[,1], pol.y = polygon[,2]
        ) == 1 &
        .point.in.polygon(point.x = mid.x - 7e-14, point.y = mid.y-7e-14,
                                      pol.x = polygon[,1], pol.y = polygon[,2]
        ) == 1 &
        .point.in.polygon(point.x = mid.x - 7e-14, point.y = mid.y+7e-14,
                                      pol.x = polygon[,1], pol.y = polygon[,2]
        ) == 1 &
        .point.in.polygon(point.x = mid.x + 7e-14, point.y = mid.y-7e-14,
                             pol.x = polygon[,1], pol.y = polygon[,2]
        ) == 1){
          for ( seg in 1:(nrow(polygon)-1) ){
            pol_seg <- polygon[seg:(seg+1),]
            if ( all(pol_seg[1,] == c(x[1],y[1])) |
                 all(pol_seg[1,] == c(x[2],y[2])) |
                 all(pol_seg[2,] == c(x[1],y[1])) |
                 all(pol_seg[2,] == c(x[2],y[2]))){
              next
            }
            X1 <- .in.range(range(x),pol_seg$x[1]) | .in.range(range(x),pol_seg$x[2])
            X2 <- .in.range(range(pol_seg$x),x[1]) | .in.range(range(pol_seg$x),x[2])
            Y1 <- .in.range(range(y),pol_seg$y[1]) | .in.range(range(y),pol_seg$y[2])
            Y2 <- .in.range(range(pol_seg$y),y[1]) | .in.range(range(pol_seg$y),y[2])

            if ( (X1 & Y1) | (X2 & Y2) | (X1 & Y2) | (X2 & Y1) ){
              b2 <- (pol_seg$y[1]-pol_seg$y[2])/(pol_seg$x[1]-pol_seg$x[2])
              a2 <- pol_seg$y[1] - b2*pol_seg$x[1]
              x.inters <- (a2 - a1) / (b1 - b2)
              if ( is.nan(x.inters) | is.infinite(x.inters) ){
                next
              }
              bad <- .in.range(range(x), x.inters) &
                .in.range(range(pol_seg$x),x.inters)
              if (bad) {
                break
              }
            }
          }
          if (!bad){
          id.min<-which.min(x)
          pairs_pts <- rbind(pairs_pts, data.frame(
            x1=x[id.min],
            y1=y[id.min],
            x2=x[-id.min],
            y2=y[-id.min],
            b=b1,
            id1=per_pts[(c(a,o)[id.min]),3],
            id2=per_pts[(c(a,o))[-id.min],3])
          )
          }
        }
      }} else{
        break
      }}
  pairs_pts
}
#' is x in range of x?
#' @noRd
.in.range <- function(range,x){
  range[1] < x & x < range[2]
}
