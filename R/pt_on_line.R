#' Find a point on a line
#'
#' This function helps to find a missing coordinate (x or y) of a point on a straight line that joins two provided points.
#' @param x1 An x coordinate of the start of a line.
#' @param x2 An x coordinate of the end of a line.
#' @param y1 An y coordinate of the start of a line.
#' @param y2 An y coordinate of the end of a line.
#' @param x3 An x coordinates, where missing y coordinates on a line must be found.
#' @param y3 An y coordinates, where missing x coordinates on a line must be found.
#' @return A missing x or y coordinates on a line that joins (x1,y1) and (x2,y2) points.
#' @author Liudas Daumantas
#' @examples
#' #A line between (3,0) & (1,-1) points
#' plot(c(3,1),c(0,-1),type='l')
#' #Finding points at x = 1, 1.5, 2, 2.5, 3
#' x <- seq(1,3,0.5)
#' y <- pt_on_line(x1 = 3,x2 = 1,y1 = 0,y2 = -1,x3 = x)
#' points(x,y,col=2,pch=19)
#' @export
pt_on_line<-function(x1,x2,y1,y2,x3=NULL,y3=NULL){
  if (is.null(y3) & is.null(x3)){
    return(print("Provide either x3 or y3"))
  } else {
    if (!is.null(y3) & !is.null(x3)){
      print("Both y3 & x3 is provided, only x3 will be used. Returning a missing y coordinate: ")
    }
    if (is.null(y3)){
      return(y1+(x3-x1)*(y2-y1)/(x2-x1))
    } else {
      return(x1+(y3-y1)/((y2-y1)/(x2-x1)))
    }
  }
}
