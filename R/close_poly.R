#' Close an open polygon
#'
#' Function closes an open polygon. If a line is provided, then it adds the line to the open ends of the polygon, else it appends the first point of the polygon at its end to close it.
#' @param open.poly A data frame of 2 columns (x,y), representing the vertices of an open polygon. E. g. an output of \emph{split_poly} function.
#' @param close.line A line to be used to close a given polygon. Must be a 2 column (x,y) data frame. The first coordinates of a line must correspond to the first point of the polygon, while the last coordinates - to the last point of the polygon. If close_line is not provided (the default), then the first point of the open polygon is used to close it.
#' @return A data frame of 2 columns (x,y), containing coordinates of a closed polygon.
#' @author Liudas Daumantas
#' @examples #A data.frame of an open polygon
#' open.pol <- data.frame(x=c(0,0,1,1),y=c(0,1,1,0))
#' plot(open.pol,type="l")
#'
#' #Closing the open polygon with its first point
#' poly <- .close_poly(open.pol)
#' plot(poly,type="l")
#'
#' #Closing the open polygon with a line
#' line <- data.frame(x=seq(0,1,length.out=10),y=c(0,.2,.1,.3,.6,.5,.3,.4,.2,0))
#' poly <- .close_poly(open.pol,line)
#' plot(poly,type="l")
#' @noRd

.close_poly<-function(open.poly, close.line=NULL){
  if (is.null(close.line)){
    return(rbind(open.poly,open.poly[1,]))
  } else {
    data.frame(x= c(open.poly[,1],rev(close.line[,1])[-1]), y=c(open.poly[,2],rev(close.line[,2])[-1]))
  }
}
