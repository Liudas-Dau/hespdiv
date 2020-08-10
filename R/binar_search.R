#' Binary search
#'
#' This function helps to find two adjacent indeces of a numeric vector, sorted in increasing order, that surround a value of interest.
#' @param X A numeric vector, sorted in increasing order.
#' @param x3 A value of interest.
#' @param l The highest possible lower index that surrounds a value of interest from left (bottom). Default value 1.
#' @param h The lowest possible higher index that surrounds a value of interest from right (top).Default value is the length of X.
#' @return A filtered data frame.
#' @note This function drops points that lie strictly outside of a given polygon. If a point lies on a relative interior of an edge of the polygon (provided that polygon is open) or point is a vertex of the polygon, then it will be retained.
#' @author Liudas Daumantas
#' @examples binar_search(X=1:100,x3=55.5)
#' @export

binar_search<-function(X,x3,l=1,h=length(X)){
  if (is.unsorted(X)){
    return(print("A numeric vector X must be sorted in increasing order"))
  }
  if (X[l]>=x3 | X[h]<=x3){
    return(print("x3 is not between provided numbers"))
  }
  recurs<-function(X,l,h,x3){
    mid_id<-round(l+(h-l)/2,0)
    if (mid_id==l|mid_id==h){
      return(c(l,h))
    }
    if (X[mid_id]>=x3){
      h<-mid_id
    } else {l<-mid_id}
    return(recurs(X,l,h,x3))
  }
  return(recurs(X,l,h,x3))
}
