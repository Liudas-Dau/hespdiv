#' Binary search
#'
#' This function helps to find two adjacent indeces of a numeric vector, sorted in increasing order, that surround a value of interest.
#' @param X A numeric vector, sorted in increasing order.
#' @param x3 A value of interest.
#' @param l The highest possible lower index that surrounds a value of interest from left (bottom). Default value 1.
#' @param h The lowest possible higher index that surrounds a value of interest from right (top).Default value is the length of X.
#' @return Two adjacent X vector indeces of values that surrounds a value of interest.
#' @note If the value of interest (\code{x3}) is equal to some value of X, than the index of this value will be returned as a higher index.
#' Function is not intended to be used in cases of duplicated values.
#' @author Liudas Daumantas
#' @examples # simple case
#' binar_search(X=1:100,x3=55.5)
#' #when x3 is equal to some value of X, index with x3 value in the returned vector will be the higher one.
#' binar_search(X=1:100,x3=50)
#' @export

binar_search<-function(X,x3,l=1,h=length(X)){
  if (is.unsorted(X[l:h])){
    return(print("A numeric vector X must be sorted in increasing order"))
  }
  if ( x3 < x[l]| x3 > x[h]){
    return(print("x3 is not between provided numbers"))
  }
  recurs<-function(X,l,h,x3){
    mid_id <- round(l+(h-l)/2,0)
    if (mid_id==l|mid_id==h){
      return(c(l,h))
    }
    if (X[mid_id]>=x3){
      h <- mid_id
    } else {
      l <- mid_id
      }
    return(recurs(X,l,h,x3))
  }
  return(recurs(X,l,h,x3))
}
