#' Remove redundant polygon points
#'
#' This function finds and removes redundant polygon points that lie on the same straight line between 2 other points
#' @param xy a data frame with polygon coordinates
#' @return a data frame with cleaned polygon from redundant points
#' @author Liudas Daumantas
#' @noRd
.clean_poly <- function(xy){
  id <- numeric()
  for (i in 1:(nrow(xy) - 2)){
    y3 <- .pt_on_line(x1 = xy[i, 1],
                      x2 = xy[i+2, 1],
                      y1 = xy[i, 2],
                      y2 = xy[i+2, 2],
                      x3 = xy[i+1, 1])
    if (xy[i+1, 2] == y3)
      id <- c(id, i+1)
  }
  if (length(id) > 0)
    return(xy[-id,])
  xy
}
