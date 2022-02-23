#' Visualise split lines.
#'
#' @description Helper function that adds the selected graphical information to
#' the existing plot.
#' @param what String that informs what graphical information to show. Can be
#' either "curves", "straight" or "both".
#' @param level Sting indicating the trace level. Can be one of the following:
#' "best", "main", "all".
#' @param when String that states the situation, where the graphical information
#' is needed.
#' @author Liudas Daumantas
#' @importFrom pracma poly_center
#' @importFrom dplyr case_when
#' @noRd
.visualise_splits <- function(what,level, when) {
  if (what != "0"){
    dplyr::case_when(
      when == "start" & what != "curves" ~ {
        if(!is.null(pnts.col)){
          plot(data$x, data$y, col=pnts.col,xlab = "x coordinate",
               ylab = "y coordinate" )
        } else {
          plot(NULL,xlim = range(data$x),ylim = range(data$y),col=0,
               xlab = "x coordinate", ylab = "y coordinate")
        }
        centras <- pracma::poly_center(margins[,1],margins[,2])
        points(centras[1],centras[2],col=2,pch=19,cex=0.65)
        points(centras[1],centras[2],col=2,cex=1.5)
        lines(rims[[1]])
        print(paste0("ID of the polygon that now will be tested is: ", testid))

        if (testid>1) {
          for (i in 2:c(testid)){
            lines(x=rims[[i]][,1],y=rims[[i]][,2],col=1,lwd=2)
          }}
        points(perim_pts[[1]],pch=19,col="black")
        0
      }
    )
  }
}
