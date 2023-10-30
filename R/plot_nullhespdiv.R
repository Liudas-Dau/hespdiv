#' Plot the results of nullhespdiv object
#'
#' @method plot nullhespdiv
#' @description  plot method for nullhespdiv class object
#' @param x A nullhespdiv class object
#' @param ... other arguments
#' @return NULL
#' @author Liudas Daumantas
#' @importFrom graphics boxplot points
#' @export
plot.nullhespdiv <- function(x, ...){
  graphics::boxplot(x[[2]], ylim = range(c(x[[1]]$performance,
                                    unlist(x[[2]]))), main = paste0("n = ",nrow(x[[2]])))
  graphics::points(1:ncol(x[[2]]), x[[1]]$performance, col=2, pch= 19, cex = 1.5)
}
