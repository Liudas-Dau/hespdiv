#' Print the results of hespdiv object
#'
#' ss
#' @method print hespdiv
#' @description  Function formats and prints the results
#' of R object of class "hespdiv". It prints rounded split.stats data frame.
#' @param x hespdiv object
#' @param ... other arguments
#' @return x
#' @author Liudas Daumantas
#' @export
print.hespdiv <- function(x, ...){
  cat("\n","Information about the splits:", "\n","\n")
  print(round(x$split.stats,2))
  if (any(names(x)=="n.m.rez")){
    cat("\n", "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  }
  invisible(x)
}
