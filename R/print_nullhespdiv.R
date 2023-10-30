#' Print the results of nullhespdiv object
#'
#' @method print nullhespdiv
#' @description  print method for nullhespdiv class object
#' @param x A nullhespdiv class object
#' @param ... other arguments
#' @return x
#' @author Liudas Daumantas
#' @export
print.nullhespdiv <- function(x, ...){
  print(x[[1]])
  invisible(x)
  }
