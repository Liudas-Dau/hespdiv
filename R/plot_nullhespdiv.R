#' Plot a nullhespdiv object
#'
#' @description Plot method for a \code{nullhespdiv} object.
#'
#' @method plot nullhespdiv
#' @param x A \code{nullhespdiv} object.
#' @param ... Additional arguments passed to \code{graphics::boxplot()}.
#'
#' @return No return value. Called for the side effect of plotting a
#'   \code{nullhespdiv} object.
#'
#' @author Liudas Daumantas
#' @family HespDiv visualization options
#' @importFrom graphics boxplot points
#' @examples
#' test_results <- nulltest(example_hespdiv)
#' plot(test_results)
#' @export
plot.nullhespdiv <- function(x, ...) {
  args <- list(...)

  if (is.null(args$ylim)) {
    args$ylim <- range(c(x[[1]]$performance, unlist(x[[2]])), na.rm = TRUE)
  }

  if (is.null(args$main)) {
    args$main <- paste0("n = ", nrow(x[[2]]))
  }

  do.call(
    graphics::boxplot,
    c(list(x = x[[2]]), args)
  )

  graphics::points(
    1:ncol(x[[2]]),
    x[[1]]$performance,
    col = 2,
    pch = 19,
    cex = 1.5
  )
}
