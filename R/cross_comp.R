#' Calculate a polygon-object cross-comparison matrix
#'
#' @description
#' Computes a cross-comparison matrix for polygon objects from a \code{hespdiv}
#' object. The matrix quantifies similarity or dissimilarity among polygon
#' objects and can be used for further analyses, such as clustering, either
#' directly or after transformation.
#'
#' @param obj A \code{hespdiv} object.
#'
#' @return
#' A numeric matrix containing pairwise comparison values among the
#' \code{hespdiv} polygon objects stored in \code{obj$poly.obj}.
#'
#' @details
#' The \code{cross_comp()} function uses the \code{compare.f} function from
#' \code{obj$call.info$Call_ARGS} to perform pairwise comparisons of
#' \code{hespdiv} polygon objects stored in \code{obj$poly.obj}. The result is
#' a cross-comparison matrix.
#'
#' @note
#' Polygon cross-comparison is currently not available for the \code{"pielou"}
#' method. It is also not supported for custom methods whose \code{compare.f}
#' function relies on variables from environments other than the function's own
#' arguments.
#'
#' @family functions for hespdiv results post-processing
#' @author Liudas Daumantas
#' @export

cross_comp <- function(obj){
  if (!inherits(obj, "hespdiv")) {
    stop("`obj` must be a `hespdiv` object.", call. = FALSE)
  }

  if (obj$call.info$METHOD$metric == "pielou") {
    # should pass original data through generalize.f
    # and compare only non-overlapping polygons
    stop(
      "Polygon cross-comparison is not yet available for the Pielou method.",
      call. = FALSE
    )
  } else {
    compare.f <- obj$call.info$Call_ARGS$compare.f
  }
  comp.mat <- matrix(NA,nrow = length(obj$poly.obj),ncol =length(obj$poly.obj) )
  for (pol.id in seq(length(obj$poly.obj))){
    comp.mat[,pol.id] <- unlist(lapply(obj$poly.obj,FUN = compare.f,
                                obj$poly.obj[[pol.id]]))
  }
  comp.mat
}

