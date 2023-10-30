#' Calculate polygon object cross-comparison matrix
#'
#' This function computes the cross-comparison matrix of hespdiv polygon objects.
#' This matrix can be used as a distance matrix, either in its original form or after
#' undergoing a transformation if necessary, for further cluster analysis.
#' The matrix provides a quantitative measure of the similarity or
#' dissimilarity between different hespdiv polygon objects, enabling the
#' exploration of spatial relationships and patterns among them.
#' @param obj A hespdiv class object.
#' @author Liudas Daumantas
#' @return A cross-comparison matrix of hespdiv polygon objects.
#' @details The \code{cross_comp} function uses the \code{compare.f} function
#' from the \code{'obj$call.info$Call_ARGS'} list to perform the
#' pairwise  comparisons of the hespdiv polygon objects from the
#' \code{'obj$poly.obj'} list, producing a cross-comparison matrix.
#' @note The polygon cross-comparison functionality is currently not available
#' for the "pielou" method. Additionally, the functionality is not supported for the
#' custom methods that like "pielou" rely on variables from other environments within the
#' 'compare.f' function.
#' @family {functions for hespdiv post-prossesing}
#' @export

cross_comp <- function(obj){
  if (!inherits(obj,"hespdiv"))
    stop("obj should have 'hespdiv' class.")
  if (obj$call.info$METHOD$metric == "pielou") { # should pass original data through generalize.f
      # and compare only non overlapping polygons.
      stop(paste0("Polygon cross-comparison is not yet available for pielou method."))
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

