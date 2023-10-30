#' Change the basal subdivision in hsa
#'
#' @description  This function allows you to select a new basal subdivision
#' from the results of hespdiv sensitivity analysis by specifying its ID. It
#' provides a convenient way to switch between different basal subdivisions.
#' You can identify a more stable subdivision alternative by plotting the
#' hespdiv sensitivity analysis results using the \code{plot_hsa} function.
#' By selecting a new basal subdivision, you can observe how it affects the
#' results of polygon object cross-comparison and the stability of hespdiv
#' clusters and polygons.
#' @param obj A \code{hsa} class object.
#' @param id An index of an alternative subdivision to be used as a new basal
#' subdivision.
#' @return The \code{hsa} class object with a new basal subdivision.
#' @author Liudas Daumantas
#' @family {functions for hespdiv sensitivity analysis}
#' @export
change_base <- function(obj, id){
  if (!inherits(obj,"hsa"))
    stop("'obj' should be of class \"hsa\" (output of 'hsa' or 'hsa_detailed' functions).")
  base <- obj[[1]][[id]][[1]]
  if (is.null(base)) stop("Selected base has no subdivisions")
  obj[[1]][[id]][[1]] <- obj[[2]]
  obj[[1]][[id]][[2]] <- obj[[2]]$call.info$Call_ARGS
  structure(list(Alternatives = obj[[1]], Basis = base),class = 'hsa')
}

