#' Remove split-lines
#'
#' @description  Function returns hespdiv object, without split-lines of
#' specified id.
#'
#' @return hespdiv object
#' @param obj hespdiv object
#' @param split.id vector of split-line ids
#' @param depend.splits logical. Remove split-lines that depend on specified
#' split-lines? If FALSE, only end-nodes of spatial dendrogram are removed.
#' @examples
#' ## DO NOT RUN:
#' ## Preparing fossil assemblage data:
#' ## loading hespdiv data package that contains completed hespdiv analysis
#' ## of fossil mammal occurrence data from US
#' # library(HDData)
#'
#' ## inspect hespdiv object
#' # plot_hespdiv(hd)
#'
#' ## identify weak splits
#' # weak_splits <- which(hd$split.stats$performance >= 0.3)
#' ## remove weak splits
#' # performance_filtered <- remove_splits(obj = hd, split.id = weak_splits)
#' # plot_hespdiv(performance_filtered)
#'
#' ## inspect split-line significance
#' # nl
#' # plot(HHDATA::nl)
#' ## identify non-significant split-lines
#' # nsig_splits <- which(nl[[1]]$quantile >= 0.05)
#' ## remove the non-significant split-line
#' # sig_filtered <- remove_splits(obj = hd, split.id = nsig_splits)
#' # plot_hespdiv(sig_filtered)
#'
#' ## In case non-significant or weak split-line has good dependent split-lines,
#' ## you may choose to retain them using depend.splits = FALSE.
#' ## Attempt to remove 2nd split-line, that has dependent split-line (ID 3)
#' # unchanged_hd <- remove_splits(obj = hd, split.id = 2, depend = FALSE)
#' # plot_hespdiv(unchanged_hd)
#'
#' ## Alternatively, if you wish to remove non-significant or weak split-lines,
#' ## even if they have good dependent split-lines, use depend.splits = TRUE.
#' ## Attempt to remove 2nd split-line, that has dependent split-line (ID 3)
#' # changed_hd <- remove_splits(obj = hd, split.id = 2, depend = TRUE)
#' # plot_hespdiv(changed_hd)
#' @export
remove_splits <- function(obj, split.id, depend.splits = TRUE) {
  all.splits <- split.id
  all.pols <- numeric()
  for (id in split.id) {
    d <- depend_splits(obj, id)
    #browser()
    if (!depend.splits & length(d$splits) != 0 & !any(d$splits %in% split.id)){
      all.splits <- all.splits[all.splits != id]
    } else {
      all.splits <- c(all.splits, d$splits)
      all.pols <- c(all.pols, d$polygons)
    }

  }
  if (length(all.splits) > 0){
    all.splits <- unique(all.splits)
    all.pols <- unique(all.pols)
    plot.id <- obj$split.stats[all.splits,"plot.id"]
    obj$poly.stats[plot.id, "has.split"] <- FALSE
    obj$poly.stats <- obj$poly.stats[-all.pols,]

    obj$split.stats <- obj$split.stats[-all.splits,]
    obj$split.lines <- obj$split.lines[-all.splits]
    obj$polygons.xy <- obj$polygons.xy[-all.pols]
    obj$poly.obj <- obj$poly.obj[-all.pols]
    obj$str.difs <- obj$str.difs[-all.pols]
  }
  obj
}
