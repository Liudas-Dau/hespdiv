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
#' \donttest{
#' if (requireNamespace("HDData")) {
#'
#'   # Inspect the hespdiv object
#'   print(plot_hespdiv(HDData::hd))
#'
#'   # Remove weak split-lines
#'   weak_splits <- which(HDData::hd$split.stats$performance >= 0.3)
#'   performance_filtered <- remove_splits(obj = HDData::hd, split.id = weak_splits)
#'   print(plot_hespdiv(performance_filtered))
#'
#'   # Remove non-significant split-lines
#'   plot(HDData::nl)
#'   nsig_splits <- which(HDData::nl[[1]]$quantile >= 0.05)
#'   sig_filtered <- remove_splits(obj = HDData::hd, split.id = nsig_splits)
#'   print(plot_hespdiv(sig_filtered))
#'
#'   # Remove only if a split-line has no dependent split-lines
#'   unchanged_hd <- remove_splits(obj = HDData::hd, split.id = 4, depend.splits = FALSE)
#'   print(plot_hespdiv(unchanged_hd))
#'
#'   # Remove the split-lines indicated as well as all other split-lines
#'   # that structurally depend on them (default behavior)
#'   changed_hd <- remove_splits(obj = HDData::hd, split.id = 4, depend.splits = TRUE)
#'   print(plot_hespdiv(changed_hd))
#' }
#' }
#' @export
remove_splits <- function(obj, split.id, depend.splits = TRUE) {
  all.splits <- split.id
  all.pols <- numeric()
  for (id in split.id) {
    d <- depend_splits(obj, id)
    #browser()
    if (!depend.splits & length(d$splits) != 0 & any(!d$splits %in% split.id)){
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
