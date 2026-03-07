#' Taxon-level leave-out contributions for split-lines
#'
#' @description
#' For each unique taxon label in \code{obj$call.info$Call_ARGS$data}
#' (for example, species or genus), remove all of its occurrences from both
#' child polygons of every split, recompute split-line performance, and record
#' both the resulting value and its difference from the original.
#'
#' @details
#' Let \eqn{P} be the split-line performance stored in
#' \code{obj$split.stats$performance}. For a focal taxon \eqn{t}, we compute
#' \eqn{P^{-t}} by removing all occurrences of \eqn{t}. We report \eqn{\Delta}
#' so that \eqn{\Delta > 0} always indicates a \strong{positive contributor}
#' (that is, removal worsens performance):
#' \itemize{
#'   \item If \code{obj$call.info$Call_ARGS$maximize = FALSE} (lower is better,
#'   for example similarity), \eqn{\Delta = P^{-t} - P}.
#'   \item If \code{obj$call.info$Call_ARGS$maximize = TRUE} (higher is better,
#'   for example distance), \eqn{\Delta = P - P^{-t}}.
#' }
#' If a taxon is absent from a split's parent polygons, elimination is a no-op
#' and \eqn{\Delta = NA}.
#'
#' @param obj A \code{hespdiv} object.
#'
#' @return A list of class \code{taxon_effect_result} with:
#' \describe{
#'   \item{\code{elim.comp.vals}}{
#'   Performance after removing each taxon; dimension
#'   \code{[n_splits x n_taxa]}.
#'   }
#'   \item{\code{delta}}{
#'   Signed contribution \eqn{\Delta} as defined above; dimension
#'   \code{[n_splits x n_taxa]}.
#'   }
#'   \item{\code{n_per_pol}}{
#'   Counts of the focal taxon per split and polygon.
#'   }
#'   \item{\code{baseline}}{
#'   The original performance vector.
#'   }
#' }
#' @family functions for hespdiv results post-processing
#' @export
taxon_effect <- function(obj) {
  if (!inherits(obj, "hespdiv")) stop("'obj' must be of class 'hespdiv'")

  dat_in_obj <- obj$call.info$Call_ARGS$data
  xy_in_obj  <- obj$call.info$Call_ARGS$xy.dat
  baseline   <- obj$split.stats$performance
  maximize   <- isTRUE(obj$call.info$Call_ARGS$maximize)

  # choose slicer + extract taxon vector
  if (is.data.frame(dat_in_obj) || is.matrix(dat_in_obj)) {
    .slicer <- .slicer.table
    taxa_vec <- as.character(dat_in_obj[[1]])
  } else if (is.list(dat_in_obj)) {
    .slicer <- .slicer.list
    taxa_vec <- as.character(unlist(dat_in_obj, use.names = FALSE))
  } else {
    .slicer <- .slicer.vect
    taxa_vec <- as.character(dat_in_obj)
  }

  group <- factor(taxa_vec)
  taxa_levels <- levels(group)
  n_taxa   <- length(taxa_levels)
  n_splits <- length(obj$split.lines)

  elim.comp.vals <- matrix(NA_real_, nrow = n_splits, ncol = n_taxa,
                           dimnames = list(rownames(obj$split.stats), taxa_levels))

  N <- n_taxa * n_splits * 2L
  n_per_pol <- data.frame(
    split.id = integer(N),
    group    = character(N),
    pol.id   = integer(N),
    n        = integer(N),
    stringsAsFactors = FALSE
  )

  for (tax_id in seq_len(n_taxa)) {
    tx <- taxa_levels[tax_id]
    idx_tax <- which(group == tx)
    xy_tax  <- xy_in_obj[idx_tax, , drop = FALSE]

    for (split.id in seq_len(n_splits)) {
      pol_ids <- which(obj$poly.stats$root.id == obj$split.stats$plot.id[split.id])
      pol_id1 <- obj$poly.stats$plot.id[pol_ids[1]]
      pol_id2 <- obj$poly.stats$plot.id[pol_ids[2]]
      child_pol1 <- obj$polygons.xy[[as.character(pol_id1)]]
      child_pol2 <- obj$polygons.xy[[as.character(pol_id2)]]
      ids1_tx <- .get_ids(child_pol1, xy_tax)
      ids2_tx <- .get_ids(child_pol2, xy_tax)

      offset_tax   <- (tax_id - 1L) * (2L * n_splits)
      offset_split <- (split.id - 1L) * 2L
      row1 <- offset_tax + offset_split + 1L
      row2 <- offset_tax + offset_split + 2L
      n_per_pol[row1, ] <- list(split.id, tx, pol_id1, length(ids1_tx))
      n_per_pol[row2, ] <- list(split.id, tx, pol_id2, length(ids2_tx))

      ids1_all <- .get_ids(child_pol1, xy_in_obj)
      ids2_all <- .get_ids(child_pol2, xy_in_obj)

      ids1_nogr <- setdiff(ids1_all, idx_tax)
      ids2_nogr <- setdiff(ids2_all, idx_tax)

      if (length(ids1_tx) == 0L && length(ids2_tx) == 0L) {
        elim.comp.vals[split.id, tax_id] <- NA
      } else {
        elim.comp.vals[split.id, tax_id] <- obj$call.info$Call_ARGS$compare.f(
          obj$call.info$Call_ARGS$generalize.f(.slicer(dat_in_obj, ids1_nogr)),
          obj$call.info$Call_ARGS$generalize.f(.slicer(dat_in_obj, ids2_nogr))
        )
      }
    }
  }

  # Base Δ as elim - baseline, then flip sign if the metric is maximised
  delta <- sweep(elim.comp.vals, 1L, baseline, FUN = "-")
  if (maximize) delta <- -delta

  structure(list(
    elim.comp.vals = elim.comp.vals,
    delta          = delta,
    n_per_pol      = n_per_pol,
    baseline       = baseline
  ), class = "taxon_effect")
}
