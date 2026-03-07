#' Investigate Group Effects on Split-Line Performance
#'
#' Re-evaluates split-line performance within each level of a grouping factor and
#' tests how much each group influences the detected split-lines in a \code{hespdiv} object.
#'
#' @description
#' For every level of \code{group}, the function:
#' \enumerate{
#'   \item subsets \code{data} and \code{xy.dat} to the focal group's observations;
#'   \item recomputes split-line performance using \code{compare.f} and \code{generalize.f};
#'   \item runs two contribution assessments:
#'     \itemize{
#'       \item \emph{group-removal}: removes the group's observations and recomputes performance;
#'       \item \emph{group-permutation (locality-block shuffle)}: shuffles intact within-locality assemblages of the focal group
#'             among the set of localities where the group occurs within the split’s parent polygons (non-group observations fixed).
#'     }
#'   \item optionally re-plots the object for that group via \code{\link{plot_hespdiv}} (with a group-specific subtitle).
#' }
#'
#' @details
#' \strong{Within-group recomputation (agreement).}
#' For each split, the group's observations are partitioned by the two child polygons and
#' performance is recomputed as \code{compare.f(generalize.f(pol1), generalize.f(pol2))}.
#' If the group's points fall on only one side, \code{within$est} is set to \code{maxdif};
#' if absent from both sides, it is \code{NA}.
#'
#' \strong{Elimination test.}
#' All observations of the focal group are removed from both polygons and performance is recomputed on the remaining data.
#' If the group's points fall on only one side, \code{elim$est} is computed normally;
#' if absent from both sides, it is \code{NA} (if desired, you could identify these cases afterwards from \code{n_per_pol}, changing to baseline performance and zero delta).
#'
#' \strong{Permutation test (locality-block shuffle).}
#' Localities are defined by identical coordinate pairs among the focal group's occurrences.
#' For each permutation, whole localities (blocks) are reassigned via a one-to-one mapping within the split’s parent polygons.
#' For speed, polygon membership of unique locality coordinates is precomputed once per split and then locality memberships are shuffled in each permutation.
#' If the group's points fall on only one side, \code{perm$est} is set to NA;
#' if absent from both sides, it is also \code{NA}.
#'
#' \strong{Baseline and deltas.}
#' The baseline vector is \code{baseline <- obj$split.stats$performance}.
#' Let \code{maximize <- obj$call.info$Call_ARGS$maximize}.
#' Deltas are defined so that \emph{positive values always indicate a positive contributor} (removal/permutation worsens performance):
#' \itemize{
#'   \item if \code{maximize == FALSE} (lower is better; e.g. similarity): \code{delta = est - baseline};
#'   \item if \code{maximize == TRUE}  (higher is better; e.g. distance):  \code{delta = baseline - est}.
#' }
#'
#' @param obj   A \code{hespdiv} object.
#' @param group A factor (or coercible to factor) giving the group label for each observation (row) in
#'              \code{obj$call.info$Call_ARGS$xy.dat}.
#' @param perm.n Integer (\eqn{>0}). Number of permutations in the locality-block permutation test. Default \code{999}.
#' @param maxdif Numeric. The performance value representing a maximal between-polygons difference for the chosen metric.
#'               If \code{NULL} and the metric is one of \code{c("pielou","morisita","sorensen","horn.morisita")},
#'               \code{maxdif} is set automatically. Otherwise, it must be provided.
#' @param plot  Logical. If \code{TRUE} (default), stores a \code{\link{plot_hespdiv}} output per group with
#'              \code{performance = TRUE} and a group-specific subtitle.
#' @param ...   Additional arguments passed to \code{\link{plot_hespdiv}}. Note that \code{obj},
#'              \code{performance}, and \code{subtitle} are set internally.
#'
#' @return A list of class \code{"group_effect"} with elements:
#' \itemize{
#'   \item \code{within = list(est, delta)} where \code{est} and \code{delta} are \code{[n_splits x n_groups]} matrices;
#'   \item \code{elim   = list(est, delta)} as above (group removed);
#'   \item \code{perm   = list(est, delta)} where each is a nested list \code{[[split]][[group]]} of numeric vectors (length \code{perm.n}).  If permutations are uninformative (e.g., one-sided/absent), the corresponding entries are \code{NA};
#'   \item \code{baseline}: the original performance vector;
#'   \item \code{n_per_pol}: \code{data.frame} with columns \code{split.id}, \code{group}, \code{pol.id}, \code{n};
#'   \item \code{plots}: list of \code{plot_hespdiv} outputs by group (or \code{NULL} if \code{plot = FALSE}).
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item Localities are defined by exact duplicate coordinate pairs among the group's observations
#'         (harmonise coordinates upstream if needed).
#'   \item The locality-block permutation preserves the \emph{number of group-present localities per polygon};
#'         only identities of the locality blocks are shuffled. The number of group occurrences per side can vary if block sizes differ.
#' }
#'
#' @family functions for hespdiv post-processing
#' @seealso \code{\link{plot_hespdiv}}
#' @importFrom stats setNames
#' @examples
#' \dontrun{
#' library(HDData)
#' group_effect(obj = hd, group = mio_mams$family, plot = FALSE)
#' }
#' @export
group_effect <- function(obj, group, perm.n = 999, maxdif = NULL, plot = TRUE, ...) {

  if (!inherits(obj, "hespdiv")) stop("'obj' must be of class 'hespdiv'")
  if (length(group) != nrow(obj$call.info$Call_ARGS$xy.dat)) {
    stop("'group' length must match the number of rows in 'xy.dat'")
  }
  if (length(perm.n) != 1L || !is.finite(perm.n) || perm.n < 1 ||
      perm.n %% 1 != 0) {
    stop("'perm.n' must be a positive integer")
  }
  perm.n <- as.integer(perm.n)

  # metric settings
  if (is.null(maxdif)) {
    supported_metrics <- c("pielou", "morisita", "sorensen", "horn.morisita")
    if (obj$call.info$METHOD$metric %in% supported_metrics) {
      if (obj$call.info$METHOD$metric == "pielou") {maxdif <- 1} else {maxdif <- 0}
    } else stop("Provide 'maxdif' value when using a custom method")
  }

  baseline <- obj$split.stats$performance
  maximize <- isTRUE(obj$call.info$Call_ARGS$maximize)

  # factorise group
  group <- as.factor(group)
  group_levels <- levels(group)
  group.n <- length(group_levels)
  n_splits <- length(obj$split.lines)

  # containers
  within_est <- matrix(NA_real_, nrow = n_splits, ncol = group.n,
                       dimnames = list(rownames(obj$split.stats), group_levels))
  elim_est   <- matrix(NA_real_, nrow = n_splits, ncol = group.n,
                       dimnames = list(rownames(obj$split.stats), group_levels))

  perm_est <- create_nested_list(list(rownames(obj$split.stats), group_levels))  # [[split]][[group]] -> numeric vector

  N <- group.n * n_splits * 2
  n_per_pol <- data.frame(
    split.id = integer(N),
    group    = character(N),
    pol.id   = integer(N),
    n        = integer(N),
    stringsAsFactors = FALSE
  )

  plots <- vector("list", length = group.n)
  names(plots) <- group_levels

  # slicer & data
  dat_in_obj <- obj$call.info$Call_ARGS$data
  xy_in_obj  <- obj$call.info$Call_ARGS$xy.dat
  if (is.data.frame(dat_in_obj) || is.matrix(dat_in_obj)) {
    .slicer <- .slicer.table
  } else if (is.list(dat_in_obj)) {
    .slicer <- .slicer.list
  } else {
    .slicer <- .slicer.vect
  }

  for (group_id in seq_len(group.n)) {
    gr <- group_levels[group_id]
    idx_group <- which(group == gr)

    data_sub <- .slicer(dat_in_obj, idx_group)
    xy_sub   <- xy_in_obj[idx_group, , drop = FALSE]

    for (split.id in seq_len(n_splits)) {

      pol_ids <- which(obj$poly.stats$root.id == obj$split.stats$plot.id[split.id])
      pol_id1 <- obj$poly.stats$plot.id[pol_ids[1]]
      pol_id2 <- obj$poly.stats$plot.id[pol_ids[2]]
      child_pol1 <- obj$polygons.xy[[as.character(pol_id1)]]
      child_pol2 <- obj$polygons.xy[[as.character(pol_id2)]]
      base <- baseline[split.id]

      # group points by child polygon (relative to xy_sub)
      split.ids1 <- .get_ids(child_pol1, xy_sub)
      split.ids2 <- .get_ids(child_pol2, xy_sub)

      # record counts
      offset_group <- (group_id - 1L) * (2L * n_splits)
      offset_split <- (split.id - 1L) * 2L
      row1 <- offset_group + offset_split + 1L
      row2 <- offset_group + offset_split + 2L
      n_per_pol[row1, ] <- list(split.id,
                                gr,
                                pol_id1,
                                length(split.ids1))
      n_per_pol[row2, ] <- list(split.id,
                                gr,
                                pol_id2,
                                length(split.ids2))

      # if absent from both sides
      if (n_per_pol[row1,"n"] == 0L && n_per_pol[row2,"n"] == 0L) {
        within_est[split.id, group_id] <- NA # cannot be defined
        elim_est  [split.id, group_id] <- NA # by definition will not change
        perm_est [[ split.id ]][[group_id]] <- NA # by definition will not change
      } else {
        # full-data indices per child polygon
        split.ids1_all <- .get_ids(child_pol1, xy_in_obj)
        split.ids2_all <- .get_ids(child_pol2, xy_in_obj)

        split.ids1_all_gr   <- split.ids1_all[split.ids1_all %in% idx_group]
        split.ids2_all_gr   <- split.ids2_all[split.ids2_all %in% idx_group]
        split.ids1_all_nogr <- split.ids1_all[!(split.ids1_all %in% idx_group)]
        split.ids2_all_nogr <- split.ids2_all[!(split.ids2_all %in% idx_group)]
        all_gr_ids_in_parent <- c(split.ids1_all_gr, split.ids2_all_gr)

        # elimination (remove group)
        elim_est[split.id, group_id] <- obj$call.info$Call_ARGS$compare.f(
          obj$call.info$Call_ARGS$generalize.f(.slicer(dat_in_obj, split.ids1_all_nogr)),
          obj$call.info$Call_ARGS$generalize.f(.slicer(dat_in_obj, split.ids2_all_nogr))
        )

        # --- Permutation: early skip if one side is empty ------------------------------

        # Check whether the focal group has occurrences on *both* sides of the split-line.
        # If one side has zero group observations, there is nothing meaningful to permute.
        has_both_sides <- (n_per_pol[row1, "n"] > 0L && n_per_pol[row2, "n"] > 0L)

        if (!has_both_sides) {
          # If the group occurs on only one side (or neither),
          # permutations are uninformative → store NA for this split/group combination.
          perm_est[[split.id]][[group_id]] <- NA

        } else {
          ## boundary-aware locality-block shuffle
          # Proceed with permutation test only if the group occupies both sides.

          # Extract coordinates of *all* occurrences of the focal group
          # that lie within the parent polygon of the current split.
          coords_agp <- xy_in_obj[all_gr_ids_in_parent, , drop = FALSE]

          # Identify unique locality coordinates
          uni_loc <- unique(coords_agp)

          # Map each observation to its unique locality index
          # These indeces allow to link locality-level info with occurrences:
          # locality_level_var[obs_to_loc] -> get locality info for each
          # all_gr_ids_in_parent
          obs_to_loc <- match(
            paste(coords_agp[[1]], coords_agp[[2]]),
            paste(uni_loc[[1]],   uni_loc[[2]])
          )
          # --- Determine polygon membership for each locality --------------------------
          # Identify which unique locality coordinates fall within the first and second
          # child polygons of this split. Localities can belong to both polygons (boundary cases).
          uni_loc_ids_in_1 <- .get_ids(child_pol1, uni_loc)
          uni_loc_ids_in_2 <- .get_ids(child_pol2, uni_loc)

          # convenience var:
          nloc <- nrow(uni_loc)

          # Create logical vectors marking polygon membership per locality.
          logi_loc_in_pol_1 <- logical(nloc)
          logi_loc_in_pol_1[uni_loc_ids_in_1] <- TRUE

          logi_loc_in_pol_2 <- logical(nloc)
          logi_loc_in_pol_2[uni_loc_ids_in_2] <- TRUE


          # Preallocate a numeric vector to store the recomputed performance
          # for each permutation iteration.
          res <- numeric(perm.n)

          # --- Main permutation loop ---------------------------------------------------
          # Each iteration permutes whole locality blocks among available positions,
          # preserving the number of blocks per polygon but randomising their identities.
          for (p in seq_len(perm.n)) {

            # Generate a random permutation of locality indices 1:nloc.
            perm_levels <- sample.int(nloc)

            # target locality index for each observation after permutation
            tgt_loc  <- perm_levels[obs_to_loc]

            # per-observation membership (after permutation)
            to_pol1 <- logi_loc_in_pol_1[tgt_loc]
            to_pol2 <- logi_loc_in_pol_2[tgt_loc]

            # observation ids for each polygon
            g_id1 <- all_gr_ids_in_parent[to_pol1]
            g_id2 <- all_gr_ids_in_parent[to_pol2]

            # Combine the non-group (fixed) observations with the permuted group observations
            # to form new assemblages for the two child polygons.
            dat1 <- .slicer(dat_in_obj, c(split.ids1_all_nogr, g_id1))
            dat2 <- .slicer(dat_in_obj, c(split.ids2_all_nogr, g_id2))

            # Recompute split-line performance for this permutation:
            #   - generalize.f() aggregates each polygon’s data
            #   - compare.f() computes between-polygons dissimilarity or distance.
            res[p] <- obj$call.info$Call_ARGS$compare.f(
              obj$call.info$Call_ARGS$generalize.f(dat1),
              obj$call.info$Call_ARGS$generalize.f(dat2)
            )
          }

          # Store the full vector of permutation results (length = nperm)
          # in the nested list corresponding to this split and group.
          perm_est[[split.id]][[group_id]] <- res
        }



        # within-group (agreement)
        if (n_per_pol[row1,"n"] == 0L || n_per_pol[row2,"n"] == 0L) {
          within_est[split.id, group_id] <- maxdif
        } else {
          dat_pol1 <- .slicer(data_sub, split.ids1)
          dat_pol2 <- .slicer(data_sub, split.ids2)
          within_est[split.id, group_id] <- obj$call.info$Call_ARGS$compare.f(
            obj$call.info$Call_ARGS$generalize.f(dat_pol1),
            obj$call.info$Call_ARGS$generalize.f(dat_pol2)
          )
        }
      }
    } # split loop

    if (plot) {
      obj_copy <- obj
      obj_copy$split.stats$performance <- within_est[, group_id]
      obj_copy$call.info$Call_ARGS$data <- data_sub
      obj_copy$call.info$Call_ARGS$xy.dat <- xy_sub
      plots[[group_id]] <- plot_hespdiv(obj_copy, performance = TRUE, subtitle = paste("Group:", gr), ...)
    } else {
      plots[[group_id]] <- NULL
    }
  } # group loop

  # ----- build deltas (consistent sign) -----
  make_delta <- function(est_mat, baseline, maximize) {
    d <- sweep(est_mat, 1L, baseline, FUN = "-")
    if (maximize) d <- -d
    d
  }
  within_delta <- -make_delta(within_est, baseline, maximize)
  elim_delta   <- make_delta(elim_est,   baseline, maximize)

  # perm deltas: transform each vector vs the split-specific baseline
  perm_delta <- create_nested_list(list(rownames(obj$split.stats), group_levels))
  for (split.id in seq_len(n_splits)) {
    base <- baseline[split.id]
    for (group_id in seq_len(group.n)) {
      v <- perm_est[[split.id]][[group_id]]
      if (is.null(v) || (length(v) == 1L && is.na(v))) {
        perm_delta[[split.id]][[group_id]] <- NA
      } else {
        d <- v - base
        if (maximize) d <- -d
        perm_delta[[split.id]][[group_id]] <- d
      }
    }
  }

  result <- list(
    within   = list(est = within_est, delta = within_delta),
    elim     = list(est = elim_est,   delta = elim_delta),
    perm     = list(est = perm_est,   delta = perm_delta),
    baseline = baseline,
    n_per_pol = n_per_pol,
    plots     = plots
  )
  class(result) <- "group_effect"
  return(result)
}

#' @noRd
create_nested_list <- function(levels) {
  if (length(levels) == 1) {
    return(setNames(replicate(length(levels[[1]]), list(), simplify = FALSE), levels[[1]]))
  }
  setNames(lapply(levels[[1]], function(x) create_nested_list(levels[-1])), levels[[1]])
}
