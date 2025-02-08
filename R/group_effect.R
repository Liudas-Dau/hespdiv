#' Investigate Group Effects on Split-Line Performance
#'
#' @description
#' This function re-evaluates the performance of split-lines within a \code{hespdiv} object,
#' separately for each subset defined by a grouping factor. Additionally, it performs permutation tests
#' to assess the influence of group membership on split-line performance.
#'
#' For each level of \code{group}, the function:
#' \enumerate{
#'   \item Subsets the data and coordinates to observations belonging to that group.
#'   \item Recomputes the performance measure using \code{compare.f} and \code{generalize.f}
#'         for each split-line in \code{obj}.
#'   \item Performs permutation tests to evaluate the robustness of group effects.
#'   \item Updates \code{obj} for that group and calls \code{plot_hespdiv()}, storing the resulting plot (with a group-specific subtitle)
#'         if \code{plot = TRUE}.
#' }
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{comp.vals}: A matrix of size \code{(number_of_splits) x (number_of_groups)},
#'   giving the recomputed performances for each split-line and each group.
#'   \item \code{elim.comp.vals}: A matrix of performance values when group-specific data is eliminated.
#'   \item \code{per_perf}: A nested list containing permutation test results for each split-line and group.
#'   \item \code{n_per_pol}: A data frame with counts of observations in each polygon for each split and group.
#'   \item \code{plots}: A list of the \code{plot_hespdiv} outputs (one for each group), or \code{NULL} if \code{plot = FALSE}.
#' }
#'
#' @param obj A \code{hespdiv} object containing all necessary components (data, XY coordinates, split-lines, etc.).
#' @param group A factor (or coercible to factor) defining the group of each observation in \code{obj$call.info$Call_ARGS$data}.
#' @param perm.n A positive integer specifying the number of permutations to perform (default is 500).
#' @param maxdif A numeric value specifying the performance measure when groups are maximally different.
#'               (Required only if \code{obj$call.info$METHOD$metric} is not one of
#'               \code{c("pielou", "morisita", "sorensen", "horn.morisita")}. Defaults to 0 for supported metrics.)
#' @param plot Logical (default \code{TRUE}). Whether to generate and store plot outputs via \code{\link{plot_hespdiv}}.
#' @param ...  Additional arguments passed to \code{\link{plot_hespdiv}}.
#'             (\code{obj}, \code{performance}, and \code{subtitle} are preset within this function.)
#'
#' @family functions for hespdiv post-processing
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Suppose we have a hespdiv object called hd_obj and a grouping factor g:
#'   result <- group_effect(hd_obj, g, perm.n = 1000, plot = TRUE)
#'   # result$comp.vals is the matrix of performance measures
#'   # result$elim.comp.vals gives performance values with group data eliminated
#'   # result$per_perf contains permutation test results for robustness evaluation
#'   # result$n_per_pol gives the counts per polygon for each split, per group
#'   # result$plots contains visualizations if plot=TRUE
#' }

group_effect <- function(obj, group, perm.n = 500, maxdif = NULL, plot = TRUE, ...) {
  if (!inherits(obj, "hespdiv")) {
    stop("'obj' must be of class 'hespdiv'")
  }
  if (length(group) != nrow(obj$call.info$Call_ARGS$xy.dat)) {
    stop("'group' length must match the number of rows in 'xy.dat'")
  }
  if (!is.numeric(perm.n) || perm.n <= 0) {
    stop("'perm.n' must be a positive integer")
  }

  # Automatic maxdif = 0 for known metrics; else require user input
  if (is.null(maxdif)) {
    supported_metrics <- c("pielou", "morisita", "sorensen", "horn.morisita")
    if (obj$call.info$METHOD$metric %in% supported_metrics) {
      maxdif <- 0
    } else {
      stop("Provide 'maxdif' value when using a custom method")
    }
  }

  # Convert to factor to ensure consistent indexing
  group <- as.factor(group)
  group_levels <- levels(group)
  group.n <- length(group_levels)
  l <- length(obj$split.lines)  # Number of split-lines

  # Initialize matrices and data structures
  comp.vals <- matrix(NA, nrow = l, ncol = group.n, dimnames = list(rownames(obj$split.stats), group_levels))
  elim.comp.vals <- matrix(NA, nrow = l, ncol = group.n, dimnames = list(rownames(obj$split.stats), group_levels))

  create_nested_list <- function(levels) {
    if (length(levels) == 1) {
      return(setNames(replicate(length(levels[[1]]), list(), simplify = FALSE), levels[[1]]))
    }
    setNames(lapply(levels[[1]], function(x) create_nested_list(levels[-1])), levels[[1]])
  }
  per_perf <- create_nested_list(list(rownames(obj$split.stats), group_levels))

  N <- group.n * l * 2
  n_per_pol <- data.frame(split.id = numeric(N), group = character(N), pol.id = numeric(N), n = numeric(N), stringsAsFactors = FALSE)
  plots <- vector("list", length = group.n)
  names(plots) <- group_levels

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
    xy_sub   <- xy_in_obj[idx_group, ]

    for (split.id in seq_len(l)) {
      pol_ids <- which(obj$poly.stats$root.id == obj$split.stats$plot.id[split.id])
      split.ids1 <- .get_ids(obj$polygons.xy[[pol_ids[1]]], xy_sub)
      split.ids2 <- .get_ids(obj$polygons.xy[[pol_ids[2]]], xy_sub)

      offset_group <- (group_id - 1) * (2 * l)
      offset_split <- (split.id - 1) * 2
      row1 <- offset_group + offset_split + 1
      row2 <- offset_group + offset_split + 2

      n_per_pol[row1, ] <- data.frame(split.id, gr, pol_ids[1], length(split.ids1), stringsAsFactors = FALSE)
      n_per_pol[row2, ] <- data.frame(split.id, gr, pol_ids[2], length(split.ids2), stringsAsFactors = FALSE)

      if (n_per_pol[row1,"n"] == 0 & n_per_pol[row2, "n"] == 0) {
        comp.vals[split.id, group_id] <- NA
        elim.comp.vals[split.id, group_id] <- obj$split.stats$performance[split.id]
        per_perf[[ split.id ]][[group_id]] <- NA
      } else {
        split.ids1_all <- .get_ids(obj$polygons.xy[[pol_ids[1]]], xy_in_obj)
        split.ids2_all <- .get_ids(obj$polygons.xy[[pol_ids[2]]], xy_in_obj)

        split.ids1_all_gr <- split.ids1_all[split.ids1_all %in% idx_group]
        split.ids2_all_gr <- split.ids2_all[split.ids2_all %in% idx_group]

        split.ids1_all_nogr <- split.ids1_all[!(split.ids1_all %in% idx_group)]
        split.ids2_all_nogr <- split.ids2_all[!(split.ids2_all %in% idx_group)]

        a_gp <- c(split.ids1_all_gr, split.ids2_all_gr)

        elim.comp.vals[split.id,group_id] <- obj$call.info$Call_ARGS$compare.f(
          obj$call.info$Call_ARGS$generalize.f(.slicer(dat_in_obj, split.ids1_all_nogr)),
          obj$call.info$Call_ARGS$generalize.f(.slicer(dat_in_obj, split.ids2_all_nogr))
        )

        n <- length(split.ids1) + length(split.ids2)
        per_perf[[split.id ]][[group_id]] <- replicate(perm.n, {
          id1 <- sample(1:n, size = sample(0:n, 1), replace = FALSE)
          id2 <- (1:n)[-id1]

          g_id1 <- a_gp[id1]
          g_id2 <- a_gp[id2]

          split.ids1_tgs <- c(split.ids1_all_nogr, g_id1)
          split.ids2_tgs <- c(split.ids2_all_nogr, g_id2)

          dat_tgs_pol1 <- .slicer(dat_in_obj, split.ids1_tgs)
          dat_tgs_pol2 <- .slicer(dat_in_obj, split.ids2_tgs)

          obj$call.info$Call_ARGS$compare.f(
            obj$call.info$Call_ARGS$generalize.f(dat_tgs_pol1),
            obj$call.info$Call_ARGS$generalize.f(dat_tgs_pol2)
          )
        })

        if ((n_per_pol[row1, "n"] > 0 | n_per_pol[row2, "n"] > 0) &&
            (n_per_pol[row1, "n"] == 0 | n_per_pol[row2, "n"] == 0)) {
          comp.vals[split.id, group_id] <- maxdif
        } else {
          dat_pol1 <- .slicer(data_sub, split.ids1)
          dat_pol2 <- .slicer(data_sub, split.ids2)

          comp.vals[split.id, group_id] <- obj$call.info$Call_ARGS$compare.f(
            obj$call.info$Call_ARGS$generalize.f(dat_pol1),
            obj$call.info$Call_ARGS$generalize.f(dat_pol2)
          )
        }
      }
    }

    if (plot) {
      obj$split.stats$performance <- comp.vals[, group_id]
      obj$call.info$Call_ARGS$data <- data_sub
      obj$call.info$Call_ARGS$xy.dat <- xy_sub

      plots[[group_id]] <- plot_hespdiv(
        obj,
        performance = TRUE,
        subtitle = paste("Group:", gr),
        ...
      )
    } else {
      plots[[group_id]] <- NULL
    }
  }

  result <- list(
    comp.vals = comp.vals,
    elim.comp.vals = elim.comp.vals,
    per_perf = per_perf,
    n_per_pol = n_per_pol,
    plots = plots
  )
  class(result) <- "group_effect_result"
  return(result)
}


