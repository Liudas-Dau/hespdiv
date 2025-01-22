#' Investigate group effect
#'
#' @description
#' Re-evaluates the performance of split-lines within a \code{hespdiv} object,
#' separately for each subset defined by a grouping factor.
#'
#' For each level of \code{group}, the function:
#' \enumerate{
#'   \item Subsets the data and coordinates to observations belonging to that group
#'   \item Recomputes the performance measure (using \code{compare.f} and \code{generalize.f})
#'         for each split-line in \code{obj}
#'   \item Updates \code{obj$split.stats$performance} for that group
#'   \item Calls \code{plot_hespdiv()}, storing the resulting plot (with a group-specific subtitle)
#' }
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \code{comp.vals}: A matrix of size \code{(number_of_splits) x (number_of_groups)},
#'   giving the recomputed performances for each split-line and each group.
#'   \item \code{plots}: A list of the \code{plot_hespdiv} outputs (one for each group).
#' }
#'
#' @param obj A \code{hespdiv} object containing all necessary components (data, XY coordinates, etc.).
#' @param group A factor (or coercible to factor) defining the group of each observation in \code{obj$call.info$Call_ARGS$data}.
#' @param maxdif A numeric value specifying the output of \code{compare.f()} applied to maximally different groups.
#' @param ...  Additional arguments passed to \code{\link{plot_hespdiv}}.
#'             (\code{obj}, \code{performance}, and \code{subtitle} are preset within this function.)
#'
#' @family functions for hespdiv post-prossesing
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Suppose we have a hespdiv object called hd_obj and a grouping factor g:
#'   result <- group_effect(hd_obj, g, n.loc = TRUE)
#'   # result$comp.vals is the matrix of performance measures
#'   # result$plots is a list of ggplot/plot objects for each group
#' }

group_effect <- function(obj,
                         group,
                         maxdif = NULL,
                         ...
) {
  # Convert to factor to ensure consistent indexing
  group <- as.factor(group)
  group_levels <- levels(group)
  group.n <- length(group_levels)

  # Number of split-lines
  l <- length(obj$split.lines)

  # Prepare a matrix to store performance for all groups
  comp.vals <- matrix(
    NA,
    nrow = l,
    ncol = group.n,
    dimnames = list(
      rownames(obj$split.stats),
      group_levels
    )
  )

  # Prepare a list to hold one plot per group
  plots <- vector("list", length = group.n)
  names(plots) <- group_levels
  # Identify the data in obj
  dat_in_obj <- obj$call.info$Call_ARGS$data
  xy_in_obj  <- obj$call.info$Call_ARGS$xy.dat

  if (obj$call.info$call.info$METHOD$metric %in% c("pielou", "morisita",
                                                   "sorensen", "horn.morisita")) {
    maxdif <- 0
  } else {
    if (is.null(maxdif)) stop("Provide 'maxdif' value when using a sutom method")
  }

  # Decide how to slice data, based on obj$call.info$Call_ARGS$data structure
  if (is.data.frame(dat_in_obj) || is.matrix(dat_in_obj)) {
    .slicer <- .slicer.table
  } else if (is.list(dat_in_obj)) {
    .slicer <- .slicer.list
  } else {
    .slicer <- .slicer.vect
  }

  # Loop over each group
  for (group_id in seq_len(group.n)) {
    # Subset indices for this group
    idx_group <- which(group == group_levels[group_id])

    # Subset the data and xy coordinates for this group
    data_sub <- .slicer(dat_in_obj, idx_group)
    xy_sub   <- xy_in_obj[idx_group,]

    # For each split-line, compute performance measure
    for (split.id in seq_len(l)) {
      pol_ids <- which(obj$poly.stats$root.id == obj$split.stats$plot.id[split.id])

      split.ids1 <- .get_ids(obj$polygons.xy[[pol_ids[1]]], xy_sub)
      split.ids2 <- .get_ids(obj$polygons.xy[[pol_ids[2]]], xy_sub)

      if ( (length(split.ids2) > 0 | length(split.ids1) > 0) &
           (length(split.ids2) == 0 | length(split.ids1) == 0)) {
        # performance is maximum if groups dot not overlap
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

    # Update obj$split.stats$performance for this group (overwrites the vector)
    # so plot_hespdiv can visualize only the current group's performance
    obj$split.stats$performance <- comp.vals[, group_id]
    obj$call.info$Call_ARGS$data <- data_sub
    obj$call.info$Call_ARGS$xy.dat <- xy_sub
    # Generate and store the plot for this group
    # The 'subtitle' differentiates each group's plot
    plots[[group_id]] <- plot_hespdiv(
      obj,
      performance = TRUE,
      subtitle = paste("Group:", group_levels[group_id]),
      ...
    )
  }

  # Optionally, if you want the final obj to contain *all* group performances,
  # you could do: obj$split.stats$performance <- comp.vals

  # Return the performance matrix and the list of plots
  return(list(
    comp.vals = comp.vals,
    plots     = plots
  ))
}
