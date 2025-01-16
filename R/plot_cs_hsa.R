#' Visualize Constrained HespDiv Sensitivity Analysis Results
#'
#' @description
#' Displays the alternative (subsampled) \code{hespdiv} subdivisions and the
#' basal (original) \code{hespdiv} subdivision on one or multiple plots,
#' illustrating how the split-lines vary across different ranks.
#'
#' @param obj An object of class \code{hsa_constrained}, typically the output
#'   of \code{\link{hsa_sample_constrained}}.
#' @param type An integer indicating the type of plot. Defaults to 1.
#'   \describe{
#'     \item{\code{type = 1}}{A single plot overlaying all alternative lines
#'       and the basal lines, colored or width-coded by rank.}
#'     \item{\code{type = 2}}{Multiple plots, each corresponding to a specific
#'       rank. Alternative lines are displayed in a user-defined color
#'       (default \code{"lightyellow3"}), with the basal line highlighted
#'       in another color.}
#'   }
#' @param col_basal Character or numeric specifying the color of basal split-lines
#'   (default \code{"gray20"}).
#' @param col_boundary Character or numeric specifying the color of the outer
#'   (first) polygon boundary (default \code{7}).
#' @param col_alternatives Character or numeric specifying the color of alternative
#'   split-lines (default \code{"lightyellow3"}).
#' @param max_lwd Numeric. The maximum line width for the highest-ranked split-line
#'   (default \code{2.5}).
#' @param min_lwd Numeric. The minimum line width for the lowest-ranked split-line
#'   (default \code{0.75}).
#' @param alpha_alt Numeric in the range \code{[0, 1]}. The transparency of
#'   alternative split-lines (default \code{0.6}).
#' @return
#'   \code{NULL}. The function is called for its side effect of generating
#'   one or more plots.
#'
#' @details
#' \itemize{
#'   \item In \code{type = 1}, the function creates a single plot showing all
#'         alternative split-lines overlaid on the first polygon boundary,
#'         plus all basal split-lines of the \code{hespdiv} basis.
#'   \item In \code{type = 2}, the function creates separate plots, each focusing
#'         on polygons of a specific rank, drawing alternative lines in the user-specified
#'         color (with transparency) and the basal line in another color or line width.
#' }
#'
#' @importFrom graphics plot lines
#' @importFrom grDevices adjustcolor
#' @family functions for hespdiv sensitivity analysis
#' @family HespDiv visualization options
#' @export
plot_cs_hsa <- function(obj,
                        type = 1,
                        col_basal = "gray20",
                        col_boundary = 7,
                        col_alternatives = "lightyellow3",
                        max_lwd = 2.5,
                        min_lwd = 0.75,
                        alpha_alt = 0.6) {
  # 1) Basic checks
  if (!inherits(obj, "hsa_constrained")) {
    stop("`obj` must be a 'hsa_constrained' object.")
  }
  if (!type %in% c(1, 2)) {
    warning("Unsupported 'type' provided; defaulting to type=1.")
    type <- 1
  }

  # 2) Retrieve the 'Basis' and the 'Alternatives'
  basis <- obj$Basis
  alt_list <- obj$Alternatives

  # ranks: vector of ranks for each polygon (matching names of alt_list)
  if (!("poly.stats" %in% names(basis))) {
    stop("The `basis` (obj$Basis) must have 'poly.stats'.")
  }
  if (!all(names(alt_list) %in% rownames(basis$poly.stats))) {
    stop("Mismatch between names(obj$Alternatives) and rownames(basis$poly.stats).")
  }

  ranks <- basis$poly.stats[names(alt_list), "rank"]
  unique_ranks <- unique(ranks)

  # pol.n: number of polygons in the analysis
  pol.n <- length(alt_list)

  # n.runs: assume each polygon has the same number of runs
  n.runs <- if (pol.n > 0) length(alt_list[[1]]) else 0

  # 3) Helper to adjust alternative color with transparency
  alt_col_transparent <- adjustcolor(col_alternatives, alpha.f = alpha_alt)

  # 4) Plot logic
  if (type == 1) {
    # ---- TYPE = 1 ----
    if (!is.null(basis$polygons.xy) && length(basis$polygons.xy) >= 1) {
      # Plot outer boundary polygon in col_boundary
      plot(basis$polygons.xy[[1]], type = 'l',
           col = col_boundary,
           main = "Constrained Sensitivity: Overlaid Subdivisions")
    } else {
      stop("basis$polygons.xy[[1]] is not available for plotting the boundary.")
    }

    # line widths for each rank, from max_lwd to min_lwd
    lw_vec <- seq(max_lwd, min_lwd, length.out = length(unique_ranks))

    # A) Plot alternative lines
    for (i in seq_len(pol.n)) {
      this_rank <- ranks[i]
      # pick appropriate line width for this rank
      lw_poly <- lw_vec[this_rank]

      for (run_idx in seq_len(n.runs)) {
        alt_obj <- alt_list[[i]][[run_idx]]
        if (!is.null(alt_obj$split.lines) && length(alt_obj$split.lines) >= 1) {
          lines(alt_obj$split.lines[[1]],
                col = alt_col_transparent,
                lwd = lw_poly)
        }
      }
    }

    # B) Plot basal lines (all lines in basis$split.lines)
    if (!is.null(basis$split.lines)) {
      for (i in seq_along(basis$split.lines)) {
        lines(basis$split.lines[[i]], lwd = 0.5, col = col_basal)
      }
    }

  } else {
    # ---- TYPE = 2 ----
    # For each unique rank, create a separate plot
    for (rk in unique_ranks) {
      # which polygons have rank == rk?
      ids <- which(ranks == rk)
      main_txt <- paste("Rank =", rk)

      # plot the boundary in col_boundary
      if (!is.null(basis$polygons.xy) && length(basis$polygons.xy) >= 1) {
        plot(basis$polygons.xy[[1]], type = 'l',
             col = col_boundary,
             main = main_txt)
      } else {
        stop("basis$polygons.xy[[1]] is not available for plotting the boundary.")
      }

      # A) Plot alt lines for polygons of this rank
      for (poly_idx in ids) {
        for (run_idx in seq_len(n.runs)) {
          alt_obj <- alt_list[[poly_idx]][[run_idx]]
          if (!is.null(alt_obj$split.lines) && length(alt_obj$split.lines) >= 1) {
            lines(alt_obj$split.lines[[1]],
                  col = alt_col_transparent, lwd = 1.2)
          }
        }
        # highlight the basal line for this polygon
        if (!is.null(basis$split.lines) && length(basis$split.lines) >= poly_idx) {
          lines(basis$split.lines[[poly_idx]], lwd = 0.5, col = "red")
        }
      }

      # B) Plot the other basal lines with default line width/col_basal
      all_lines <- seq_along(basis$split.lines)
      others <- setdiff(all_lines, ids)
      for (ln_idx in others) {
        lines(basis$split.lines[[ln_idx]], lwd = 0.5, col = col_basal)
      }
    }
  }

  invisible(NULL)
}
