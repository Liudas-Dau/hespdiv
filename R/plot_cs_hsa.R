#' Visualize Constrained HespDiv Sensitivity Analysis Results
#'
#' @description
#' Displays the alternative (subsampled) \code{hespdiv} subdivisions and the
#' basal (original) \code{hespdiv} subdivision on one or multiple plots,
#' illustrating how the split-lines vary across different ranks. Additionally,
#' for each alternative split-line (which is defined by a start and end coordinate),
#' the function aggregates identical endpoints and overlays their counts on the plot.
#'
#' @param obj An object of class \code{hsa_constrained}, the output
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
#' @param rank Integer. Optional. When \code{type = 2}, if provided the function
#'   will only generate plots for the specified rank. If \code{NULL} (default),
#'   plots for all unique ranks will be produced.
#' @param col_basal Character or numeric specifying the color of basal split-lines
#'   (default \code{"gray20"}).
#' @param main Character. Title for the plot(s).
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
#'         If a specific \code{rank} is provided, only that rank is plotted.
#'   \item In both cases, after drawing the alternative split-lines the function
#'         aggregates their endpoints (start and end points) and overlays the count
#'         at each unique coordinate using \code{text()}.
#' }
#'
#' @importFrom graphics plot lines text
#' @importFrom grDevices adjustcolor
#' @family functions for hespdiv sensitivity analysis
#' @family HespDiv visualization options
#' @export
plot_cs_hsa <- function(obj,
                        type = 1,
                        rank = NULL,
                        col_basal = "gray20",
                        main,
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
    warning("Unsupported 'type' provided; defaulting to type = 1.")
    type <- 1
  }

  # 2) Retrieve the 'Basis' and the 'Alternatives'
  basis <- obj$Basis
  alt_list <- obj$Alternatives

  if (!("poly.stats" %in% names(basis))) {
    stop("The `basis` (obj$Basis) must have 'poly.stats'.")
  }
  if (!all(names(alt_list) %in% rownames(basis$poly.stats))) {
    stop("Mismatch between names(obj$Alternatives) and rownames(basis$poly.stats).")
  }

  ranks <- basis$poly.stats[names(alt_list), "rank"]
  unique_ranks <- unique(ranks)

  # If a specific rank is provided (for type==2) then use only that rank
  if (type == 2 && !is.null(rank)) {
    if (!rank %in% unique_ranks) {
      warning(paste("Provided rank", rank, "is not found. Nothing will be plotted."))
      return(invisible(NULL))
    }
    unique_ranks <- rank
  }

  pol.n <- length(alt_list)
  n.runs <- if (pol.n > 0) length(alt_list[[1]]) else 0

  alt_col_transparent <- adjustcolor(col_alternatives, alpha.f = alpha_alt)

  # Helper function to overlay counts for unique endpoints
  overlay_density_endpoints <- function(alt_coords) {
    if (length(alt_coords) == 0) return()
    # Each element in alt_coords is a matrix with two rows (start, end)
    # Combine all endpoints into one matrix (each row = one endpoint)
    all_points <- do.call(rbind, alt_coords)
    # Create a unique key for each point by concatenating its coordinates
    keys <- apply(all_points, 1, function(row) paste(row, collapse = ","))
    counts <- table(keys)
    unique_keys <- names(counts)
    for (key in unique_keys) {
      coords <- as.numeric(strsplit(key, split = ",")[[1]])
      text(coords[1], coords[2], labels = counts[[key]], col = "blue", cex = 0.8)
    }
  }

  # 3) Plot logic
  if (type == 1) {
    # ---- TYPE = 1 ----
    if (!is.null(basis$polygons.xy) && length(basis$polygons.xy) >= 1) {
      plot(basis$polygons.xy[[1]], type = 'l', col = col_boundary, main = main)
      points(obj$Basis$call.info$Call_ARGS$xy.dat,pch=19,cex = 0.25)
    } else {
      stop("basis$polygons.xy[[1]] is not available for plotting the boundary.")
    }

    # alt_coords will store the endpoints of each alternative split-line
    alt_coords <- list()

    lw_vec <- seq(max_lwd, min_lwd, length.out = length(unique_ranks))

    # A) Plot alternative lines and collect endpoints
    for (i in seq_len(pol.n)) {
      this_rank <- ranks[i]
      lw_poly <- if (length(lw_vec) > 1) lw_vec[this_rank] else lw_vec[1]
      for (run_idx in seq_len(n.runs)) {
        alt_obj <- alt_list[[i]][[run_idx]]
        if (!is.null(alt_obj$split.lines) && length(alt_obj$split.lines) >= 1) {
          line_coords <- alt_obj$split.lines[[1]]
          lines(line_coords, col = alt_col_transparent, lwd = lw_poly)
          # Record the start and end points of this split-line
          endpoints <- line_coords[c(1, nrow(line_coords)), , drop = FALSE]
          alt_coords[[length(alt_coords) + 1]] <- endpoints
        }
      }
    }

    # B) Plot basal lines
    if (!is.null(basis$split.lines)) {
      for (i in seq_along(basis$split.lines)) {
        lines(basis$split.lines[[i]], lwd = 0.5, col = col_basal)
      }
    }

    # C) Overlay the density (number of split-lines) at each unique endpoint
    overlay_density_endpoints(alt_coords)

  } else {
    # ---- TYPE = 2 ----
    for (rk in unique_ranks) {
      alt_coords <- list()
      ids <- which(ranks == rk)
      main_txt <- paste("Rank =", rk)

      if (!is.null(basis$polygons.xy) && length(basis$polygons.xy) >= 1) {
        plot(basis$polygons.xy[[1]], type = 'l', col = col_boundary, main = main_txt)
        points(obj$Basis$call.info$Call_ARGS$xy.dat,pch=19,cex = 0.25)
      } else {
        stop("basis$polygons.xy[[1]] is not available for plotting the boundary.")
      }

      # A) Plot alternative lines for polygons of this rank and collect endpoints
      for (poly_idx in ids) {
        for (run_idx in seq_len(n.runs)) {
          alt_obj <- alt_list[[poly_idx]][[run_idx]]
          if (!is.null(alt_obj$split.lines) && length(alt_obj$split.lines) >= 1) {
            line_coords <- alt_obj$split.lines[[1]]
            lines(line_coords, col = alt_col_transparent, lwd = 1.2)
            endpoints <- line_coords[c(1, nrow(line_coords)), , drop = FALSE]
            alt_coords[[length(alt_coords) + 1]] <- endpoints
          }
        }
        if (!is.null(basis$split.lines) && length(basis$split.lines) >= poly_idx) {
          lines(basis$split.lines[[poly_idx]], lwd = 0.5, col = "red")
        }
      }

      all_lines <- seq_along(basis$split.lines)
      others <- setdiff(all_lines, ids)
      for (ln_idx in others) {
        lines(basis$split.lines[[ln_idx]], lwd = 0.5, col = col_basal)
      }

      overlay_density_endpoints(alt_coords)
    }
  }

  invisible(NULL)
}
