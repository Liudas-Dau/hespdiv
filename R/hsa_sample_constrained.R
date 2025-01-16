#' Constrained HespDiv Sensitivity Analysis by Subsampling
#'
#' @description
#' Conduct a constrained sensitivity analysis on a \code{hespdiv} object by
#' repeatedly subsampling observations within each polygon. Each subsample
#' is used to call \code{hespdiv} with recursion disabled (i.e., single-split only).
#'
#' @details
#' For each polygon in the \code{hespdiv} object, this function draws \code{subsample_factor}
#' of the data (by default 70%), creating multiple random subsamples (\code{n.runs}).
#' These are processed in chunks (as given by \code{chunks}), each chunk being
#' parallelized to manage memory usage. A higher \code{chunks} value reduces the
#' simultaneous load on memory but can slightly increase total computation time.
#'
#' @param obj A \code{hespdiv} object, containing the original data and polygons.
#' @param n.runs Integer. The number of subsampling runs to perform (default: 100).
#' @param subsample_factor Numeric proportion of data to subsample within each polygon (0 to 1].
#'   For example, 0.7 means 70% of the data in each polygon are retained.
#' @param RAM Integer. Approximate amount of RAM in GB to guide how many parallel
#'   workers to use. The function uses up to 80\% of your available CPU cores
#'   but also caps the number of workers at \code{RAM}.
#' @param chunks Integer. Controls how many chunks the \code{n.runs} are split into.
#'   Increasing this value reduces the number of runs processed simultaneously in
#'   parallel (helpful if you run out of memory). By default, \code{chunks=1}, meaning
#'   all runs are processed at once.
#'
#' @return A \code{hsa_constrained} class object, which is a list with two elements:
#' \itemize{
#'   \item \strong{Alternatives}: A named list corresponding to each polygon,
#'         where each entry is another list of \code{hespdiv} results for each subsample run.
#'   \item \strong{Basis}: The original \code{hespdiv} object (\code{obj}).
#' }
#'
#' @family functions for hespdiv sensitivity analysis
#' @family functions for hespdiv post-prossesing
#'
#' @seealso
#' \code{\link{hespdiv}} for details on the main function.
#' \code{\link{hsa}} for alternative hespdiv sensitivity analysis
#' \code{\link{future.apply}} for parallel processing methods.
#'
#' @importFrom future.apply future_lapply
#' @importFrom future plan multisession availableCores
#' @export
hsa_sample_constrained <- function(obj,
                                   n.runs = 100,
                                   subsample_factor = 0.7,
                                   RAM = 8,
                                   chunks = 1) {
  # ---- Input Validation ----
  if (!is.numeric(n.runs) || length(n.runs) != 1 || n.runs < 1) {
    stop("'n.runs' must be a positive integer.")
  }
  if (!is.numeric(subsample_factor) || subsample_factor <= 0 || subsample_factor > 1) {
    stop("'subsample_factor' must be in the range (0, 1].")
  }
  if (!inherits(obj, "hespdiv")) {
    stop("'obj' must be of class 'hespdiv'.")
  }

  # ---- Extract Core Information ----
  call_args <- obj$call.info$Call_ARGS
  dat_in_obj <- call_args$data

  # Decide on .slicer function
  if (is.data.frame(dat_in_obj) || is.matrix(dat_in_obj)) {
    .slicer <- .slicer.table
  } else if (is.list(dat_in_obj)) {
    .slicer <- .slicer.list
  } else {
    .slicer <- .slicer.vect
  }

  # ---- Prepare Output Structure ----
  split_ids <- obj$split.stats$plot.id
  hespdivs <- vector("list", length(split_ids))
  names(hespdivs) <- split_ids

  # ---- Parallel Setup ----
  num_cores <- future::availableCores()
  safe_cores <- max(1, floor(num_cores * 0.8))  # Use ~80% of available CPU
  safe_workers <- max(1, min(safe_cores, RAM)) # Also limit by declared RAM
  message("Using ", safe_workers, " parallel workers.")

  # Configure future
  future::plan(future::multisession, workers = safe_workers)
  options(future.globals.maxSize = 1024^3 * 2)  # 2 GB limit for serialized objects

  # ---- Main Loop: Iterate Over Polygons ----
  for (i in seq_along(split_ids)) {
    pol.id <- split_ids[i]

    # Subset data to this polygon
    polygon <- obj$polygons.xy[[as.character(pol.id)]]
    data_for_poly <- obj$poly.obj[[as.character(pol.id)]]
    xy.dat <- call_args$xy.dat[.get_ids(polygon, call_args$xy.dat), ]

    l <- nrow(xy.dat)

    # Generate random subsets (indices) for each run
    indices_list <- lapply(seq_len(n.runs), function(.) {
      sample.int(l, size = round(l * subsample_factor), replace = FALSE)
    })

    # Split the runs into chunks
    run_indices <- split(seq_len(n.runs), ceiling(seq_len(n.runs) / chunks))

    # Process each chunk in parallel
    result_list <- lapply(run_indices, function(runs_for_this_chunk) {
      future_lapply(runs_for_this_chunk, function(run) {
        # Subset data for this run
        l_xy_dat <- xy.dat[indices_list[[run]], , drop = FALSE]
        l_data   <- .slicer(data_for_poly, indices_list[[run]])

        # Call hespdiv with single-split recursion
        hespdiv(
          data         = l_data,
          xy.dat       = l_xy_dat,
          n.split.pts  = call_args$n.split.pts,
          same.n.split = call_args$same.n.split,  # usually TRUE
          method       = NULL,
          generalize.f = call_args$generalize.f,
          compare.f    = call_args$compare.f,
          maximize     = call_args$maximize,
          N.crit       = call_args$N.crit,
          N.rel.crit   = call_args$N.rel.crit,
          N.loc.crit   = call_args$N.loc.crit,
          N.loc.rel.crit = call_args$N.loc.rel.crit,
          S.crit       = call_args$S.crit,
          S.rel.crit   = call_args$S.rel.crit,
          Q.crit       = call_args$Q.crit,
          c.splits     = call_args$c.splits,
          c.Q.crit     = call_args$c.Q.crit,
          c.crit.improv= call_args$c.crit.improv,
          c.X.knots    = call_args$c.X.knots,
          c.Y.knots    = call_args$c.Y.knots,
          c.max.iter.no= call_args$c.max.iter.no,
          c.fast.optim = call_args$c.fast.optim,
          c.corr.term  = call_args$c.corr.term,
          study.pol    = polygon,
          use.chull    = FALSE,
          tracing      = NULL,
          pnts.col     = 1,
          display      = FALSE,
          pacific.region = call_args$pacific.region,
          .do_recurse  = FALSE
        )
      })
    })

    # Flatten the results from chunks into a single list for this polygon
    hespdivs[[i]] <- unlist(result_list, recursive = FALSE)
  }

  # ---- Return 'hsa_constrained' object ----
  structure(
    list(Alternatives = hespdivs, Basis = obj),
    class = "hsa_constrained"
  )
}
