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
#' These are processed in chunks (as given by \code{chunk_size}) and runs in each chunk being
#' parallelized to manage memory usage.
#'
#' @param obj A \code{hespdiv} object.
#' @param n.runs Integer. The number of subsampling runs to perform (default: 100).
#' @param subsample_factor Numeric proportion of data to subsample within each polygon (0 to 1].
#'   For example, 0.7 means 70% of the data in each polygon are retained.
#' @param RAM Integer. Approximate amount of RAM in GB to guide how many parallel
#'   workers to use. The function uses up to 80\% of your available CPU cores
#'   but also caps the number of workers at \code{RAM}.
#' @param chunk_size Integer. Determines the number of hespdiv runs to be processed in parallel.
#' @param load_prop Numeric value (0,1]. Specifies the proportion of available
#'   CPU cores or RAM to be used for setting up parallel workers. For example,
#'   \code{load_prop = 0.8} uses 80% of the available resources. If both
#'   \code{load_prop} and \code{RAM} are provided, the number of workers will
#'   be the minimum based on the constraints imposed by both. Defaults to
#'   \code{0.8} if not provided.
#' @param workers a number of parrallel workers.
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
#' \CRANpkg{future.apply}
#'
#' @importFrom future.apply future_lapply
#' @importFrom future plan multisession availableCores
#' @importFrom pracma polyarea
#' @export
hsa_sample_constrained <- function(obj,
                                   n.runs = 100,
                                   subsample_factor = 0.7,
                                   RAM = NULL,
                                   load_prop = NULL,
                                   chunk_size = 8,
                                   workers = NULL) {
  # ---- Extract and Validate Basic Args ----
  if (!is.numeric(subsample_factor) || subsample_factor <= 0 || subsample_factor > 1) {
    stop("'subsample_factor' must be in the range (0, 1].")
  }
  if (!inherits(obj, "hespdiv")) {
    stop("'obj' must be of class 'hespdiv'.")
  }
  if (!is.numeric(n.runs) || n.runs < 1 || n.runs %% 1 || length(n.runs) != 1) stop("'n.runs' must be a positive integer.")
  if (!is.numeric(chunk_size) || chunk_size < 1 || chunk_size %% 1 || length(n.runs) != 1) stop("'chunks' must be a positive integer.")

  # Extract call_args AFTER we confirm obj is a hespdiv
  call_args <- obj$call.info$Call_ARGS
  if (!call_args$same.n.split) {
    stop("function currently works only when 'same.n.split' is TRUE.")
  }

  # ---- Extract Core Information ----
  dat_in_obj <- call_args$data

  # The bounding polygon's area
  big_poly_area <- pracma::polyarea(x = obj$polygons.xy[[1]][,1],
                                    y = obj$polygons.xy[[1]][,2])
  S_crit_abs <- call_args$S.crit * abs(big_poly_area)

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
  hespdivs  <- vector("list", length(split_ids))
  names(hespdivs) <- split_ids

  # ---- Parallel Setup ----
  if (is.null(workers)) {
    num_cores <- future::availableCores()

    # if load_prop is NULL, default to 0.8 (or 1, or any sensible default)
    if (is.null(load_prop)) {
      load_prop <- 0.8
    }
    if (!is.numeric(load_prop) || load_prop <= 0 || load_prop > 1) {
      stop("'load_prop' must be a numeric in (0,1].")
    }
    safe_cores <- max(1, floor(num_cores * load_prop))

    # if RAM is NULL, treat it as infinite
    if (is.null(RAM)) {
      RAM <- Inf
    } else {
      if (!is.numeric(RAM) || length(RAM) != 1 || RAM < 1) {
        stop("'RAM' must be a positive number if not NULL.")
      }
    }

    safe_workers <- max(1, min(safe_cores, RAM))
  } else {
    if (!is.numeric(workers) || workers < 1) {
      stop("'workers' must be a positive number of workers.")
    }
    safe_workers <- workers
  }
  message("Using ", safe_workers, " parallel workers.")

  # Configure future
  future::plan(future::multisession, workers = safe_workers)
  options(future.globals.maxSize = 1024^3 * 2)  # 2 GB limit for serialized objects

  # ---- Main Loop: Iterate Over Polygons ----
  for (i in seq_along(split_ids)) {
    pol.id <- split_ids[i]

    # Subset data to this polygon
    polygon <- obj$polygons.xy[[as.character(pol.id)]]
    ids <- .get_ids(polygon, call_args$xy.dat)
    xy.dat <- call_args$xy.dat[ids, ]
    data_for_poly <- .slicer(dat_in_obj, ids)

    l <- nrow(xy.dat)

    # Generate random subsets (indices) for each run
    indices_list <- lapply(seq_len(n.runs), function(.) {
      sample.int(l, size = round(l * subsample_factor), replace = FALSE)
    })

    # Split the runs into chunks of a given size
    run_indices <- split(seq_len(n.runs), ceiling(seq_len(n.runs) / chunk_size))

    # Process each chunk in parallel
    result_list <- lapply(run_indices, function(runs_for_this_chunk) {
      future_lapply(runs_for_this_chunk, function(run) {
        # Subset data for this run
        l_xy_dat <- xy.dat[indices_list[[run]], , drop = FALSE]
        l_data   <- .slicer(data_for_poly, indices_list[[run]])

        # Calculate S.crit for this polygon area
        poly_area <- abs(pracma::polyarea(x = polygon[,1],
                                          y = polygon[,2]))
        S_crit_local <- S_crit_abs / poly_area

        # Call hespdiv with single-split recursion
        hespdiv(
          data         = l_data,
          xy.dat       = l_xy_dat,
          n.split.pts  = call_args$n.split.pts,
          same.n.split = call_args$same.n.split,  # must be TRUE
          method       = NULL,
          generalize.f = call_args$generalize.f,
          compare.f    = call_args$compare.f,
          maximize     = call_args$maximize,
          N.crit       = call_args$N.crit,
          N.rel.crit   = call_args$N.rel.crit,
          N.loc.crit   = call_args$N.loc.crit,
          N.loc.rel.crit = call_args$N.loc.rel.crit,
          S.crit       = S_crit_local,  # local area-based
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

    # Flatten the results from chunks into a single list for this polygons
    hespdivs[[i]] <- unlist(result_list, recursive = FALSE)
  }

  # ---- Return 'hsa_constrained' object ----
  structure(
    list(Alternatives = hespdivs, Basis = obj),
    class = "hsa_constrained"
  )
}
