#' HespDiv Sensitivity Analysis
#'
#' This function performs sensitivity analysis of the \code{hespdiv} method.
#' Starting from a reference \code{hespdiv} object, it generates a specified
#' number of alternative \code{hespdiv} calls by randomly sampling new values
#' for selected arguments from user-provided sets.
#'
#' @param obj A \code{hespdiv} class object.
#'
#' @param n.runs Integer. The number of alternative \code{hespdiv} calls to evaluate.
#'
#' @param parallel Logical. If \code{TRUE}, runs are evaluated in parallel using
#'   \pkg{future.apply}. If \code{FALSE} (default), runs are evaluated sequentially.
#'
#' @param workers Integer. Number of parallel workers (CPU cores) to use when
#'   \code{parallel = TRUE}. At most this many \code{hespdiv} calls will be executed
#'   simultaneously. If \code{NULL}, the number of workers is determined
#'   automatically based on \code{RAM} and \code{load_prop}.
#'
#' @param RAM Integer. Approximate amount of available RAM (in GB) used to limit
#'   the number of parallel workers when \code{parallel = TRUE}.
#'
#' @param load_prop Numeric in (0, 1]. Proportion of available CPU cores to use
#'   when determining the number of parallel workers automatically. Defaults to 0.8.
#'
#' @param chunk_size Integer. Number of runs submitted per batch (chunk) to the
#'   parallel backend when \code{parallel = TRUE}. Chunking does not change how
#'   many runs execute simultaneously (controlled by \code{workers}). Smaller
#'   values reduce the number of queued futures at once (often lowering RAM
#'   overhead). Defaults to \code{workers*2}.
#'
#' @param future_seed Logical. Passed to \code{future.apply::future_lapply()}.
#'   If \code{TRUE}, RNG is managed by \pkg{future} to provide parallel-safe,
#'   reproducible random streams.
#'
#' @param data.paired Logical. Controls whether alternative values of \code{data}
#'   are paired with corresponding elements of \code{xy.dat}.
#'
#' @param display Logical. Controls the value of the \code{display} argument
#'   in each \code{hespdiv} call.
#'
#' @param images.path Character or \code{NULL}. Path to an existing directory
#'   where PNG images of displayed results will be saved. If \code{NULL}
#'   (default), images are not saved.
#'
#' @param pnts.col Value passed to the \code{pnts.col} argument of each
#'   \code{hespdiv} call. Can be provided as atomic vector,  vector or list of
#'   vectors.
#'
#' @param data A list of data objects (matrices, data frames, vectors, lists,
#'   or other supported data structures) used as alternative inputs for
#'   sensitivity analysis.
#'
#' @param xy.dat,study.pol Lists of data frames with two columns: \code{x} and \code{y}.
#'
#' @param same.n.split,c.fast.optim,use.chull,c.splits Logical vectors specifying
#'   alternative values for corresponding \code{hespdiv} arguments.
#'
#' @param n.split.pts,c.max.iter.no,N.crit,N.loc.crit,c.X.knots,c.Y.knots
#'   Integer vectors specifying alternative values for corresponding arguments.
#'
#' @param N.rel.crit,N.loc.rel.crit,S.crit,S.rel.crit Numeric vectors with values
#'   between 0 and 1.
#'
#' @param Q.crit,c.Q.crit,c.crit.improv Numeric vectors specifying alternative
#'   threshold or improvement criteria.
#'
#' @param c.corr.term Numeric vector with values between 0.01 and 0.2.
#'
#' @param generalize.f,compare.f Lists of functions defining custom similarity
#'   or generalization methods.
#'
#' @param maximize Logical vector of the same length as \code{compare.f},
#'   indicating whether the corresponding similarity metric should be maximized.
#'
#' @param method Character vector specifying predefined similarity metrics.
#'
#' @param .run.id Integer. Runs with indices less than or equal to this value
#'   will be skipped. This can be used to resume an interrupted analysis.
#'
#' @return An object of class \code{hsa}, containing:
#' \describe{
#'   \item{\bold{Alternatives}}{A list of alternative \code{hespdiv} results
#'   produced during the sensitivity analysis.}
#'   \item{\bold{Basis}}{The reference \code{hespdiv} object whose arguments
#'   were perturbed.}
#' }
#'
#' @details
#' \subsection{Difference Between \code{"hsa"} and other sensitivity analysis functions}{
#' The \code{hsa_detailed} function evaluates all combinations of provided
#' argument values, resulting in dense sampling of the parameter space at
#' substantial computational cost. In contrast, \code{hsa} samples the parameter
#' space stochastically and is generally more suitable for exploratory or
#' large-scale sensitivity analyses.
#'
#' In \code{hsa_detailed}, alternative argument values are provided as lists,
#' whereas in \code{hsa} they are supplied as vectors or lists depending on
#' the argument.
#'
#'
#' The \code{hsa_sample_constrained} function performs non-recursive \code{hespdiv}
#' runs for each split-line produced based on different data subsamples. Thus,
#' \code{hsa} is more general, as it allows to inspect sensitivity to other arguments.
#' }
#'
#' \subsection{Paired Arguments}{
#' When \code{data.paired = TRUE}, the same index is used to sample elements
#' of \code{data} and \code{xy.dat}, allowing sensitivity analysis across
#' datasets of differing size or composition. When \code{FALSE}, data and
#' coordinates are sampled independently, enabling analyses based on noise
#' addition or spatial shuffling.
#'
#' Arguments defining custom methods (\code{compare.f}, \code{generalize.f},
#' \code{maximize}) are always treated as paired and must therefore have
#' equal lengths.
#' }
#'
#' @note
#' If a particular run produces a warning or an error, the corresponding list
#' element will contain two components. In case of a warning, these are the
#' resulting \code{hespdiv} object and the warning message. In case of an error,
#' they are the arguments used for the call and the error message.
#'
#' @family functions for hespdiv sensitivity analysis
#' @family functions for hespdiv results post-processing
#'
#' @export
hsa <- function(obj,
                n.runs = 100,
                data.paired = TRUE,
                display = FALSE,
                images.path = NULL,
                pnts.col = 1,
                data = NULL,
                xy.dat = NULL,
                same.n.split = NULL,
                n.split.pts = NULL,
                N.crit = NULL,
                N.rel.crit = NULL,
                N.loc.crit = NULL,
                N.loc.rel.crit = NULL,
                S.crit = NULL,
                S.rel.crit = NULL,
                Q.crit = NULL,
                c.splits = NULL,
                c.Q.crit = NULL,
                c.crit.improv = NULL,
                c.X.knots = NULL,
                c.Y.knots = NULL,
                c.max.iter.no = NULL,
                c.fast.optim = NULL,
                c.corr.term = NULL,
                study.pol = NULL,
                use.chull = NULL,
                generalize.f = NULL,
                maximize = NULL,
                method = NULL,
                compare.f = NULL,
                .run.id = NULL,
                parallel = FALSE,
                RAM = NULL,
                load_prop = 0.8,
                chunk_size = workers*2,
                workers = NULL,
                future_seed = TRUE) {
  eval_one <- function(i, d) {

    if (!is.null(images.path)) {
      fn <- file.path(images.path, sprintf("hsa_run_%04d.png", i))
      png(fn)
      on.exit(dev.off(), add = TRUE)
    }

    call_args <- quote(
      hespdiv(
        data = data,
        xy.dat = xy.dat,
        same.n.split = same.n.split,
        n.split.pts = n.split.pts,
        generalize.f = generalize.f,
        maximize = maximize,
        compare.f = compare.f,
        method = method,
        N.crit = N.crit,
        N.rel.crit = N.rel.crit,
        N.loc.crit = N.loc.crit,
        N.loc.rel.crit = N.loc.rel.crit,
        S.crit = S.crit,
        S.rel.crit = S.rel.crit,
        Q.crit = Q.crit,
        c.splits = c.splits,
        c.Q.crit = c.Q.crit,
        c.crit.improv = c.crit.improv,
        c.X.knots = c.X.knots,
        c.Y.knots = c.Y.knots,
        c.max.iter.no = c.max.iter.no,
        c.fast.optim = c.fast.optim,
        c.corr.term = c.corr.term,
        use.chull = use.chull,
        study.pol = study.pol,
        tracing = NULL,
        pnts.col = pnts.col,
        display = display,
        pacific.region = pacific.region
      )
    )

    env <- list2env(c(d, list(display = display, pacific.region = pacific.region)),
                    parent = parent.frame())

    out <- tryCatch(
      eval(call_args, envir = env),
      warning = function(w) {
        list(Subdivision = suppressWarnings(eval(call_args, envir = env)),
             warning = w)
      },
      error = function(e) {
        list(Arguments = as.list(d), error = e)
      }
    )

    list(Subdivison = out, Arguments = as.list(d))
  }

  draw_one <- function(i) {
    p.id <- sample.int(dat_len, 1)

    local_use_m <- use_m
    if (m.compete) local_use_m <- sample(0:1, 1, prob = c(length(maximize), length(method)))

    if (!local_use_m) cm.id <- sample.int(len_cm, 1) else cm.id <- NA_integer_

    v1 <- data[[p.id]]
    v2 <- if (data.paired) xy.dat[[p.id]] else sample(xy.dat, 1)[[1]]
    v3 <- same.n.split[sample.int(length(same.n.split), 1)]
    v4 <- n.split.pts[sample.int(length(n.split.pts), 1)]
    v5 <- .ifelse(!local_use_m, generalize.f[[cm.id]], NULL)
    v6 <- .ifelse(!local_use_m, maximize[cm.id], NULL)
    v7 <- .ifelse(!local_use_m, compare.f[[cm.id]], NULL)
    v8 <- .ifelse(local_use_m, sample(method, 1), NULL)
    v9 <- N.crit[sample.int(length(N.crit), 1)]
    v10 <- N.rel.crit[sample.int(length(N.rel.crit), 1)]
    v11 <- N.loc.crit[sample.int(length(N.loc.crit), 1)]
    v12 <- N.loc.rel.crit[sample.int(length(N.loc.rel.crit), 1)]
    v13 <- S.crit[sample.int(length(S.crit), 1)]
    v14 <- S.rel.crit[sample.int(length(S.rel.crit), 1)]
    v15 <- .ifelse(local_use_m,
                   ifelse(v8 == "pielou",
                          (Q.crit[Q.crit < 1])[
                            sample.int(length(Q.crit[Q.crit < 1]), 1)
                          ],
                          (Q.crit[Q.crit > 0])[
                            sample.int(length(Q.crit[Q.crit > 0]), 1)
                          ]),
                   NULL)
    v16 <- c.splits[sample.int(length(c.splits), 1)]
    is.maxim <- ifelse(local_use_m, v8 == "pielou", v6)

    if (is.null(v15)) {
      v17 <- NULL
    } else {
      v17 <- .ifelse(v16,
                     ifelse(is.maxim, (c.Q.crit[c.Q.crit <= v15])[
                       sample.int(length(c.Q.crit[c.Q.crit <= v15]), 1)
                     ],
                     (c.Q.crit[c.Q.crit >= v15])[
                       sample.int(length(c.Q.crit[c.Q.crit >= v15]), 1)
                     ]),
                     NULL)
    }

    v18 <- .ifelse(v16, c.crit.improv[sample.int(length(c.crit.improv), 1)], NULL)
    v19 <- .ifelse(v16, c.X.knots[sample.int(length(c.X.knots), 1)], NULL)
    v20 <- .ifelse(v16, c.Y.knots[sample.int(length(c.Y.knots), 1)], NULL)
    v21 <- .ifelse(v16, c.max.iter.no[sample.int(length(c.max.iter.no), 1)], NULL)
    v22 <- .ifelse(v16, c.fast.optim[sample.int(length(c.fast.optim), 1)], NULL)
    v23 <- .ifelse(v16, c.corr.term[sample.int(length(c.corr.term), 1)], NULL)
    v24 <- use.chull[sample.int(length(use.chull), 1)]
    v25 <- .ifelse(v24, obj$call$Call_ARGS$study.pol, sample(study.pol, 1)[[1]])
    v26 <- .ifelse(display, .ifelse(is.pnts, pnts.col[[p.id]], pnts.col), NULL)

    list(data=v1, xy.dat=v2, same.n.split=v3, n.split.pts=v4, generalize.f=v5,
         maximize=v6, compare.f=v7, method=v8, N.crit=v9, N.rel.crit=v10,
         N.loc.crit=v11, N.loc.rel.crit=v12, S.crit=v13, S.rel.crit=v14,
         Q.crit=v15, c.splits=v16, c.Q.crit=v17, c.crit.improv=v18,
         c.X.knots=v19, c.Y.knots=v20, c.max.iter.no=v21, c.fast.optim=v22,
         c.corr.term=v23, use.chull=v24, study.pol=v25, pnts.col=v26)
  }

  if (!inherits(obj, "hespdiv")) stop("obj should have 'hespdiv' class.")

  .check_hsa_arg_structures(
    data = data,
    xy.dat = xy.dat,
    study.pol = study.pol,
    same.n.split = same.n.split,
    c.fast.optim = c.fast.optim,
    use.chull = use.chull,
    c.splits = c.splits,
    n.split.pts = n.split.pts,
    c.max.iter.no = c.max.iter.no,
    N.crit = N.crit,
    N.loc.crit = N.loc.crit,
    c.X.knots = c.X.knots,
    c.Y.knots = c.Y.knots,
    N.rel.crit = N.rel.crit,
    N.loc.rel.crit = N.loc.rel.crit,
    S.crit = S.crit,
    S.rel.crit = S.rel.crit,
    Q.crit = Q.crit,
    c.Q.crit = c.Q.crit,
    c.crit.improv = c.crit.improv,
    c.corr.term = c.corr.term,
    generalize.f = generalize.f,
    compare.f = compare.f,
    maximize = maximize,
    method = method
  )

  if (is.null(.run.id)) .run.id <- -1

  if (!is.numeric(n.runs) || n.runs < 1 || n.runs %% 1 || length(n.runs) != 1)
    stop("'n.runs' must be a positive integer.")
  if (!is.numeric(chunk_size) || chunk_size < 1 || chunk_size %% 1 || length(chunk_size) != 1)
    stop("'chunk_size' must be a positive integer.")

  if (!is.null(images.path) && !display) stop("Attemptig to save empty images.")
  if (!is.null(images.path) && !dir.exists(images.path))
    stop("'images.path' must be an existing directory.")

  # ---  pairing checks ---
  if (data.paired) {
    if (length(data) != length(xy.dat)) {
      stop("If 'data' is paired with 'xy.dat', then their lists must be of equal length.")
    }
  } else {
    if (!is.null(xy.dat)) xy.len <- length(xy.dat)
  }

  # basal data
  if ((!list(obj$call.info$Call_ARGS$xy.dat) %in% xy.dat |
       !list(obj$call.info$Call_ARGS$data) %in% data)) {
    if (data.paired) {
      xy.dat <- c(list(obj$call.info$Call_ARGS$xy.dat), xy.dat)
      data   <- c(list(obj$call.info$Call_ARGS$data), data)
    } else {
      xy.dat <- .set("xy.dat", obj)
      data   <- .set("data", obj)
    }
  }

  same.n.split <- .set("same.n.split", obj)
  n.split.pts  <- .set("n.split.pts", obj)

  # Custom method bookkeeping
  is.custom <- (!is.null(maximize) | obj$call.info$METHOD$method.type == "custom")

  if (obj$call.info$METHOD$method.type == "custom" &
      (!list(obj$call.info$Call_ARGS$generalize.f) %in% generalize.f |
       !list(obj$call.info$Call_ARGS$compare.f) %in% compare.f |
       !list(obj$call.info$Call_ARGS$maximize) %in% maximize)) {
    generalize.f <- c(list(obj$call.info$Call_ARGS$generalize.f), generalize.f)
    compare.f    <- c(list(obj$call.info$Call_ARGS$compare.f), compare.f)
    maximize     <- c(obj$call.info$Call_ARGS$maximize, maximize)
  }

  if (!is.null(method)) {
    method <- sapply(method, .arg_check, name = "metric",
                     NAMES = names(.get_methods()[["biozonation"]]))
    names(method) <- NULL
  }
  method    <- .set("method", obj)
  is.method <- !is.null(method)

  # Pull defaults/augment vectors
  N.crit       <- .set("N.crit", obj)
  N.rel.crit   <- .set("N.rel.crit", obj)
  N.loc.crit   <- .set("N.loc.crit", obj)
  N.loc.rel.crit <- .set("N.loc.rel.crit", obj)
  S.crit       <- .set("S.crit", obj)
  S.rel.crit   <- .set("S.rel.crit", obj)
  Q.crit       <- .set("Q.crit", obj)
  c.Q.crit     <- .set("c.Q.crit", obj)
  c.splits     <- .set("c.splits", obj)
  c.crit.improv <- .set("c.crit.improv", obj)
  c.X.knots    <- .set("c.X.knots", obj)
  c.Y.knots    <- .set("c.Y.knots", obj)
  c.max.iter.no <- .set("c.max.iter.no", obj)
  c.fast.optim <- .set("c.fast.optim", obj)
  c.corr.term  <- .set("c.corr.term", obj)
  use.chull    <- .set("use.chull", obj)

  if (!is.null(study.pol)) {
    if (!is.null(obj$call$Call_ARGS$study.pol)) {
      if (!list(obj$call$Call_ARGS$study.pol) %in% study.pol) {
        study.pol <- c(list(obj$call$Call_ARGS$study.pol), study.pol)
      }
    }
  } else {
    study.pol <- list(obj$call$Call_ARGS$study.pol)
  }

  pnts.col <- .set("pnts.col", obj)
  is.pnts  <- is.list(pnts.col) & length(pnts.col) > 1L
  pacific.region <- obj$call.info$Call_ARGS$pacific.region

  if (any(!use.chull) & is.null(study.pol) & is.null(obj$call$Call_ARGS$study.pol)) {
    stop("If any value of 'use.chull' is FALSE, 'study.pol' must be provided when missing in obj.")
  }

  m.compete <- ifelse(is.custom & is.method, TRUE, FALSE)
  if (!m.compete) use_m <- ifelse(is.method, TRUE, FALSE)

  if (is.method) {
    if ("pielou" %in% method) {
      if (!any(Q.crit < 1))    Q.crit   <- c(-Inf, Q.crit)
      if (!any(c.Q.crit < 1))  c.Q.crit <- c(-Inf, c.Q.crit)
    } else {
      if (!any(Q.crit > 0))    Q.crit   <- c(+Inf, Q.crit)
      if (!any(c.Q.crit > 0))  c.Q.crit <- c(+Inf, c.Q.crit)
    }
  }

  # ---- Serial or parallel execution shares the same draw logic ----
  len_cm  <- length(compare.f)
  dat_len <- length(data)
  run_ids <- seq_len(n.runs)
  active_runs <- run_ids[run_ids > .run.id]

  # Pre-size results (important for resume/skip)
  hes.res <- vector("list", length = n.runs)
  names(hes.res) <- as.character(run_ids)


  # =========================
  # SERIAL PATH
  # =========================
  if (!parallel) {
    for (i in active_runs) {
      message(paste0("Evaluating call: ", i))
      d <- draw_one(i)
      hes.res[[i]] <- eval_one(i, d)
    }

    check.warns(hes.res, .message = FALSE)
    check.errs(hes.res, .message = FALSE)
    return(structure(list(Alternatives = hes.res, Basis = obj), class = "hsa"))
  }

  # =========================
  # PARALLEL PATH
  # =========================
  # Decide workers
  if (is.null(workers)) {
    num_cores <- future::availableCores()
    if (is.null(load_prop)) load_prop <- 0.8
    if (!is.numeric(load_prop) || load_prop <= 0 || load_prop > 1)
      stop("'load_prop' must be a numeric in (0,1].")
    safe_cores <- max(1, floor(num_cores * load_prop))

    if (is.null(RAM)) {
      RAM <- Inf
    } else {
      if (!is.numeric(RAM) || length(RAM) != 1 || RAM < 1)
        stop("'RAM' must be a positive number if not NULL.")
    }
    safe_workers <- max(1, min(safe_cores, RAM))
  } else {
    if (!is.numeric(workers) || length(workers) != 1 || workers < 1)
      stop("'workers' must be a positive number.")
    safe_workers <- workers
  }
  message("Using ", safe_workers, " parallel workers.")

  future::plan(future::multisession, workers = safe_workers)
  options(future.globals.maxSize = 1024^3 * 2)

  # Precompute draws for all runs so workers only evaluate
  draws <- vector("list", length = n.runs)
  for (i in active_runs) draws[[i]] <- draw_one(i)

  # Chunk submission (queue-size control)
  run_chunks <- split(active_runs, ceiling(seq_along(active_runs) / chunk_size))

  res_active <- unlist(
    lapply(run_chunks, function(chunk) {
      future.apply::future_lapply(
        chunk,
        function(i) eval_one(i, draws[[i]]),
        future.seed = future_seed
      )
    }),
    recursive = FALSE
  )

  # Write results back (preserve indexing)
  hes.res[active_runs] <- res_active

  check.warns(hes.res, .message = FALSE)
  check.errs(hes.res, .message = FALSE)

  structure(list(Alternatives = hes.res, Basis = obj), class = "hsa")
}


#' @noRd
.ifelse <- function(test,yes,no){
  if (test){
    yes
  } else  {
    no
  }
}
#' @noRd
check.warns <- function(hes.res, .message = TRUE){
  ids <- which(unlist(lapply(hes.res,function(o)length(o)==2)))

  if (length(ids) > 0) {
    ids <- ids[which(unlist(lapply(hes.res[ids],function(o) names(o)[2] ==
                                     "warning")))]
    if (length(ids) > 0) {
      message(paste0("Warnings were detected in calls: ",
                     paste(as.character(ids),collapse = ", ")))
    } else {
      message("No warnings detected.")
      return(NULL)
    }
  }else {
    message("No warnings detected.")
    return(NULL)
  }
  if (.message) sapply(ids, function(o) hes.res[[o]][[2]]$message)
}

#' @noRd
check.errs <- function(hes.res, .message = TRUE){
  ids <- which(unlist(lapply(hes.res,function(o)length(o)==2)))

  if (length(ids) > 0) {
    ids <- ids[which(unlist(lapply(hes.res[ids],function(o) names(o)[2] ==
                                     "error")))]
    if (length(ids) > 0) {
      message(paste0("Errors were detected in calls: ",
                     paste(as.character(ids),collapse = ", ")))
    } else {
      message("No errors detected.")
      return(NULL)
    }
  }else {
    message("No errors detected.")
    return(NULL)
  }
  if (.message) sapply(ids, function(o) hes.res[[o]][[2]]$message)
}
#' @noRd
.set <- function(name, obi){
  x <- get(name,envir = parent.frame())
  if (is.null(x)){
    x <- eval(parse( text = paste0("obi$call.info$Call_ARGS$",name)))
  } else {
    if (is.list(x)){
      y <- list(eval(parse( text = paste0("obi$call.info$Call_ARGS$",name))))
    } else {
      y <- eval(parse( text = paste0("obi$call.info$Call_ARGS$",name)))
    }
    if (!y %in% x) {
      x <- c(y, x)
    }
  }
  x
}

#' @noRd
.check_hsa_arg_structures <- function(
    data = NULL,
    xy.dat = NULL, study.pol = NULL,
    same.n.split = NULL, c.fast.optim = NULL, use.chull = NULL, c.splits = NULL,
    n.split.pts = NULL, c.max.iter.no = NULL, N.crit = NULL, N.loc.crit = NULL,
    c.X.knots = NULL, c.Y.knots = NULL,
    N.rel.crit = NULL, N.loc.rel.crit = NULL, S.crit = NULL, S.rel.crit = NULL,
    Q.crit = NULL, c.Q.crit = NULL, c.crit.improv = NULL,
    c.corr.term = NULL,
    generalize.f = NULL, compare.f = NULL,
    maximize = NULL,
    method = NULL
) {
  .bad <- function(name, msg) stop(paste0("'", name, "': ", msg), call. = FALSE)

  .is_intish <- function(x) is.numeric(x) && all(is.finite(x)) && all(x %% 1 == 0)

  .chk_list <- function(x, name) {
    if (is.null(x)) return(invisible(TRUE))
    if (!is.list(x)) .bad(name, "must be a list (or NULL).")
    invisible(TRUE)
  }

  .chk_vec_type <- function(x, name, type = c("logical","numeric","integer","character")) {
    if (is.null(x)) return(invisible(TRUE))
    type <- match.arg(type)
    ok <- switch(
      type,
      logical   = is.logical(x),
      numeric   = is.numeric(x),
      integer   = .is_intish(x),
      character = is.character(x)
    )
    if (!ok) .bad(name, paste0("must be a ", type, " vector (or NULL)."))
    invisible(TRUE)
  }

  .chk_range <- function(x, name, lo, hi) {
    if (is.null(x)) return(invisible(TRUE))
    if (!is.numeric(x)) .bad(name, "must be a numeric vector (or NULL).")
    if (any(x < lo | x > hi, na.rm = TRUE)) {
      .bad(name, paste0("values must be in [", lo, ", ", hi, "]."))
    }
    invisible(TRUE)
  }

  .chk_xy_list <- function(x, name) {
    if (is.null(x)) return(invisible(TRUE))
    if (!is.list(x)) .bad(name, "must be a list of data.frames (or NULL).")
    for (i in seq_along(x)) {
      d <- x[[i]]
      if (!is.data.frame(d)) .bad(name, paste0("element ", i, " must be a data.frame."))
      if (ncol(d) != 2) .bad(name, paste0("element ", i, " must have exactly 2 columns: x and y."))
      if (!all(c("x","y") %in% names(d))) .bad(name, paste0("element ", i, " must have columns named 'x' and 'y'."))
      if (!is.numeric(d[["x"]]) || !is.numeric(d[["y"]])) .bad(name, paste0("element ", i, " columns 'x' and 'y' must be numeric."))
    }
    invisible(TRUE)
  }

  .chk_fun_list <- function(x, name) {
    if (is.null(x)) return(invisible(TRUE))
    if (!is.list(x)) .bad(name, "must be a list of functions (or NULL).")
    if (any(!vapply(x, is.function, logical(1)))) .bad(name, "must contain only functions.")
    invisible(TRUE)
  }

  # data: list
  .chk_list(data, "data")

  # xy.dat, study.pol: lists of 2-col data.frames (x,y)
  .chk_xy_list(xy.dat, "xy.dat")
  .chk_xy_list(study.pol, "study.pol")

  # logical vectors
  .chk_vec_type(same.n.split, "same.n.split", "logical")
  .chk_vec_type(c.fast.optim, "c.fast.optim", "logical")
  .chk_vec_type(use.chull,    "use.chull",    "logical")
  .chk_vec_type(c.splits,     "c.splits",     "logical")

  # integer vectors (integer-like numeric)
  .chk_vec_type(n.split.pts,   "n.split.pts",   "integer")
  .chk_vec_type(c.max.iter.no, "c.max.iter.no", "integer")
  .chk_vec_type(N.crit,        "N.crit",        "integer")
  .chk_vec_type(N.loc.crit,    "N.loc.crit",    "integer")
  .chk_vec_type(c.X.knots,     "c.X.knots",     "integer")
  .chk_vec_type(c.Y.knots,     "c.Y.knots",     "integer")

  # numeric vectors with bounds
  .chk_range(N.rel.crit,     "N.rel.crit",     0, 1)
  .chk_range(N.loc.rel.crit, "N.loc.rel.crit", 0, 1)
  .chk_range(S.crit,         "S.crit",         0, 1)
  .chk_range(S.rel.crit,     "S.rel.crit",     0, 1)

  # numeric vectors (no bounds stated)
  .chk_vec_type(Q.crit,        "Q.crit",        "numeric")
  .chk_vec_type(c.Q.crit,      "c.Q.crit",      "numeric")
  .chk_vec_type(c.crit.improv, "c.crit.improv", "numeric")

  # c.corr.term bounds
  .chk_range(c.corr.term, "c.corr.term", 0.01, 0.2)

  # lists of functions
  .chk_fun_list(generalize.f, "generalize.f")
  .chk_fun_list(compare.f,    "compare.f")

  # maximize: logical vector; if compare.f also provided, lengths must match
  .chk_vec_type(maximize, "maximize", "logical")
  if (!is.null(maximize) && !is.null(compare.f) && length(maximize) != length(compare.f)) {
    .bad("maximize", paste0("must have the same length as 'compare.f' (", length(compare.f), ")."))
  }

  # method: character vector
  .chk_vec_type(method, "method", "character")

  invisible(TRUE)
}
