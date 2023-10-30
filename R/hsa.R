#' HespDiv Sensitivity Analysis
#'
#' This function is the main function that performs HespDiv sensitivity analysis.
#' The function takes a provided hespdiv object as a base and produces a
#' desired number of alternative versions by randomly sampling new values for
#' its arguments from a provided selection of values for each argument.
#'
#' @param obj A 'hespdiv' class object.
#' @param n.runs An integer representing how many alternative hespdiv calls to evaluate.
#' @param data.paired A Boolean controlling whether the provided alternative
#' values of 'data' are paired with 'xy.dat'.
#' @param display A Boolean indicating the value of the "display" argument in
#' each hespdiv call.
#' @param images.path  A path to an existing directory where PNG images of the
#' displayed results will be saved. If NULL (default), images won't be saved.
#' @param pnts.col the value of the "pnts.col" argument in each hespdiv call.
#' @param data A list containing matrices, time-series, lists, data frames,
#' vectors, or other data structures.
#' @param xy.dat,study.pol Lists of data frames with two columns: 'x' and 'y'.
#' @param same.n.split,c.fast.optim,use.chull,c.splits A Boolean value (if used, should be different from the
#' one in the basal hespdiv call).
#' @param n.split.pts,c.max.iter.no,N.crit,N.loc.crit,c.X.knots,c.Y.knots Numeric integer vectors.
#' @param N.rel.crit,N.loc.rel.crit,S.crit,S.rel.crit Numeric vectors with values
#' between 0 and 1.
#' @param Q.crit,c.Q.crit,c.crit.improv Numeric vectors
#' @param c.corr.term  A numeric vector with values between 0.01 and 0.2.
#' @param generalize.f,compare.f Lists of functions.
#' @param maximize A logical vector of the same length as 'compare.f' list.
#' @param method A character vector.
#' @param .run.id An integer. Runs up to this id will be skipped.
#' @return A 'hsa' class object. It is a list of two items:
#' \describe{
#' \item{\bold{Alternatives:}}{ A list containing the produced
#' alternative hespdiv objects. }
#' \item{\bold{Basis}}{ The basal hespdiv object whose call was modified to produce
#' alternative subdivisions.}
#'   }
#' @details
#' \subsection{Difference Between "hsa" And "hsa_detailed"}{
#' The major difference between "hsa_detailed" and "hsa" is that the former produces
#' all possible hespdiv calls from combinations of the provided hespdiv arguments.
#' Therefore, it samples a much smaller segment of the parameter space but more
#' densely, requiring much more computation time. Although such behavior may
#' be desired in some cases, the "hsa" function is generally more suitable for
#' performing hespdiv sensitivity analysis.
#'
#' Additionally, alternative values for hespdiv arguments in the "hsa_detailed"
#' function are provided in lists, whereas in the "hsa" function, they are
#' provided in vectors or lists (depending on the argument).
#'}
#' \subsection{Paired Arguments}{
#' The argument 'data.paired = TRUE' means that in each produced hespdiv
#' call, the index that samples lists of 'xy.dat' and 'data' will be the same.
#' This way, a different number of observations or entirely different data sets
#' can be used for sensitivity analysis. If 'data.paired = FALSE', then the number
#' of observations must be kept the same in all provided elements of 'data' and
#' 'xy.dat' lists. This option allows to add noise to the data or coordinates,
#' or shuffle them entirely in order to test the significance of the detected spatial
#' structure.
#'
#' Arguments that define custom methods ('compare.f', 'generalize.f' and
#' 'maximize') are paired by default. Therefore, their lengths must also be
#' the same.
#' }
#' @note If a particular call produces a warning or error, a list of length 2
#' will be returned for that call. If a warning was produced, the first element
#' of the list will hold the created hespdiv object, and the second element
#' will contain the warning message. In case of an error, the first element will
#'  be a list of arguments used to produce the call, and the second element will
#'  contain the error message.
#' @family {functions for hespdiv sensitivity analysis}
#' @family {functions for hespdiv post-prossesing}
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
                .run.id = NULL){
  if (!inherits(obj,"hespdiv"))
    stop("obj should have 'hespdiv' class.")
  if (is.null(.run.id)){
    .run.id <- -1
  }

  if (!is.null(images.path) & !display){
    stop("Attemptig to save empty images.")
  }
  if (data.paired){
    if(length(data)!= length(xy.dat)){
      stop(paste0("If 'data' is paired with 'xy.dat',",
                  " then their lists must be of equal length."))
    }
  } else {
    if (!is.null(xy.dat)){
      xy.len <- length(xy.dat)
    }
  }

  if((!list(obj$call.info$Call_ARGS$xy.dat) %in% xy.dat |
      !list(obj$call.info$Call_ARGS$data) %in% data)){
    if (data.paired){
      xy.dat <- c(list(obj$call.info$Call_ARGS$xy.dat), xy.dat)
      data <- c(list(obj$call.info$Call_ARGS$data), data)
    } else {
      xy.dat <- .set("xy.dat",obj)
      data <- .set("data",obj)
    }
  }
  same.n.split <- .set("same.n.split",obj)
  n.split.pts <- .set("n.split.pts",obj)
  if (!is.null(maximize) | obj$call.info$METHOD$method.type == "custom"){
    is.custom <- TRUE
  } else {
    is.custom <- FALSE
  }
  if (obj$call.info$METHOD$method.type == "custom" &
      (!list(obj$call.info$Call_ARGS$generalize.f) %in% generalize.f |
       !list(obj$call.info$Call_ARGS$compare.f) %in% compare.f |
       !list(obj$call.info$Call_ARGS$maximize) %in% maximize) ){
    generalize.f <- c(list(obj$call.info$Call_ARGS$generalize.f), generalize.f)
    compare.f <- c(list(obj$call.info$Call_ARGS$compare.f), compare.f)
    maximize <- c(obj$call.info$Call_ARGS$maximize, maximize)
  }


  if(!is.null(method)){
    method <- sapply(method,.arg_check,name = "metric", NAMES = names(
      .get_methods()[['biozonation']]))
    names(method) <- NULL
  }
  method <- .set("method",obj)
  is.method <- !is.null(method)

  N.crit <- .set("N.crit",obj)
  N.rel.crit <- .set("N.rel.crit",obj)
  N.loc.crit <- .set("N.loc.crit",obj)
  N.loc.rel.crit <- .set("N.loc.rel.crit",obj)
  S.crit <- .set("S.crit",obj)
  S.rel.crit <- .set("S.rel.crit",obj)
  Q.crit <- .set("Q.crit",obj)
  c.Q.crit <- .set("c.Q.crit",obj)
  c.splits <- .set("c.splits",obj)
  c.crit.improv <- .set("c.crit.improv",obj)
  c.X.knots <- .set("c.X.knots",obj)
  c.Y.knots <- .set("c.Y.knots",obj)
  c.max.iter.no <- .set("c.max.iter.no",obj)
  c.fast.optim <- .set("c.fast.optim",obj)
  c.corr.term <- .set("c.corr.term",obj)
  use.chull <- .set("use.chull",obj)
  if (!is.null(study.pol)){
    if (!is.null(obj$call$Call_ARGS$study.pol)){
      if (!list(obj$call$Call_ARGS$study.pol) %in% study.pol){
        study.pol <- c(list(obj$call$Call_ARGS$study.pol), study.pol)
      }
    }
  } else {
    study.pol <- list(obj$call$Call_ARGS$study.pol)
  }
  pnts.col <- .set("pnts.col",obj)

  if (ifelse(is.null(maximize), TRUE, all(maximize)) &
      !any(c("morisita","sorensen","horn.morisita")
           %in% method)){
    if (all(max(c.Q.crit) > Q.crit)){
      stop(paste0("Maximum value of 'c.Q.crit' is higher than all 'Q.crit'"))
    }
  } else {
    if (ifelse(is.null(maximize), TRUE, all(!maximize)) & !"pielou" %in% method){
      if (all(min(c.Q.crit) < Q.crit)){
        stop(paste0("Minimum value of 'c.Q.crit' is higher than all 'Q.crit'",))
      }
    } else {
      if (all(min(c.Q.crit) < Q.crit) | all(max(Q.crit) > c.Q.crit)){
        if (all(min(c.Q.crit) < Q.crit))
          stop(paste0("Minimum value of 'c.Q.crit' is lower than all 'Q.crit'."))
        if (all(max(Q.crit) > c.Q.crit))
          stop(paste0("Maximum value of 'Q.crit' is higher than all 'Q.crit'."))
      }
      if (all(max(c.Q.crit) > Q.crit) | all(min(Q.crit) < c.Q.crit)){
        if (all(max(c.Q.crit) > Q.crit))
          stop(paste0("Maximum value of 'c.Q.crit' is higher than all 'Q.crit'."))
        if (all(min(Q.crit) < c.Q.crit))
          stop(paste0("Minimum value of 'Q.crit' is lower than all 'c.Q.crit'."))
      }
    }
  }


  is.pnts <- is.list(pnts.col) & length(pnts.col) > 1L

  pacific.region <- obj$call.info$Call_ARGS$pacific.region

  if (any(!use.chull) & is.null(study.pol) & is.null(obj$call$Call_ARGS$study.pol)){
    stop(paste0("If any value of 'use.chull' is FALSE, 'study.pol' must be",
                " provided, when its missing in obj."))
  }

  m.compete <- ifelse(is.custom & is.method, TRUE, FALSE)
  if (!m.compete){
    use_m <- ifelse(is.method, TRUE, FALSE)
  }
  if (is.method){
    if ("pielou" %in% method){
      if (!any(Q.crit < 1)){
        Q.crit <- c(-Inf, Q.crit)
      }
      if (!any(c.Q.crit < 1)){
        c.Q.crit <- c(-Inf, c.Q.crit)
      }
    } else {
      if (!any(Q.crit > 0)){
        Q.crit <- c(+Inf, Q.crit)
      }
      if (!any(c.Q.crit > 0)){
        c.Q.crit <- c(+Inf, c.Q.crit)
      }
    }
  }

  len_m <- length(generalize.f)
  len_cm <- length(compare.f)
  dat_len <- length(data)
  hes.res <- vector(mode = "list", length = n.runs)
  for (i in 1:n.runs){
    message(paste0("Evaluating call: ", i))
    p.id <- sample(1:dat_len,1)

    if (m.compete){ use_m <- sample(0:1,1,prob = c(length(maximize),
                                                   length(method)))}
    if (!use_m) {
      cm.id <- sample(1:len_cm,1)
    }
    if (!is.null(images.path)) png(paste0(images.path,"//",n.runs))
    v1 <- data[[p.id]]
    v2 <- ifelse(data.paired, xy.dat[p.id], sample(xy.dat,1))[[1]]
    v3 <- sample(same.n.split,1)
    v4 <- sample(n.split.pts,1)
    v5 <- .ifelse(!use_m, generalize.f[[cm.id]], NULL)
    v6 <- .ifelse(!use_m, maximize[cm.id], NULL)
    v7 <- .ifelse(!use_m, compare.f[[cm.id]], NULL)
    v8 <- .ifelse(use_m, sample(method,1),NULL)
    v9 <- sample(N.crit,1)
    v10 <- sample(N.rel.crit,1)
    v11 <- sample(N.loc.crit,1)
    v12 <- sample(N.loc.rel.crit,1)
    v13 <- sample(S.crit,1)
    v14 <- sample(S.rel.crit,1)
    v15 <- .ifelse(use_m, ifelse(v8 == "pielou", sample(Q.crit[Q.crit < 1],1),
                         sample(Q.crit[Q.crit > 0],1)), NULL )
    v16 <- sample(c.splits,1)
    is.maxim <- ifelse(use_m, v8 == "pielou", v6)
    if (is.null(v15)){
      v17 <- NULL
    } else {
      v17 <- .ifelse(v16, ifelse(is.maxim, sample(c.Q.crit[c.Q.crit <= v15],1),
                                 sample(c.Q.crit[c.Q.crit >= v15],1)), NULL)
    }
    v18 <- .ifelse(v16, sample(c.crit.improv,1), NULL)
    v19 <- .ifelse(v16, sample(c.X.knots,1), NULL)
    v20 <- .ifelse(v16, sample(c.Y.knots,1), NULL)
    v21 <- .ifelse(v16, sample(c.max.iter.no,1), NULL)
    v22 <- .ifelse(v16, sample(c.fast.optim,1), NULL)
    v23 <- .ifelse(v16, sample(c.corr.term,1), NULL)
    v24 <- sample(0:1, 1, prob = c(length(study.pol),ifelse(any(use.chull),1,0)))
    v25 <- .ifelse(v24, obj$call$Call_ARGS$study.pol, sample(study.pol,1))[[1]]
    v26 <- .ifelse(display, .ifelse(is.pnts, pnts.col[p.id],
                                    list(pnts.col))[[1]],NULL)
    if (i < .run.id) {
      next
    }
    call_args <- expression(f1(data = v1,
                               xy.dat = v2,
                               same.n.split = v3,
                               n.split.pts = v4,
                               generalize.f = v5,
                               maximize = v6,
                               compare.f = v7,
                               method = v8,
                               N.crit = v9,
                               N.rel.crit = v10,
                               N.loc.crit = v11,
                               N.loc.rel.crit = v12,
                               S.crit = v13,
                               S.rel.crit = v14,
                               Q.crit = v15,
                               c.splits = v16,
                               c.Q.crit = v17,
                               c.crit.improv = v18,
                               c.X.knots = v19,
                               c.Y.knots = v20,
                               c.max.iter.no = v21,
                               c.fast.optim = v22,
                               c.corr.term = v23,
                               use.chull = v24,
                               study.pol = v25,
                               tracing = NULL,
                               pnts.col = v26,
                               display = display,
                               pacific.region = pacific.region
    ))
    f1 <- hespdiv
    hes.res[[i]] <- list(
      Subdivison = tryCatch(eval(call_args),
                            error = function(cond){
                              f1 <- list
                              message(paste(cond,"\n"))
                              return(list(Arguments = eval(call_args),
                                          error = cond))
                            },
                            warning = function(cond){
                              message(paste(cond,"\n"))
                              return(list(Subdivision = eval(call_args),
                                          warning = cond))
                            }
      ),
      Arguments = {f1 <- list; eval(call_args)}
    )
    if (!is.null(images.path)) dev.off()
  }
  if (length(hes.res) != n.runs){
    hes.res <- c(hes.res,list(NULL))
  }
  names(hes.res) <- 1:n.runs
  check.warns(hes.res, .message = FALSE)
  check.errs(hes.res, .message = FALSE)
  structure(list(Alternatives = hes.res, Basis = obj),class = 'hsa')
}
#' @noRd
.use_args <- function(obj ,display = TRUE,tracing = NULL, args = NULL){
  if (is.null(args)){
    args <- obj$call.info$Call_ARGS
  }
  hespdiv(data = args$data,
          n.split.pts = args$n.split.pts,
          generalize.f = args$generalize.f,
          maximize = args$maximize,
          method = args$method,
          same.n.split = args$same.n.split,
          compare.f = args$compare.f,
          N.crit = args$N.crit,
          N.rel.crit = args$N.rel.crit,
          N.loc.crit = args$N.loc.crit,
          N.loc.rel.crit = args$N.loc.rel.crit,
          S.crit = args$S.crit,
          S.rel.crit = args$S.rel.crit,
          Q.crit = args$Q.crit,
          c.splits = args$c.splits,
          c.Q.crit = args$c.Q.crit,
          c.crit.improv = args$c.crit.improv,
          c.X.knots = args$c.X.knots,
          c.Y.knots = args$c.Y.knots,
          xy.dat = args$xy.dat,
          c.max.iter.no = args$c.max.iter.no,
          c.fast.optim = args$c.fast.optim,
          c.corr.term = args$c.corr.term,
          study.pol = args$study.pol,
          use.chull = args$use.chull,
          tracing = tracing,
          pnts.col = args$pnts.col,
          display = display,
          pacific.region = args$pacific.region)
}
#' @noRd
.use <- function(x) ifelse(is.null(x), FALSE, TRUE)
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
