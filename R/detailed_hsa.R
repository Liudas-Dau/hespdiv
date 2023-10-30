#' Detailed HespDiv Sensitivity Analysis
#'
#' This function is one of the two that perform HespDiv sensitivity analysis.
#' It creates and evaluates alternative hespdiv calls, according to the desired
#' changes in method, data and other subdivision criteria arguments. As a result,
#' it returns alternative hespdiv objects that can be directly compared with
#' the original hespdiv object and with each other using \code{plot.hsa} and
#' \code{hsa_quant} functions.
#'
#' @param obj An object of hespdiv class. The base object whose call will be
#' modified to produce alternative hespdiv objects.
#' @param comb.args A Boolean value. Do you want to combine the provided
#' argument values to make alternative hespdiv calls? If not, then at once only one
#' argument will be modified, trying all provided values for it one by one.
#' @param pick.n.args A numeric vector that controls how many arguments would
#' you like to change at once in hespdiv runs. Multiple values allowed.
#' @param comb.type A character determining how combinations of argument values
#' are selected. Possible values: "all", "random", or "handpicked".
#' @param n.combs An integer controlling how many argument value combinations
#' should be randomly selected from all possible combinations when comb.type
#' is "random".
#' @param paired Boolean. Are the provided hespdiv arguments 'data' and 'xy.dat'
#' paired?
#' @param display A Boolean value. The value of the "display" argument in each
#' hespdiv call.
#' @param images.path  A path to an existing directory where PNG images of the
#' displayed results will be saved. If NULL (default), images won't be saved.
#' @param pnts.col The value of the "pnts.col" argument in each hespdiv call.
#' @param data A list containing matrices, time-series, lists, data frames,
#' vectors, or other data structures.
#' @param xy.dat,study.pol Lists of data frames with two columns: 'x' and 'y'.
#' @param same.n.split,c.fast.optim,use.chull,c.splits Lists with a Boolean value (if used, should be different from the
#' one in the basal hespdiv call).
#' @param n.split.pts,c.max.iter.no,N.crit,N.loc.crit,c.X.knots,c.Y.knots Lists with integer values.
#' @param N.rel.crit,N.loc.rel.crit,S.crit,S.rel.crit Lists with values between 0 and 1.
#' @param Q.crit,c.Q.crit,c.crit.improv Lists of numeric values.
#' @param c.corr.term  A list of numeric values between 0.01 and 0.2.
#' @param generalize.f,compare.f Lists of functions.
#' @param maximize A list of logical values of the same length as 'compare.f' list.
#' @param method A list of character values.
#' @return A 'hsa' class object. It is a list of two items:
#' \describe{
#' \item{\bold{Alternatives:}}{ A list containing the produced
#' alternative hespdiv objects. }
#' \item{\bold{Basis}}{ The basal hespdiv object whose call was modified to produce
#' alternative subdivisions.}
#'   }
#' @importFrom utils combn
#' @importFrom grDevices dev.off png
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
#' \subsection{Internally Set Default Argument Values}{
#' When comb.args is TRUE, the default value of comb.type is "all".
#'
#' When comb.args is TRUE and pick.n.args is NULL (default), the value of
#' pick.n.args will be changed to a vector 1:N, where N is the maximum possible
#' value of pick.n.args. The maximum possible value for pick.n.args depends on
#' the hespdiv arguments provided. Each hespdiv argument that
#' influences the results is counted as one, except for "data" and "xy.dat"
#' when paired is TRUE, and all four arguments ("method", "compare.f",
#' "generalize.f", and "maximize") that define the subdivision method, as the
#' pair/group of them is counted as one. Therefore, N can vary from 1 (single
#' argument provided) to 22 (all arguments provided and paired is FALSE). If
#' comb.args is FALSE, then pick.n.args should be NULL. Using pick.n.args = 1
#' is the same as setting comb.type to FALSE.
#'}
#'\subsection{Paired Arguments}{
#' If paired is TRUE, the "data" and "xy.dat" elements with the same index are
#' treated as one value of the same argument. Therefore, the provided lists of
#' "data" and "xy.dat" should be of the same length. Pairing of "data" and
#' "xy.dat" can be useful, for example, when you want to re-run hespdiv after
#' adding or removing some observations (these changes should be made in both
#' "xy.dat" and "data") to test how hespdiv results are influenced by some
#' extra observations or the number of observations in general. When paired
#' is FALSE, the number of observations in "data" and "xy.dat" must be the same
#' as it was in the call of the base hespdiv object. This option allows you to
#' re-run hespdiv after adding some noise to the object features (via changes
#' in "data") or coordinates (via changes in "xy.dat") to test how hespdiv
#' results are influenced by the data itself or localization.
#'
#' By default, arguments determining the custom method ("compare.f",
#' "generalize.f", "maximize") are paired, similar to how "data" and "xy.dat"
#' are paired when paired is TRUE. Thus, the lists of "compare.f",
#' "generalize.f", and "maximize" should be of the same length.
#' }
#' @family {functions for hespdiv sensitivity analysis}
#' @family {functions for hespdiv post-prossesing}
#' @note Use "pnts.col" of length >1 only when the number of observations does
#' not change.
#'
#' If a particular call produced a warning or error, then a list of
#' length 2 will be returned for that call. If a warning was produced, then the
#' first element of the list will hold the created hespdiv object, and the
#' second element will contain the warning message. In the case of an error,
#' the first element will be a list of arguments used to produce the call,
#' and the second element will contain the error message.
#' @author Liudas Daumantas
#' @export
hsa_detailed <- function(
    obj,
    comb.args = TRUE,
    pick.n.args = NULL,
    comb.type = NULL,
    n.combs = NULL,
    display = TRUE,
    images.path = NULL,
    paired = NULL,
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
    compare.f = NULL
) {
  use.chull <- .to.list(use.chull)
  c.splits <- .to.list(c.splits)
  same.n.split <- .to.list(same.n.split)
  c.fast.optim <- .to.list(c.fast.optim)
  if (is.null(display))
    display <- FALSE
  if (!is.null(images.path)){
    if (!display)
      stop(paste0("If you wish to save images, then set 'display' = TRUE"))
  }
  if (obj$call.info$METHOD$method.type == "preset"){
    obj$call.info$Call_ARGS[which(names(obj$call.info$Call_ARGS)
                                  == "generalize.f")] <- list(NULL)
    obj$call.info$Call_ARGS[which(names(obj$call.info$Call_ARGS)
                                  == "compare.f")] <- list(NULL)
  }
  if (comb.args & is.null(comb.type)){
    comb.type <- "all"
  }
  if (!comb.args & !is.null(comb.type)){
    stop(paste0("Check 'comb.args' and 'comb.type' arguments.",
                "\nArgument 'comb.args' is FALSE, but 'comb.type' is not NULL."
    ))
  }
  if (comb.args) {
    comb.type <- .arg_check(name = "comb.type", given = comb.type, NAMES =
                              c("all", "random", "handpicked"))
    if (comb.type == "random" & (is.null(n.combs) |
                                 length(n.combs) > 1)) {
      stop(paste0("Check 'n.combs' and 'n.combs' arguments.",
                  "\nWhen 'comb.type' is set to 'random'",
                  ", 'n.combs' should be a numeric integer."))

    }
  }
  if (!is.null(n.combs)){
    if (!comb.args)
      stop(paste0("Check 'comb.args' and 'n.combs' arguments.",
                  "\nArgument 'comb.args' is FALSE, but 'n.combs' is not NULL."
      ))
    if (comb.type != "random")
      stop(paste0("Check 'comb.type' and 'n.combs' arguments.",
                  "\nIf 'n.combs' is not NULL, 'comb.type' should be set to",
                  " 'random'."
      ))
  }

  prov.args.names <- ls()[!ls() %in% c("obj", "comb.args","pick.n.args","paired",
                                       "display","images.path", "pnts.col", "comb.type",
                                       "n.combs")]
  c.pars <- any(c(!is.null(c.Q.crit),!is.null(c.crit.improv),!is.null(c.X.knots),
                  !is.null(c.Y.knots),!is.null(c.max.iter.no),!is.null(c.fast.optim),
                  !is.null(c.corr.term)))
  if (is.null(c.splits) & !obj$call.info$Call_ARGS$c.splits & c.pars){
    stop(paste0("You provide curve generation arguments (starting with 'c.'), ",
                "but 'c.splits' in the base hespdiv object is FALSE and ",
                "you do not provide its alternative value",
                " (c.splits = list(TRUE))."))
  }
  prov.args.l <- sapply(prov.args.names,get,environment())
  ind <- which(!sapply(prov.args.l,is.null))
  if(length(ind) == 0) {
    stop("All provided hespdiv arguments for sensitivity analysis are NULL")
  }

  hes.names <- (names(obj$call.info$Call_ARGS)[
    !names(obj$call.info$Call_ARGS) %in%
      c("tracing", "display", "pnts.col","pacific.region")])[ind]
  hes.mods <- prov.args.l[hes.names]

  if ("data" %in% hes.names | "xy.dat" %in% hes.names){

    if (is.null(paired) & "data" %in% hes.names & "xy.dat" %in% hes.names){
      stop(paste0("Argument 'paired' is NULL, but 'data' and 'xy.dat' are",
                  " provided. \nPlease set 'paired' to either FALSE or TRUE."))
    }
    if ("data" %in% hes.names){
      if (any(!unlist(lapply(lapply(hes.mods[["data"]],class),identical,
                             class(obj$call.info$Call_ARGS$data)))) |
          any(!unlist(lapply(lapply(hes.mods[["data"]],mode),identical,
                             mode(obj$call.info$Call_ARGS$data))))) {
        message(paste0("Class or mode of the provided 'data' is different from",
                       " 'data' in 'obj'.",
                       " \nAre you sure you wish to continue? \nAnswer Y / N."))
        answer <- readline(prompt = "")
        answer <- .answer_check(given = answer,NAMES = c("Yes", "No"))
        if (answer == 'No'){
          return(NULL)
        }
      }
      dat_lens <- unlist(lapply(hes.mods[["data"]],function(o){
        if (is.data.frame(o) | is.matrix(o)){
          nrow(o)
        } else {length(o)}
      }))
    }
    if ("xy.dat" %in% hes.names){
      xy_lens <- unlist(lapply(hes.mods[["xy.dat"]], nrow))
    }
    if (paired){
      if (!("data" %in% hes.names & "xy.dat" %in% hes.names)){
        stop("Argument 'paired' = TRUE, but 'data' or 'xy.dat' is not provided.")
      }
      if (length(hes.mods[["data"]]) != length(hes.mods[["xy.dat"]])){
        stop(paste0("Argument 'paired' = TRUE, but 'data' and 'xy.dat' are not",
                    " of the same lengths"))
      }
      if(!all(dat_lens == xy_lens)){
        stop(paste0("Argument 'paired' = TRUE, but some pairs of 'data' and",
                    " 'xy.dat' have different number of observations"))
      }

    } else {
      if ("data" %in% hes.names) {
        if (!all(dat_lens == nrow(obj$call.info$Call_ARGS$xy.dat))){
          stop(paste0("Length of some elements of the provided 'data' is not the",
                      " same as the number of coordinates in 'obj'."))
        }
      }
      if ("xy.dat" %in% hes.names){
        if (!all(xy_lens == nrow(obj$call.info$Call_ARGS$xy.dat))){
          stop(paste0("Length of some elements of the provided 'xy.dat' is not the",
                      " same as the number of data points in 'obj'."))
        }
      }
    }
  } else{

    if (!is.null(paired)){
      stop(paste0("Argument 'paired' is not NULL, but both 'data' and 'xy.dat'",
                  " were not provided. \nPlease check the input."))
    }

  }

  if (!all(unlist(lapply(hes.mods, is.list)) &
           unlist(lapply(hes.mods, class)) == "list" )){
    stop("Arguments provided for hespdiv sensitivity analysis should put into lists.")
  }
  if ("study.pol" %in% hes.names){

    if ("study.pol" %in% hes.names & obj$call.info$Call_ARGS$use.chull &
        !"use.chull" %in% hes.names){
      stop(paste0("If 'use.chull' = TRUE, changes in 'study.pol' can not change",
                  " the subdivision results.","\nIf you want to test the impact of",
                  "'study.pol on subdivisions, please provide 'use.chull' = FALSE."))
    }
    if (any(unlist(lapply(hes.mods[["study.pol"]], is.null)))){
      stop(paste0('"study.pol = list(NULL, ...)" is invalid, since NULL polygon',
                  " values cannot change the results."))
    }
  }
  if ("method" %in% hes.names){
    hes.mods[["method"]] <- lapply(hes.mods[["method"]],.arg_check,name =
                                     "metric", NAMES = names(
                                       .get_methods()[['biozonation']]))
  }
  if (is.null(paired)){
    paired <- FALSE
  }

  if(paired){
    org_dataxy <- obj$call.info$Call_ARGS[c("data", "xy.dat")]
    names(org_dataxy) <- NULL
    given_dataxy <- vector(mode = "list", length = length(hes.mods[["data"]]))
    for (i in seq(length(given_dataxy))) {
      given_dataxy[[i]] <-  list(hes.mods[["data"]][[i]],
                                 hes.mods[["xy.dat"]][[i]])
      if (identical(given_dataxy[[i]],org_dataxy)){
        stop(paste0("Element ",i," of the paired data and xy.data datasets ",
                    "is not new (the same as in base hespdiv object)."))
      }
    }
    if (any(duplicated(given_dataxy))){
      stop(paste0("There are duplicated values in paired data and xy.data
                  arguments. Remove duplicated pairs.'"))
    }
  }

  cond <- c("compare.f", "generalize.f", "maximize") %in% hes.names
  if (any(cond) & !all(cond)){
    stop(paste0("If you provide custom method, all required arguments",
                ' must be provided:"compare.f", "generalize.f", "maximize".',
                "\nMissing argument(s): ",
                paste(c("compare.f", "generalize.f", "maximize")[!cond],
                      collapse = ", ")))
  }
  if (all(cond)){
    if (length(unique(c(
      length(hes.mods[["compare.f"]]),
      length(hes.mods[["generalize.f"]]),
      length(hes.mods[["maximize"]])))) != 1)
      stop(paste0("Custom method arguments are paired, therefore, ",
                  "they must have the same number of elements."))
    org_method <-  obj$call.info$Call_ARGS[c("compare.f", "generalize.f",
                                             "maximize")]
    names(org_method) <- NULL
    given_methods <- vector(mode = "list",
                            length = length(hes.mods[["compare.f"]]))
    for (i in seq(length(given_methods))) {
      given_methods[[i]] <- list(hes.mods[["compare.f"]][[i]],
                                 hes.mods[["generalize.f"]][[i]],
                                 hes.mods[["maximize"]][[i]])
      if (identical(given_methods[[i]], org_method) &
          !is.null(org_method[[1]])){
        stop(paste0("Element ",i," of the paired custom method arguments ",
                    "is not new (the same as in base hespdiv object)."))
      }
    }
    if (any(duplicated(given_methods))){
      stop(paste0("There are duplicated values in the paired custom method",
                  " arguments. Remove duplicated custom methods.'"))
    }
  }


  for (var.n in hes.names){
    if ((!var.n %in%  c("data", "xy.dat") | !paired) & !var.n %in%
        c("compare.f", "generalize.f", "maximize")){
      if (any(
        unlist(
          lapply(
            hes.mods[[var.n]],
            identical,
            obj$call.info$Call_ARGS[[var.n]])
        )
      )) {
        stop(paste0("All provided argument values for sensitivity analysis",
                    " should be new.","\nAt least one value of the '", var.n,
                    "' arg. is the same as in the basal hespdiv object \n"))
      }
      if(any(duplicated(hes.mods[[var.n]]))){
        stop(paste0("There are duplicated values in provided '", var.n,
                    "' argument. Remove duplicated values"))
      }
    }
    if (!var.n %in% c("data", "xy.dat", "study.pol")){
      if(any(unlist(lapply(hes.mods[[var.n]], length)) != 1)){
        stop(paste0(
          "Some of '", var.n, "' elements have length >1, but should",
          " be atomic."))
      }

    }


  }

  base <- paste0("hespdiv(data = obj$call.info$Call_ARGS$data,
                n.split.pts = obj$call.info$Call_ARGS$n.split.pts,
                generalize.f = obj$call.info$Call_ARGS$generalize.f,
                maximize = obj$call.info$Call_ARGS$maximize,
                method = obj$call.info$Call_ARGS$method,
                same.n.split = obj$call.info$Call_ARGS$same.n.split,
                compare.f = obj$call.info$Call_ARGS$compare.f,
                N.crit = obj$call.info$Call_ARGS$N.crit,
                N.rel.crit = obj$call.info$Call_ARGS$N.rel.crit,
                N.loc.crit = obj$call.info$Call_ARGS$N.loc.crit,
                N.loc.rel.crit = obj$call.info$Call_ARGS$N.loc.rel.crit,
                S.crit = obj$call.info$Call_ARGS$S.crit,
                S.rel.crit = obj$call.info$Call_ARGS$S.rel.crit,
                Q.crit = obj$call.info$Call_ARGS$Q.crit,
                c.splits = obj$call.info$Call_ARGS$c.splits,
                c.Q.crit = obj$call.info$Call_ARGS$c.Q.crit,
                c.crit.improv = obj$call.info$Call_ARGS$c.crit.improv,
                c.X.knots = obj$call.info$Call_ARGS$c.X.knots,
                c.Y.knots = obj$call.info$Call_ARGS$c.Y.knots,
                xy.dat = obj$call.info$Call_ARGS$xy.dat,
                c.max.iter.no = obj$call.info$Call_ARGS$c.max.iter.no,
                c.fast.optim = obj$call.info$Call_ARGS$c.fast.optim,
                c.corr.term = obj$call.info$Call_ARGS$c.corr.term,
                study.pol = obj$call.info$Call_ARGS$study.pol,
                use.chull = obj$call.info$Call_ARGS$use.chull,
                tracing = NULL,
                pnts.col = ", pnts.col,",
                display = ", display,",
                pacific.region = obj$call.info$Call_ARGS$pacific.region)")

  if (any(c("data","xy.dat", "study.pol", "compare.f") %in% hes.names)){
    if ("compare.f" %in% hes.names){

      methods_ids <- which(hes.names %in% c("compare.f","generalize.f",
                                            "maximize"))
      name.list <- hes.mods
      name.list$method <- c(name.list$method,
                            paste0("Custom_method_",
                                   1:length(hes.mods[["compare.f"]])))
      name.list <- name.list[-methods_ids]
    } else {
      name.list <- hes.mods
    }
    if (any(c("data","xy.dat") %in% hes.names)) {
      if (paired) {
        name.list$pdataxy <- paste0("Dataset_",
                                    1:length(hes.mods[["data"]]))
        name.list <- name.list[-which(names(name.list) %in% c("data","xy.dat"))]
        if ("study.pol" %in% hes.names){
          name.list[["study.pol"]] <- as.list(paste0("study.pol_",
                                                     1:length(hes.mods[["study.pol"]])))
        }
      } else {
        df_ids <- which(hes.names %in% c("data","xy.dat","study.pol"))
        for (i in df_ids){
          name.part <- names(hes.mods[i])
          name.list[[i]] <- as.list(paste0(name.part,"_",
                                           1:length(hes.mods[[i]])))
        }
      }
    }
    vec.par.vals <- lapply(name.list,unlist)
  } else {
    vec.par.vals <- lapply(hes.mods,unlist)
  }


  if (!comb.args){
    model.names <- unlist(lapply(1:length(vec.par.vals),
                                 function(i) paste0(names(vec.par.vals[i]),
                                                    " = ",
                                                    vec.par.vals[[i]] )))
  } else {
    N <- ifelse(!"compare.f" %in% hes.names,length(hes.names),ifelse(
      !"method" %in% hes.names,length(hes.names)-2,length(hes.names)-3))
    N <- ifelse(paired, N-1, N)
    argum <- lapply(1:N, function(i) paste0(names(vec.par.vals[i]), " = ",
                                            vec.par.vals[[i]] ))
    if (is.null(pick.n.args)){
      pick.n.args <- 1:N
    } else {
      if (max(pick.n.args) > N){
        stop(paste0("Some provided values of 'pick.n.args' are higher then the",
                    " possible maximum: ",N))
      }
    }
    arg.combs <- vector( mode = "list", length = length(pick.n.args))
    for ( i in seq(length(pick.n.args))){
      arg.combs[[i]] <- Reduce(rbind,
                               lapply(combn(1:N, m = sort(pick.n.args)[i],
                                            simplify = FALSE),
                                      function(id) expand.grid(argum[id])))
    }
    model.names <- unlist(lapply(arg.combs, apply, 1, paste, collapse = " & "))
    if (comb.type == "random") {
      if (n.combs >= length(model.names)){
        message(paste0(
          "The required number of argument"," combinations ('n.combs' = ",
          n.combs,") \nis higher than the number of available argument",
          " combinations (",length(model.names),").\n",
          "Do you wish to procced and test all the available combinations?",
          "\n Choose: Yes or No"))
        answer <- readline(prompt = "")
        answer <- .answer_check(given = answer,NAMES = c("Yes", "No"))
        if (answer == 'No'){
          return(NULL)
        } else {
          comb.type <- "all"
          n.combs <- NULL
        }
      } else {
        model.names <- sample(model.names, n.combs, replace = FALSE)
      }
    } else {
      if (comb.type == "handpicked"){
        print(model.names)
        message(paste0(
          "\nFrom the vector printed above choose argument combinations ",
          "you would like to test. \nThen type a call that creates a vector of",
          " their indeces (i.e. 'c(43:55, 67:68, 101)')"))
        cond <- TRUE
        while(cond) {
          modinds <- readline(prompt = "")
          modinds <- tryCatch({eval(parse(text = modinds))},
                              error = function(cond) {
                                message(cond)
                                message("\nRetype the call.")
                                return(NULL)
                              },
                              warning = function(cond) {
                                message(cond)
                                message("\nRetype the call.")
                                return(NULL)
                              })
          if (!is.null(modinds)){
            if (is.numeric(modinds) & length(modinds)>0 & is.vector(modinds) &
                !is.list(modinds) & all(!is.na(modinds)) &
                all(!is.nan(modinds))){
              selection <- tryCatch({model.names[modinds]},
                                    error = function(cond) {
                                      message("Use of provided indeces produced error.")
                                      message(cond)
                                      message("\nRetype the call.")
                                      return(NULL)},
                                    warning = function(cond) {
                                      message("Use of provided indeces issued a warning.")
                                      message(cond)
                                      message("\nRetype the call.")
                                      return(NULL)
                                    }
              )
              if (!is.null(selection)){
                if (any(is.nan(selection) | is.na(selection))) {
                  message("Use of provided indeces produced NA or NaN.")
                  print(selection)
                  message("Retype the call.")
                } else {
                  cond <- FALSE
                }
              }
            } else {
              message("Call is invalid. The output is not a numeric integer ",
                      "vector. Its values is: ")
              print(modinds)
              message("Retype the call.")
            }
          }
        }

        model.names <- selection

      }
    }
  }
  model.names <- .clean_names(model.names, obj, hes.names, c.splits, c.pars)
  message("Changes to be made in the basal hespdiv call: ")
  print(model.names)

  l <- length(model.names)
  mods <- vector(mode = "list", length = l)
  names(mods) <- model.names

  for (mod.id in 1:l){
    message(paste0("Calling hespdiv alternative: [",mod.id,"] ",
                   model.names[mod.id]))
    filt.names <- unlist(
      strsplit(gsub("= | &", "",  model.names[mod.id]), split = " ",
               perl = TRUE))
    f.id <- 1:length(filt.names)
    arg.names <- filt.names[f.id %% 2 == 1]
    arg.vals <- filt.names[f.id %% 2 == 0]
    mod.base <- base
    if (any(grepl("Custom_method_.*", arg.vals, fixed = F))){
      del.id <- which(grepl("Custom_method_", arg.vals, fixed = TRUE))
      vals.id <- as.numeric(strsplit(arg.vals[del.id], split = "_")[[1]][3])
      arg.names <- c(arg.names[-del.id],"compare.f","generalize.f",
                     "maximize")
      arg.vals <- c(arg.vals[-del.id], rep(vals.id,3))
      mod.base <- sub('obj\\$call.info\\$Call_ARGS\\$method',"NULL",
                      mod.base)
    } else {
      if (any(arg.names == "method")) {
        args <- c("compare.f","generalize.f", "maximize")
        for (arg in args)
          mod.base <- sub(paste0('obj\\$call.info\\$Call_ARGS\\$',arg),
                          "NULL",
                          mod.base)
      }
    }
    if ("pdataxy" %in% arg.names) {
      del.id <- which(arg.names == "pdataxy")
      vals.id <- as.numeric(strsplit(arg.vals[del.id], split = "_")[[1]][2])
      arg.names <- c(arg.names[-del.id],"data","xy.dat")
      arg.vals <- c(arg.vals[-del.id], rep(vals.id,2))
    }
    for (arg.id in 1:length(arg.names)){
      if (arg.names[arg.id] == "study.pol" |
          (arg.names[arg.id] %in%  c("data", "xy.dat") & !paired ) ){
        val.id <- as.numeric(strsplit(arg.vals[arg.id], split = "_")[[1]][2])
        arg.name <- arg.names[arg.id]
      } else {
        arg.name <- arg.names[arg.id]
        val.id <- ifelse(!arg.name %in% c("compare.f","generalize.f",
                                          "maximize", "data", "xy.dat"),
                         yes =
                           which(as.character(unlist(hes.mods[arg.name])) ==
                                   arg.vals[arg.id]),
                         no = as.numeric(arg.vals[arg.id]))
      }
      mod.base <- sub(paste0('obj\\$call.info\\$Call_ARGS\\$',
                             arg.name),
                      paste0("hes.mods[['",arg.name,"']][[",
                             val.id,"]]"),
                      mod.base)
    }
    if (!is.null(images.path))
      png(paste0(images.path,"\\",paste0(gsub(pattern = "&", replacement = "_",
                                              x = gsub(pattern = " ",
                                                       replacement = "",
                                                       x =  model.names[mod.id]))
                                         ,"_",mod.id,".png")))
    mods[[mod.id]] <- list(Subdivision = tryCatch({eval(parse(text = mod.base))},
                                                  error = function(cond) {
                                                    base::message(paste0(cond,"\n"))
                                                    return(cond)
                                                  },
                                                  warning = function(cond) {
                                                    base::message(paste0(cond,"\n"))
                                                    return(list(Subdivision =
                                                                  eval(parse(text = mod.base)),
                                                                Warning = cond))
                                                  }))#, Changes = model.names[mod.id]) removed so that hsa and hsa_detailed would be similar
    if (!is.null(images.path))
      dev.off()
  }
  structure(list(Alternatives = mods, Basis = obj), class = "hsa")
}
#' @noRd
.to.list <- function(x){
  if(!is.null(x)){
    if (is.list(x)){
      x <- unlist(x)
    }
    if (length(x) > 1){
      stop("Alternative Boolean hespdiv arguments should only contain one value")
    }
    as.list(x)
  }
}
#' @noRd
.clean_names <- function(model.names, obj, hes.names, c.splits, c.pars){
  if("study.pol" %in% hes.names){
    ids <- which(
      (grepl("study.pol", model.names, fixed = TRUE) &
         grepl("use.chull = TRUE", model.names, fixed = TRUE)) |
        (grepl("study.pol", model.names, fixed = TRUE) &
           obj$call.info$Call_ARGS$use.chull &
           !grepl("use.chull = FALSE", model.names, fixed = TRUE))
    )
    if (length(ids) > 0)
      model.names <- model.names[-ids]
  }
  if (c.pars & !is.null(c.splits)){
    if (c.splits[[1]]){
      ids <- which(
        grepl("c.splits", model.names, fixed = TRUE) |
          !grepl("^c\\.| c\\.", model.names)
      )
      if (length(ids) > 0)
        model.names <- model.names[ids]
    } else {
      ids <- which(
        (!grepl("c.splits", model.names, fixed = TRUE) &
           grepl("^c\\.| c\\.", model.names)) |
          (grepl("c.splits", model.names, fixed = TRUE) &
             !grepl("^c\\.| c\\.",
                    sub(pattern = "c.splits", replacement = "",
                        x = model.names))) |
          !grepl("^c\\.| c\\.", model.names)
      )
      if (length(ids) > 0)
        model.names <- model.names[ids]
    }
  }
  model.names
}

#' @noRd
.answer_check <- function(given,NAMES){
  matched.i <- pmatch(tolower(given), tolower(NAMES))
  while(is.na(matched.i)){
    cat("Invalid input: ", paste0('"', given,'".'),
        paste0("\nPlease select viable option: "),
        paste(NAMES,collapse = " or ",sep = "'"))
    given <- readline(prompt = "")
    matched.i <- pmatch(tolower(given), tolower(NAMES))
  }
  NAMES[matched.i]
}
