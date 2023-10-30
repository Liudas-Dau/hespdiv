#' Print the results of hespdiv object
#'
#' @method print hespdiv
#' @description  Function formats and prints the rounded split.stats data frame.
#'  from "hespdiv" class R object.
#' @param x A hespdiv class object.
#' @param ... other arguments
#' @return x
#' @author Liudas Daumantas
#' @export
print.hespdiv <- function(x, ...){
  dat <- x$split.stats
  if (!inherits(x,"hespdiv"))
    stop("x should have 'hespdiv' class.")
  if (x$call.info$METHOD$method.type == "custom"){
    type <- "Custom"
    cat(paste0("\n","Information about the split-lines:", "\n\n",type,
               " method was used.\n\n"))
  } else {
    type <- "Preset"
    if (x$call.info$METHOD$method == "biozonation"){
      method <- "Spatial biozonation"
      if (x$call.info$METHOD$metric == "sorensen") {
        metric <- paste0("So",
        "rensen-Dice coefficient")
        variant <- x$call.info$METHOD$variant
        if (variant == '1'){
          variant <- "ref. Sorensen, T. A. (1948)"
        }
      } else {
        if (x$call.info$METHOD$metric == "morisita"){
          metric <- paste0("Morisita Overlap index")
          variant <- x$call.info$METHOD$variant
          if (variant == '1'){
            variant <- "ref. Morisita, M. (1959)"
          }
        } else {
          if (x$call.info$METHOD$metric == "pielou"){
            metric <- paste0("Proportion of mean reduction in Pielou entropy")
            variant <- x$call.info$METHOD$variant
            if (variant == '1'){
              variant <- paste0("(1 - \u0394Pe / Pe.base)(",
              "ref. in hespdiv documentation)")
            }
          } else {
            if (x$call.info$METHOD$metric == "horn.morisita"){
              metric <- paste0("Morisita Overlap index")
              variant <- x$call.info$METHOD$variant
              if (variant == '2'){
                variant <- "Horn modification, ref. Horn (1966)"
              }
            }
          }
        }
      }
    }

  if (x$call.info$METHOD$similarity){
    similarity <- "Similarity"
  } else {
    similarity <- "Distance"
  }
  if (variant == "" | is.null(variant)){
    ending <- ""
  } else {
    ending <- paste0(", ",variant)
  }

  metric.row <- paste0("Metric: ", metric," (",
                       similarity,")", ending, ".")
  nm<- nchar(metric.row)

  cat(paste0("\n","Information about the split-lines:", "\n\n",type,
             " method was used.\n",
             "Method: ",method,
             ".\n", metric.row,"\n",paste0(rep("-",nm),collapse = ""),"|\n"))
  }
  if (is.null(dat)){
    cat("No Split-lines were found.", "\n")
  } else {
    id <- which(colnames(dat) == "is.curve")
    dat[,-id] <- round(dat[,-id],2)
    print(dat)
    if (any(names(dat)=="n.m.rez")){
      cat("\n", "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    }
  }
  invisible(x)
}
