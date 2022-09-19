#' Calculate polygon cross-comparison matrix
#'
#' This function calculates cross-comparison matrix of polygons identified by
#' hespdiv. This matrix can be transformed into distance matrix to be used in
#' cluster analysis.
#' @param obj hespdiv object
#' @author Liudas Daumantas
#' @export

cross_comp <- function(obj){
  if (obj$call.info$METHOD$method.type == "preset") {
    method <- obj$call.info$METHOD$metric
    if (method == "pielou") {
      compare.f <- function(eveness1,eveness2) {
        base.eveness <- poly.obj[[testid]]
        (1 - mean(c(eveness1, eveness2)) / base.eveness) * 100 # percent change in eveness
      }

    } else {
      if (method == "sorensen"){
        compare.f <- function(uniq_tax1,uniq_tax2) {
          sum <- length(uniq_tax1) + length(uniq_tax2)
          int_2x <- length(which(duplicated(c(uniq_tax1,uniq_tax2))))*2
          if (length(int_2x)!=0){
            int_2x/sum
          } else {
            0
          }
        }
      } else {
        if (method == "morisita"){
          compare.f <- function(x,y) {
            all_sp <- unique(c(x,y))
            x_f <- factor(x,levels = all_sp)
            y_f <- factor(y,levels = all_sp)
            (2*sum(table(x_f) * table(y_f)))/
              (length(x) * length(y) *
                 ((sum(table(x_f)*(table(x_f)-1)) /
                     (length(x)* (length(x)-1))) +
                    (sum(table(y_f)*(table(y_f)-1)) /
                       (length(y)* (length(y)-1)))))
          }
        } else {
          if (method == "horn.morisita"){
            compare.f <- function(x,y) {
              all_sp <- unique(c(x,y))
              x_f <- factor(x,levels = all_sp)
              y_f <- factor(y,levels = all_sp)
              (2*sum(table(x_f) * table(y_f))) /
                (length(x) * length(y) *
                   ( sum(table(x_f)^2) / (length(x)^2)  +
                       sum(table(y_f)^2) / (length(y)^2)))
            }
          }
        }
      }
    }
  } else {
    compare.f <- obj$call.info$Call_ARGS$compare.f
  }
  comp.mat <- matrix(NA,nrow = length(obj$poly.obj),ncol =length(obj$poly.obj) )
  for (pol.id in seq(length(obj$poly.obj))){
    comp.mat[,pol.id] <- unlist(lapply(obj$poly.obj,FUN = compare.f,
                                obj$poly.obj[[pol.id]]))
  }
  comp.mat
}

