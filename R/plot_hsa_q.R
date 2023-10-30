#' Plot the stability of hespdiv clusters
#'
#' @description  This function visualizes the stability of basal subdivision
#' clusters obtained from \code{hsa_quant}
#' @param obj The output of \code{hsa_quant} (an object of class \code{hsa_quant}).
#' @param hist A Boolean value. If FALSE, EPDFs obtained with \code{plot(density(jac.sim))}
#' will be displayed instead of histograms.
#' @return None
#' @details The stability of each basal cluster is revealed by the distribution
#' of Jaccard similarity values with the 'analogues' clusters found in alternative
#' hespdiv subdivisions. For example, a unimodal distribution with a peak at
#' high similarity values (>0.8) indicates that the basal hespdiv cluster is
#' stable, even if the polygon boundaries are not. This situation may arise when
#' there is indeed a spatial structure within the data, but there are also wide
#' gaps between sampled regions (or more generally when there is limited spatial
#' data coverage). A unimodal distribution with a peak at medium values
#' (0.4-0.6) and a tail to higher values could also indicate a more persistent
#' spatial structure. On the other hand, a single peak at low values (<0.4)
#' indicates low cluster stability (e.g., bioregion does not exist). Finally,
#' uniform, bimodal, or other more complex distributions may indicate that the
#' stability and existence of the corresponding basal cluster depend on the
#' parameters used in alternative hespdiv calls.
#' @family {functions for hespdiv sensitivity analysis}
#' @family {HespDiv visualization options}
#' @family {functions to evaluate hesdpiv cluster stability}
#' @author Liudas Daumantas
#' @importFrom graphics plot hist layout
#' @importFrom stats density
#' @export
plot_hsa_q <- function(obj, hist = FALSE){
  if (!inherits(obj,"hsa_quant"))
    stop("'obj' should be of class \"hsa_quant\" (output of 'hsa_quant' function).")

  jac.sim <- obj[[2]]
  l <- ncol(jac.sim)
  if (!hist){
    n <- 5
    while(!l%%n == 0){
      n <- n-1
    }
    graphics::layout(1:n)
    for (i in 1:l)
      graphics::plot(stats::density(jac.sim[,i]),
                     main = paste0("Stability of the HespDiv cluster - ",
                                   colnames(jac.sim)[i]), xlim= c(0,1),
                     xlab = "Jaccard similarity with the 'analogue' clusters")
    graphics::layout(1)
  } else {
    n <- 5
    while(!l%%n == 0){
      n <- n-1
    }
    graphics::layout(1:n)
    for (i in 1:l)
      graphics::hist(jac.sim[,i],
                     main = paste0("Stability of the HespDiv cluster - ",
                                   colnames(jac.sim)[i]),xlim= c(0,1),
                     xlab = "Jaccard similarity with the 'analogue' clusters")
    graphics::layout(1)
  }
}
