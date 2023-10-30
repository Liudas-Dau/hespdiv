#' Quantify the stability of hespdiv clusters
#'
#' @description This function evaluates the stability of the basal subdivision
#' clusters using hespdiv sensitivity analysis results. It does so by
#' calculating Jaccard similarities between the observations of basal
#' subdivision clusters and the observations of alternative subdivision
#' clusters. For each basal cluster, the function identifies the most similar
#' ('analog') cluster within each alternative subdivision. The stability of
#' each basal cluster can be assessed by examining the distribution of
#' similarity values with their corresponding 'analog' clusters. If a highly
#' similar cluster reappears in multiple alternative subdivisions, it
#' indicates that the basal cluster is stable.
#' @param obj An object of class \code{hsa}.
#' @param probs A numeric vector of probabilities with values in the range
#' [0, 1]. This argument is used to calculate quantiles of Jaccard similarity
#' values.
#' @return The function returns a list containing three data frames, each providing specific information for each basal hespdiv cluster:
#'  \describe{
#'  \enumerate{
#'  \item{\bold{jaccard.quantiles:} This data frame contains the quantiles of Jaccard similarities between the basal cluster and the 'analog' clusters from alternative subdivisions.}
#'  \item{\bold{jaccard.similarity:} This data frame provides the Jaccard similarity values between the basal cluster and the 'analog' cluster from each alternative subdivision.}
#'  \item{\bold{analog.clusters:} This data frame includes the IDs of the hespdiv polygons that produced the 'analog' clusters in each alternative subdivision.}
#'  }}
#' @details
#' If a basal subdivision cluster obtains a distribution of high similarity values,
#' then it is considered a 'stable' and 'existing' cluster. On the other hand,
#' low 'analog' cluster similarity values may signal that a basal cluster is an
#' artifact of hespdiv computation.
#'
#' The more technical description of how \code{hsa_quant} works:
#' \describe{
#' \item{Obtaining alternative hespdiv clusters:}{ The function filters the \code{xy.dat} coordinates of the basal subdivision using all the polygons of alternative subdivisions, obtaining alternative hespdiv clusters.}
#' \item{Quantifying Jaccard similarity:}{ The function measures the Jaccard overlap index between the observations of the basal subdivision clusters and the observations of the alternative clusters. }
#' \item{Identification of 'analog' clusters and value assingments:}{ Each basal hespdiv cluster from each alternative subdivision is assigned the ID of the cluster that produced the maximum Jaccard similarity value, along with the corresponding similarity value.}
#' }
#' The purpose of the \code{hsa_quant} function is to address situations where hespdiv
#' polygons, despite having different geometry and location, may filter nearly
#' identical sets of observations, leading to similar hespdiv clusters. This can occur when the spatial
#' coverage of observations is incomplete and irregular, or when the boundaries
#' between hespdiv polygons are expected to be open, soft, or fuzzy, such as in
#' the case of boundaries between bioregions. In such cases, visual hespdiv
#' sensitivity analysis alone may show irregular and non-converging
#' distributions of split-lines. However, \code{hsa_quant} can reveal that these
#' irregular polygons are based on nearly identical clusters of observations,
#' indicating a strong spatial structure within the analyzed data. Conversely,
#' if the observations within these clusters significantly differ, it indicates
#' that the basal clusters are specific to the hespdiv parameters used and
#' likely lack ontological meaning.
#'
#' Thus, by analyzing the similarity values between clusters of observations, \code{hsa_quant}
#' facilitates the assessment of the stability and reliability of basal
#' subdivision clusters, aiding in evaluating their significance.
#'
#' @note You can use the hsa_quant function to track the evolution of hespdiv
#' subdivisions over time by providing correctly formatted input. For
#' instance, you can obtain the basal subdivision for time bin 1 using the
#' hespdiv function. Then, using the hsa function, you can specify the paired
#' xy.dat and data from time bin 2. The resulting hsa object can be inputted
#' into \code{hsa_quant} The \code{hsa_quant} result will then provide insights into
#' extinctions, speciations, fusions, and splits of hespdiv polygons/clusters
#' that occur between time bin 1 and 2. This allows for the analysis of changes
#' and dynamics in hespdiv subdivisions over time.
#' 1 and 2.
#' @family {functions for hespdiv sensitivity analysis}
#' @family {functions to evaluate hesdpiv cluster stability}
#' @author Liudas Daumantas
#' @importFrom stats quantile density
#' @importFrom graphics plot hist
#' @export

hsa_quant <- function(obj, probs = c(0.05,0.5,0.95)){
  if (!inherits(obj,"hsa"))
    stop("'obj' should be of class \"hsa\" (output of 'hsa' or 'hsa_detailed' functions).")

  subs <- lapply(obj$Alternatives, function(o) o[[1]])
  # Check if all alternative subdivisions are NULL
  if (all(il <- sapply(subs, is.null) | sapply(subs,class) != "hespdiv"))
    stop("All alternative subdivisions are NULL or contained errors")
  # Remove NULL subdivisions or subdivisions with
  if (any(il))
    subs <- subs[-which(il)]
  connections <- matrix(NA, nrow = length(subs),ncol = length(obj$Basis$polygons.xy) -1 )
  colnames(connections) <- names(obj$Basis$polygons.xy[-1])
  rownames(connections) <- names(subs)
  jac.sim <- matrix(NA, nrow = length(subs),ncol = length(obj$Basis$polygons.xy) -1 )
  colnames(jac.sim) <- names(obj$Basis$polygons.xy[-1])
  rownames(jac.sim) <- names(subs)

    for (basis_pol in 2:length(obj$Basis$polygons.xy)){
      b.ids <- .get_ids(obj$Basis$polygons.xy[[basis_pol]],
                                  obj$Basis$call.info$Call_ARGS$xy.dat)
      l <- lapply(subs, function(o, bas.ids) {
        lapply(o$polygons.xy[-1],function(e, dat, bs.ids) {
          .jaccard.sim(.get_ids(e, dat), bs.ids)
        }, dat = obj$Basis$call.info$Call_ARGS$xy.dat, bs.ids = bas.ids)
      }, b.ids)

      jac.sim[,basis_pol-1] <- unlist(lapply(l,function(o) max(unlist(o))))
      connections[,basis_pol-1] <- unlist(lapply(l,function(o) as.numeric(names(
        which.max(unlist(o))))))
    }
    # code to be considered, if fusion/split of polygons would be disallowed
    # (e.g. alternative subdivision has cluster that is the best analog
    # for more than 1 basal subdivision cluster) :
    # bas.ids <- sapply(de$Basis$polygons.xy[-1], function(o,ids){
    #   .get_ids(o, de$Basis$call.info$Call_ARGS$xy.dat)})
    # for (alt.subs in 2:length(de)){
    #   conn <- matrix(NA, nrow = length(de[[alt.subs]]$polygons.xy)-1,
    #                  ncol = length(de$Basis$polygons.xy) -1 )
    #   colnames(conn) <- names(de$Basis$polygons.xy[-1])
    #   rownames(conn) <- names(de[[alt.subs]]$polygons.xy[-1])
    #   jac.sim <- matrix(NA, nrow = length(de)-1,ncol = length(de$Basis$polygons.xy) -1 )
    #   colnames(jac.sim) <- names(de$Basis$polygons.xy[-1])
    #   rownames(jac.sim) <- names(de[-1])
    #   for (alt.pols in 2:length(de[[alt.subs]]$polygons.xy)){
    #     alt.ids <- .get_ids(de[[alt.subs]]$polygons.xy[[
    #       alt.pols]], de$Basis$call.info$Call_ARGS$xy.dat)
    #     conn[alt.pols-1,] <- sapply(bas.ids, function(o,ids){
    #       .jaccard.sim(o, ids)
    #     },alt.ids)
    #   }
    #   if (any(dup.id <- duplicated(apply(conn,2,which.max)))){
    #     dup.pols <- unique(apply(conn,2,which.max)[dup.id])
    #     for (dup.pol in 1:length(dup.pols)){
    #
    #     }
    #   }
    # }

  structure(list(jaccard.quantiles = apply(jac.sim,2, stats::quantile,
                                 probs = probs),
    jaccard.similarity = jac.sim, analog.clusters = connections),
    class = "hsa_quant")
}

#' @noRd
.jaccard.sim <- function(x,y){
  inter <- length(intersect(x,y))
  inter / (length(x) + length(y) - inter)
}


