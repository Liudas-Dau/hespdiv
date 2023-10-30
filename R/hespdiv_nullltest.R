#' Test significance of spatial structure in data
#'
#' @description  Function n times shuffles the data and each time
#' checks what is the performance of each split-line dividing the new data subsets.
#'
#' @return A "nullhespdiv" class object. It is a list containing two elements:
#' \itemize{
#' \item{1)}{ A data.frame which reports the obtained mean and standard deviation of split-line
#' performances after data shuffling. Also, it compares these performances with
#' the ones estimated with hespdiv analysis by providing a quantiles and z-scores.}
#' \item{2)}{ A data frame with all obtained comparison values for each split-line after data shuffling.}
#' }
#' @param obj An object of class \code{hsa}
#' @param n Number of permutations
#' @export
nulltest <- function(obj, n){
  if (!inherits(obj,"hespdiv"))
    stop("obj should have 'hespdiv' class.")
  coords <- obj$call.info$Call_ARGS$xy.dat
  data <- obj$call.info$Call_ARGS$data
  N <- nrow(coords)
  l <- length(obj$split.lines)
  comp.vals <- vector(mode = "list", length = n)
  p.vals <- vector(mode = "list", length = n)
  if (obj$call.info$Call_ARGS$maximize){
    pal <- function(x,criteria){ x > criteria}
  } else {
    pal <- function(x,criteria){ x < criteria}
  }
  for (i in 1:n){
    comp.vals[[i]] <- numeric(l)
    p.vals[[i]] <- numeric(l)
  }
  if (is.data.frame(data) | is.matrix(data)){
    .slicer <- .slicer.table
  } else {
    if (is.list(data)) {
      .slicer <- .slicer.list
    } else {
        .slicer <- .slicer.vect
    }
  }
  for ( i in 1:n) {
    ids <- sample(1:N,N,replace = FALSE)
    co <- coords[ids,]
    dat <- .slicer(data, ids)

    for (split.id in 1 : l){
      pol_ids <- which(obj$poly.stats$root.id == obj$split.stats$plot.id[split.id])

      split.ids1 <- .get_ids(obj$polygons.xy[[pol_ids[1]]], co)
      split.ids2 <- .get_ids(obj$polygons.xy[[pol_ids[2]]], co)
      dat_pol1 <- .slicer(data, split.ids1)
      dat_pol2 <- .slicer(data, split.ids2)
      comp.vals[[i]][split.id] <- obj$call.info$Call_ARGS$compare.f(
        obj$call.info$Call_ARGS$generalize.f(dat_pol1),
        obj$call.info$Call_ARGS$generalize.f(dat_pol2))
      p.vals[[i]][split.id] <- pal(comp.vals[[i]][split.id],
                                     obj$split.stats$performance[split.id])
    }

  }

  p.vals <- Reduce("+", p.vals)/n
  comp.vals <- t(as.data.frame(comp.vals))
  rownames(comp.vals) <- 1:n
  stats <- data.frame(ID = 1:l, performance = obj$split.stats$performance,
                      mean.random = apply(comp.vals,2,mean),
                      sd.random = apply(comp.vals,2,sd),
                      quantile = p.vals)
  stats$z.score.random <- (stats$performance - stats$mean.random)/stats$sd
  return(print.nullhespdiv(structure(list(stats, comp.vals), class = "nullhespdiv")))

}


