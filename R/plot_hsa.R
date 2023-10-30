#' Visualize the stability of hespdiv polygons
#'
#' @description  The function uses hespdiv sensitivity analysis results to visually
#' demonstrate the stability of the basal hespdiv subdivision. This is achieved
#' by displaying both alternative and basal hespdiv subdivisions on the same plot.
#' @param obj An object of class \code{hsa}
#' @param alpha The alpha value for transparency of split lines. Default is 0.6.
#' @param split.col Color of alternative subdivision split-lines. Default is "gray20".
#' @param pnts.col The color of data points. Default is NULL.
#' @param pol.col The color of polygons. Default is "7".
#' @param type An integer indicating the type of plot. Default is 1.
#' @param max.lwd The maximum line width for split-lines. Default is 3.
#' @param min.lwd The minimum line width for split-lines. Default is 0.5.
#' @param basal.col The color of basal subdivision split-lines.
#' @param split.col.seed A seed for generating random colors for split lines. Default is NULL.
#' @param seperated Boolean. When \code{type} is >= 3, open a new graphical device for each rank?
#' @param newplot Create a plot in new device?
#' @return NULL
#' @details The type parameter determines the type of plot to be generated:
#' \describe{
#' \item{type = 1:}{
#' Basic plot - Displays the alternative and basal hespdiv subdivisions on the
#' same plot without split-line ranks or titles.}
#' \item{type = 2:}{
#'  Plot with split-line ranks - Includes split-lines ranks in the plot.
#'  Each split-line is assigned a different line width based on its rank.}
#'\item{type = 3:}{
#'  Plot with separate ranks - Generates multiple plots, each representing
#'  split-line ranks up to a certain value.}
#'\item{type = 4:}{
#'  Plot with separate and isolated ranks - Similar to mode 3 but isolates the
#'  split-line ranks. Generates multiple plots, each representing a
#'  specific split-line rank.}
#'}
#' @author Liudas Daumantas
#' @details If the alternative subdivisions spatially converge, but the basal
#' subdivision lies far from the zone of convergence, then you can use
#' the \code{change_base} function to select a more representative alternative
#' subdivision to be used as the basal subdivision. However, you should check
#' if the arguments used in that subdivision makes sense.
#' @family {functions for hespdiv sensitivity analysis}
#' @family {HespDiv visualization options}
#' @note 'newplot' allows to correctly render legend in types 2 and 3, and lines
#' in general when drawing in an active device (use 'broom' otherwise to delete
#' devices).
#' @importFrom graphics plot lines points
#' @importFrom scales alpha
#' @importFrom grDevices dev.new
#' @export
plot_hsa <- function(obj, alpha = 0.6, split.col = "gray20", pnts.col = NULL,
                     pol.col = "7", type = 1, basal.col = 2,
                     max.lwd = 3, min.lwd = 0.5, split.col.seed = NULL,
                     newplot = TRUE, seperated = TRUE){
  if (!inherits(obj,"hsa"))
    stop("'obj' should be of class \"hsa\" (output of 'hsa' or 'hsa_detailed' functions).")
  if ( length(pol.col) > 1){
    if ( length(pol.col) != length(subs)){
      stop(paste0("Length of colors for study area polygons is not equal to",
                  "the length of subdivisions"))
    }}
  if (newplot & !seperated & type == 3) {
    stop("Change seperated to TRUE or newplot to FALSE to make all plots visible.")
  }
  # Combine the basal and alternative subdivisions
  subs <- c(list(obj$Basis),lapply(obj$Alternatives, function(o) o[[1]]))
  # Check if all alternative subdivisions are NULL
  if (all(il <- sapply(subs[-1], is.null) | sapply(subs[-1],class) != "hespdiv"))
    stop("All alternative subdivisions are NULL or contained errors/warnings")
  # Remove NULL subdivisions or subdivisions with
  if (any(il))
    subs <- subs[-(which(il) + 1)]
  # Determine the plot settings based on the type parameter
  if (type == 1){
    ranks <- separate <- FALSE
    title <- ""
  } else {
    if (newplot) {
      dev.new()
      if (type == 2 | type == 3){
        op <- graphics::par(no.readonly = TRUE)
        graphics::par(mar=c(5, 4, 4, 8)+0.1, xpd = TRUE)
      }
    }
    if (type == 2){
      ranks <- TRUE
      separate <- FALSE
      title <- ""
    } else {
      if (type == 3){
        ranks <- separate <- TRUE
        isolate.ranks <- FALSE
      } else {
        if (type == 4){
          ranks <- separate <- isolate.ranks <- TRUE
        }
      }
    }
  }
  # Plot the initial plots if not separate
  if (!separate){
    .initial_plots(subs = subs, pnts.col = pnts.col,max.lwd, title, basal.col)
  }
 # Plot the split lines
  if (!ranks) {
    for (i in 1:length(subs)){
      for (split in 1:length(subs[[i]][["split.lines"]]))
        graphics::lines(x=subs[[i]][["split.lines"]][[split]],
                        col=scales::alpha(split.col[1],alpha),
                        lwd=0.5)
    }
  } else {
    all_ranks <- range(sapply(subs,function(o) range(o[["split.stats"]]$rank)))
    if (is.null(split.col.seed)){
    if (length(split.col) < all_ranks[2] ){
      if (length(split.col) != 1)
        warning(paste0("More split ranks than provided split.col.",
                     "\nGenerating random colors"))
      split.col <- .generate_cols(all_ranks[2], sample(1:9999,1))
    }} else {
      split.col <- .generate_cols(all_ranks[2], split.col.seed)
    }
    if (separate){
      if (!isolate.ranks){
      for (rangas in all_ranks[1]:all_ranks[2]){
        if (seperated) {
          if (rangas !=1)
            grDevices::dev.new()
          if (!(rangas == 1 & newplot)){
          op <- graphics::par(no.readonly = TRUE)
          graphics::par(mar=c(5, 4, 4, 8)+0.1, xpd = TRUE)}
        }
        title <- paste0("Split-line Rank - up to ", rangas)
        .initial_plots(subs, pnts.col,max.lwd, title, basal.col)
        if (any(viz <- sapply(subs,function(o,rangas) any(
          o[["split.stats"]]$rank <= rangas), rangas ))){
          for (i in (1:length(subs))[viz]){
            for (split in 1:length(subs[[i]][["split.lines"]])){
              if (subs[[i]][["split.stats"]]$rank[split] <= rangas)
                graphics::lines(x=subs[[i]][["split.lines"]][[split]],
                                col=scales::alpha(split.col[
                                  subs[[i]][["split.stats"]]$rank[split]
                                ],alpha),
                                lwd=seq(max.lwd,min.lwd,length.out = all_ranks[2])[
                                  subs[[i]][["split.stats"]]$rank[split]])
            }
          }
        }
        graphics::legend("right", inset=c(-0.2,0),
                         legend=c(1:all_ranks[2]),
                         col = scales::alpha(split.col[1:all_ranks[2]],alpha),
                         lty = 1,
                         cex = 0.7,
                         lwd = seq(max.lwd,min.lwd,length.out = all_ranks[2]),
                         title="Rank")
        .final_plots(subs, pol.col, max.lwd)
      }} else {
        for (rangas in all_ranks[1]:all_ranks[2]){
          if (seperated)
            grDevices::dev.new()
          title <- paste0("Split-line Rank - ", rangas)
          .initial_plots(subs, pnts.col,max.lwd, title, basal.col)
          if (any(viz <- sapply(subs,function(o,rangas) any(
            o[["split.stats"]]$rank == rangas), rangas ))){
            for (i in (1:length(subs))[viz]){
              for (split in 1:length(subs[[i]][["split.lines"]])){
                if (subs[[i]][["split.stats"]]$rank[split] == rangas)
                  graphics::lines(x=subs[[i]][["split.lines"]][[split]],
                                  col=scales::alpha(split.col[
                                    subs[[i]][["split.stats"]]$rank[split]
                                  ],alpha),
                                  lwd=seq(max.lwd,min.lwd,length.out = all_ranks[2])[
                                    subs[[i]][["split.stats"]]$rank[split]])
              }
            }
          }
          .final_plots(subs, pol.col, max.lwd)
        }
      }
    } else{
      for (i in 1:length(subs)){
        for (split in 1:length(subs[[i]][["split.lines"]])){
          graphics::lines(subs[[i]][["split.lines"]][[split]],
                          col=scales::alpha(split.col[
                            subs[[i]][["split.stats"]]$rank[split]
                          ],alpha),
                          lwd=seq(max.lwd,min.lwd,length.out = all_ranks[2])[
                            subs[[i]][["split.stats"]]$rank[split]])
        }

      }
      graphics::legend("right", inset=c(-0.2,0),
             legend=c(1:all_ranks[2]),
             col = scales::alpha(split.col[1:all_ranks[2]],alpha),
             lty = 1,
             cex = 0.7,
             lwd = seq(max.lwd,min.lwd,length.out = all_ranks[2]),
             title="Rank")}
    }
  if (!separate){
    .final_plots(subs, pol.col, max.lwd)
  }
  if (type == 2 & newplot){
    graphics::par(op)
  }
}
#' @noRd
.initial_plots <- function(subs, pnts.col,max.lwd, title, basal.col){
  graphics::plot(NULL,
                 ylim = range(sapply(subs,
                                     function(o) range(c(o$polygons.xy[[1]]$y,
                                                         o$call.info$Call_ARGS$study.pol$y)))),
                 xlim =  range(sapply(subs,
                                      function(o) range(c(o$polygons.xy[[1]]$x,
                                                          o$call.info$Call_ARGS$study.pol$x)))),
                 col = 0,
                 xlab = "x coordinate", ylab = "y coordinate", main = title)
  if (!is.null(pnts.col)){
    graphics::points(subs[[1]]$call.info$Call_ARGS$xy.dat,col=pnts.col,
                     pch = 19, cex = 0.5)
  }
  if (!is.null(subs[[1]]$call.info$Call_ARGS$study.pol)){
    if (!identical(subs[[1]]$call.info$Call_ARGS$study.pol, subs[[1]]$polygons.xy[[1]])){
      graphics::lines(subs[[1]]$call.info$Call_ARGS$study.pol, col = "lightyellow3")
    }
  }
  for (split in 1:length(subs[[1]][["split.lines"]])){
    graphics::lines(x=subs[[1]][["split.lines"]][[split]],
                    col = basal.col,
                    lwd = max.lwd)
  }
}
#' @noRd
.final_plots <- function(subs, pol.col, max.lwd){
  if ( length(pol.col) > 1){
    for (i in 1:length(subs)){
      lines(subs[[i]]$polygons.xy[[1]],
            col = pol.col[i],lwd = c(max.lwd,rep(1,length(subs)-1))[i])
    }
  } else {
    for (i in 1:length(subs)){
      lines(subs[[i]]$polygons.xy[[1]],
            col = pol.col[1],lwd = c(max.lwd,rep(1,length(subs)-1))[i])
    }
  }
}


