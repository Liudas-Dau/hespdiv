#' Schematic plot of hespdiv polygons
#'
#' Produces a schematic visualization of subdivided territory, in which
#' the location of each polygon is made clear via their centroinds, their ID
#' label, and punctuated lines that join polygon centroids with the split-lines
#' that produced that polygons.
#' @param obj hespdiv object
#' @param segment Boolean. Display the segments joining the centroids with
#' the split-lines?
#' @param id Boolean. Display the IDs of polygons?
#' @param seed Integer. Randomization seed to produce colors. Optimize by trial
#' and error, if too confusing colors appear.
#' @export
poly_scheme <- function(obj,segment = TRUE, id = TRUE, seed = 1){

  split.stats <- obj$split.stats
  ord <- order(split.stats[,"performance"], decreasing = FALSE)
  split.stats <- split.stats[ord,]
  split.lines <- lapply(ord,function(id){obj$split.lines[[id]]})
  df <- Reduce(rbind,split.lines)
  npt.in.split <- as.numeric(lapply(split.lines,nrow))
  df$group <- factor(rep(1:length(split.lines),times=npt.in.split))

  base <- ggplot2::ggplot(obj$polygons.xy[[1]],aes(x,y),xlab = '', ylab = '') +
    ggplot2::geom_path(data= obj$polygons.xy[[1]],aes(x,y), size=0.5,
              lineend = "round",linejoin = "round",color = 1)
  scale_id <- 1
  color <- .generate_cols(nrow(split.stats), seed)
  df$color <- rep(color, times=npt.in.split)
  base<-base + ggplot2::geom_path(data = df, aes(x,y),group=df$group, color=df$color,size = 1)+
    ggplot2::theme_set(ggplot2::theme_void())  +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),panel.background =
                     ggplot2::element_blank())
  centrai <- as.data.frame(Reduce(rbind,lapply(obj$polygons.xy, FUN = function(o){
    pracma::poly_center(x = o[,1],y = o[,2])})))
  names(centrai) <- c("x1","y1")

  base <- base + ggplot2::geom_point(data=centrai,aes(x1,y1),shape = 20,size=3,
                                     color = c(1,
                                               color[
                                                 match(obj$poly.stats$root.id[-1],
                                                       split.stats$plot.id)]))
  if (segment){
    id.cord <- which.min(apply(obj$polygons.xy[[1]],1, FUN = function(o){
      sqrt(sum((o - unlist(centrai[1,]))^2))
    }))
    split.mid.p <- as.data.frame(rbind(unlist(obj$polygons.xy[[1]][id.cord,]),
                                       Reduce(rbind,lapply(obj$split.lines,function(o){if (nrow(o)==2){
                                         if (o[1,1]-o[2,1]!=0){
                                           c(mean(o[,1]),.pt_on_line(x1 = o[1,1],x2 = o[2,1], y1 = o[1,2],y2 = o[2,2],
                                                                     x3 = mean(o[,1])))} else {
                                                                       c(o[1,1],mean(o[,2]))
                                                                     }
                                       } else {
                                         c(o[round(nrow(o)/2,0),1],o[round(nrow(o)/2,0),2])
                                       }}))))
    names(split.mid.p) <- c("x2","y2")
    segments <- cbind(centrai,
                      split.mid.p[c(1,match(obj$poly.stats$root.id[-1],
                                            obj$split.stats$plot.id)+1),])
    base <- base + ggplot2::geom_segment(data=segments, linetype = 3,
                                aes(x=x1,y=y1,xend = x2,yend =y2),
                                color=c(1,color[match(obj$poly.stats$root.id[-1],
                                                      split.stats$plot.id)]))}
  if (id){
    base <- base + ggplot2::geom_text(data=centrai,
                             aes(x=x1,y=y1,label=obj$poly.stats$plot.id),
                             nudge_y = 1,
                             color=c(1,color[match(obj$poly.stats$root.id[-1],
                                                   split.stats$plot.id)]))}
  base
}

