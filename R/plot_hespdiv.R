#' Plot the results of hespdiv object
#'
#' @description  Function prints the results of hespdiv object.
#' @param obj hespdiv object
#' @param xy.dat data.frame. coordinates of locations of data points used in
#' hespdiv.
#' @param type character. Either "width" or "color" (default "color").
#' Determines whether quality of split-lines is expressed by line width or
#' color.
#' @param n.loc logical. Would you like to visualize the number of observations
#' at each location? Only possible, when there are localities with more than
#' one observation. If type is 'color', then no. of obs. are expressed via point
#' sizes. Otherwise, they are expressed by color in logarithmic scale.
#' @param title character. Metric title in legend. The default is built
#' according to the obj$call.info$METHOD$metric.
#' @param seed integer. Seed used to randomized the colors of the split-lines.
#' Only meaningful, when argument type is 'width'.
#' Try setting different value, if colors of parallel split-lines or nearby
#' labels look too similar or to increase the general appeal of the graph.
#' @param pnts.col color of data points.
#' @param ... other arguments
#' @return ggplot object
#' @importFrom ggplot2 aes geom_point geom_path guides guide_legend
#' @author Liudas Daumantas
#' @export
ggplot_hespdiv <- function(obj, xy.dat, type = "color",n.loc = FALSE,
                          title = NULL,  pnts.col =
                           NULL,seed = 10,...){

  type <- .arg_check("type",type,c("width","color"))

  if (n.loc){
    xy_df <- xy.dat
    xy_df$n <- 1
    uni.loc.n <- aggregate(n~., data = xy_df ,FUN = sum)
    uni.loc.n <- uni.loc.n[order(uni.loc.n$n,decreasing = FALSE),]
    if (nrow(uni.loc.n) == nrow(xy_df)){
      stop("All observations are from unique locations.",
           " Please use n.loc = FALSE.")
    }
    if (type == "color") {
      if (is.null(pnts.col)) {
        pnts.col <- rep(1, nrow(uni.loc.n))
      } else {
        if (length(pnts.col) != nrow(uni.loc.n))
          stop("Length of pnts.col is not equal to the number of unique ",
               "locations.")}
    } else {
      if (!is.null(pnts.col)){
        stop("Conflincting arguments: pnts.col is not null, when n.loc is TRUE",
             " and type is 'width'.", "\nCannot add two same",
             " type (color & color, size & size) aesthetics.")
      }
    }
  } else {
    if (is.null(pnts.col)) {
      pnts.col <- rep(1, nrow(xy.dat))
    } else {
      if (length(pnts.col) != nrow(xy.dat))
        stop("Length of pnts.col is not equal to the number of observations.")
    }
  }

  split.stats <- obj$split.stats
  maximize <- obj$call.info$METHOD$maximize
  if (is.null(title)){
    if (obj$call.info$METHOD$method.type == "custom"){
      title <- obj$call.info$METHOD$metric
    } else {
      if (obj$call.info$METHOD$metric == "sorensen"){
        title <- paste0("S",rawToChar(as.raw(184)),"rensen-Dice\ncoefficient")
      } else {
        if (obj$call.info$METHOD$metric == "morisita"){
          title <- paste0("Morisita\n","overlap")
        } else {
          if (obj$call.info$METHOD$metric == "pielou"){
            title <- paste0("Pielou\nentropy\nreduction\n(%)")
          } else {
            if (obj$call.info$METHOD$metric == "horn.morisita"){
              title <- paste0("Morisita\n","overlap\n(Horn)")
            }
          }
        }
      }
    }
  }
  if (!maximize){
    split.stats[,"performance"] <- -obj$split.stats[,"performance"]
  }
  ord <- order(split.stats[,"performance"], decreasing = FALSE)
  split.stats <- split.stats[ord,]
  split.lines <- lapply(ord,function(id){obj$split.lines[[id]]})
  df <- Reduce(rbind,split.lines)
  npt.in.split <- as.numeric(lapply(split.lines,nrow))
  if (type == "width"){
    size <- rep(split.stats[,"performance"], times = npt.in.split)
  } else {
    color <- rep(split.stats[,"performance"], times = npt.in.split)
  }
  df$group <- factor(rep(1:length(split.lines),times=npt.in.split))

  base <- ggplot2::ggplot(obj$polygons.xy[[1]],aes(x,y),xlab = 'x', ylab = 'y') +
    geom_path(data= obj$polygons.xy[[1]],aes(x,y), size=0.5,
              lineend = "round",linejoin = "round",color = 1)

  if (n.loc){
    scale_id <- 2
    if (type == "color"){
      base <- base + geom_point(data=uni.loc.n,aes(x=x,y=y,size=n),
                                pch =rep(1,nrow(uni.loc.n)),color = pnts.col) +
        guides(size=guide_legend(title=paste0("Number of", "\nobservations" ,
       "\nin a location"),title.hjust = 0.5,label.position = "left",
       label.hjust = 1))+
        ggplot2::scale_size_area(max_size = 8)
    } else {
      base <- base + geom_point(data=uni.loc.n,aes(x=x,y=y,color=n),
                                pch =rep(19,nrow(uni.loc.n)),size =2) +
        viridis::scale_color_viridis(guide ="colourbar",trans = "log") +
        guides(color = ggplot2::guide_colourbar(
          title = paste0("Number of", "\nobservations" ,
                         "\nin a location"),label.position = "left",
          label.hjust = 1, title.hjust = 0.5,title.vjust = 1))

      base$scales$scales[[1]]$limits <-  range(log(uni.loc.n$n))
      base$scales$scales[[1]]$labels <- round(exp(as.numeric(na.omit(
        base$scales$scales[[1]]$get_breaks() ))))

    }

  } else {
    scale_id <- 1
  base <- base + geom_point(
    pch=16,color=pnts.col) +
    geom_path(data= obj$polygons.xy[[1]],aes(x,y), size=.5,
                           lineend = "round",linejoin = "round")
  }
  if (type == "width"){
    color <- .generate_cols(nrow(split.stats), seed)
    df$color <- rep(color, times=npt.in.split)
    df$size <- size
    base<-base + geom_path(data = df, aes(x,y,group=group,
                                          size = size),
                           color=df$color) +
      ggplot2::scale_size_continuous(range = c(0.5,2))
    size.l <- seq(0.5,2,0.5)
    base <- base + guides(size = guide_legend(override.aes =
                                                list(size = size.l))) +
      guides(size=guide_legend(title=title, title.hjust = 0.5))
  } else {
    df$color <- color
    base<-base + geom_path(data = df, aes(x,y,group=group, color=color),size = 2) +
      viridis::scale_color_viridis(guide ="colourbar") +
    guides(color = ggplot2::guide_colourbar(title = title, title.hjust = 0.5,
                                   label.position = "left",label.hjust = 1))
  }
  if (!maximize) {
      base$scales$scales[[scale_id]]$limits <-  range(split.stats[,"performance"])
      base$scales$scales[[scale_id]]$breaks <-
        as.numeric(na.omit( base$scales$scales[[scale_id]]$get_breaks() ))
      base$scales$scales[[scale_id]]$labels <-
        -base$scales$scales[[scale_id]]$breaks
}
  if (type == "color") {
    if (maximize){
      base$scales$scales[[scale_id]]$limits <-  range(split.stats[,"performance"])
    }
    self <- base$scales$scales[[scale_id]]
    x <- self$rescale(self$oob(split.stats[,"performance"], range = self$limits), self$limits)
    color <- self$palette(x)
  }
  base <- base +
    ggplot2::theme_set(ggplot2::theme_bw())  +
    ggplot2::theme(legend.key= ggplot2::element_rect(fill = "white"),
          legend.title = ggplot2::element_text(hjust = 0.5),
          panel.grid = ggplot2::element_blank(),panel.background =
            ggplot2::element_rect(colour = "black", size=0.5,fill = "white"))

  mid.pt <- data.frame(x=numeric(),y=numeric())
  for (a in 1:length(obj$split.lines)){
    if(nrow(split.lines[[a]])!=2){
      mid.pt <- rbind(mid.pt,split.lines[[a]][
        round(length(split.lines[[a]]$x)/2,0),])} else {
          mid.pt<-rbind(mid.pt,data.frame(x=mean(split.lines[[a]]$x),
                                          y=mean(split.lines[[a]]$y)))
        }
  }
  if (is.null(split.stats$p.val)){
    base<-base + ggrepel::geom_label_repel(data=mid.pt, aes(x,y),alpha=rep(3/5, nrow(mid.pt)),
                                  label = paste0(ord,") ",
                                                round(obj$split.stats[ord,"performance"],2)),
                                  fill  = color, size = 4,
                                  direction="both",fontface='bold',
                                  colour  = rep(1, nrow(mid.pt)))
  } else {
    base<-base + ggrepel::geom_label_repel(data=mid.pt, aes(x,y),alpha=rep((3/5), (nrow(mid.pt))),
                                  label = paste(ord,") p = ",
                                                split.stats$p.val,
                                                "\n  Div. qual. = ",
                                                round(split.stats$delta.E,2)),
                                  fill = color,
                                  size = 3.5,direction="both",fontface='bold')
  }
  base
}
