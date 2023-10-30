#' Remove polygons from rgl device
#'
#' This function allows you to interactively select and remove unwanted polygons
#' from a 3D plot created with the \code{blok3d} function.
#'
#' @param obj The hespdiv object used to create the currently active rgl
#' device with the \code{blok3d} function.
#' @param height A character value that indicates the height co-ordinate.
#' @family {function for hespdiv visualization in 3D}
#' @author Liudas Daumantas
#' @importFrom pracma poly_center
#' @importFrom rgl identify3d rgl.ids rgl.pop
#' @export
polypop <- function(obj,height){
  height <- .arg_check("height",height, c("mean","sd","best","z.score",
                                          "str.best", "str.z.score","rank"))
  poly.stats <- obj$poly.stats
  if (height == "best"){
    poly.stats$best <- ifelse(is.na(obj$poly.stats$is.curve),
                              obj$poly.stats$str.best,
                              ifelse(obj$poly.stats$is.curve,
                                     obj$poly.stats$crv.best,
                                     obj$poly.stats$str.best))
  } else {
    if (height == "z.score" ){
      poly.stats$z.score <- ifelse(is.na(obj$poly.stats$is.curve),
                                   obj$poly.stats$str.z.score,
                                   ifelse(obj$poly.stats$is.curve,
                                          obj$poly.stats$crv.z.score,
                                          obj$poly.stats$str.z.score))
    }
  }
  if (height == "rank") {
    ZZ <- data.frame(zmin = obj$poly.stats$rank-1, zmax = obj$poly.stats$rank)
    del.id <- numeric()
  } else {
    ZZ <- .Zcoords(poly.stats,height)
    del.id <- which(ZZ[,1]==ZZ[,2])
  }
  if (length(del.id) != 0){
    true.ids <- poly.stats$plot.id[-del.id]
  } else {
    true.ids <- poly.stats$plot.id
  }
  basic.id<-5
  OIDS <- seq(length(true.ids))

  while(length(OIDS)>1){
    centrai <- data.frame(x=rep(NaN,length(OIDS)),
                          y=rep(NaN,length(OIDS)))
    for (i in seq(length(OIDS))) {
      centrai[i,] <- pracma::poly_center(
        obj$polygons.xy[[true.ids[OIDS[i]]]][,1],
        obj$polygons.xy[[true.ids[OIDS[i]]]][,2])
    }

    XO <- centrai[,1]
    YO <- centrai[,2]

    ZO <- ZZ[true.ids[OIDS],2]
    LABS <- poly.stats[true.ids[OIDS],"plot.id"]
    cat("Select the centers of polygons you wish to remove.\n")
    pts <- rgl::identify3d(x=XO, y=YO, z=ZO, labels = LABS)
    if(length(pts)==1) {
      id.start<-basic.id+((pts-1)*5)
      OIDS<-OIDS[-which(OIDS==OIDS[pts])]
      id.finito<-id.start+4
      ids<-id.finito:id.start
      ids<-rgl::rgl.ids( type = "shapes", subscene = NA )[ids,1]
      rgl::rgl.pop(type = "lights",id=ids)
    } else {
      if (length(pts>1)){
        pts<-sort(pts,T)
        for (i in 1:length(pts)) {
          id.start<-basic.id+((pts[i]-1)*5)
          OIDS<-OIDS[-which(OIDS==OIDS[pts[i]])]
          id.finito<-id.start+4
          ids<-id.finito:id.start
          ids<-rgl::rgl.ids( type = "shapes", subscene = NA )[ids,1]
          rgl::rgl.pop(type = "lights",id=ids)
        }
      } else {cat("Nothing was selected\n")}
    }

  }
  cat("Nothing to remove - there is only one polygon.\n")
}
