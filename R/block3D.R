#' Draw hespdiv polygons in 3D space
#'
#' This function vizualise polygons obtained by hespdiv function. Polygon height
#' (z coordinate) encodes chosen information from the poly.stats data frame.
#' @param obj hespdiv object
#' @param height character (default 'mean'). What information from poly.stats
#' data frame do you wish to encode as the height of the polygons? Options:
#' "mean", "sd", "best", "z.score", "str.best", "str.z.score". See details
#' for explanations.
#' @param color.seed integer. Controls the color of polygons. Change it to a
#' different number if you want to get a different set of colors.
#' @param lines logical. Do you want split-lines to be displayed over the top of
#' the polygons?
#' @param obs logical. Do you want observations to be displayed over the top of
#' the polygons?
#' @param pnts.col vector. Color codes to be used for displaying observations.
#' @author Liudas Daumantas
#' @export
blok3d <- function(obj,height = "mean", color.seed=1, lines=TRUE, pnts.col = NULL,
                   obs = TRUE) {
  height <- as.vector(sapply(height, .arg_check, name = "height", NAMES = c("mean","sd","best","z.score",
                                                 "str.best", "str.z.score")))
  for (height in height){
    rgl::open3d()
    poly.stats <- obj$poly.stats
    if (obj$call.info$Call_ARGS$c.splits){
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
    } else {
      poly.stats$best <- obj$poly.stats$str.best
      poly.stats$z.score <- obj$poly.stats$str.z.score
    }

    if (is.null(obj$call.info$Call_ARGS$xy.dat)){
      xy.dat <- obj$call.info$Call_ARGS$data[,c('x','y')]
    } else {
      xy.dat <- obj$call.info$Call_ARGS$xy.dat
    }
    ZZ <- .Zcoords(poly.stats,height)
    del.id <- which(ZZ[,1]==ZZ[,2])
    Zoff <- mean(apply(ZZ[-del.id,], 1, function(x){x[2]-x[1]}))/5
    if (height %in% c('z.score','str.z.score')){
      highest.z <- min(ZZ,na.rm = T)
    } else {
      highest.z <- max(ZZ,na.rm = T)
    }
    if (is.null(pnts.col)){
      pnts.col <- 1
    }
    # ZZ - is data.frame of z coordinates of blocks. Col 1 is z cord of the bottom
    # and col 2 is z cord of the ceiling. Rows correspond to polygons in
    # obj$poly.stats


    #NUPIESIAM RIBAS
    #sukuriam erdve, kurioj piesim blokus - min, max koordinates padedam, nematomi taskiukai plius
    #rgl.material(alpha = 0)
    #spalvu palete sukuriam
    n <- nrow(ZZ)-length(del.id)
    palete <- .generate_cols(n, seed=color.seed)


    # jei uzduoti duomenys, nupaisom erdves pavirsiuje duomenis, spalvu paletes paskutiniai elementai - rusiu spalvos
    if (obs){
      alfa <- 1
    } else {
      alfa <- 0
    }


    rgl::plot3d(x = xy.dat$x, y = xy.dat$y, z = highest.z,
                col = pnts.col , alpha = alfa, add = F,
                ylim = range(obj$polygons.xy$`1`$y),
                xlim = range(obj$polygons.xy$`1`$x),
                zlim = range(unlist(c(ZZ)),
                             na.rm = T),
                zlab = height, xlab = "x", ylab = "y")
#
    #paisom kiekvieno bloko pavirsius
    for (i in seq(n)) {
      .draw_poly(obj,(seq(nrow(ZZ))[-del.id])[i],color=palete[i], ZZ, Zoff)
    }
    if (lines==TRUE){
      for (i in seq(nrow(ZZ))){
        rgl::rgl.linestrips(x=obj$polygons.xy[[i]]$x,y=obj$polygons.xy[[i]]$y,
                            z=highest.z,col=1)
      }
    }
  }
}


#' Iterative function that provides the id of the closest polygon parent
#' @noRd
.iterator<-function(poly.id,roots,id.list){
  id.list<-c(id.list,roots[poly.id])
  if (roots[poly.id]==0) {return(id.list)
  } else {
    poly.id <- roots[poly.id]
    return(.iterator(poly.id,roots,id.list))
  }
}
#' Return the ids of parent polygons for a given polygon id.
#' @noRd
.collector<-function(poly.id,roots){
  id.list <- roots[poly.id]
  if (id.list == 0) {return(0)
  } else {
    poly.id <- roots[poly.id]
  }
  rezas <- .iterator(poly.id,roots,id.list)
  return(rezas[-length(rezas)])
}
#' Return bottom z coordinate for a polygon of given poly.id
#' @noRd
.zmin <- function(poly.id, poly.stats, height){
  roots <- poly.stats$root.id #[[3]][,5]-1, turi prasideti nuo 1
  ids <- .collector(poly.id = poly.id,roots = roots) # ids surenka iki duoto polio
  zmin <- sum(poly.stats[ids,height])
  return(zmin)
}
#' Obtain top and bottom Z coordinates for each polygon that has non-NA
#' performance scores.
#' @noRd
.Zcoords<-function(poly.stats,height){
  len <- nrow(poly.stats)
  zmini <- numeric(len) # bottom z coord. for each poly.
  zmaxi <- numeric(len) # ceiling z coord. for each poly.
  for (poly.id in 1:len){
    zmini[poly.id] <- .zmin(poly.id,poly.stats,height)
    zmaxi[poly.id] <- zmini[poly.id]+
      ifelse(is.na(poly.stats[poly.id,height]),0,
             poly.stats[poly.id,height]) # [[3]]
  }
  return(data.frame(zmini,zmaxi))
}
# Add a block of polygon to the rgl device
#' @noRd
.draw_poly<-function(obj,i,color,ZZ, Zoff){
  #bloko koordinates
  x <- obj$polygons.xy[[i]]$x
  y <- obj$polygons.xy[[i]]$y

  zmini <- ZZ[i,1]
  zmaxi <- ZZ[i,2]
  #virsus
  rgl::polygon3d(x, y, z=rep(zmaxi,length(x)), col = color,add=T,fill = F)
  #apacia
  rgl::polygon3d(x, y, z=rep(zmini,length(x)), col =color, add=T, fill =F)
  #sonai
  xmat <- matrix(NA, 2, length(x))
  ymat <- matrix(NA, 2, length(x))
  zmat <- matrix(NA, 2, length(x))
  ZZZ<-c(zmini,zmaxi)
  for (a in 0:1) {
    xmat[a+1,] <- x
    ymat[a+1,] <- y
    zmat[a+1,] <- ZZZ[a+1]
  }
  rgl::persp3d(x=xmat,y=ymat,z=zmat,col=color,add = T)
  #ID

  centras <- pracma::poly_center(x,y)
  rgl::rgl.points(x=centras[1],y=centras[2],z=zmaxi,col=1,pch=19,labels=i)
  rgl::rgl.texts(x=centras[1],y=centras[2],z=zmaxi+Zoff, text = c(i),
                 justify="left")

}

