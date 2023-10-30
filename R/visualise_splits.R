#' Visualise split lines.
#'
#' @description Helper functions that adds the selected graphical information to
#' the existing plot. After second dot (e.g. '.start') is defined the situation
#' when the function is called.
#' @param what String that informs what graphical information to show. Can be
#' either "curve", "straight" or "both".
#' @author Liudas Daumantas
#' @importFrom grDevices dev.set dev.list dev.set dev.new
#' @importFrom graphics lines points
#' @noRd
.visualise_splits.start <- function(what,
                                    pnts.col, xy.dat, rims, perim_pts,
                                    testid,study.pol.viz) {
  if (!is.null(what)){

    plot(NULL,xlim = range(rims[[testid]]$x),ylim = range(rims[[testid]]$y),col=0,
         xlab = "x coordinate", ylab = "y coordinate")
    if (!is.null(study.pol.viz)){
      graphics::lines(study.pol.viz[,1:2], lwd=2)
    }
    graphics::lines(rims[[1]])
    graphics::points(perim_pts[[1]],pch=17,col="darkgreen",cex = 1)
    if(!is.null(pnts.col)){
      graphics::points(xy.dat$x, xy.dat$y, col=pnts.col,pch=19,cex=0.5)
    }
    centras <- pracma::poly_center(rims[[testid]][,1],rims[[testid]][,2])
    points(centras[1],centras[2],col=2,pch=19,cex=0.65)
    points(centras[1],centras[2],col=2,cex=1.5)

    if (testid>1) {
      for (i in 2:c(testid)){
        graphics::lines(x=rims[[i]][,1],y=rims[[i]][,2],col=1,lwd=2)
      }}
    readline(prompt = cat(paste0('\nGoing to test polygon No.: '
                                 ,testid,
                                 "\n\nPress enter to continue...\n")))
  }}
.visualise_splits.end <- function(pnts.col, xy.dat, rims, study.pol.viz) {
  if(!is.null(pnts.col)){
    plot(xy.dat$x, xy.dat$y, col=pnts.col,xlab = "x coordinate",
         ylab = "y coordinate" ,pch=19, ylim = range(rims[[1]]$y),
         xlim = range(rims[[1]]$x))
  } else {
    plot(NULL, ylim = range(rims[[1]]$y),
         xlim = range(rims[[1]]$x), col=0,
         xlab = "x coordinate", ylab = "y coordinate")
  }
  if (!is.null(study.pol.viz)){
    graphics::lines(study.pol.viz[,1:2], lwd=2)
  }
  lines(rims[[1]])

  if (length(rims)>1) {
    for (i in 2:length(rims)){
      graphics::lines(x=rims[[i]][,1],y=rims[[i]][,2],col=21,lwd=2)
    }}
}
.visualise_splits.best_straight <- function(what,
                                            pairs_pts, maxid, maxdif) {
  if (!is.null(what)){
    graphics::lines(pairs_pts[maxid,c(1,3)],pairs_pts[maxid,c(2,4)],lwd=2)
    readline(prompt = cat(paste0('\nThe best straight split-line',
                                 ' is displayed.\nIt\'s quality is: ',
                                 round(maxdif,2),
                                 "\n\nPress enter to continue...\n")))
  }
}
.visualise_splits.curve_start <- function(what,
                                          rot.dat.cords, pnts.col, rot.poli, AE,
                                          c.X.knots, split.line.x, c.Y.knots,
                                          knot.y.matrix) {
  if (!is.null(what)){
    if (what != 'straight') {
      grDevices::dev.new()
      if(!is.null(pnts.col)){
        plot(rot.dat.cords$x, rot.dat.cords$y, col=pnts.col,
             xlab = "shifted & rotated X coordinate",
             ylab = "shifted & rotated Y coordinate" ,pch=19,
             xlim = range(rot.poli$x),ylim = range(rot.poli$y),cex = 0.5)
      } else {
        plot(NULL,xlim = range(rot.poli$x),ylim = range(rot.poli$y),col=0,
             xlab = "shifted & rotated X coordinate",
             ylab = "shifted & rotated Y coordinate")
      }
      graphics::lines(rot.poli)

      graphics::lines(c(0,AE),c(0,0),col = "lightblue")
      graphics::points(split.line.x[c(-1, - length(split.line.x))],
             rep(0,length(split.line.x)-2),
             pch=8,cex=2,col="blue")
      for (i in 1:c.X.knots){
        graphics::points(rep(split.line.x[i],c.Y.knots-2),
               knot.y.matrix[c(-1,-nrow(knot.y.matrix)),i])
      }
    }
  }
}

.visualise_splits.no_best_curve <- function(what){
  if (!is.null(what)){
    if (what != 'straight'){
      invisible(grDevices::dev.set(grDevices::dev.list()[length(grDevices::dev.list())]))
      readline(prompt = cat(paste0('\nNo better nonlinear split-line',
                                   ' was found.',
                                   "\nNow the window with nonlinear",
                                   " split-lines will close.",
                                   "\n\nPress enter to continue...\n")))
      grDevices::dev.off(grDevices::dev.list()[length(grDevices::dev.list())])
    }
  }
}

.visualise_splits.best_curve <- function(what,
                                         curve.final, SS) {
  if (!is.null(what)){
    if (what != 'straight'){
      invisible(grDevices::dev.set(grDevices::dev.list()[length(grDevices::dev.list())]))
      graphics::lines(curve.final,lwd=2)
      readline(prompt = cat(paste0('\nThe best curvial split-line',
                                   ' is displayed.\nIt\'s quality is: ',
                                   round(SS[[1]],2),
                                   "\nNow the window with nonlinear",
                                   " split-lines will close.",
                                   "\n\nPress enter to continue...\n")))
      grDevices::dev.off(grDevices::dev.list()[length(grDevices::dev.list())])
    }
  }
}
.visualise_splits.good_curve <- function(what, level,
                                         curve, SS, best.old.curve) {
  if (!is.null(what)){
    if (level != 'best' & what != 'straight'){
      invisible(dev.set(dev.list()[length(dev.list())]))
      graphics::lines(best.old.curve,lwd=3, col = 'gray70')
      graphics::lines(curve,lwd=2, col = '#56B4E9')
      readline(prompt = cat(paste0('\nThe displayed curvial split-line',
                                   ' is so far the best.',
                                   '\nIt\'s quality is: ',
                                   round(SS[[1]],2),
                                   "\n\nPress enter to continue...\n")))


    }}}
.visualise_splits.try_int_curve <- function(what, level,
                                            curve, x, y) {
  if (!is.null(what)){
    if (what != "straight" & level == "all"){
      invisible(grDevices::dev.set(grDevices::dev.list()[length(grDevices::dev.list())]))
      graphics::lines(curve, col = '#E69F00')
      graphics::points(x,y,col="#E69F00", pch = 8,cex= 2 )
      readline(prompt = cat(paste0("\nPotentially better curve was found by using",
                                   " knot, interpolated from other knot performances",
                                   ".\n\nPress enter to",
                                   " test it...\n")))
    }}}
.visualise_splits.good_int_curve <- function(what, level,
                                         curve, SS, best.old.curve,
                                         x, y) {
  if (!is.null(what)){
    if (level != 'best' & what != 'straight'){
      invisible(grDevices::dev.set(grDevices::dev.list()[length(grDevices::dev.list())]))
      graphics::lines(best.old.curve,lwd=3, col = 'gray70')
      graphics::lines(curve,lwd=2, col = '#56B4E9')
      if (level == "main"){
        graphics::points(x,y,col="#E69F00", pch = 8,cex= 2 )
      }
      readline(prompt = cat(paste0(
        '\nThe curve produced by using interpolated, potentially better knot',
                                   ' was indeed the best so far.',
                                   '\nIt\'s quality is: ',
                                   round(SS[[1]],2),
                                   "\n\nPress enter to continue...\n")))


    }}
}
.visualise_splits.bad_int_curve <- function(what, level,
                                            curve, x, y, message) {
  if (!is.null(what)){
    if (level == 'all' & what != 'straight'){
      invisible(grDevices::dev.set(grDevices::dev.list()[length(grDevices::dev.list())]))
      graphics::lines(curve, col = '#999999')
      graphics::points(x, y, col="gray70", pch = 8,cex = 2 )
      readline(prompt = cat(paste0("\nThe interpolated knot did not produce",
                                   ' a better curve.',
                                   ' \nReason: ',message,
                                   "\n\nPress enter to continue...\n")))

    }}
}
.visualise_splits.selected_knot <- function(what, level,
                                            x, y, old.knot.y) {
  if (!is.null(what)){
    if (level != 'best' & what != 'straight'){

      invisible(grDevices::dev.set(grDevices::dev.list()[length(grDevices::dev.list())]))
      graphics::points(x,y,col="#D55E00", pch = 8,cex= 2 )
      graphics::points(x,old.knot.y,col="gray70", pch = 8,cex= 2 )
      cat(paste0('\nSelected knot is displayed.\n'))

    }}}

.visualise_splits.interpol_knot<- function(what = trace.object,
                                level = trace.level,
                                x , y, old.knot.y){
  if (!is.null(what)){
    if (level != 'best' & what != 'straight'){
      invisible(grDevices::dev.set(grDevices::dev.list()[length(grDevices::dev.list())]))
      graphics::points(x,y,col="purple", pch = 8,cex= 2 )
      graphics::points(x,old.knot.y,col="gray70", pch = 8,cex= 2 )
      cat(paste0('\nInterpolated knot is displayed.\n'))

    }}}


.visualise_splits.good_straight <- function(what, level,
                                            pairs_pts, Skirtumas, i, maxid) {
  if (!is.null(what)){
    if (level != 'best' & what != 'curve'){
      graphics::lines(x= pairs_pts[maxid,c(1,3)], y = pairs_pts[maxid,c(2,4)],lwd=3,
                      col = 'gray70')
      graphics::lines(x= pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],lwd=2,
                      col = '#56B4E9')
      readline(prompt = cat(paste0('\nThe displayed straigth split-line',
                                   ' is so far the best.',
                                   '\nIt\'s quality is: ',
                                   round(Skirtumas,2),
                                   "\n\nPress enter to continue...\n")))

    }}}
.visualise_splits.try_straight <- function(what, level, pairs_pts, i) {
  if (!is.null(what)){
    if (what != "curve" & level == "all"){
      graphics::lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
            col = "#E69F00")
      readline(prompt = cat(paste0("\nTesting straight split-line",
                                   " No. ", i,".\n\nPress enter to",
                                   " continue...\n")))
    }}}
.visualise_splits.bad_straight <- function(what, level,
                                           pairs_pts, message, i) {
  if (!is.null(what)){
    if (what != "curve" & level == "all"){

      readline(prompt = cat(paste0("\nSplit-line was not selected.",
                                   '\nReason: ',message,
                                   "\n\nPress enter to continue...\n")))
      graphics::lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
            col = "#999999")
    }}}
.visualise_splits.try_curve <- function(what, level,
                                        counter, curve) {
  if (!is.null(what)){
    if (what != "straight" & level == "all"){
      invisible(grDevices::dev.set(grDevices::dev.list()[length(grDevices::dev.list())]))
      graphics::lines(curve, col = '#E69F00')
      readline(prompt = cat(paste0("\nTesting curvial split-line",
                                   " No. ", counter,
                                   ".\n\nPress enter to",
                                   " continue...\n")))
    }}}
.visualise_splits.bad_curve <- function(what, level,
                                        curve, message) {
  if (!is.null(what)){
    if (what != "straight" & level == "all"){

      readline(prompt = cat(paste0("\nSplit-line was not selected.",
                                   ' \nReason: ',message,
                                   "\n\nPress enter to continue...\n")))
      invisible(grDevices::dev.set(grDevices::dev.list()[length(grDevices::dev.list())]))
      graphics::lines(curve, col = '#999999')
    }}}
.visualise_splits.best_split <- function(what,
                                         best.splitl, maxdif) {
  if (!is.null(what)){
    graphics::lines(best.splitl,col = "gold", lwd=2)
    readline(prompt = cat(paste0('\nThe best split-line',
                                 ' that will be used to',
                                 ' divide polygon',
                                 ' is displayed.',
                                 '\nIt\'s quality is: ',
                                 round(maxdif,2),
                                 "\n\nPress enter to continue...\n")))
  }
}

