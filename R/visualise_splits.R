#' Visualise split lines.
#'
#' @description Helper functions that adds the selected graphical information to
#' the existing plot. After second dot (e.g. '.start') is defined the situation
#' when the function is called.
#' @param what String that informs what graphical information to show. Can be
#' either "curve", "straight" or "both".
#' @author Liudas Daumantas
#' @importFrom pracma poly_center
#' @noRd
.visualise_splits.start <- function(what,
                                    pnts.col, data, margins, rims, testid) {
  if (!is.null(what)){
    if(!is.null(pnts.col)){
      plot(data$x, data$y, col=pnts.col,xlab = "x coordinate",
           ylab = "y coordinate" ,pch=19)
    } else {
      plot(NULL,xlim = range(data$x),ylim = range(data$y),col=0,
           xlab = "x coordinate", ylab = "y coordinate")
    }
    centras <- pracma::poly_center(margins[,1],margins[,2])
    points(centras[1],centras[2],col=2,pch=19,cex=0.65)
    points(centras[1],centras[2],col=2,cex=1.5)
    lines(rims[[1]])

    if (testid>1) {
      for (i in 2:c(testid)){
        lines(x=rims[[i]][,1],y=rims[[i]][,2],col=1,lwd=2)
      }}
    readline(prompt = cat(paste0('\nGoing to test polygon No.: '
                                 ,testid,
                                 "\n\nPress enter to continue...\n")))
  }}
.visualise_splits.best_straight <- function(what,
                                            pairs_pts, maxid, maxdif) {
  if (!is.null(what)){
        lines(pairs_pts[maxid,c(1,3)],pairs_pts[maxid,c(2,4)],lwd=2)
        readline(prompt = cat(paste0('\nThe best straight split-line',
                                     ' is displayed.\nIt\'s quality is: ',
                                     round(maxdif,2),
                                     "\n\nPress enter to continue...\n")))
  }
}
.visualise_splits.curve_start <- function(what,
                                          rot.data, pnts.col, rot.poli, AE,
                                          c.X.knots, split.line.x, c.Y.knots,
                                          knot.y.matrix) {
  if (!is.null(what)){
    if (what != 'straight') {
      dev.new()
      if(!is.null(pnts.col)){
        plot(rot.data$x, rot.data$y, col=pnts.col,xlab = "shifted & rotated X coordinate",
             ylab = "shifted & rotated Y coordinate" ,pch=19,
             xlim = range(rot.poli$x),ylim = range(rot.poli$y))
      } else {
        plot(NULL,xlim = range(rot.poli$x),ylim = range(rot.poli$y),col=0,
             xlab = "shifted & rotated X coordinate", ylab = "shifted & rotated Y coordinate")
      }
      lines(rot.poli)
      lines(c(0,AE),c(0,0),col = "lightblue")
      for (i in 1:c.X.knots){
        points(rep(split.line.x[i],c.Y.knots),knot.y.matrix[,i])
      }
    }
  }
}
.visualise_splits.best_curve <- function(what,
                                         curve.final, SS) {
  if (!is.null(what)){
    if (what != 'straight'){
      invisible(dev.set(dev.list()[length(dev.list())]))
      lines(curve.final,lwd=2)
      readline(prompt = cat(paste0('\nThe best curvial split-line',
                                   ' is displayed.\nIt\'s quality is: ',
                                   round(SS[[1]],2),
                                   "\nNow the window with nonlinear",
                                   " split-lines will close.",
                                   "\n\nPress enter to continue...\n")))
      dev.off(dev.list()[length(dev.list())])
    }
  }
}
.visualise_splits.good_curve <- function(what, level,
                                         curve, SS) {
  if (!is.null(what)){
    if (level != 'best' & what != 'straight'){
      invisible(dev.set(dev.list()[length(dev.list())]))
      lines(curve,lwd=2, col = 'gray70')
      readline(prompt = cat(paste0('\nThe displayed curvial split-line',
                                   ' is so far the best.',
                                   '\nIt\'s quality is: ',
                                   round(SS[[1]],2),
                                   "\n\nPress enter to continue...\n")))
    }}}
.visualise_splits.good_straight <- function(what, level,
                                            pairs_pts, Skirtumas, i) {
  if (!is.null(what)){
    if (level != 'best' & what != 'curves'){
      lines(x= pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],lwd=2,
            col = 'gray70')
      readline(prompt = cat(paste0('\nThe displayed straigth split-line',
                                   ' is so far the best.',
                                   '\nIt\'s quality is: ',
                                   round(Skirtumas,2),
                                   "\n\nPress enter to continue...\n")))
    }}}
.visualise_splits.try_straight <- function(what, level, pairs_pts, i) {
  if (!is.null(what)){
    if (what != "curves" & level == "all"){
      lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
            col = "orange")
      readline(prompt = cat(paste0("\nTesting straight split-line",
                                   " No. ", i,".\n\nPress enter to",
                                   " continue...\n")))
    }}}
.visualise_splits.bad_straight <- function(what, level,
                                           pairs_pts, message, i) {
  if (!is.null(what)){
    if (what != "curves" & level == "all"){
      lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
            col = "gray70")
      readline(prompt = cat(paste0("\nSplit-line was not selected.",
                                   '\nReason: ',message,
                                   "\n\nPress enter to continue...\n")))
    }}}
.visualise_splits.try_curve <- function(what, level,
                                        counter, curve) {
  if (!is.null(what)){
    if (what != "straight" & level == "all"){
      invisible(dev.set(dev.list()[length(dev.list())]))
      lines(curve, col = 'orange')
      readline(prompt = cat(paste0("\nTesting curvial split-line",
                                   " No. ", counter,
                                   ".\n\nPress enter to",
                                   " continue...\n")))
    }}}
.visualise_splits.bad_curve <- function(what, level,
                                        curve, message) {
  if (!is.null(what)){
    if (what != "straight" & level == "all"){
      invisible(dev.set(dev.list()[length(dev.list())]))
      lines(curve, col = 'gray70')
      readline(prompt = cat(paste0("\nSplit-line was not selected.",
                                   ' \nReason: ',message,
                                   "\n\nPress enter to continue...\n")))
    }}}
.visualise_splits.best_split <- function(what,
                                         best.splitl, maxdif) {
  if (!is.null(what)){
    lines(best.splitl,col = "gold", lwd=2)
    readline(prompt = cat(paste0('\nThe best split-line',
                                 ' that will be used to',
                                 ' divide polygon',
                                 ' is displayed.',
                                 '\nIt\'s quality is: ',
                                 round(maxdif,2),
                                 "\n\nPress enter to continue...\n")))
  }
}

