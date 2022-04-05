#' Visualise split lines.
#'
#' @description Helper function that adds the selected graphical information to
#' the existing plot.
#' @param what String that informs what graphical information to show. Can be
#' either "curves", "straight" or "both".
#' @param level Sting indicating the trace level. Can be one of the following:
#' "best", "main", "all".
#' @param when String that states the situation, where the graphical information
#' is needed.
#' @author Liudas Daumantas
#' @importFrom pracma poly_center
#' @importFrom dplyr case_when
#' @noRd
.visualise_splits <- function(what,level, when) {
  if (!is.null(what)){
    if( when == "start"  ) {
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
      } else {
      if (when == "best.straight" & what != "curves"){
        lines(pairs_pts[maxid,c(1,3)],pairs_pts[maxid,c(2,4)],lwd=2)
        readline(prompt = cat(paste0('\nThe best straight split-line',
                                     ' is displayed.\nIt\'s quality is: ',
                                     round(maxdif,2),
                                     "\n\nPress enter to continue...\n")))
      } else {
        if (when == "best.curve" & what != 'straight'){
          lines(curve.final,lwd=2)
          readline(prompt = cat(paste0('\nThe best curvial split-line',
                                       ' is displayed.\nIt\'s quality is: ',
                                       round(SS[[1]],2),
                                       "\n\nPress enter to continue...\n")))
        } else {
          if (when == 'good.curve' & level != 'best' & what != 'straight'){
            lines(curve,lwd=2, col = 'gray70')
            readline(prompt = cat(paste0('\nThe displayed curvial split-line',
                                         ' is so far the best.',
                                         '\nIt\'s quality is: ',
                                         round(SS[[1]],2),
                                         "\n\nPress enter to continue...\n")))
          } else {
            if (when == 'good.straight' & level != 'best' &
                what != 'curves'){
              lines(x= pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],lwd=2,
                    col = 'gray70')
              readline(prompt = cat(paste0('\nThe displayed straigth split-line',
                                           ' is so far the best.',
                                           '\nIt\'s quality is: ',
                                           round(Skirtumas,2),
                                           "\n\nPress enter to continue...\n")))
            } else {
              if (when == 'try.straight' & what != "curves" & level == "all"){
                lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
                      col = "orange")
                readline(prompt = cat(paste0("\nTesting straight split-line",
                                             " No. ", i,".\n\nPress enter to",
                                             " continue...\n")))
              } else {
                if (when == 'bad.straight'& what != "curves" & level == "all"){
                  lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
                        col = "gray70")
                  readline(prompt = cat(paste0("\nSplit-line was not selected.",
                                               '\nReason: ',message,
                                               "\n\nPress enter to continue...\n")))
                } else {
                  if (when == 'try.curve' & what != "straight" & level == "all"){
                    lines(curve, col = 'orange')
                    readline(prompt = cat(paste0("\nTesting curvial split-line",
                                                 " No. ", counter,
                                                 ".\n\nPress enter to",
                                                 " continue...\n")))
                  } else {
                    if (when == 'bad.curve' & what != "straight" & level == "all"){
                      lines(curve, col = 'gray70')
                      readline(prompt = cat(paste0("\nSplit-line was not selected.",
                                                   ' \nReason: ',message,
                                                   "\n\nPress enter to continue...\n")))
                    } else {
                      if (when == "best"){
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
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


