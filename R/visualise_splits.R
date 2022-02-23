#' Visualise split lines.
#'
#' @description Helper function that adds the selected graphical information to
#' the existing plot.
#' @param what String that informs what graphical information to add. It can be
#' one of the following: 0, best, best.straight, best.curves, main,
#' main.straight, main.curves, all, all.straight, all.curves.
#' @param when String that states the situation, where the graphical information
#' is needed.
#' @author Liudas Daumantas
#' @noRd
.visualise_splits <- function(what, when) {
  if (what != "0"){

    if (when == "start"){

      if(!is.null(pnts.col)){
        plot(data$x, data$y, col=pnts.col )
      } else {
        plot(0,0,xlim = range(data$x),ylim = range(data$y),col=0)
      }
      lines(study.pol)
    }

    if (when == "1"){

    }
    if (when == "1"){

    }
    if (when == "1"){

    }
    if (when == "1"){

    }
  }
}
