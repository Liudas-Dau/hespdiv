#' Pairing connectable polygon points
#'
#' This function pairs points that are located on a perimeter of a polygon and can be joined by a straight line that is strictly contained within the polygon.
#' @param per_pts a data frame of 4 columns, first item of list output of \code{perimeter_pts} function.
#' @param polygon A data frame containing coordinates of a polygon. The polygon can be either closed or open.
#' @param check.geometry If convex hull is used as a polygon and only straight split-lines are generated, they cannot cross polygon boundary. Therefore, no need to check for this.
#' @return A data frame with 7 columns. Each row represents a connectable pair of points located on a perimeter of a polygon.
#' \itemize{
#'   \item \code{x1} - x coordinate of a point that has lower x coordinate.
#'   \item \code{y1} - y coordinate of a point that has lower x coordinate.
#'   \item \code{x2} - x coordinate of a point that has higher x coordinate.
#'   \item \code{y2} - y coordinate of a point that has higher x coordinate.
#'   \item \code{b} - slope of a straight line connecting the paired points.
#'   \item \code{id1} - relative location of a point that has lower x coordinate. This ID is derived from \code{per_pts} third column (\code{ID}).
#'   \item \code{id2} - relative location of a point that has higher x coordinate. This ID is derived from \code{per_pts} third column (\code{ID}).
#' }
#' @note If both points have the same x coordinate, a point with lower ID will be treated as the one with lower x coordinate.
#' @author Liudas Daumantas
#' @importFrom utils combn
#' @examples #Creating data.frame of a polygon
#' poly<-data.frame(X=c(3.38,3.30,1.70,0.78,-0.06,-2.30,-2.94,-3.97,-1.61,
#'                     -0.39,0.68,1.28,1.60,3.38),
#'                 Y=c(-0.12,-0.31,-2.73,-3.22,0,-2.19,-1.62,0.94,3.10,
#'                 3.00,2.91,2.49,2.20,-0.12))
#' plot(poly,type='o')
#' #Generating 10 points on a polygon perimeter
#' a <- .perimeter_pts(poly,n.pts = 15)
#' #location of points
#' points(a[[1]][,-3],col=2,pch=19)
#' #ID of points
#' text(x=a[[1]][,1],y=a[[1]][,2],a[[1]][,3],-0.3)
#' #pairing points
#' b <- .pair_pts(a[[1]],poly)
#' #all arrows point to the right and all lines are within the polygon.
#' arrows(x0=b[,1],y0=b[,2],x1=b[,3],y1=b[,4],length = 0.15,col=2)
#' @noRd
.pair_pts <- function(per_pts, polygon, check.geometry = TRUE) {

  # Fast path:
  # For a convex polygon, a straight line connecting two perimeter points
  # cannot cross outside the polygon.
  if (!check.geometry) {
    message("No geometry checking")
    pairs <- utils::combn(seq_len(nrow(per_pts)), 2)

    # Remove pairs whose two points are on the same polygon segment
    keep <- per_pts[pairs[1, ], 4] != per_pts[pairs[2, ], 4]
    pairs <- pairs[, keep, drop = FALSE]

    if (ncol(pairs) == 0) {
      return(data.frame())
    }

    i1 <- pairs[1, ]
    i2 <- pairs[2, ]

    xa <- per_pts[i1, 1]
    ya <- per_pts[i1, 2]
    xb <- per_pts[i2, 1]
    yb <- per_pts[i2, 2]

    # x1 is the endpoint with the smaller x coordinate
    first <- xa <= xb

    x1 <- ifelse(first, xa, xb)
    y1 <- ifelse(first, ya, yb)
    x2 <- ifelse(first, xb, xa)
    y2 <- ifelse(first, yb, ya)

    b <- (ya - yb) / (xa - xb)

    id1 <- ifelse(first, per_pts[i1, 3], per_pts[i2, 3])
    id2 <- ifelse(first, per_pts[i2, 3], per_pts[i1, 3])

    return(
      data.frame(
        x1 = x1,
        y1 = y1,
        x2 = x2,
        y2 = y2,
        b = b,
        id1 = id1,
        id2 = id2
      )
    )
  }


  # Full geometry checking for non-convex polygons
  pairs_pts <- list()
  pair_id <- 0L

  for (a in seq_len(nrow(per_pts) - 1L)) {

    candidates <- seq.int(a + 1L, nrow(per_pts))

    # No need to test points lying on the same polygon segment
    candidates <- candidates[
      per_pts[candidates, 4] != per_pts[a, 4]
    ]

    for (o in candidates) {

      x <- per_pts[c(a, o), 1]
      y <- per_pts[c(a, o), 2]

      bad <- FALSE

      b1 <- (y[1] - y[2]) / (x[1] - x[2])
      a1 <- y[1] - b1 * x[1]

      # Midpoint can be calculated directly
      mid.x <- mean(x)
      mid.y <- mean(y)

      if (
        .point.in.polygon(
          point.x = mid.x + 7e-14,
          point.y = mid.y + 7e-14,
          pol.x = polygon[, 1],
          pol.y = polygon[, 2]
        ) == 1 &&

        .point.in.polygon(
          point.x = mid.x - 7e-14,
          point.y = mid.y - 7e-14,
          pol.x = polygon[, 1],
          pol.y = polygon[, 2]
        ) == 1 &&

        .point.in.polygon(
          point.x = mid.x - 7e-14,
          point.y = mid.y + 7e-14,
          pol.x = polygon[, 1],
          pol.y = polygon[, 2]
        ) == 1 &&

        .point.in.polygon(
          point.x = mid.x + 7e-14,
          point.y = mid.y - 7e-14,
          pol.x = polygon[, 1],
          pol.y = polygon[, 2]
        ) == 1
      ) {

        range_x <- range(x)
        range_y <- range(y)

        for (seg in seq_len(nrow(polygon) - 1L)) {

          pol_seg <- polygon[seg:(seg + 1L), ]

          if (
            all(pol_seg[1, ] == c(x[1], y[1])) ||
            all(pol_seg[1, ] == c(x[2], y[2])) ||
            all(pol_seg[2, ] == c(x[1], y[1])) ||
            all(pol_seg[2, ] == c(x[2], y[2]))
          ) {
            next
          }

          X1 <- .in.range(range_x, pol_seg$x[1]) ||
            .in.range(range_x, pol_seg$x[2])

          X2 <- .in.range(range(pol_seg$x), x[1]) ||
            .in.range(range(pol_seg$x), x[2])

          Y1 <- .in.range(range_y, pol_seg$y[1]) ||
            .in.range(range_y, pol_seg$y[2])

          Y2 <- .in.range(range(pol_seg$y), y[1]) ||
            .in.range(range(pol_seg$y), y[2])

          if ((X1 && Y1) || (X2 && Y2) ||
              (X1 && Y2) || (X2 && Y1)) {

            b2 <- (pol_seg$y[1] - pol_seg$y[2]) /
              (pol_seg$x[1] - pol_seg$x[2])

            a2 <- pol_seg$y[1] - b2 * pol_seg$x[1]

            x.inters <- (a2 - a1) / (b1 - b2)

            if (is.nan(x.inters) || is.infinite(x.inters)) {
              next
            }

            bad <- .in.range(range_x, x.inters) &&
              .in.range(range(pol_seg$x), x.inters)

            if (bad) {
              break
            }
          }
        }

        if (!bad) {

          id.min <- which.min(x)

          pair_id <- pair_id + 1L

          pairs_pts[[pair_id]] <- data.frame(
            x1 = x[id.min],
            y1 = y[id.min],
            x2 = x[-id.min],
            y2 = y[-id.min],
            b = b1,
            id1 = per_pts[c(a, o)[id.min], 3],
            id2 = per_pts[c(a, o)[-id.min], 3]
          )
        }
      }
    }
  }

  if (length(pairs_pts) == 0) {
    return(data.frame())
  }

  do.call(rbind, pairs_pts)
}


.in.range <- function(range, x) {
  range[1] < x && x < range[2]
}
#' is x in range of x?
#' @noRd
.in.range <- function(range,x){
  range[1] < x & x < range[2]
}
