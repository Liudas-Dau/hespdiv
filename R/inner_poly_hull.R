#' Create data frame of x and y intervals of selected polygon segments
#' @param rule logic vector indicating which polygon segments are selected
#' @param cut_poly data.frame of a polygon
#' @author Liudas Daumantas
#' @noRd
.create_ints <- function(rule,cut_poly){
  data.frame(x1 = cut_poly[which(rule),1],
             y1 = cut_poly[which(rule),2],
             x2 = cut_poly[which(rule)+1,1],
             y2 = cut_poly[which(rule)+1,2])
}

#' Indicate whether the provided polygon is above zero
#' @param lowest_up_intervals data frame of filtered polygon segments are the
#' nearest to the y axis.
#' @param polygon data frame of a polygon
#' @author Liudas Daumantas
#' @importFrom  sp point.in.polygon
#' @noRd
.is_pol_up <- function(lowest_up_intervals,polygon) {
  sp::point.in.polygon(point.x = lowest_up_intervals[2,1],
                   point.y = lowest_up_intervals[2,2]-lowest_up_intervals[2,2]/2,
                   pol.x = polygon[,1],
                   pol.y = polygon[,2]) != 0

}

#' Convert polygon segments data.frame into points data.frame
#' @param inters data frame of x and y intervals of polygon segments
#' @noRd
.ints_to_pnts <- function(inters) {
  idx<-cbind(rep(1:nrow(inters),each=2),rep(c(1,3),nrow(inters)))
  idy<-cbind(rep(1:nrow(inters),each=2),rep(c(2,4),nrow(inters)))
  pol <- round(data.frame(x=inters[idx],y=inters[idy]),6)
  if (any(duplicated(pol))) {
    pol <- pol[-which(duplicated(pol)),]
  }
  if (pol[1,2] != 0){
    pol <- rbind(data.frame(x=0,y=0),pol)
  }
  if (pol[nrow(pol),2] != 0) {
    pol <- rbind(pol,data.frame(x=pol$x[nrow(pol)],y=0))
  }
  return(pol)
}


#' Insert given rows to a given data frame at given locations
#' @param df data frame object where rows will be inserted
#' @param newrows data frame of rows to be inserted
#' @param r numeric vector of indeces that show where the rows should be placed
#' @author Liudas Daumantas
#' @noRd
.insert_filter_Rows <- function(df, newrows, r) {
  df <- rbind(df,newrows)
  df <- df[order(c(1:(nrow(df)-nrow(newrows)),r+0.5)),]
  if (dynGet(x = "x3",ifnotfound = stop("ERROR"),minframe = 1,
             inherits = T) == 0) {
    df = df[df$x >= 0,] } else {
      df = df[df$x <= dynGet(x = "x3",ifnotfound = stop("ERROR"),
                             minframe = 1,inherits = T),]
    }
  row.names(df) <- 1:nrow(df)
  return(df)
}

#' Instead of indicated polygon vertexes, add new vertexes
#' with interpolated y coordinates
#'
#' @description Function removes indicated polygon vertexes. Instead of the
#' removed vertexes it adds new vertexes at provided x coordinates
#' for each mismatching pair of adjacent indicators of removal.
#' @param changes_log logical indicator, showing which vertexes are to be
#' removed from the polygon. Adjacent mismatches (TRUE FALSE and FALSE TRUE)
#' indicate where in the modified polygon vertexes should be added.
#' @param x3 an x coordinate at which y coordinates will be interpolated for
#' added vertexes.
#' @param polygon A data frame containing coordinates of a polygon vertexes.
#' @return  A data frame containing coordinates of vertexes of a modified
#'  polygon
#' @author Liudas Daumantas
#' @noRd
.find_y <- function(changes_log,x3,polygon){
  change_id0 <- which( changes_log[ -1 ] - changes_log[ -length(polygon$x) ]
                       !=0 )
  if (change_id0[1] == 1 & changes_log[1] ) {
    change_id0 <- change_id0[-1]
  } else {
    if (change_id0[ length(change_id0) ] == length(polygon$x) - 1 &
        changes_log[length(polygon$x) ] ) {
      change_id0 <- change_id0[ -length(change_id0) ]
    }
  }

  new_pts <- data.frame(x = rep(x3, length(change_id0)),
                        y = pt_on_line(x1 = polygon$x[change_id0],
                                       x2 = polygon$x[ change_id0 + 1 ],
                                       y1 = polygon$y[change_id0],
                                       y2 = polygon$y[ change_id0 + 1 ],
                                       x3 = rep(x3,length(change_id0)) )
  )


  return( .insert_filter_Rows(polygon, new_pts, change_id0) )
}

#' Remove parts of a polygon that are outside of x range
#' of a split line
#'
#' @description Function removes polygon verteces that are outside of x range
#' of a split line, as well as interpolates and adds vertexes where polygon
#' segments from deleted vertexes enters the allowed x range.
#' @param polygon A data frame containing coordinates of a polygon vertexes, an
#' output of .split_poly function that is rotated so that the split line is
#' horizontal, starting at x=0, y0 and ending at some positive x, y = 0.
#' The provided polygon can be either closed or open.
#' @return A data frame containing the coordinates of vertexes of a provided
#' polygon that are inside an x range of a split line.
#' @author Liudas Daumantas
#' @noRd
.cut_by_x_margins <- function(polygon) {
  pol1 <- .find_y(polygon$x <= 0, x3= 0,polygon)
  return (.find_y(pol1$x >= pol1$x[ length(pol1$x) ],
                  x3 = pol1$x[ length(pol1$x) ], pol1) )
}

#' Filter segments of a polygon that are not OK (out of allowed y range,
#' turns left or are not the nearest to the Y axis)
#'
#' @param intervals a data frame of x and y intervals of a polygon segments.
#' These segments are pre-filtered by the allowed x ranges.
#' @return  A data frame containing coordinates of vertexes of a modified
#'  polygon
#' @author Liudas Daumantas
#' @noRd
.segment_filter <- function(intervals) {
  if (nrow(intervals)>1){
    intervals <- round(intervals,6)
    int_x_ranges <- cbind(t(apply(intervals[,c(1,3)],1,range)),1:nrow(intervals))

    overlap_log <-  apply(int_x_ranges[,1:2], 1,function(o) {
      (o[1] >= int_x_ranges[,1] & o[1] < int_x_ranges[,2] ) |
        ( o[2] > int_x_ranges[,1] & o[2] <= int_x_ranges[,2] ) |
        (o[1] <= int_x_ranges[,1] & o[2] >= int_x_ranges[,2])
    })
    diag(overlap_log) <- FALSE
    noverlap_int_log <- apply(!overlap_log,2,all)
    if (any(noverlap_int_log)) {
      overlap_int_ranges <- int_x_ranges[-which(noverlap_int_log),]
    } else {
      overlap_int_ranges <- int_x_ranges
    }

    if (nrow(overlap_int_ranges)==0) {
      return(.ints_to_pnts(intervals))
    }

    x_end_points <- sort(unique(c(overlap_int_ranges[,1:2])))
    x_ranges <- cbind(x_end_points[-length(x_end_points)],x_end_points[-1])
    x_ranges_int_id <-  apply(x_ranges, 1,function(o) {
      ids <- which((o[1] >= overlap_int_ranges[,1] &
                      o[1] < overlap_int_ranges[,2] ) |
                     ( o[2] > overlap_int_ranges[,1] &
                         o[2] <= overlap_int_ranges[,2] ) |
                     (o[1] <= overlap_int_ranges[,1] &
                        o[2] >= overlap_int_ranges[,2]))
      overlap_int_ranges[ids[which.min(abs(pt_on_line(
        x1 = intervals$x1[overlap_int_ranges[ids,3]],
        x2 = intervals$x2[overlap_int_ranges[ids,3]],
        y1 = intervals$y1[overlap_int_ranges[ids,3]],
        y2 = intervals$y2[overlap_int_ranges[ids,3]],
        x3 = rep(o[1]+(o[2]-o[1])/2,length(ids)))))],3]
    })
    if (is.list(x_ranges_int_id)) { # then this side should not be the one
      # needed. Could be improved by returning indicator and moving on to
      # the next side. If both sides gets indicator, then abort
      return(NULL)
    }
    if (any(noverlap_int_log)) {
      x_ranges_int_id <- c(x_ranges_int_id,which(noverlap_int_log))


      x_ranges_int_id_order <- order(sapply(x_ranges_int_id,
                                            function(o){
                                              min(intervals[
                                                o,c(1,3)])}))
      x_ranges_int_id <- x_ranges_int_id[x_ranges_int_id_order]
      x_ranges <- rbind(x_ranges,t(apply(intervals[which(noverlap_int_log),
                                                   c(1,3)]
                                         ,1,range)))[x_ranges_int_id_order,]

    }

    intervals_selected <- t(apply(intervals[x_ranges_int_id,],1,
                                  function(o){
                                    id <- which.min(c(o[1],o[3]))
                                    if (id == 1 ) o else o[c(3,4,1,2)]
                                  }))

    #intervals_selected <- intervals_selected[paste(x_ranges_int_id),]
    inters <- data.frame(x1 = x_ranges[,1],y1 = pt_on_line(
      x1 = intervals_selected[,1],
      x2 = intervals_selected[,3],
      y1 = intervals_selected[,2],
      y2 = intervals_selected[,4],
      x3 = x_ranges[,1]), x2 = x_ranges[,2],
      y2 = pt_on_line(
        x1 = intervals_selected[,1],
        x2 = intervals_selected[,3],
        y1 = intervals_selected[,2],
        y2 = intervals_selected[,4],
        x3 = x_ranges[,2]))

  } else{
    inters <- intervals
  }
  in.pol <- .ints_to_pnts(inters)

  return(in.pol)
}


#' Remove left-turning segments of a polygon divided with horizontal split line
#'
#' @description This function from a polygon derived with .split_poly
#' function extracts the maximum area part that is directly above or below
#' split line (depending on the which side of divided polygon is provided) and
#' that has no left-turning segments.
#' @param polygon A data frame containing coordinates of a polygon, output of
#' .split_poly function that is rotated so that the split line is horizontal,
#' starting at x=0, y0 and ending at some positive x, y = 0.
#' The provided polygon can be either closed or open.
#' @return A list of two values:
#' [[1]] A data frame containing the coordinates of the extracted polygon part.
#' The returned polygon is open and it is not strictly inside a provided
#' polygon, since it inherits vertical & right-turning segments that are
#' directly above or below the split line.
#' [[2]] -1 or +1 integer indicating the side of the extracted polygon in
#' relation to the split line (+1 above, -1 below)
#' @author Liudas Daumantas
#' @examples #Creating a data.frame of an irregular polygon
#' poly<-data.frame(x=c(0,-1,12,0,5,12,12,7,-5,
#' 11,-3,10,-10,-10,15,15,-1,13,13,-7,-7,0),
#' y=c(0,5,10,4,4,3,-3,-10,-11,-4,-4,0,-2,20,20,-20,-20,-15,15,15,2,0))
#' plot(poly,type='o',main='Horizontal Split Line')
#' #Horizontal Split Line
#' lines(poly[c(1,12),1],poly[c(1,12),2],lty='dotted')
#' #Dividing a polygon
#' upper_half <- split_poly(polygon = poly, split_ids = c(1,12), min_id = 1,
#' trivial_side = F, poli_side = T)
#' lower_half <- split_poly(polygon = poly, split_ids = c(1,12), min_id = 1,
#' trivial_side = F, poli_side = F)
#' lines(upper_half,col=4)
#' lines(lower_half,col=2)
#' #Filtering left turning segments
#' upper_half_lf <- .inner_poly_hull(upper_half)
#' lower_half_lf <- .inner_poly_hull(lower_half)
#' lines(upper_half_lf[[1]],col=4,lwd=3)
#' lines(lower_half_lf[[1]],col=2,lwd=3)
#' @noRd
.inner_poly_hull<-function(polygon){

  if ( all(polygon[-1,1]-polygon[-nrow(polygon),1]>=0 ) ){
    return(
      list(polygon,sign(polygon[2,2]))
    )
  }
  v_cutted_poly <- .cut_by_x_margins(polygon)

  x_shift <- v_cutted_poly[,1][-1] -
    v_cutted_poly[,1][-nrow(v_cutted_poly)]

  if (all(x_shift>=0)){
    return(
      list(v_cutted_poly,sign(v_cutted_poly[2,2]))
    )
  }
  is_x_shift <- x_shift != 0

  if (all(v_cutted_poly$y <=0) | all(v_cutted_poly$y >=0) ) {

    side <- sign(v_cutted_poly$y[2])
    int <- .create_ints(is_x_shift,v_cutted_poly)

    return(
      list(.segment_filter(int),side)
    )

  } else {

    upper_int <- .create_ints(is_x_shift & v_cutted_poly[,2][-1] >=0 &
                                v_cutted_poly[,2][-nrow(v_cutted_poly)] >=0
                              ,v_cutted_poly)

    lowest.up.intervals <- .segment_filter(upper_int)

    if (is.null(lowest.up.intervals)){
      side <- FALSE
    } else {
      side <- .is_pol_up(lowest.up.intervals,polygon)
    }

    if (side) {

      return(
        list(lowest.up.intervals, 1)
      )

    } else {

      lower_int <- .create_ints(is_x_shift & v_cutted_poly[,2][-1] <=0 &
                                  v_cutted_poly[,2][-nrow(v_cutted_poly)] <=0,
                                v_cutted_poly)
      return(
        list(.segment_filter(lower_int),-1)
      )
    }
  }
}