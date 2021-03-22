#' Split a polygon
#'
#' This function splits a polygon in two halfs by a split line that links two provided polygon vertices and returns one half of this polygon.
#'
#' @param polygon A data frame containing coordinates of a polygon that will be split. The polygon can be either closed or open.
#' @param split_ids A vector containing 2 integers. These integers should be indeces of polygon data frame rows that contain coordinates
#' of vertices that are used to split a polygon.
#' @param min_id (Optional) Integer (1 or 2) showing which of the polygon vertices indicated by the \code{split_ids} has lower x coordinate
#' (or lower y coordinate, if split line is vertical).
#' @param trivial_side A logical argument indicating whether it matters which side of the divided polygon will be returned.
#' @param poli_side A logical argument. When \code{trivial_side} is \code{FALSE}, \code{poli_side = TRUE} returns the upper half
#' (or the left half, if split line is vertical) of a divided polygon and \code{poli_side = FALSE} returns the lower (the right) half.
#' @return A data frame containing coordinates of one half of a polygon divided by the split line. The resulting polygon is open.
#' @details  When \code{trivial_side = TRUE}, computations are much faster. The duration of computations also depend on the complexity of
#' polygons - they are longer when none of the sides of the divided polygon has all vertices on one side of the split line. The computations
#' are also faster, when split line is horizontal or vertical. The use of \code{min_id} only marginally reduces the amount of calculations.
#' @author Liudas Daumantas
#' @examples #Creating a data.frame of an irregular polygon
#' poly<-data.frame(X=c(0,-1,12,0,5,12,12,7,-5,
#' 11,-3,10,-10,-10,15,15,-1,13,13,-7,-7,0),
#' Y=c(0,5,10,4,4,3,-3,-10,-11,-4,-4,0,-2,20,20,-20,-20,-15,15,15,2,0))
#' plot(poly,type='o',main='Horizontal Split Line')
#' #Horizontal Split Line
#' lines(poly[c(1,12),1],poly[c(1,12),2],lty='dotted')
#' #Dividing a polygon
#' upper_half <- split_poly(polygon = poly, split_ids = c(1,12), min_id = 1, trivial_side = FALSE, poli_side = TRUE)
#' lower_half <- split_poly(polygon = poly, split_ids = c(1,12), min_id = 1, trivial_side = FALSE, poli_side = FALSE)
#' lines(upper_half,col=4)
#' lines(lower_half,col=2)
#'
#'
#' plot(poly,type='o',main='Vertical Split Line')
#' #Vertical Split Line
#' lines(poly[c(1,4),1],poly[c(1,4),2],lty='dotted')
#' #Dividing a polygon
#' left_half <- split_poly(polygon = poly, split_ids = c(1,4), min_id = 1, trivial_side = FALSE, poli_side = TRUE)
#' right_half <- split_poly(polygon = poly, split_ids = c(1,4), min_id = 1, trivial_side = FALSE, poli_side = FALSE)
#' lines(left_half,col=4)
#' lines(right_half,col=2)
#' print(FALSE)
#' plot(poly,type='o',main='Inclined Split Line')
#' #Inclined Split Line
#' lines(poly[c(1,6),1],poly[c(1,6),2],lty='dotted')
#' #Dividing a polygon
#' upper_half <- split_poly(polygon = poly, split_ids = c(1,6), min_id = 1, trivial_side = FALSE, poli_side = TRUE)
#' lower_half <- split_poly(polygon = poly, split_ids = c(1,6), min_id = 1, trivial_side = FALSE, poli_side = FALSE)
#' lines(upper_half,col=4)
#' lines(lower_half,col=2)
#' @export

split_poly<-function(polygon,split_ids,min_id=NULL, trivial_side=TRUE, poli_side){
  if (is.null(min_id)){
    if (polygon[split_ids[1],1]==polygon[split_ids[2],1]){
      min_id<-which.min(c(polygon[split_ids[1],2],polygon[split_ids[2],2]))
    } else {
      min_id<-which.min(c(polygon[split_ids[1],1],polygon[split_ids[2],1]))
    }
  }
  min_x_id<-split_ids[min_id]
  max_x_id<-split_ids[-min_id]
  ids<-min_x_id:max_x_id
  if (all(polygon[1,]!=polygon[nrow(polygon),])){
    polygon<-close_poly(polygon)
  }
  if (trivial_side==T){
    case <- TRUE
  } else {
    if (polygon[min_x_id,2]==polygon[max_x_id,2]){
      zero <- polygon[min_x_id,2]
      if ((all(polygon[ids,2]>=zero) | all(polygon[ids,2]<=zero)) | (all(polygon[-ids,2]>=zero) | all(polygon[-ids,2]<=zero))) {
        if (all(polygon[ids,2]>=zero) | all(polygon[ids,2]<=zero)){
          case<-all(polygon[ids,2]>=zero)
        } else {
          case<-all(polygon[-ids,2]<=zero)
        }
      } else {
        pol<-polygon[ids,]
        direction<-numeric()
        side<-numeric()
        mid_x<-polygon[max_x_id,1]/2
        i=1
        while (i < nrow(pol)){
          if (pol[i,1]>polygon[min_x_id,1] & pol[i,1] < polygon[max_x_id,1] & i>1){
            start=i-1
            i=i+1
            while (pol[i,1]>polygon[min_x_id,1] & pol[i,1] < polygon[max_x_id,1] & i < nrow(pol)){
              i=i+1
            }
            end=i
            i=i-1
            if (pol[start,1]<=polygon[min_x_id,1] & pol[end,1]>=polygon[max_x_id,1]) {
              direction <- c(direction,1)
            } else {
              if (pol[start,1]>=polygon[max_x_id,1] & pol[end,1]<=polygon[min_x_id,1]){
                direction <- c(direction,-1)
              } else {direction <- c(direction,0)}
            }

            side <- c(side,pol[i,2]>zero)
            i=end
          } else {
            if ((pol[i,1]<=polygon[min_x_id,1] | pol[i+1,1]<=polygon[min_x_id,1]) &
                (pol[i,1]>=polygon[max_x_id,1] | pol[i+1,1]>=polygon[max_x_id,1])){
              side <- c(side,pt_on_line(x1=pol[i,1],x2=pol[i+1,1],y1=pol[i,2],y2=pol[i+1,2],x3=mid_x)>zero)
              if (pol[i,1]<pol[i+1,1]){
                direction<-c(direction, 1)
              } else{
                direction<-c(direction, -1)
              }
              i=i+1
            } else{
              i=i+1
            }
          }
        }
        case<-sum(direction[side==1])==1
      }
    } else {
      if (polygon[min_x_id,1]==polygon[max_x_id,1]){
        zero <- polygon[min_x_id,1]
        if ((all(polygon[ids,1]>=zero) | all(polygon[ids,1]<=zero)) | (all(polygon[-ids,1]>=zero) | all(polygon[-ids,1]<=zero))){
          if (all(polygon[ids,1]>=zero) | all(polygon[ids,1]<=zero)){
            case<-all(polygon[ids,1]<=zero)
          } else {
            case<-all(polygon[-ids,1]>=zero)
          }
        } else {
          pol<-polygon[ids,]
          direction<-numeric()
          side<-numeric()
          mid_y<-max(polygon[c(max_x_id,min_x_id),2])/2
          i=1
          while (i < nrow(pol)){
            if (pol[i,2]>polygon[min_x_id,2] & pol[i,2] < polygon[max_x_id,2] & i>1){
              start=i-1
              i=i+1
              while (pol[i,2]>polygon[min_x_id,2] & pol[i,2] < polygon[max_x_id,2]){
                i=i+1
              }
              end=i
              i=i-1
              if (pol[start,2]<=polygon[min_x_id,2] & pol[end,2]>=polygon[max_x_id,2]) {
                direction <- c(direction,1)
              } else {
                if (pol[start,2]>=polygon[max_x_id,2] & pol[end,2]<=polygon[min_x_id,2]){
                  direction <- c(direction,-1)
                } else {direction <- c(direction,0)}
              }

              side <- c(side,pol[i,1]<zero)
              i=end
            } else {
              if ((pol[i,2]<=polygon[min_x_id,2] | pol[i+1,2]<=polygon[min_x_id,2]) &
                  (pol[i,2]>=polygon[max_x_id,2] | pol[i+1,2]>=polygon[max_x_id,2])){
                side <- c(side,pt_on_line(x1=pol[i,1],x2=pol[i+1,1],y1=pol[i,2],y2=pol[i+1,2],y3=mid_y)<zero)
                if (pol[i,2]<pol[i+1,2]){
                  direction<-c(direction, 1)
                } else{
                  direction<-c(direction, -1)
                }
                i=i+1
              } else{
                i=i+1
              }
            }
          }
          case<-sum(direction[side==1])==1
        }
      }  else {
        b<-(polygon[max_x_id,2]-polygon[min_x_id,2])/(polygon[max_x_id,1]-polygon[min_x_id,1])
        b1<-tan((atan(b)*180/pi+90)*pi/180)
        y_pts_on_line<-polygon[min_x_id,2]+b*polygon[ids,1]
        y2_pts_on_line<-polygon[min_x_id,2]+b*polygon[-ids,1]
        if ((all(polygon[ids,2]>=y_pts_on_line) | all(polygon[ids,2]<=y_pts_on_line)) |
            (all(polygon[-ids,2]>=y2_pts_on_line) | all(polygon[-ids,2]<=y2_pts_on_line)) ) {
          if (all(polygon[ids,2]>=y_pts_on_line) | all(polygon[ids,2]<=y_pts_on_line)) {
            case <- all(polygon[ids,2]>=y_pts_on_line)
          } else {
            case <- all(polygon[-ids,2]<=y2_pts_on_line)
          }
        } else {
          pol<-polygon[ids,]
          direction<-numeric()
          side<-numeric()
          mid_x<-polygon[min_x_id,1]+(polygon[max_x_id,1]-polygon[min_x_id,1])/2
          mid_y<-polygon[min_x_id,2]+b*mid_x
          i=1
          y_b1_pts<-polygon[min_x_id,2]+b1*pol[,1]
          y_b2_pts<-polygon[max_x_id,2]+b1*pol[,1]
          while (i < nrow(pol)){
            if (pol[i,2]>y_b1_pts[i] & pol[i,2] < y_b2_pts[i] & i>1){
              start=i-1
              i=i+1
              while (pol[i,2]>y_b1_pts[i] & pol[i,2] < y_b2_pts[i] & i < nrow(pol)){
                i=i+1
              }
              end=i
              i=i-1
              if (pol[start,2]<=y_b1_pts[i] & pol[end,2]>=y_b2_pts[i]) {
                direction <- c(direction,1)
              } else {
                if (pol[start,2]>=y_b2_pts[i] & pol[end,2]<=y_b1_pts[i]){
                  direction <- c(direction,-1)
                } else {direction <- c(direction,0)}
              }

              side <- c(side,pol[i,2]>y_pts_on_line[i])
              i=end
            } else {
              if ((pol[i,2]<=y_b1_pts[i] | pol[i+1,2]<=y_b1_pts[i+1]) &
                  (pol[i,2]>=y_b2_pts[i] | pol[i+1,2]>=y_b2_pts[i+1])){
                side <- c(side,pt_on_line(x1=pol[i,1],x2=pol[i+1,1],y1=pol[i,2],y2=pol[i+1,2],x3=mid_x)>mid_y)
                if (pol[i,2]<=y_b1_pts[i]){
                  direction<-c(direction, 1)
                } else{
                  direction<-c(direction, -1)
                }
                i=i+1
              } else{
                i=i+1
              }
            }
          }
          case<-sum(direction[side==1])==1
        }
      }
    }
  }
  if (case){
    if (poli_side==T){
      return(polygon[ids,])
    } else {
      polygon<-cbind(polygon,1:nrow(polygon))
      polygon<-polygon[-ids[-c(1,length(ids))],]
      polygon<-polygon[-nrow(polygon),]
    }
  } else {
    if (poli_side==F){
      return(polygon[ids,])
    } else {
      polygon<-cbind(polygon,1:nrow(polygon))
      polygon<-polygon[-ids[-c(1,length(ids))],]
      polygon<-polygon[-nrow(polygon),]
    }
  }
  if (min_x_id<max_x_id){
    return(polygon[match(c(min_x_id:1,polygon[nrow(polygon),3]:max_x_id),polygon[,3]),1:2])
  } else {
    return(polygon[match(c(min_x_id:polygon[nrow(polygon),3],1:max_x_id),polygon[,3]),1:2])
  }
}
