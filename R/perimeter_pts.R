#' Generate a specified number of regularly spaced points on a perimeter of a polygon
#'
#' This function generates regularly spaced points on a perimeter of a polygon. Their number is specified by arguments n.pts or dst.pts.
#' @param polygon A data frame of 2 columns (x,y) that contain coordinates of polygon vertices. Both, closed and open polygons are accepted.
#' @param n.pts A number of regularly spaced points to be generated on a perimeter of a polygon. If n.pts is not specified, then it is calculated according to argument dst.pts.
#' @param dst.pts A distance along a polygon perimeter between adjecent points to be generated on a perimeter of a polygon. If dist is not specified, then it is calculated according to argument n.pts.
#' @return A list of 3 elements. coords data.frame provides coordinates of generated points (X,Y - columns 1, 2) and their ID (column 3).
#' This ID reflects the relative location of a point along a perimeter of a polygon in relation to other generated points and polygon vertices.
#' segment.no indicates the ID of a polygon segment on which a generated point is located. It helps to indentify points located on the same
#' polygon segment. full.poly data.frame contains coordinates of the provided polygon vertices and generated points. coords[,"ID"] can be used to
#' extract rows of generated points.
#'
#' @note If both, n.pts and dst.pts, are specified, then points are generated according to n.pts.
#' @author Liudas Daumantas
#' @examples #Creating data.frame of a polygon
#' poly<- data.frame(X=c(3.38,3.30,1.70,0.78,-0.06,-2.30,-2.94,-3.97,-1.61,-0.39,0.68,1.28,1.60,3.38),
#'                   Y=c(-0.12,-0.31,-2.73,-3.22,-3.29,-2.19,-1.62,0.94,3.10,3.00,2.91,2.49,2.20,-0.12))
#' plot(poly,type='o')
#' #Generating 10 points on a polygon perimeter
#' a<-perimeter_pts(poly,n.pts = 10)
#' #location of points
#' points(a[[1]][,-3],col=2,pch=19)
#' #ID of points
#' text(x=a[[1]][,1],y=a[[1]][,2],a[[1]][,3],-0.3)
#' @export
#'

perimeter_pts<-function (polygon,n.pts=NULL,dst.pts=NULL){
  if (all(polygon[1,]!=polygon[nrow(polygon),])){
    polygon<-close_poly(polygon)
  }
  x<-polygon[,1]
  y<-polygon[,2]
  perimetras<-0
  for( i in 1:c(length(x)-1)){
    perimetras<-perimetras+sqrt((x[i]-x[i+1])^2+(y[i]-y[i+1])^2)
  }
  if (!is.null(n.pts)){
    dst.pts<-perimetras/n.pts
  } else {
    if (!is.null(dst.pts)) {
    n.pts<-round(perimetras/dst.pts,0)
  } else {
    return(print("Specify either n.pts or dst.pts"))
    }
  }
  praeiti<-0
  gaubx_pt<-x[1]
  gauby_pt<-y[1]
  ID<-1
  i=2
  segment.no<-i
  x.polio<-x[1]
  y.polio<-y[1]
  lastx<-gaubx_pt
  lasty<-gauby_pt
  while (length(gaubx_pt)<n.pts){
    #linijinis atstumas tarp dvieju poligono tasku
    m<-sqrt((lastx-x[i])^2+(lasty-y[i])^2)
    if (m+praeiti>=dst.pts){ #jei nueito atstumo palei poligono krastine ilgis ilgesnis nei turetu buti,
      #pasiziurim kokia dali atstumo tarp analizuojamu dvieju tasku  sudaro, sudaro atstumo palei poligono
      #krastine trukumas - si dalis - objektas kof
      kof<-(dst.pts-praeiti)/m
      if (lastx > x[i]) { #jei poligonas suka i kaire, X koordinate yra apibreziama viena formule:
        X<-lastx-(lastx-x[i])*kof
      } else { # jei i desine - kita
        X<-lastx+(x[i]-lastx)*kof
      }
      if (lasty > y[i]){ #tas pats galioja ir Y koordinatei
        Y<-lasty-(lasty-y[i])*kof
      } else {
        Y<-lasty+(y[i]-lasty)*kof
      }
      gaubx_pt<-c(gaubx_pt,X)
      gauby_pt<-c(gauby_pt,Y)
      x.polio<-c(x.polio,X)
      y.polio<-c(y.polio,Y)
      ID<-c(ID,length(x.polio))
      lasty<-Y
      lastx<-X
      segment.no<-c(segment.no,i)
      praeiti<-0
    } else {
      praeiti<-praeiti+m
      x.polio<-c(x.polio,x[i])
      y.polio<-c(y.polio,y[i])
      lasty<-y[i]
      lastx<-x[i]
      i<-i+1
    }
  }

  coords<-data.frame(gaubx_pt,gauby_pt,ID)
  full.poly<-rbind(data.frame(x.polio,y.polio),data.frame(x.polio=x[i:length(x)],y.polio=y[i:length(x)]))
  return(list(coords,segment.no,full.poly))
}
