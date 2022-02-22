#' Find the curve that divides the polygon best (main function)
#'
#' @description function prepares the data to start the searching procedure of the curve that provides the best spatial separation.
#' As curve are generated using splines, matrix of knot coordinates are prepared. Also, polygon is rotated so that split line would be horizontal and at Y = 0, and start at X = 0.
#' @param poly.x a vector of x coordinates of a polygon perimeter points
#' @param poly.y a vector of y coordinates of a polygon perimeter points
#' @param min.x.id index of split line vertex in poly.x and poly.y objects that has lower x coordinate
#' @param max.x.id index of split line vertex in poly.x and poly.y objects that has lower y coordinate
#' @param b slope of a split line
#' @param data data frame of data being analized
#' @param knot.density.X number of spline knots along the split line
#' @param knot.density.Y number of spline knots orthogonal to the split line
#' @param N.condminimum minimum number of fossils required to establish subdivision of a plot
#' @param S.cond minimum area required to establish subdivision of a plot
#' @param n.curve.iter number of curve iterations
#' @param correction.term term that defines how much the a problematic spline
#' will be corrected (in terms of proportion of polygon width where spline
#' intersects the polygon boundary) if the spline is not contained within the
#' plot. Small values recommended (default is 0.05).
#' @return A list of two elements: 1) curve in shape of a spline that produces the best data separation; 2) quality of the division
#' @author Liudas Daumantas
#' @importFrom DescTools Rotate
#' @noRd
.curvial_split<-function(poly.x,poly.y,min.x.id,max.x.id,b,
                         data,knot.density.X=20,knot.density.Y=20,
                         N.cond,S.cond,n.curve.iter,
                         correction.term){
  #nustatau, kuris padalinimo linijos id yra kairej, kuris desinej
  # length of a split line
  AE<-sqrt(sum((c(poly.x[min.x.id],
                  poly.y[min.x.id])-c(poly.x[max.x.id],poly.y[max.x.id]))^2))
  #randu kampa, kuriuo reikia paversti poligona, kad padalinimo linija butu
  # horizontali
  teta<-atan(b)
  #paverciu poligona ir duomenis
  rot.pol.cords <- Rotate(x=poly.x,y=poly.y,mx=poly.x[min.x.id],
                          my=poly.y[min.x.id],theta=-teta)
  rot.dat.cords <- Rotate(x=data$x,y=data$y,mx=poly.x[min.x.id],
                          my=poly.y[min.x.id],theta=-teta)
  #pastumiu duomenis ir poligona, kad padalinimo linijos kairinis taskas butu
  # koordinaciu sistemos pradzioje
  rot.data<-data.frame(data[,1],rot.dat.cords$x-poly.x[min.x.id],
                       rot.dat.cords$y-poly.y[min.x.id])
  #sukuriu duomenu masyva pasukto poligono
  rot.poli<-data.frame(x= rot.pol.cords$x-poly.x[min.x.id],
                       y= rot.pol.cords$y-poly.y[min.x.id])
  #padalinimo linijos pabaigoje y turetu buti 0, bet del paklaidu jis gali buti
  # nelygus 0 ir filtracija gali blogai ivykti - nunulinu
  rot.poli[max.x.id,2]<-0

  #"atidarau" poligona
  # rot.poli<-rot.poli[-nrow(rot.poli),]
  #nufiltruoju virsutine ir apatine poligono dali, t
  #askai tampa organizuoti palei poligono kreive is kaires i desine
  rot.poli.up<-.close_poly(.split_poly(polygon = rot.poli, min_id = 1,
                                       split_ids = c(min.x.id,max.x.id),
                                       trivial_side = TRUE,poli_side = TRUE))
  rot.poli.do<-.close_poly(.split_poly(polygon = rot.poli, min_id = 1,
                                       split_ids = c(min.x.id,max.x.id),
                                       trivial_side = TRUE,poli_side = FALSE))

  #randu vidinio pataisyto poligono apatine ir virsutine dali
  bottom.inner.poli<-.inner_poly_hull(polygon = rot.poli.do)
  upper.inner.poli<-.inner_poly_hull(polygon = rot.poli.up)
  if( bottom.inner.poli[[2]] == 1 ){
    change <- bottom.inner.poli
    bottom.inner.poli<-upper.inner.poli
    upper.inner.poli<-change
    change<-rot.poli.do
    rot.poli.do<-rot.poli.up
    rot.poli.up<-change
  }
  #pasirupinam splino knot'ais pradiniais

  split.line.x<-seq(0,AE,length.out = knot.density.X)
  split.line.x[knot.density.X]<-AE
  #ant pataisyto poligono reikia rasti Y koordinates duotom x koordinatem,
  # atitinkanciom padalinimo linijos atkarpu galus
  up.y<-numeric(knot.density.X)
  do.y<-numeric(knot.density.X)
  for(i in 1:length(split.line.x)){
    up.y[i]<-.y.online(x=upper.inner.poli$x,y=upper.inner.poli$yp,
                       x3=split.line.x[i])
    do.y[i]<-.y.online(x=bottom.inner.poli$x,y=bottom.inner.poli$yp,
                       x3=split.line.x[i])
  }
  #randam ploti pataisyto poligono
  range<-up.y-do.y
  #nustatom i kiek daliu dalinsim polygona pagal Y koordinate
  B<-seq(0,1,length.out = knot.density.Y)
  #sugeneruojam matrica su testuojamu spline knotu y koordinatem.
  #Kiekvienas stulpelis atitinka skirtinga X verte ant padalinimo linijos
  #spit.line.x, kiekviena eilute atitinka skirtinga poligono pjuvio ties x
  # koordinate dali - virsutines eilutes apatine poligono dalis
  knot.y.matrix<-matrix(B,knot.density.Y,1)%*%matrix(range,1,knot.density.X)+
    matrix(rep(do.y,each=knot.density.Y),knot.density.Y,knot.density.X)
  #randam geriausia padalinimo kreive
  best..curvi_split<-.curvi_split(
    knot.y.matrix,split.line.x,Xup=upper.inner.poli$x,
    Xdown=bottom.inner.poli$x,Yup=upper.inner.poli$y,Ydown=bottom.inner.poli$y,
    N.cond,S.cond,knot.density.Y,knot.density.X,rot.poli.up,rot.poli.do,
    rot.data,n.curve.iter =n.curve.iter,correction.term = correction.term
    )

# Up --> Yup, Down --> Ydown; xp --> x; yp --> y.  visur pakeist
  best.curve<-data.frame(x=best..curvi_split[[1]]$x+poly.x[min.x.id],
                         y=best..curvi_split[[1]]$y+poly.y[min.x.id])
  best.curve<-Rotate(x=best.curve$x,y=best.curve$y,mx=poly.x[min.x.id],
                     my=poly.y[min.x.id],theta=teta)
  return(list(best.curve,best..curvi_split[[2]]))
}

#' Find the curve that divides the polygon best (recursive function)
#'
#' @description At this stage of an algorithm, polygon is rotated so that the split line is horizontal and positioned along X axis (at Y = 0), and split line starts at X = 0.
#' Function iterates trough knot matrix several times (number of times is defined by n.curve.iter argument) and tries different splines to seperate data in space.
#' Iteration through knots starts from left to right along the split line. For each position along the split line several positions orthogonal to the split line are tried.
#' The best orthogonal positions are estimated from a spline model (Split Quality ~ Y coordinate of a knot for a given X coordinate)
#' So the best Y knot coordinate is provided for each knot along X. Thus, knots along the split line are fixated, but along Y are not. When one iteration is complete, the next
#' one starts in different direction (e.g. from right to left along the split line).
#'
#' @param knot.y.matrix a matrix with Y coordinates in columns for each knot along the split line
#' @param split.line.x a vector of x coordinates for each knot along the split line
#' @param Xup X coordinates of upper part (above split line) of standartized polygon
#' @param Xdown X coordinates of lower part (below split line) of standartized polygon
#' @param Yup Y coordinates of upper part (above split line) of standartized polygon
#' @param Ydown Y coordinates of lower part (below split line) of standartized polygon
#' @param N.cond minimum number of fossils required to establish subdivision of a plot
#' @param S.cond minimum area required to establish subdivision of a plot
#' @param knot.density.Y number of spline knots orthogonal to the split line
#' @param knot.density.X number of spline knots along the split line
#' @param rot.poli.up rotated, but not standartized, full upper polygon
#' @param rot.poli.do rotated, but not standartized, full lower polygon
#' @param rot.data rotated data that is analyzed
#' @param n.curve.iter number of curve iterations
#' @param correction.term term that defines how much the a problematic spline
#' will be corrected (in terms of proportion of polygon width where spline
#' intersects the polygon boundary) if the spline is not contained within the
#' plot. Small values recommended (default is 0.05).
#' @return A list of two elements: 1) rotated (not suitable for the original polygon) curve in shape of a spline that produces the best data separation; 2) quality of the division
#' @author Liudas Daumantas
#' @noRd
.curvi_split<-function(knot.y.matrix,split.line.x,Xup,Xdown,Yup,Ydown,N.cond,
                      S.cond,knot.density.Y,knot.density.X,rot.poli.up,
                      rot.poli.do,rot.data,n.curve.iter=n.curve.iter,
                      correction.term){
  best.y.knots<-numeric(knot.density.X)
  SSk<-numeric(knot.density.Y-2)
  SSks<-numeric(knot.density.X-2)
  knot.y.matrix<-knot.y.matrix[-c(1,knot.density.Y),-c(1,knot.density.X)]
  for (it in seq(n.curve.iter)){
    if(it%%2==1){
      if(it==1){
        seq.of.x.knots <- 1:(knot.density.X-2)
      } else {
        seq.of.x.knots <- 2:(knot.density.X-2)
      }} else {
        seq.of.x.knots <- (knot.density.X-3):1
      }
    for (l in seq.of.x.knots) {
      for (k in 1:(knot.density.Y-2)) {
        #isbandom vis kita Y koordinate knotui su duota x koordinate
        best.y.knots[1+l]<-knot.y.matrix[k,l]
        #sugeneruojam splina per knotus
        curve<-spline(split.line.x,best.y.knots,n=1000)
        #iteratyviai darom, kol nebemeta klaidos:
        #patikrinam ar poligono viduj kreive, gaunam pataisos tasku koordinates
        #skeliam split.line.x ir best.y.knots i tiek daliu, kiek reikia ir
        # reikiamose vietose pridedam koordinates
        #sukuriam nauja split.line.x best.y.knots varianta
        #kartojam splina
        curve<-.spline_corrections(curve,Xup,Xdown,Yup,Ydown,best.y.knots,
                                 split.line.x,correction.term=correction.term)
        #ivertinam poligono padalinimo su sugeneruota kreive kokybe
        SS<-.curve_quality(curve=curve,rot.poli.up, rot.poli.do, rot.data,
                          N.cond,S.cond)
        #fiksuojam verte
        SSk[k]<-SS
      }
      #sugeneruojam spline, kurio X asyje yra knotu Y koordinates, o Y asyje
      #padalinimo kreives kokybe
      #skirta interpoliuojant aproksimuoti tarpiniu Y koordinaciu vertes
      proj<-spline(knot.y.matrix[,l],SSk,n=1000)
      SSks[l]<-proj$y[which.max(proj$y)]
      #issaugom, kaip geriausia Y koordinate X koordinates knotui ta, kuri turi
      # didziausia kokybe
      best.y.knots[1+l]<-proj$x[which.max(proj$y)]
    }
  }
  #siame etape kiekvienam X koordinates knotui turime po geriausia Y koordinate
  #sugeneruojam per siuos knotus kreive
  curve.final<-spline(split.line.x,best.y.knots,n=1000)
  #kreive gali islisti is uz polygono (pvz. kokybes dideli iverciai buvo
  # gauti pataisius duotus knotus)
  #taigi, kreive dar karta bandoma pataisyti, jei reikia
  curve.final<-.spline_corrections(curve.final,Xup,Xdown,Yup,Ydown,best.y.knots,
                                 split.line.x,correction.term = correction.term)
  #ivertinam galutines padalinimo kreives kokybe
  lines(curve.final,col=3)
  SS<-.curve_quality(curve.final,rot.poli.up, rot.poli.do, rot.data, N.cond,
                    S.cond)
  #grazinam padalinimo kreive, jos kokybes iverti ir visu kitu lokaliai
  # geriausiu kreiviu ivercius - SSk
  #SSk tik tam, kad pasizeti, ar tikrai grizta pati geriausia kreive
  if (any(SSks>SS)){
    warning("Intermediate curves performed better than the final curve")
  }
  return(list(curve.final,SS))
}

#' Evaluate the quality of curve in terms of data spatial separation
#'
#' @description function forms two polygons from a provided curve (upper and lower),
#' filters data using each polygon and compares them.
#'
#' @param curve a matrix with Y coordinates in columns for each knot along the split line
#' @param rot.poli.up rotated, but not standartized, full upper polygon
#' @param rot.poli.do rotated, but not standartized, full lower polygon
#' @param rot.data rotated data that is analyzed
#' @param N.cond minimum number of fossils required to establish subdivision of a plot
#' @param S.cond minimum area required to establish subdivision of a plot
#' @return A list of two elements: 1) rotated (not suitable for the original polygon) curve in shape of a spline that produces the best data separation; 2) quality of the division
#' @author Liudas Daumantas
#' @importFrom sp point.in.polygon
#' @importFrom pracma polyarea
#' @noRd
.curve_quality<-function(curve,rot.poli.up,rot.poli.do,rot.data,N.cond,S.cond){
  #sudarom poligonus padalintus kreive
  I.poli<-data.frame(x=c(rot.poli.up$xp,rev(curve$x)[-1]),y=c(rot.poli.up$yp,rev(curve$y)[-1]))
  II.poli<-data.frame(x=c(rot.poli.do$xp,rev(curve$x)[-1]),y=c(rot.poli.do$yp,rev(curve$y)[-1]))
  #atrenkam taskus patenkancius i skirtingas poligono dalis, padalintas kreive
  I.poli.data<-rot.data[which(sp::point.in.polygon(point.x = rot.data$X.rot,point.y =rot.data$Y.rot ,pol.x = I.poli$x,pol.y =I.poli$y)!=0),]
  II.poli.data<-rot.data[which(sp::point.in.polygon(point.x = rot.data$X.rot,point.y = rot.data$Y.rot,pol.x =II.poli$x ,pol.y =II.poli$y)!=0),]
  #atliekam plotu palyginima, t.y. padalinimo gerumo ivertinima, jei kreive netenkina minimaliu kriteriju - padalinimo kokybe nuline
  S1<-abs(pracma:polyarea(I.poli$x,I.poli$y))
  S2<-abs(pracma:polyarea(II.poli$x,II.poli$y))
  if(min(nrow(I.poli.data),nrow(II.poli.data))>N.cond&min(S1,S2)>S.cond){
    SS<-.dif_fun(I.poli.data[,1],II.poli.data[,1])
  } else{
    SS<-0
  }
  SS
}

#' Correct the curve
#'
#' @description function corrects the curve if it is not contained within the polygon.
#' It does so by finding intervals of a curve that cross the boundaries of a polygon,
#' then it adds additional knots at the middle of these intervals inside the polygon.
#' Finally, curve is generated again using new knots, but it may still have some outlying intervals.
#' Thus, the procedure is repeated recursively until no outlying intervals are left.
#' @param curve a curve (output of spline function from stats package)
#' @param Xup X coordinates of upper part (above split line) of standartized polygon
#' @param Xdown X coordinates of lower part (below split line) of standartized polygon
#' @param Yup Y coordinates of upper part (above split line) of standartized polygon
#' @param Ydown Y coordinates of lower part (below split line) of standartized polygon
#' @param best.y.knots a vector of Y coordinates of knots of a spline that were used to produce the curve.
#' Each Y coordinate is dedicated for corresponding elements of split.line.x, that is X coordinates of spline knots.
#' @param split.line.x a vector of x coordinates for each knot along the split line
#' @param correction.term term that defines how much the a problematic spline
#' will be corrected (in terms of proportion of polygon width where spline
#' intersects the polygon boundary) if the spline is not contained within the
#' plot. Small values recommended (default is 0.05).
#' @return a curve (output of a spline function from stats package)
#' @author Liudas Daumantas
#' @noRd
.spline_corrections<-function(curve,Xup,Xdown,Yup,Ydown,best.y.knots,
                             split.line.x,correction.term){

  corrected.coords<-.y_corrections(spline.x = curve$x,spline.y =curve$y ,
                                  Xup = Xup,Xdown = Xdown,Yup = Up,
                                  Ydown = Ydown,
                                  correction.term = correction.term)
  if (length(dim(corrected.coords))>1){
    ind<-numeric()
    k<-0
    corrected.split.line.x<-split.line.x
    corrected.best.y.knots<-best.y.knots
    for (a in seq(length(corrected.coords[,1]))){
      for (i in seq(length(split.line.x))){
        if(split.line.x[i]>=corrected.coords[a,1]){
          ind<-i-1
          if(ind==0){
            if(k==0) {
              corrected.split.line.x<-c(corrected.coords[a,1],
                                        corrected.split.line.x)
              corrected.best.y.knots<-c(corrected.coords[a,2],
                                        corrected.best.y.knots)
            } else {
              corrected.split.line.x<-c(corrected.split.line.x[1:k],
                                        corrected.coords[a,1],
                                        corrected.split.line.x[
                                          c((k+1):
                                              length(corrected.split.line.x))])
              corrected.best.y.knots<-c(corrected.best.y.knots[1:k],
                                        corrected.coords[a,2],
                                        corrected.best.y.knots[
                                          c((k+1):
                                              length(corrected.best.y.knots))])
            }} else {
              if (ind==length(split.line.x)){
                corrected.split.line.x<-c(corrected.split.line.x,
                                          corrected.coords[a,1])
                corrected.best.y.knots<-c(corrected.best.y.knots,
                                          corrected.coords[a,2])
              } else {
                corrected.split.line.x<-c(corrected.split.line.x[c(1:(ind+k))],
                                          corrected.coords[a,1],
                                          corrected.split.line.x[
                                            c((ind+k+1):
                                                length(corrected.split.line.x))
                                            ])
                corrected.best.y.knots<-c(corrected.best.y.knots[c(1:(ind+k))],
                                          corrected.coords[a,2],
                                          corrected.best.y.knots[
                                            c((ind+k+1):
                                                length(corrected.best.y.knots))
                                            ])
              }}
          k<-k+1
          break
        }
      }
    }
    curve<-spline(corrected.split.line.x,corrected.best.y.knots,n=1000)
    best.y.knots<-corrected.best.y.knots
    split.line.x<-corrected.split.line.x
    return(.spline_corrections(curve =curve ,Xup =Xup ,Xdown =Xdown ,Yup=Yup,
                              Ydown = Ydown,best.y.knots=best.y.knots,
                              split.line.x = split.line.x,
                              correction.term = correction.term) )
  } else {
    return(curve)
  }
}



#' Find coordinates for new knots
#'
#' @description function checks if there are curve intervals that cross the boundary of a polygon.
#' If present, for each interval it calculates the middle X coordinate.
#' Then, Y coordinates are produced for each X coordinate in a way so that each produced Y are within the polygon.
#' The Y coordinates are calculated as follows: 1) width of the polygon is calculated
#' at X coordinate (distance in Y coord. between upper and lower boundaries of the
#' standartized polygon at X coord.); 2) the correction.term is then multiplied by the calculated width
#' to get the correction size in Y units; 3) correction is then added to the Y of lower polygon boundary
#' if curve crosses the boundary there, or is subtracted from the upper boundary Y if else.
#'
#' @param spline.x x coordinates of a curve
#' @param spline.y y coordinates of a curve
#' @param Xup X coordinates of upper part (above split line) of standartized polygon
#' @param Xdown X coordinates of lower part (below split line) of standartized polygon
#' @param Yup Y coordinates of upper part (above split line) of standartized polygon
#' @param Ydown Y coordinates of lower part (below split line) of standartized polygon
#' @param best.y.knots a vector of Y coordinates of knots of a spline that were used to produce the curve.
#' Each Y coordinate is dedicated for corresponding elements of split.line.x, that is X coordinates of spline knots.
#' @param correction.term term that defines how much the a problematic spline
#' will be corrected (in terms of proportion of polygon width where spline
#' intersects the polygon boundary) if the spline is not contained within the
#' plot. Small values recommended (default is 0.05).
#' @return If there are curve intervals that cross the boundary of a polygon,
#' then function returns a data frame of x and y coordinates that should be combined with coordinates of other knots used to produce the curve.
#' If no such intervals are present function returns "FALSE".
#' @author Liudas Daumantas
#' @importFrom sp point.in.polygon
#' @noRd
.y_corrections<-function(spline.x,spline.y,Xup,Xdown,Yup,Ydown,
                        correction.term){
  #tikrinam, ar kreive polygone, nuimam po viena taska krastini, nes jie ant poligono ribos
  curve.in.polygon<-sp::point.in.polygon(
    point.x = spline.x,
    point.y = spline.y,
    pol.x = c(Xup,rev(Xdown)),
    pol.y = c(Up,rev(Down))
    )[-c(1,length(spline.x))]

  #jei kreive ne polygone
  #v2 jei neatsivelgti,kuris taskas labiausiai nutoles nuo poligono
  if (any(curve.in.polygon!=1)) {
    #pasiziurim, kurie kreives taskai nepatenka i poligona
    ind<-which(curve.in.polygon!=1)
    #patikrinam, kiek yra ind verciu, jei viena, reiska, tik vienas taskas iseina is poligono, isisaugom taska
    if (length(ind)==1){
      vid.x<-spline.x[ind+1]
      vid.y<-spline.y[ind+1]
    } else { # JEI NE - gali buti kelios kreives atkarpos, nepatenkancios i poligona
      #surandam centrus intervalu, kurie nepatenka i poligona
      min.id<-ind[1]
      vid.x<-numeric()
      vid.y<-numeric()
      for (i in 1:(length(ind)-1)){
        #jei randam luzi indeksuose, reiskia priejom naujo intervalo pradzia,
        #fiksuojam buvusio intervalo viduri ir atnaujiname intervalo pradzia
        if(ind[i+1]-ind[i]>1){
          x3<-spline.x[min.id+1]+(spline.x[ind[i]+1]-spline.x[min.id+1])/2
          vid.x<-c(vid.x,x3)
          #y3 ant splino
          y3.spline<-.y.online(x=spline.x,y=spline.y,x3=x3)
          vid.y<-c(vid.y,y3.spline)
          #fiksuojam naujo intervalo pradzia
          min.id<-ind[i+1]}
      }
      # jei naujo intervalo neprieita, reiskia yra vienas intervalas - issisaugom jo centra
      if (min.id==ind[1]){
        if(ind[1]==1){
          vid.x<-spline.x[1]+(spline.x[ind[length(ind)]+1]-spline.x[1])/2
        } else {
          vid.x<-spline.x[ind[1]+1]+(spline.x[ind[length(ind)]+1]-spline.x[ind[1]+1])/2
        }
        vid.y<-.y.online(x=spline.x,y=spline.y,x3=vid.x)
      } else {# jei rasta intervalu, tai ufiksuoti ju viduriai, bet paskutinio intervalo vidurys neuzfiksuotas
        x3<-spline.x[min.id+1]+(spline.x[ind[length(ind)]+1]-spline.x[min.id+1])/2
        vid.x<-c(vid.x,x3)
        #y3 ant splino
        y3.spline<-.y.online(x=spline.x,y=spline.y,x3=x3)
        vid.y<-c(vid.y,y3.spline)
      }
    }
    #siame etape turim tasku koordinates, kuriuos reikia "nustumti" i poligona ir pakartoti splina
    #poligonas skaidytas i dvi dalis, zemiau 0 ir virs 0.
    #surandam y3 ant poligono virsaus ir apacios, x3 vietoje
    y3.up<-numeric()
    y3.down<-numeric()
    RANGES<-numeric()
    corrected.y<-numeric()
    for (i in 1:length(vid.x)){
      y3.up<-c(y3.up,y.online(x=Xup,y=Yup,x3=vid.x[i]))
      y3.down<-c(y3.down,y.online(x3=vid.x[i],x=Xdown,y=Down))
      #surandam atstuma tarp y3 up ir down
      RANGES<-c(RANGES,y3.up[i]-y3.down[i])
      #paziurim, per kuria puse kreives atkarpa iseina uz poligono (per virsu ar ne?)
      #pataisom y3 taip, kad jis atsidurtu poligone, 5% atstumu nuo krasto poligono
      if (vid.y[i]>=0){
        corrected.y<-c(corrected.y,y3.up[i]-correction.term*RANGES[i])
      } else {
        corrected.y<-c(corrected.y,y3.down[i]+correction.term*RANGES[i])
      }
    }
    return(data.frame(vid.x,corrected.y))
  } else {
    return(FALSE)}
}
