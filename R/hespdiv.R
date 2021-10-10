#' Hierarchically subdivide spatial data
#'
#' This is the main function of HespDiv package that performs hierarchical
#' spatial data subdivision. It recursively divides the data in space using
#' random split-lines and evaluates their performance by estimating how well
#' have they separated the data in space in terms of the magnitude of difference
#' in emergent data qualities of interest. The best performing straight
#' split-line is then used as an axis around which various shape curves are
#' created and tested in a similar manner. Thus, each recursive iteration can
#' produce either a straight or curvi-linear split-line that divides the study
#' area and data into two parts. Each part can be further subdivided in
#' the subsequent iterations to ultimately produce a hierarchical subdivision of
#' space and data.
#'
#' Two main functions must be provided to estimate the separation of data in
#' space. First one (argument = \code{generalize.f}) is needed to calculate some
#' emergent data quality (eg. some model, summary statistic, etc.).
#' The second one (argument = \code{compare.f}) defines how the difference between
#' emergent data qualities estimated from different areas should be
#' quantified (e.g. prediction error, change in model structure, absolute
#' difference in statistic, etc).
#'
#' In some sense, data generalization functions similar to the distance
#' calculation method and comparison function - to the linkage function of
#' cluster analysis. The difference is, that the top-down approach used here
#' allows to quantify the distance between groups of data in ways that require
#' certain amount of data points. As bottom-up approaches proceed from
#' comparing and first grouping single data points, they cannot use more
#' emergent data qualities (whose estimation requires certain amount of points
#' distributed in space) to calculate distance between clusters.
#'
#' @param data a data frame with columns containing the variables analyzed and
#' rows - observations, potentially from different locations.
#' @param xy.data a data frame with coordinate information (x,y columns in this
#' order). \code{xy.data} rows correspond to rows (observations) of the
#' \code{data} data frame.
#' @param n.split.pts number of points that are used in creation of split-lines
#' since these points serve as endings / origins of straight, as well as
#' curvi-linear split-lines. Thus, the bigger this number, the more split-lines
#' will be created and tested. Higher values of this parameter greatly
#' increase the computation time dedicated to the search of straight
#' split-lines, but increase the fit to the data.
#' @param generalize.f a function used to estimate some emergent data quality.
#' As an input it  should use a subset of \code{data}, though it can also use
#' the values of other \code{hespidv} arguments. As an output, it should produce
#' an R object, that is recognized and used by \code{compare.f} function to
#' estimate the difference between two data groups separated in space by a
#' split-line.
#' @param compare.f a function used to quantify the difference between two data
#' groups separated in space by a split-line. This function should have
#' arguments \code{plot.1} and \code{plot.2}, since these are the outputs of
#' \code{generalize.f} function applied to the data sets from different areas
#' that were separated by a split-line. Arguments of \code{hespidv} function
#' can also be used as arguments in \code{compare.f} function. The output of
#' \code{compare.f} should be a single numerical value that represents the
#' difference between two data sets (\code{plot.1} and \code{plot.2} objects).
#' The higher this number, the greater the difference will be recognized.
#' @param method A pre-set combination of \code{generalize.f} and
#' \code{compare.f} that serve some distinct purpose. Available methods:
#'  "Pielou_biozonation" (distinguishes paleoprovinces by reductions in Pielou
#'  entropy after data sets of two paleoprovinces are divided).
#' @param N.crit Algorithm stopping criteria - number of observations.
#' Minimum number of observations (rows in data) that should
#' be present in areas separated by a split-line in order to establish the
#' split-line. Default is 1.
#' @param S.crit Algorithm stopping criteria - size of plots.
#' Minimum area in squared coordinate units that plots separated by a
#' split-line should have in order to establish the split-line. Default is 0.
#' @param P.crit Algorithm stopping criteria - magnitude of difference.
#' Minimum difference as estimated by \code{compare.f} function that plots
#' should have in order to establish the split-line. Default is -Inf.
#' @param c.splits Logical (default TRUE).
#' Should curvi-linear split-lines be estimated?
#' @param c.axis.knots Curve parameter. The number of columns in the net of
#' spline knots. These columns are distributed regularly along to the straight
#' split-line. This parameter controls wiggliness (wave length) of the
#' curvi-linear split-lines. Higher values allow  wigglier curves, thus
#' inceasing the fit to the data, but increases the optimization time of
#' curvi-linear split-lines. Default value is 5.
#' @param c.ort.knots Curve parameter. The number of rows in the net of
#' spline knots. These rows are distributed regularly orthogonal to the straight
#' split-line. This parameter controls wiggliness (resolution of tested wave
#' amplitudes) of the curvi-linear split-lines. Higher values allow higher
#' variety of wave amplitudes to be tested, when optimizing the shape of
#' curvi-linear split-lines. Thus higher values increase the fit to the data at
#' the cost of optimization time. Default value is 10.
#' @param c.iter.no number of times the algorithm iterates through the net of
#' spline knots (default 2). Odd number iterations iterate through the net
#' of spline knots along the split-line in eastward direction, while even number
#' iterations - in westward direction. Thus, it is recommended to
#' set an even number for this parameter, in order to keep the balance between
#' eastward and westward iteration biases. Higher values should increase the
#' fit to the data, at the cost of optimization time.
#' @param c.corr.term The term that defines the correction size of problematic
#' curvi-linear split-lines which intersect the boundary of the polygon.
#' Possible values are between 0 and 1, though small values are recommended
#' (default is 0.05). These values define how much the interval of a generated
#' spline, that crosses the boundary of polygon, should be shifted away from the
#' boundary, inside the polygon, in direction orthogonal to the straight
#' split-line, in  terms of proportion of polygon width where spline intersects
#' the polygon boundary.
#' @param n.m.test Logical (default is FALSE). Should the established
#' split-lines be tested with null models? These test are made by counting the
#' proportion of how many times the established boundaries worked better with
#' the same data that were randomly shifted in space. The used spatial
#' randomization method (toroidal shift) quite well preserves the
#' spatial relationship between points (this relationship can be driven by
#' ecological, geological, oceanographic or other natural physical processes
#' and laws that produce predictable spatial changes in environmental and
#' ecological features), but changes their over-all configuration and
#' density distribution (this configuration is idiosyncratic feature since it
#' can depend on a particular instance of environment (e.g. the same
#' environemnt in a different region) or spatial perspective of the same area).
#' Thus, it allows to check whether the same entities gowevern by the same
#' physical laws would be clustered similarly given different instances of
#' "worlds". If considerable proportion of runs (e. g. > %5) produces better
#' boundary performance scores, then it may be either that clustering of data
#' qualities is quite poor or that spatial autocorrelation structure caused by
#' all these physical proccesses is responsible for the observed clusterization.
#' TWO models (toroidal and totally random).
#' @param n.m.N number of spatial simulations in null models. Default is
#' 1000.
#' @param n.m.seed randomization seed that is set before analysis (default 1).
#' The use of the same seed allows to obtain the same stochastic process
#' simulation results (e. g. p values calculated from null model simulations).
#' @param study.pol A polygon of study area (optional). It should be data.frame
#' with two columns containing coordinates (x and y, respectively) of vertixes
#' of a polygon that encompasses the locations of \code{xy.data}. If not
#' provided (default is NULL), convex hull of \code{xy.data} will be used a
#' study area polygon.
#' @param trace.level Integer from 0 to 7, indicates how much information should
#' the algorithm communicate during computations. 0 (default) - no algorithm
#' tracing; 1 - only selected split-splits are reported; 2 - best
#' intermediate straight split-lines are reported; 3 - best intermediate
#' curvi-linear split-lines are reported; 4 - best intermediate split-lines are
#' reported; 5 - all straight split lines are reported; 6 - all curvi-linear
#' split lines are reported; 7 - all split-lines are reported.
#' @param pnts.col Color of data points. Argument is when \code{trace.level} > 0
#' @param inherit.f A function whose output is saved internally
#' as \code{inher_dat} variable. \code{inher_dat} can be used as an argument in
#' \code{generalize.f} or \code{compare.f} functions. Thus, \code{inherit.f} can
#' be used to extract some information from parent plots and data, and pass it
#' to other recursive iterations to be used when dividing offspring plots and
#' data.
#' @param root.heritage The assumed \code{inher_dat} variable, used in the first
#' iteration.
#'
#' @return A list of 2 elements:
#' \describe{
#'   \item{\code{per_pts}}{A data frame of 4 columns, providing the information about the generated points on a perimeter of a polygon. This data frame is used as an input in \code{\link{pair_pts}} function.}
#'   \itemize{
#'   \item \code{x} - X coordinates of generated points.
#'   \item \code{y} - Y coordinates of generated points.
#'   \item \code{ID} - An ID that reflects the relative location of a point along a perimeter of a polygon in relation to other generated points and polygon vertices.
#'   \item \code{segment.no} = A vector indicating the ID of a polygon segment on which a generated point is located. It helps to indentify points located on the same
#' polygon segment.
#'   }
#'   \item{\code{full.poly}}{ A data frame that contains coordinates of the provided polygon vertices and generated points. \code{coords[,"ID"]} can be used to
#' extract rows of generated points.This data frame is used as an input in \code{\link{curvial.split}} function.}
#' }
#' @note If both, n.pts and dst.pts, are specified, then points are generated according to n.pts.
#' @author Liudas Daumantas
#' @importFrom grDevices chull
#' @export

hespdiv<-function(data,polygon=NULL,method=NA,variation=NA,metric=NA,criteria=NA,
                           C.cond=0,E.cond=0,N.cond=0,S.cond=0,divisions=NULL,lim=NULL,
                           knot.density.X=knot.density.X,knot.density.Y=knot.density.Y,curve.iterations,
                           correction.term=0.05,null.models=T,seed.t=round(runif(1,0,9999),0),test.n) {
  ##### pirmas data stulpelis turi buti rusis - faktorius, antras X, trecias Y., dependencies = dplyr
  #S.cond pateikti kaip proporcija ploto
  require(dplyr)
  require(pracma)
  library(gstat)
  library(sp)
  library(spatstat)
  #duomenys pirminiai
  {
    if(names(xy.data)!= c("x","y")){
      stop("xy.data should contain columns named \"x\" and \"y\"")
    }

    ch <- chull(xy.data$x, xy.data$y)
    ids <- c(ch,ch[1])
    x <- xy.data$x[ids]
    y <- xy.data$y[ids]
    origins.chull <- data.frame(X=x,Y=y)

    if (trace.level > 0 ) {
      if(!is.null(pnts.col)){
        plot(xy.data$x, xy.data$y, col=pnts.col )
      } else { plot(xy.data$x,xy.data$y) }

      lines(x,y)
      lines(origins.chull)
    }
    rims <- list(origins.chull)
    blokai <- data.frame(performance=0,iteration=0,saknis=0)
    split.reliability <- numeric()
    split.performance <- numeric()
    split.reliability2 <- numeric()
    n.splits <- numeric()
    ave.split.abE <- numeric()
  }
  S.cond<-abs(polyarea(x,y))*S.cond
  splits<-numeric()
  checks<-list()
  original.quality<-numeric()
  iteration<-1
  #motinine rekursyvine funkcija
  spatial_div<-function(samp.dat, root=2){
    #testuojamas plotas
    testid<-length(rims)
    margins<-rims[[testid]]
    original.qual<- -2*entropija(samp.dat[,1])
    iteration<<- iteration +1
    #grafikas pirminis nupaisomas
    {
      plot(data$X,data$Y,col=data$rusis)
      points(samp.dat$X,samp.dat$Y,pch=19)
      centras<-poly_center(margins[,1],margins[,1])
      points(centras[1],centras[2],col=3,pch=19)
      lines(rims[[1]])
    }
    print(c("testid: ", testid))
    # padalinimai savo ribose nupaisomi
    if (testid>1) {
      for (i in 2:c(testid)){
        print(c("brai?om r?m?", i-1))
        lines(x=rims[[i]][,1],y=rims[[i]][,2],col=1,lwd=2)
      }}
    #testavimui pjuviai paruosiami
    perim_pts<-.perimeter_pts(polygon = margins,n.pts = divisions)
    points(perim_pts[[2]],pch=19,col="purple")
    pairs_pts<-.pair_pts(perim_pts[[1]],polygon = margins)
    maxdif<- original.qual
    print(maxdif)
    any.split<-numeric()
    maxid<-0
    if (nrow(pairs_pts)!=0){
      points(pairs_pts[,c(1:2)],col="yellow",pch=19)
      points(pairs_pts[,c(3:4)],col="red",pch=19)
      #pjaustymo ir testavimo ciklas
      {
        for (i in 1:nrow(pairs_pts)){
          print('testuojamas padalinimas Nr.:')
          print(i)

          virs <- .close_poly(.split_poly(polygon = perim_pts[[2]], min_id = 1,
                              split_ids = pairs_pts[i,c(6:7)],
                              trivial_side = TRUE,poli_side = TRUE))
          po <- .close_poly(.split_poly(polygon = perim_pts[[2]], min_id = 1,
                            split_ids = pairs_pts[i,c(6:7)],
                            trivial_side = TRUE,poli_side = FALSE))

          # padalinami duomenys i dvi dalis pagal pjuvio koordinates
          Puses <- list(.get_data(po,samp.dat),.get_data(virs,samp.dat))
          if (all(c(length(Puses[[1]]$rusis),length(Puses[[2]]$rusis))>N.cond)){
            SpjuvioI<-abs(polyarea(x=virs[,1],y=virs[,2]))
            SpjuvioII<-abs(polyarea(x=po[,1],y=po[,2]))
            if (SpjuvioI>S.cond&SpjuvioII>S.cond){
              #nupiesiam padalinima ir paskaiciuojam kokybe
              Skirtumas<-alfa(Puses[[1]]$rusis,Puses[[2]]$rusis)
              any.split<-c(any.split,Skirtumas)
              #Paskaiciuojam plotus padalintu bloku
              if (Skirtumas > maxdif){
                #Jei padalinimas patenkina minimalias saligas ir yra geresnis nei pries tai - pasizymim ir issisaugom ji
                maxdif<-Skirtumas
                maxid<-i
                print(c('max Skirtumas=',Skirtumas))
              }
            }
          }
        }
      }
      if(length(any.split)>0){
        performance<-mean(any.split)
        #jei performance 0 ir visu padalinimu performance 0, tai padalinimo nera - idealiai atskirtas plotas
        if(all(any.split==0)){
          maxid<-0
        }
      } else {
        performance<-original.qual # negalejom ivertinti ne vieno padalinimo, taigi performance lygu max.
      }

      blokai<<-rbind(blokai,data.frame(performance=1-performance/original.qual,
                                      iteration=iteration,saknis=root))
      print(c("blokai: ", blokai))
      print(c("performance: ", performance))
    }
    # duomenu saugojimas
    #Jei rastas tinkamas padalinimas - ieskom geriausios padalinimo kreives,
    #issaugom duomenis ir ziurim ar galima skaidyti toliau
    if (maxid>0){
      print(perim_pts)
      print(pairs_pts)
      best.curve<-curvial.split(poly.x=perim_pts[[2]]$x.poly,poly.y=perim_pts[[2]]$y.poly,
                                min.x.id = pairs_pts[maxid,6],max.x.id = pairs_pts[maxid,7],b=pairs_pts[maxid,5],samp.dat,knot.density.X=knot.density.X,
                                knot.density.Y=knot.density.Y,N.cond,S.cond,iteracija=curve.iterations,correction.term=correction.term,original.qual)
      lines(best.curve[[1]],col=2,lwd=3)
      if ((1-(max(best.curve[[2]],maxdif)/original.qual))<C.cond ){
        maxid<-0}}
    if (maxid>0){
      if(best.curve[[2]]>maxdif) {
        splits<<-do.call(c,list(splits,list(data.frame(x=best.curve[[1]]$x,y=best.curve[[1]]$y))))
        split.performance<<-do.call(c,list(split.performance,best.curve[[2]]))
        split.reliability<<-do.call(c,list(split.reliability,(best.curve[[2]]-performance)/sd(any.split)))
        #dalinam duomenis padalinimo kreive
        #reikia sukurti poligonus du ir nufiltruoti duomenis  - galima padaryti geriau
        virsus.h<-filter.reorg.poly(polygon = data.frame(xp=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],yp=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),min.x.id =pairs_pts[maxid,6],
                                    max.x.id =pairs_pts[maxid,7],poli.side = T,b= pairs_pts[maxid,5])
        O.poli<-close.poly(split.poly = virsus.h,split.line.x = best.curve[[1]]$x, split.line.y = best.curve[[1]]$y)
        O<-get.data(O.poli,samp.dat)
        lines(O.poli,col=5,lwd=4)

        apacia.h<-filter.reorg.poly(polygon = data.frame(xp=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],yp=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),min.x.id =pairs_pts[maxid,6],max.x.id =pairs_pts[maxid,7],poli.side = F,b=pairs_pts[maxid,5])
        OO.poli<-close.poly(split.poly = apacia.h,split.line.x = best.curve[[1]]$x, split.line.y = best.curve[[1]]$y)
        lines(OO.poli,col=5,lwd=4)

        OO<-get.data(OO.poli,samp.dat)

        #issaugom duomenis padalinimo
        ribs<-list(O.poli,OO.poli)
      } else {
        #issaugom duomenis padalinimo
        splits<<-do.call(c,list(splits,list(data.frame(x=as.numeric(c(pairs_pts[maxid,c(1,3)])),y=as.numeric(c(pairs_pts[maxid,c(2,4)]))))))
        split.performance<<-do.call(c,list(split.performance,maxdif))
        split.reliability<<-do.call(c,list(split.reliability,(maxdif-performance)/sd(any.split)))
        #padalinam duomenis pagal pjuvi
        virs<-close.poly(split.line.x = as.numeric(pairs_pts[maxid,c(1,3)]),split.line.y = as.numeric(pairs_pts[maxid,c(2,4)]),
                         split.poly =filter.reorg.poly(polygon = data.frame(xp=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],yp=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
                                                       min.x.id =pairs_pts[maxid,6],max.x.id =pairs_pts[maxid,7],poli.side = T,b=pairs_pts[maxid,5]))

        po<-close.poly(split.line.x = as.numeric(pairs_pts[maxid,c(1,3)]),split.line.y = as.numeric(pairs_pts[maxid,c(2,4)]),
                       split.poly =filter.reorg.poly(polygon = data.frame(xp=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],yp=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
                                                     min.x.id =pairs_pts[maxid,6],max.x.id =pairs_pts[maxid,7],poli.side = F,b=pairs_pts[maxid,5]))

        O<-get.data(virs,samp.dat)
        OO<-get.data(po,samp.dat)

        #paryskinam atrinkta pjuvi
        print(c('total max Skirtumas =', maxdif))

        lines(pairs_pts[maxid,c(1,3)],pairs_pts[maxid,c(2,4)],col=2)

        # sukuriam koordinates ribines zemu koordinaciu plotui ir
        # bandom skaidyti ji. Jei pavyksta suskaidyti issaugom updatintus rezus, jei ne istrinam gogolis
        # ir ribines koordinates neegzistuojancio pjuvio
        ribs<-list(virs,po)
      }
      #maisom duomenis ir vertinam aptiktu erdviniu strukturu patimuma
      if (null.models){
        set.seed(seed.t)
        test.samp.dat<-sim.testdat(samp.dat,test.n)
        pseudo.kokybe<-numeric(test.n)
        for (a in 1:test.n){
          I.dat<-get.data(ribs[[1]],test.samp.dat[[a]])
          II.dat<-get.data(ribs[[2]],test.samp.dat[[a]])
          pseudo.kokybe[a]<-alfa(I.dat[,1],II.dat[,1])
        }
        checks<<-do.call(c,list(checks,list(pseudo.kokybe)))
        split.reliability2<<-do.call(c,list(split.reliability2,sum(last(split.performance)<pseudo.kokybe)/test.n))
      }
      n.splits<<-do.call(c,list(n.splits,length(any.split)))
      ave.split.abE<<-do.call(c,list(ave.split.abE,performance))
      original.quality<<-do.call(c,list(original.quality,original.qual))
      #
      rims<<-do.call(c,list(rims,ribs[1]))
      lines(ribs[[1]],col="purple")
      spatial_div(O,root = iteration)
      print(paste("griztam i", testid, "padalinima [po mazu koord bloko]", sep=" "))


      # Skaidom antra bloka
      #jei egzistuoja gogolis (updatinti duomenys) tuomet sukuriam ribines koordinates ir lipdom prie
      #gogolis masyvo. Duotu koordinaciu ribose ir bandom ieskoti pjuvio bei toliau updatinti duomenis.
      #Jei pavyksta rasti pjuvi, updatinam duomenis, jei ne trinam null bobolis ir priklijuotas koo-
      #rdinates.

        rims<<-do.call(c,list(rims,ribs[2]))
        lines(ribs[[2]],col=2)
        spatial_div(OO,root=iteration)
        print(paste("griztam i", testid, "padalinima [po aukstu koord bloko (gogolis exists)]", sep=" "))

    } else{
      #Jei tinkamo padalino nerasta, grizta tuscias masyvas
      if (testid>1){
        print("tinkamo padalinimo nerasta, grizta tuscias masyvas NULL")
      } else{
        print("nebuvo skaldomu bloku")
}
    }
  }
  environment(spatial_div) <- environment()

  spatial_div(data,root=2)


  if (null.models==F){
    split.reliability2<-rep(NaN,length(n.splits))
  }

  Signif <- symnum(split.reliability2, corr = FALSE, na = FALSE,
                   cutpoints = c(0 ,0.001,0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))


  rezas <- structure(list(
    split.lines = splits,
    boundaries = rims,
    block.stats = blokai[-1,],
    split.stats = data.frame(
      n.splits = n.splits,
      z.score = round(split.reliability,2),
      ave.abE = round(-ave.split.abE,2),
      split.abE = round(-split.performance,2),
      parent.2E = round(-original.quality,2),
      delta.E = round(split.performance -  original.quality,2),
      p_value = split.reliability2,
      signif. = format(Signif)
    ),
    null.m.st= checks
  ),
  class = "spdiv"
  )
  print.spdiv <- function(x){
    cat("\n","Information about splits:", "\n","\n")
    print(x[[4]])
    cat("\n", "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    invisible(x)
  }
  return(print.spdiv(rezas))
}


