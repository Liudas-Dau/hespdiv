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
#' Two main functions must be provided to estimate the separation of data in
#' space. First one (argument = \code{generalize.f}) is needed to calculate some
#' emergent data quality (eg. some model, summary statistic, etc.).
#' The second one (argument = \code{compare.f}) defines how the difference between
#' emergent data qualities estimated from different areas should be
#' quantified (e.g. prediction error, change in model structure, absolute
#' difference in statistic, etc).
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
#' rows - observations, potentially from different locations. \code{data} must
#' contain columns named "x" and "y" that contain coordinate information.
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
#' Minimum area expressed as a proportion of original study area
#' (provided polygon or estimated as convex hull of \code{xy_dat}) that plots
#' separated by a split-line should have so that the split-line could be
#' established. Default is 0.
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
#' with two columns 'x' and 'y' containing coordinates of vertexes
#' of a polygon that encompasses the locations of \code{data}. If not
#' provided (default is NULL), convex hull of \code{data} will be used a
#' study area polygon.
#' @param trace.level Integer from 0 to 7, indicates how much information should
#' the algorithm communicate during computations. 0 (default) - no algorithm
#' tracing; 1 - only selected split-splits are reported; 2 - best
#' intermediate straight split-lines are reported; 3 - best intermediate
#' curvi-linear split-lines are reported; 4 - best intermediate split-lines are
#' reported; 5 - all straight split lines are reported; 6 - all curvi-linear
#' split lines are reported; 7 - all split-lines are reported.
#' @param pnts.col Color of data points, default is 1. Argument is used when
#' \code{trace.level} > 0. If set to NULL, data points will not be displayed.
#' @param inherit.f A function whose output is saved internally
#' as \code{inher.dat} variable. \code{inher.dat} can be used as an argument in
#' \code{generalize.f} or \code{compare.f} functions. Thus, \code{inherit.f} can
#' be used to extract some information from parent plots and data, and pass it
#' to other recursive iterations to be used when dividing offspring plots and
#' data.
#' @param root.heritage The assumed \code{inher.dat} variable, used in the first
#' iteration, where there are no inherited information from parent plots.
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
#' @importFrom pracma polyarea
#' @export

hespdiv<-function(data,polygon=NULL,method=NA,variation=NA,metric=NA,criteria=NA,
                           C.cond=0,E.cond=0,N.cond=0,S.cond=0,divisions=NULL,lim=NULL,
                           knot.density.X=knot.density.X,knot.density.Y=knot.density.Y,curve.iterations,
                           correction.term=0.05,null.models=T,seed.t=round(runif(1,0,9999),0),test.n) {

  if (method == "Pielou_biozonation") {
    if (ncol(data) != 3){
      stop("There should be one column in data besides 'x' and 'y' columns that
           records which taxa is present at a given location. Taxa should be
           coded numerically. If multiple taxa is present at the same location,
           multiple rows should be dedicated for the same location in data.")
    }
    generalize.f <- function(plot.dat){
      plot.dat <- plot.dat[,-c("x","y")]
      p <- table(as.numeric(paste(x)))/sum(table(as.numeric(paste(x))))
      H <- -sum(log(p)*p)
      ifelse(H==0,0,H/log(length(p)))
    }
    compare.f <- function(plot1,plot2) {
      samp.dat <- samp.dat[,-c("x","y")]
      base.eveness <- -2 * generalize.f(samp.dat[,1])
      base.eveness + plot1 + plot2
    }
  }


  ##### pirmas data stulpelis turi buti rusis - faktorius, antras X, trecias Y., dependencies = dplyr
  #S.cond pateikti kaip proporcija ploto
  require(dplyr)
  library(gstat)
  library(sp)
  library(spatstat)
  #duomenys pirminiai
  {
    if(any(names(data)== "x") & any(names(data)== "y")){
      stop("data should contain columns named \"x\" and \"y\" that contain
           coordinate information")
    }

    ch <- chull(data$x, data$y)
    ids <- c(ch,ch[1])
    x <- data$x[ids]
    y <- data$y[ids]
    origins.chull <- data.frame(X=x,Y=y)

    if (trace.level > 0 ) {
      if(!is.null(pnts.col)){
        plot(data$x, data$y, col=pnts.col )
      } else {
        plot(0,0,xlim = range(data$x),ylim = range(data$y),col=0)
        }

      lines(x,y)
      lines(origins.chull)
    }
    rims <- list(origins.chull)
    blocks <- data.frame(performance=numeric(0),iteration=numeric(0),
                         root=numeric(0))
    split.reliability <- numeric()
    split.performance <- numeric()
    split.reliability2 <- numeric()
    n.splits <- numeric()
    ave.split.abE <- numeric()
  }
  S.cond <- abs(polyarea(x,y)) * S.crit
  splits <- numeric()
  checks <- list()
  original.quality <- numeric()
  iteration <- 1

  environment(.spatial_div) <- environment()

  .spatial_div(data,root=2)


  if (null.models==F){
    split.reliability2 <- rep(NaN,length(n.splits))
  }

  Signif <- symnum(split.reliability2, corr = FALSE, na = FALSE,
                   cutpoints = c(0 ,0.001,0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))


  rezas <- structure(list(
    split.lines = splits,
    boundaries = rims,
    block.stats = blocks[-1,],
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
  class = "hespdiv"
  )

  return(print.hespdiv(rezas))
}
#' Main recursive hespdiv inner helper function
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
#' @importFrom pracma poly_center
#' @noRd

.spatial_div <- function(samp.dat,samp_xy, root=2){
  #testuojamas plotas
  testid <- length(rims)
  margins <- rims[[testid]]
  original.qual <- -2*entropija(samp.dat[,1]) #### THIS IS P(Q).crit, use inherit.f or even better: extract information from intermediate results
  iteration <<- iteration +1
  perim_pts <- .perimeter_pts(polygon = margins,n.pts = divisions)
  #grafikas pirminis nupaisomas
  {
    if (trace.level > 0 ) {
      if(!is.null(pnts.col)){
        plot(data$x, data$y, col=pnts.col )
      } else {
        plot(0,0,xlim = range(data$x),ylim = range(data$y),col=0)
        }

      lines(x,y)
      lines(origins.chull)
      points(samp_xy$x,samp_xy$y,pch=19)
      centras <- poly_center(margins[,1],margins[,1])
      points(centras[1],centras[2],col=3,pch=19)
      lines(rims[[1]])
      print(c("testid: ", testid))
      # padalinimai savo ribose nupaisomi
      if (testid>1) {
        for (i in 2:c(testid)){
          print(c("polygons drawed:", i-1))
          lines(x=rims[[i]][,1],y=rims[[i]][,2],col=1,lwd=2)
        }}
      points(perim_pts[[2]],pch=19,col="purple")
    }
  }
  #testavimui pjuviai paruosiami



  pairs_pts <- .pair_pts(perim_pts[[1]],polygon = margins)
  maxdif <- original.qual # first split minimum quality. P.crit
  print(maxdif)
  any.split <- numeric()
  maxid <- 0

  if (nrow(pairs_pts)!=0){
    if (trace.level > 0 ) {
      points(pairs_pts[,c(1:2)],col="yellow",pch=19)
      points(pairs_pts[,c(3:4)],col="red",pch=19)
    }
    #pjaustymo ir testavimo ciklas
    {
      for (i in 1:nrow(pairs_pts)){
        print('testuojamas padalinimas Nr.:')
        print(i)

        virs <- .close_poly(
          open.poly =
            .split_poly(
              polygon = perim_pts[[2]],
              min_id = 1,
              split_ids = pairs_pts[i,c(6:7)],
              trivial_side = TRUE,
              poli_side = TRUE
            ))
        po <- .close_poly(
          open.poly =
            .split_poly(
              polygon = perim_pts[[2]],
              min_id = 1,
              split_ids = pairs_pts[i,6:7],
              trivial_side = TRUE,
              poli_side = FALSE
            ))

        # padalinami duomenys i dvi dalis pagal pjuvio koordinates
        Puses <- list(.get_data(po,samp.dat),.get_data(virs,samp.dat))

        if (all(c(nrow(Puses[[1]]),
                  nrow(Puses[[2]]))>N.cond)){
          if (S.crit > 0){
          SpjuvioI <- abs(polyarea(x=virs[,1],y=virs[,2]))
          SpjuvioII <- abs(polyarea(x=po[,1],y=po[,2]))
          cond <- SpjuvioI > S.cond & SpjuvioII > S.cond
          } else {
            cond <- TRUE
          }
          if (cond) {
            #nupiesiam padalinima ir paskaiciuojam kokybe
            environment(generalize.f) <- environment()
            environment(compare.f) <- environment()
            Skirtumas <- .dif_fun()
            any.split <- c(any.split,Skirtumas)
            #Paskaiciuojam plotus padalintu bloku
            if (Skirtumas > maxdif){
              #Jei padalinimas patenkina minimalias saligas ir yra geresnis nei pries tai - pasizymim ir issisaugom ji
              maxdif <- Skirtumas
              maxid <- i
              print(c('max Skirtumas=',Skirtumas))
            }
          }
        }
      }
    }
    if(length(any.split) > 0){

      performance <- mean(any.split) # mean spatial heterogeneity
      # sd(any.split) -> anisotropy of heterogeneity
      # if all are equal, then complete randomness: no anisotropy - no splits present
      if( all(any.split==any.split[1]) & length(any.split) > 1  ){
        warning(paste(c(
          "All tested splits are of equal quality. No anisotropy - no spilts in"
          , iteration, "iteration.")))
        maxid<-0
      }
    } else {
      performance <- original.qual # negalejom ivertinti ne vieno padalinimo, taigi performance lygu max.
      # P.crit
    }

    blocks <<- rbind(blocks,data.frame(
      performance = 1 - performance/original.qual, # perdatyr i tiesiog santyki. Ty. how many times better
      # than null? DEL UNIVERSALUMO, REIKTU P.CRIT NAUDOTI. BET KA DARYT TUO ATVEJU, JEI KAIP CIA, NAUDOJAM
      # LOKALU P.CRIT IS PAVELDETU DUOMENU ARBA IS PASKAICIUOTU SU INHERIT FUNKCIJA?
      # GAL TADA GERIAU TIESIOG VIDUTINI PALIKTI? TUOMET PERSKAICIUOTI GALIMA.
      # MINUSAS TAS, JOG SU ATASKAITA DAR REIKIA ZAISTI.
      # TAIP PAT MINUSAS TAME, JOG ATASKAITOS DUOMENU GALI NEPAKAKTI, JEI KITAIP
      # LOKALUS P.CRIT ISSISKAICIUOJA.
      # PIRMA REIK ATSAKYT I KLAUSIMA KAM BUTENT TOKS RODIKLIS REIKALINGAS
      # KAIP IR LOGISKA BUTU SANTYKI NAUDOT SU LOKALIU NULL MODELIU, JEI
      # LOKALUS NULL MODELIS VERTINAMAS, NES TAS ATSPINDES OBJEKTYVIAU PERFORMANCE
      # NEGU ABSOLIUTUS, KURIO KOKYBE IS TIESU PRIKLAUSO NUO NULL MODELIO
      # BET TUOMET GAL REIKETU PATEIKTI: VIDURKI; SANTYKI SU ABSOLIUCIU NULL;
      # SANTYKI SU LOKALIU NULL?
      # KAIP PAIMTI TA LOKALU NULL? NEV NE VISI GI NAUDOS... IR TAIP SUDETINGAS NAUDOJIMAS
      iteration = iteration,
      root = root
      ))
    print(c("blocks: ", blocks))
    print(c("performance: ", performance))
  # duomenu saugojimas
  #Jei rastas tinkamas padalinimas - ieskom geriausios padalinimo kreives,
  #issaugom duomenis ir ziurim ar galima skaidyti toliau
  if (maxid>0) {
    print(perim_pts)
    print(pairs_pts)

    best.curve <- .curvial_split(

      poly.x = perim_pts[[2]]$x.poly,
      poly.y = perim_pts[[2]]$y.poly,
      min.x.id = pairs_pts[maxid,6],
      max.x.id = pairs_pts[maxid,7],
      b = pairs_pts[maxid,5],
      data = samp.dat,
      knot.density.X = knot.density.X,
      knot.density.Y = knot.density.Y,
      N.cond =N.cond,
      S.cond = S.cond,
      iteracija = curve.iterations,
      n.curve.iter = correction.term

      )
    if (trace.level > 2) {
      lines(best.curve[[1]],col=2,lwd=3)
    }
    if ((1-(max(best.curve[[2]],maxdif)/original.qual)) < P.crit ){ # vel santykinis base line. Be to,
      # galima gi reikalaut, kad atotrukis nuo base line butu tam tikro dydzio. Dar vienas P.crit
      # argumentas reikalingas?
      maxid<-0}
  }

  if (maxid>0){
    if(best.curve[[2]] > maxdif) {
      best.splitl <- data.frame(x = best.curve[[1]]$x, y = best.curve[[1]]$y)
      maxdif <- best.curve[[2]]
    } else {
      best.splitl <- data.frame(x=as.numeric(c(pairs_pts[maxid,c(1,3)])),
                                y=as.numeric(c(pairs_pts[maxid,c(2,4)])))
    }
    splits <<- do.call(c,list(splits,list(best.splitl)))
    split.performance <<- do.call(c,list(split.performance,maxdif ))
    split.reliability <<- do.call(c,list(split.reliability,
                                         (maxdif - performance)/sd(any.split))) # z-score of performance
    #dalinam duomenis padalinimo kreive
    #reikia sukurti poligonus du ir nufiltruoti duomenis  - galima padaryti geriau
    up.pol <- .close_poly(
      open.poly =
        .split_poly(
          polygon = data.frame(x=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],
                               y=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
          split_ids = pairs_pts[maxid,6:7],
          min_id = 1,
          trivial_side = TRUE,
          poli_side = TRUE
        ),
      close.line = best.splitl
    )
    up.dat <- .get_data(up.pol,samp.dat)

    do.pol <- .close_poly(
      open.poly =
        .split_poly(
          polygon = data.frame(x=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],
                               y=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
          split_ids = pairs_pts[maxid,6:7],
          min_id = 1,
          trivial_side = TRUE,
          poli_side = FALSE
        ),
      close.line = best.splitl
    )
    do.dat <- .get_data(do.pol,samp.dat)

    if (trace.level > 0) {
      lines(up.pol,col=5,lwd=4)
      lines(do.pol,col=5,lwd=4)
    }

    #issaugom duomenis padalinimo
    ribs <- list(up.pol,do.pol)
    print(c('total max Skirtumas =', maxdif))

    #maisom duomenis ir vertinam aptiktu erdviniu strukturu patimuma
    if (n.m.test){

      set.seed(n.m.seed)
      test.samp.dat <- sim.testdat(samp.dat,n.m.N) # reikia pakurti dar sita funkcija
      pseudo.kokybe <- numeric(n.m.N)

      for (a in 1:n.m.N){
        I.dat <- .get_data(ribs[[1]], test.samp.dat[[a]])
        II.dat <- .get_data(ribs[[2]], test.samp.dat[[a]])
        pseudo.kokybe[a] <- alfa(I.dat[,1], II.dat[,1])
      }
      checks <<- do.call(c,list(checks,list(pseudo.kokybe)))
      split.reliability2 <<- do.call(c,list(
        split.reliability2,
        sum(split.performance[length(split.performance)]<pseudo.kokybe)/n.m.N
        )) # kvantilis
      # kvantilio reiksme - empirine p verte
      # jei gausinis skirstinys,tai:
      # 1 - qnorm(sum(last(split.performance),mean(pseudo.kokybe),sd(pseudo.kokybe)) duotu teorine p-verte
    }
  # updatinam ka reikia in hespdiv env.
    n.splits <<- do.call(c,list(n.splits,length(any.split)))
    ave.split.abE <<- do.call(c,list(ave.split.abE,performance))
    original.quality <<- do.call(c,list(original.quality,original.qual))
    #
    rims <<- do.call(c,list(rims,ribs[1]))
    # where to next?
    if (trace.level > 1) {
      lines(ribs[[1]],col="purple")
    }
    .spatial_div(up.dat, root = iteration)
    print(paste("griztam i", testid, "padalinima [po mazu koord bloko]", sep=" "))


    # Skaidom antra bloka
    #jei egzistuoja gogolis (updatinti duomenys) tuomet sukuriam ribines koordinates ir lipdom prie
    #gogolis masyvo. Duotu koordinaciu ribose ir bandom ieskoti pjuvio bei toliau updatinti duomenis.
    #Jei pavyksta rasti pjuvi, updatinam duomenis, jei ne trinam null bobolis ir priklijuotas koo-
    #rdinates.

    rims <<- do.call(c,list(rims,ribs[2]))

    if (trace.level > 1) {
      lines(ribs[[2]],col="purple")
    }
    .spatial_div(do.dat, root = iteration)
    print(paste("griztam i", testid, "padalinima [po aukstu koord bloko (gogolis exists)]", sep=" "))

  }} else{
    #Jei tinkamo padalino nerasta, grizta tuscias masyvas
    if (testid>1){
      print("tinkamo padalinimo nerasta, grizta tuscias masyvas NULL")
    } else{
      print("nebuvo skaldomu bloku")
    }
  }
}

.dif_fun <- function() {
  compare.f( generalize.f(Puses[[1]]), generalize.f(Puses[[2]]) )
}

print.hespdiv <- function(x){
  cat("\n","Information about splits:", "\n","\n")
  print(x[[4]])
  cat("\n", "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  invisible(x)
}
