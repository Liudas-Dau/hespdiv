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
#' space and data. \cr Two main functions must be provided to estimate the
#' separation of data in
#' space. First one (argument = \code{generalize.f}) is needed to calculate some
#' emergent data quality (eg. some model, summary statistic, etc.).
#' The second one (argument = \code{compare.f}) defines how the difference between
#' emergent data qualities estimated from different areas should be
#' quantified (e.g. prediction error, change in model structure, absolute
#' difference in statistic, etc). \cr In some sense, data generalization
#' functions similar to the distance
#' calculation method and comparison function - to the linkage function of
#' cluster analysis. The difference is, that the top-down approach used here
#' allows to quantify the distance between groups of data in ways that require
#' certain amount of data points. As bottom-up approaches proceed from
#' comparing and first grouping single data points, they cannot use more
#' emergent data qualities (whose estimation requires certain amount of points
#' distributed in space) to calculate distance between clusters.
#' @param data a data frame with columns containing the variables analyzed and
#' rows - observations, potentially from different locations. \code{data} must
#' have columns named "x" and "y" that contain coordinate information.
#' @param n.split.pts number of points that are used in creation of split-lines
#' since these points serve as endings / origins of straight, as well as
#' curvi-linear split-lines. Thus, the bigger this number, the more split-lines
#' will be created and tested. Higher values of this parameter greatly
#' increase the computation time dedicated to the search of straight
#' split-lines, but increase the fit to the data.
#' @param generalize.f a function used to estimate some emergent data quality.
#' It must have one argument to which a spatially filtered subset of \code{data}
#' could be assigned. However, \code{generalize.f} can access the
#' \code{hespdiv} environment, thus inside this function other variables from
#' \code{hespdiv} environment could be used, without requiring them as
#' arguments (see the list of free variables available in \code{generalize.f}
#' - list). As an output, \code{generalize.f} should produce
#' an R object, that is recognized and used by \code{compare.f} function to
#' estimate the quality of the split-lines (difference between two data groups
#' separated in space by a split-line)
#' @param compare.f a function used to quantify the difference between two data
#' groups separated in space by a split-line. The estimated difference
#' represents the quality of a split-line. This function must have two
#' arguments to which the outputs of \code{generalize.f} could be assigned.
#' These outputs are the emergent data qualities estimated at opposite sides of
#' a split-line. Thus, \code{compare.f} must define how these data qualities
#' should be compared to produce a difference measure that represents the
#' quality of a split line. Similar to \code{generalize.f}, \code{compare.f}
#' can access the \code{hespdiv} environment. Thus, inside \code{compare.f}
#' function other variables from \code{hespdiv} environment could be used,
#' without requiring them as arguments (see the list of free variables
#' available in \code{generalize.f} - list).
#' @param method A pre-set combination of \code{generalize.f} and
#' \code{compare.f} functions that serve some distinct purpose.
#' Available methods:
#'  "Pielou_biozonation" - distinguishes paleoprovinces by maximum reductions in
#'  Pielou entropy that are observed in the occurrence data of fossil taxa, when
#'  the split-line manages to correctly separate different paleoprovinces.
#' @param N.crit Subdivision stopping criteria - number of observations.
#' Minimum number of observations (rows in data) that should
#' be present in areas separated by a split-line in order to establish the
#' split-line. Default is 1.
#' @param S.crit Subdivision stopping criteria - size of plots.
#' Minimum area expressed as a proportion of original study area
#' (provided polygon or estimated as convex hull of \code{xy_dat}) that plots
#' separated by a split-line should have so that the split-line could be
#' established. Default is 0.
#' @param lower.Q.crit Subdivision stopping criteria - lower limit of split-line
#' quality applied to the straight split-lines. This is a minimum difference as
#' estimated by \code{compare.f} function that separated
#' plots should exhibit, so that a straight split-line would be accepted. If
#' the best straight split line has lower quality than \code{upper.Q.crit}, but
#' passes this limit, then curvi-linear split-lines will be generated with
#' expectation that they will improve the quality of a split above
#'  \code{upper.Q.crit}.
#' If \code{c.splits} is FALSE, then \code{lower.Q.crit} is set equal to
#' \code{upper.Q.crit}
#' Default is -Inf.
#' @param upper.Q.crit Subdivision stopping criteria - upper limit of split-line
#' quality applied to the final split-line. This is a minimum difference as
#' estimated by \code{compare.f} function that separated
#' plots should exhibit, so that a subdivision of a plot using the
#' best split-line would be established. Default is -Inf.
#' @param c.splits Logical (default TRUE).
#' Should curvi-linear split-lines be estimated?
#' @param c.axis.knots Curve parameter. The number of columns in the net of
#' spline knots. These columns are distributed regularly along to the straight
#' split-line. This parameter controls wiggliness (wave length) of the
#' curvi-linear split-lines. Higher values allow  wigglier curves, thus
#' increasing the fit to the data, but increases the optimization time of
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
#' @param n.m.keep logical (default FALSE). Do you wish to keep null model
#' simulations?
#' @param study.pol A polygon of study area (optional). It should be data.frame
#' with two columns 'x' and 'y' containing coordinates of vertexes
#' of a polygon that encompasses the locations of \code{data}. If not
#' provided (default is NULL), convex hull of \code{data} will be used a
#' study area polygon.
#' @param trace.level Integer from 0 to 7, indicates how much information should
#' the algorithm communicate during computations. 0 (default) - no algorithm
#' tracing; 1 - only the best split-lines are reported; 2 - best
#' intermediate straight split-lines are reported; 3 - best intermediate
#' curvi-linear split-lines are reported; 4 - best intermediate split-lines are
#' reported; 5 - all straight split lines are reported; 6 - all curvi-linear
#' split lines are reported; 7 - all split-lines are reported.
#' @param pnts.col Color of data points, default is 1. Argument is used when
#' \code{trace.level} > 0. If set to NULL, data points will not be displayed.
#' @return hespdiv class object, a list of at least 5 elements (see details):
#' \describe{
#'   \item{\code{split.lines}}{ a list containing data frames of
#'   split-lines coordinates}
#'   \item{\code{poly.stats}}{ a data frame containing information about
#'   polygons established by the split-lines. Columns:}
#'   \itemize{
#'   \item \code{mean.dif} - mean quality of straight split-lines that were
#'   produced and tested inside a polygon. Can be interpreted as a spatial
#'   heterogeneity of emergent data quality.
#'   \item \code{sd.dif} - standard deviation of quality of straight split-lines
#'   that were produced and tested inside a polygon. Can be interpreted as the
#'   extent of anisotropy in spatial heterogeneity of emergent data quality.
#'   \item \code{str.z.score} - z-score of the best straight split-line quality
#'   produced in a polygon. Indicates, how outstanding the produced straight
#'   split-line is, when compared to other tested straight split-lines.
#'   \item \code{iteration} - ID of \code{hespdiv} iteration, in which a
#'   polygon was analyzed. Can be considered as the ID of a polygon.
#'   \item \code{root} - the ID of \code{hespdiv} iteration, which produced the
#'   split-line, isolating a polygon. Can be considered as the ID of
#'   a polygon's parent-polygon.
#'   }
#'   \item{\code{polygons.xy}}{a list, containing data frames of polygons,
#'   produced by split-lines. Who corresponds to who??????????????}
#'   \item{\code{poly.obj}}{a list of \code{generalize.f} outputs of each
#'   polygon. Names of \code{poly.obj} elements correspond
#'   to the \code{iteration} column of \code{poly.stats}}
#'   \item{\code{split.stats}}{a data frame containing information about the
#'   produced split-lines. Its contents depend on the choice of method and
#'   the use of null models (see details)}
#'   }
#' @note Nothing to note yet.
#' @details Nothing to detail yet (split.stats, null model results).
#' @author Liudas Daumantas
#' @importFrom grDevices chull
#' @importFrom pracma polyarea
#' @export

hespdiv<-function(data, n.split.pts = 15 ,generalize.f = NULL,
                  compare.f = NULL, method = "Pielou_biozonation", N.crit = 0,
                  S.crit = 0, lower.Q.crit = -Inf, upper.Q.crit = -Inf,
                  c.splits = TRUE, c.axis.knots = 5, c.ort.knots = 10,
                  c.iter.no = 2, c.corr.term = 0.05, n.m.test = FALSE,
                  n.m.N = 1000, n.m.seed = 1,  n.m.keep = FALSE,
                  study.pol = NULL, trace.level = 0, pnts.col = 1){

    if (c.splits == FALSE & upper.Q.crit != lower.Q.crit) {
    print(paste("Since 'c.splits' is FALSE, 'lower.Q.crit' is set equal to
          'upper.Q.crit'"))
      lower.Q.crit <- upper.Q.crit
  }
  if (method == "Pielou_biozonation") {
    if (ncol(data) != 3){
      stop("There should be one column in data besides 'x' and 'y' columns that
           records which taxa is present at a given location. Taxa should be
           coded numerically. If multiple taxa is present at the same location,
           multiple rows should be dedicated for the same location in data.")
    }

generalize.f <- function(plot.dat){
      x <- table(as.numeric(paste(subset(plot.dat, select = -c(x, y))[,1])))
      p <- x/sum(x)
      H <- -sum(log(p)*p)
    }
    compare.f <- function(eveness1,eveness2) {
      base.eveness <- poly.obj[[testid]]
      (1 - mean(eveness1, eveness2) / base.eveness) * 100 # percent change in eveness
    }
  }


  ##### pirmas data stulpelis turi buti rusis - faktorius, antras X, trecias Y., dependencies = dplyr
  #S.cond pateikti kaip proporcija ploto
  #require(dplyr)
  #library(gstat)
  #library(sp)
  #library(spatstat)
  #duomenys pirminiai
  {
    if( all(names(data) != "x") | all(names(data) != "y") ){
      stop("data should contain columns named \"x\" and \"y\" that contain
           coordinate information")
    }
    if (is.null(study.pol)){

      ch <- chull(data$x, data$y)
      ids <- c(ch,ch[1])
      x <- data$x[ids]
      y <- data$y[ids]
      study.pol <- data.frame(x=x,y=y)
    }
    if (trace.level > 0 ) {
      if(!is.null(pnts.col)){
        plot(data$x, data$y, col=pnts.col )
      } else {
        plot(0,0,xlim = range(data$x),ylim = range(data$y),col=0)
      }

      lines(study.pol)
    }
    rims <- list(study.pol)
    poly.info <- data.frame(mean.dif = numeric(), # mean spatial heterogeneity irrespective of split-line position
                         sd.dif = numeric(), # anizotropy of heterogeneity based on straight split-lines
                         str.z.score = numeric(), # level of outstandingness. Are there other competetive candidate splits?
                         iteration=numeric(),
                         root=numeric() )
    poly.obj <- list()
    plot.id <- numeric()
    split.z.score <- numeric()
    split.quality <- numeric()

    if (n.m.test){
      p.val1 <- numeric()
      p.val2 <- numeric()
      sim1.difs <- numeric()
      sim2.difs <- numeric()
      if (n.m.keep){
        n.m.sims1 <- list()
        n.m.sims2 <- list()
      }
    }

    n.splits <- numeric()
    mean.difs <- numeric()
  }
  S.cond <- abs(pracma::polyarea(x,y)) * S.crit
  splits <- numeric()

  iteration <- 1

  e <- environment()
  environment(.spatial_div) <- e


  .spatial_div(data,root=2)


  names(poly.obj) <- poly.info$iteration
  names(rims) <- poly.info$iteration # bad line?



  if(method == "Pielou_biozonation"){
    parent.E <- unlist(poly.obj[paste(plot.id)])
    result <- structure(list(
      split.lines = splits,
      polygons.xy = rims,
      poly.stats = poly.info,
      poly.obj = poly.obj,
      split.stats = data.frame(
        plot.id = plot.id,
        n.splits = n.splits,
        z.score = split.z.score,
        mean.p.red = mean.difs,
        split.p.red = split.quality,
        parent.E = parent.E,
        delta.E = -parent.E * split.quality/100
      )
    ),
    class = "hespdiv"
    )
  } else {
    result <- structure(list(
      split.lines = splits,
      polygons.xy = rims,
      poly.stats = poly.info,
      poly.obj = poly.obj,
      split.stats = data.frame(
        plot.id = plot.id,
        n.splits = n.splits,
        z.score = split.z.score,
        mean.dif = mean.difs,
        split.quality = split.quality
      )
    ),
    class = "hespdiv"
    )
  }

  if (n.m.test){
    Signif1 <- symnum(p.val1, corr = FALSE, na = FALSE,
                      cutpoints = c(0 ,0.001,0.01, 0.05, 0.1, 1),
                      symbols = c("***", "**", "*", ".", " "))
    Signif2 <- symnum(p.val2, corr = FALSE, na = FALSE,
                      cutpoints = c(0 ,0.001,0.01, 0.05, 0.1, 1),
                      symbols = c("***", "**", "*", ".", " "))

    result$split.stats <- cbind(result$split.stats,
                               data.frame(p.val1 = p.val1,
                                          signif.1 = format(Signif1),
                                          p.val2 = p.val2,
                                          signif.2 = format(Signif2)))
    result <- do.call(c,list(result,
                            list(n.m.rez =
                                   list(sim1.difs = sim1.difs, sim2.difs))))

    if (n.m.keep){
      result <- do.call(c,list(result,
                               list(n.m.sim = list(n.m.sims1,n.m.sims2))))
    }
  }


  return(print.hespdiv(result,n.m.test))
}
#' Main recursive hespdiv inner helper function
#'
#' @description  During each recursive iteration this function fits straight and
#' curvi-linear split-lines, then it uses the better of the two to separate
#' data in space, saves information about the polygon and split-lines,
#' and calls itself again using the extracted data samples.
#' @param samp.dat a spatial subset of \code{data} that lies entirely within
#' some polygon, produced using split-line subdivisions of original study area.
#' At first iteration, while no subdivisions are established,
#' \code{samp.dat = data}.
#' @param root An id of recursive iteration that produced the \code{samp.dat}.
#' An id of parent-polygon.
#' @return No return. Function updates variables in \code{hespdiv} environment.
#' @author Liudas Daumantas
#' @note Function inherits the environment of \code{hespdiv} function.
#' @importFrom pracma poly_center
#' @noRd

.spatial_div <- function(samp.dat, root=2){
  #testuojamas plotas
  testid <- length(rims)
  margins <- rims[[testid]]

  assign(x = "iteration",value = iteration +1, envir = e)
  assign(x = "poly.obj" ,
         value = do.call(c,list(poly.obj,list(generalize.f(samp.dat)))),
         envir = e)

  perim_pts <- .perimeter_pts(polygon = margins,n.pts = n.split.pts)
  #grafikas pirminis nupaisomas
  {
    if (trace.level > 0 ) {
      if(!is.null(pnts.col)){
        plot(data$x, data$y, col=pnts.col )
      } else {
        plot(0,0,xlim = range(data$x),ylim = range(data$y),col=0)
        }

      centras <- pracma::poly_center(margins[,1],margins[,2])
      points(centras[1],centras[2],col=3,pch=19)
      lines(rims[[1]])
      print(paste0("Polygon tested No. : ", testid))
      # padalinimai savo ribose nupaisomi
      if (testid>1) {
        for (i in 2:c(testid)){
          print(paste0("polygons drawed: ", i-1))
          lines(x=rims[[i]][,1],y=rims[[i]][,2],col=1,lwd=2)
        }}
      points(perim_pts[[1]],pch=19,col="purple")
    }
  }
  #testavimui pjuviai paruosiami

  environment(.dif_fun) <- e

  pairs_pts <- .pair_pts(perim_pts[[1]],polygon = margins)
  maxdif <- lower.Q.crit # first split minimum quality. P.crit
  any.split <- numeric()
  maxid <- 0

  if (nrow(pairs_pts)!=0){
    if (trace.level > 0 ) {
      points(pairs_pts[,c(1:2)], col = "green", pch=19)
      points(pairs_pts[,c(3:4)], col = "green", pch=19)
    }
    #pjaustymo ir testavimo ciklas
    {
      for (i in 1:nrow(pairs_pts)){
        virs <- .close_poly(
          open.poly =
            .split_poly(
              polygon = perim_pts[[2]],
              min_id = 1,
              split_ids = as.numeric(pairs_pts[i,c(6:7)]),
              trivial_side = TRUE,
              poli_side = TRUE
            ))
        po <- .close_poly(
          open.poly =
            .split_poly(
              polygon = perim_pts[[2]],
              min_id = 1,
              split_ids = as.numeric(pairs_pts[i,6:7]),
              trivial_side = TRUE,
              poli_side = FALSE
            ))

        # padalinami duomenys i dvi dalis pagal pjuvio koordinates
        Puses <- list(.get_data(po,samp.dat),.get_data(virs,samp.dat))
        if (any(trace.level == c(5,7)) ) {
          readline(prompt = paste0('Going to test straight split-line No. : ',i))
          points(pairs_pts[i,c(3,4)],col = 4, pch = 19, cex = 1.5)
          lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
                col = "yellow", pch = 19)
          # filtravimas spalvu su get.dat - imanoma, bet gal nereikia?
          points(Puses[[1]]$x,Puses[[1]]$y,col="darkgreen",pch=19)
          readline(prompt = "Press enter, to see points from other polygon")
          points(Puses[[2]]$x,Puses[[2]]$y,col=6,pch=19)
        }
        if (all(c(nrow(Puses[[1]]),
                  nrow(Puses[[2]]))>N.crit)){
          if (S.crit > 0){
            SpjuvioI <- abs(pracma::polyarea(x=virs[,1],y=virs[,2]))
            SpjuvioII <- abs(pracma::polyarea(x=po[,1],y=po[,2]))
            cond <- SpjuvioI > S.cond & SpjuvioII > S.cond
          } else {
            cond <- TRUE
          }
          if (cond) {
            #nupiesiam padalinima ir paskaiciuojam kokybe
            environment(.dif_fun) <- environment()
            Skirtumas <- .dif_fun(Puses[[1]],Puses[[2]])
            any.split <- c(any.split,Skirtumas)
            #Paskaiciuojam plotus padalintu bloku
            if (Skirtumas > maxdif){
              if (any(trace.level == c(2,4,5,7)) ) {
                lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
                      col = "blue", pch = 19)
                print(paste0('Difference obtained between polygons: ',round(Skirtumas,2)))
                print(paste0('Old best difference: ',round(maxdif,2)))
                readline(prompt = "The blue line is currently the
                         best straight split-line. Press enter to continue...")
              }

              #Jei padalinimas patenkina minimalias saligas ir yra geresnis nei pries tai - pasizymim ir issisaugom ji
              maxdif <- Skirtumas
              maxid <- i
            } else {
              if (any(trace.level == c(5,7)) ) {
                lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
                      col = "gray60", pch = 19)
                points(pairs_pts[i,c(3,4)],col = "gray60", pch = 19, cex = 1.5)
                print("Estimated difference was too small.")
                print(paste0("Required: ",round(maxdif,2)))
                print(paste0("Was obatained: ",round(Skirtumas,2)))
              }
            }
          } else {
            if (any(trace.level == c(5,7)) ){
              print("The obtained polygons were too small.")
              lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
                    col = "gray60", pch = 19)
              points(pairs_pts[i,c(3,4)],col = "gray60", pch = 19, cex = 1.5)
              print(paste0("Required area: ",round(S.cond,2)))
              print(paste0("S_pol1: ", round(SpjuvioI,2)))
              print(paste0("S_pol2: ", round(SpjuvioII,2)))
            }
          }
        } else {
          if (any(trace.level == c(5,7)) ){
            lines(x = pairs_pts[i,c(1,3)], y = pairs_pts[i,c(2,4)],
                  col = "gray60", pch = 19)
            points(pairs_pts[i,c(3,4)],col = "gray60", pch = 19, cex = 1.5)
            print("Not enough data points in at least one of the polygons.")
            print(paste0("N required: ", N.crit))
            print(paste0("N1: ",nrow(Puses[[1]])))
            print(paste0("N2: ",nrow(Puses[[2]])))
          }
        }
      }
    }
    if(length(any.split) > 0){

      mean.dif <- mean(any.split) # mean spatial heterogeneity
      sd.dif <- sd(any.split) # NA if 1 split
      # sd(any.split) -> anisotropy of heterogeneity
      # if all are equal, then complete randomness: no anisotropy - no splits present
      if( all(any.split==any.split[1]) & length(any.split) > 1  ){
        warning(paste(c(
          "All tested splits are of equal quality. No anisotropy - no spilts in"
          , iteration, "iteration.")))
        maxid<-0
      }
    } else {
      mean.dif <- NA # negalejom ivertinti ne vieno padalinimo, taigi performance lygu max.
      sd.dif <- NA # P.crit
    }
    assign(x = "poly.info" ,value = rbind(poly.info, data.frame(
      mean.dif = mean.dif, # NA if 0
      sd.dif = sd.dif, # NA if 1 or 0 split
      str.z.score = (maxdif - mean.dif) / sd.dif, # NA if 1 or 0 split
      iteration = iteration,
      root = root
    ))
    ,envir = e)

    print(paste0("poly.info: ", poly.info))
    print(paste0("mean quality of straight splits: ", mean.dif))
    print(paste0("anysotropy of quality of straight splits: ", sd.dif))
    print((maxdif - mean.dif) / sd.dif)
    # duomenu saugojimas
    #Jei rastas tinkamas padalinimas - ieskom geriausios padalinimo kreives,
    #issaugom duomenis ir ziurim ar galima skaidyti toliau
    if (maxid>0) {
      if ( c.splits) {

        best.curve <- .curvial_split(

          poly.x = perim_pts[[2]]$x.poly,
          poly.y = perim_pts[[2]]$y.poly,
      min.x.id = pairs_pts[maxid,6],
      max.x.id = pairs_pts[maxid,7],
      b = pairs_pts[maxid,5],
      data = samp.dat,
      knot.density.X = knot.density.X,
      knot.density.Y = knot.density.Y,
      N.cond = N.crit,
      S.cond = S.cond,
      n.curve.iter = curve.iterations,
      correction.term = c.corr.term

      )
    if (trace.level > 2) {
      lines(best.curve[[1]],col=2,lwd=3)
    }
    if ( max(best.curve[[2]],maxdif) < upper.Q.crit ){ # vel santykinis base line. Be to,
      # galima gi reikalaut, kad atotrukis nuo base line butu tam tikro dydzio. Dar vienas P.crit
      # argumentas reikalingas?
      maxid <- 0
    }
    } else {
      if (maxdif < upper.Q.crit){
        maxdif <- 0
      }
    }
  }

  if (maxid > 0){ # save the split and perform new splits if TRUE
    curve.best <- FALSE
    if (c.splits) {
      if(best.curve[[2]] > maxdif) {
        curve.best <- TRUE
      }
    }
    if(curve.best) {
      best.splitl <- data.frame(x = best.curve[[1]]$x, y = best.curve[[1]]$y)
      maxdif <- best.curve[[2]]
    } else {
      best.splitl <- data.frame(x=as.numeric(c(pairs_pts[maxid,c(1,3)])),
                                y=as.numeric(c(pairs_pts[maxid,c(2,4)])))
    }
    assign(x = "splits" ,value = do.call(c,list(splits,list(best.splitl))),
           envir = e)
    assign(x = "split.quality" ,
           value = do.call(c,list(split.quality,maxdif )),
           envir = e)
    assign(x = "split.z.score" ,
           value = do.call(c,list(split.z.score,
                                  (maxdif - mean.dif)/sd.dif)),
           envir = e)
    # z-score of performance
    #dalinam duomenis padalinimo kreive
    #reikia sukurti poligonus du ir nufiltruoti duomenis  - galima padaryti geriau
    up.pol <- .close_poly(
      open.poly =
        .split_poly(
          polygon = data.frame(x=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],
                               y=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
          split_ids = as.numeric(pairs_pts[maxid,6:7]),
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
          split_ids = as.numeric(pairs_pts[maxid,6:7]),
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
    print(paste0('total max Skirtumas =', maxdif))

    #maisom duomenis ir vertinam aptiktu erdviniu strukturu patimuma
    if (n.m.test){

      set.seed(n.m.seed)
      test1 <- .sp.n.m(samp.dat,ribs,n.m.N,n.m.keep,type = 1)
      test2 <- .sp.n.m(samp.dat,ribs,n.m.N,n.m.keep,type = 2)

      if (n.m.keep){
        assign(x = "n.m.sims1",value = do.call(c,list(n.m.sims1,list(test1[[2]]))),
               envir = e)
        assign(x = "n.m.sims2",value = do.call(c,list(n.m.sims2,list(test2[[2]]))),
               envir = e)
      }
      assign(x = "sim2.difs",value = do.call(c,list(sim2.difs,list(test2[[1]]))),
                                          envir = e)
      assign(x = "sim1.difs",value = do.call(c,list(sim1.difs,list(test1[[1]]))),
             envir = e)

      assign(x = "p.val1", value =
               do.call(c,list(
                 p.val1,
        sum(maxdif < test1[[1]])/n.m.N
        )),
        envir = e) # kvantilis
      # kvantilio reiksme - empirine p verte
      # jei gausinis skirstinys,tai:
      # 1 - qnorm(sum(last(split.quality),mean(pseudo.kokybe),sd(pseudo.kokybe)) duotu teorine p-verte
      assign(x = "p.val2", value =
               do.call(c,list(
                 p.val,
                 sum(maxdif<test2[[1]])/n.m.N
               )),
             envir = e)
      }
  # updatinam ka reikia in hespdiv env.
    assign(x = "n.splits",value = do.call(c,list(n.splits,length(any.split))),
           envir = e)

    assign(x = "mean.difs",
           value = do.call(c,list(mean.dif,mean.dif)),
           envir = e)
    # kaip su situo rims blet?
    assign(x = "rims" ,value = do.call(c,list(rims,ribs[1])) ,envir = e)
    assign(x = "plot.id",value = do.call(c,list(plot.id,iteration)),
           envir = e)

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
    assign(x = "rims" ,value = do.call(c,list(rims,ribs[2])) ,envir = e)

    if (trace.level > 1) {
      lines(ribs[[2]],col="purple")
    }
    .spatial_div(do.dat, root = iteration)
    print(paste("griztam i", testid, "padalinima [po aukstu koord bloko (gogolis exists)]", sep=" "))

  }} else{
    #Jei tinkamo padalino nerasta, grizta tuscias masyvas
    if (testid>1){
      if (trace.level > 0)
      print("There were no suitable points on the polygon perimeter
            to generate straight-split lines")
    } else{
      if (trace.level > 0)
      print("There were no suitable points on the polygon perimeter
            to generate straight-split lines")
    }
  }
}

#' Calculate difference between two data sets
#'
#' @description  Function combines generalize.f and compare.f functions, to
#' estimate the difference between two data sets obtained from opposite sides
#' of a split-line. This difference reflects the ability of a split-line to
#' spatially separate data. Thus, it is considered to reflect the quality of a
#' split-line.
#' @param dat1 a spatial subset of \code{samp.dat} that is located at opposite
#' side of a split-line than \code{dat2} is.
#' \code{samp.dat = data}.
#' @param dat2 a spatial subset of \code{samp.dat} that is located at opposite
#' side of a split-line than \code{dat1} is.
#' @return Numeric value, expressing the difference between two data sets and
#' the quality of a split-line.
#' @author Liudas Daumantas
#' @note Function inherits the environment of \code{hespdiv} function. Also,
#' it forces this environment to inherited by generalize.f and compare.f
#' functions.
#' @noRd
.dif_fun <- function(dat1,dat2) {
  environment(generalize.f) <- environment()
  environment(compare.f) <- environment()
  compare.f( generalize.f(dat1), generalize.f(dat2) )
}

#' Print the results of hespdiv object
#'
#' @description  Function formats and prints the results
#' of R object of class "hespdiv". It prints rounded split.stats data frame.
#' @param x hespdiv object
#' @param n.m.test logical - were splits tested using null model simulations?
#' @return x
#' @author Liudas Daumantas
#' @noRd
print.hespdiv <- function(x, n.m.test){
  cat("\n","Information about the splits:", "\n","\n")
  print(round(x$split.stats,2))
  if (n.m.test){
    cat("\n", "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  }
  invisible(x)
}

#' Test the split-line with null model simulations
#'
#' @description  Function simulates spatial null models of the provided data and
#' then checks the performance of the established split-line fitted to the
#' simulated data.
#' @param data a samp.dat from .spatial_div function.
#' @param ribs a list of two data frames of polygons, established by the
#' fitted, "best" split-line.
#' @param n integer - number of data simulations to perform.
#' @param n.m.keep logical - should the produced simulations be kept
#' @param type type of null model simulations:
#' 1 - completely random: toroidal rotation of individual data points (the only
#' things maintained is polygon shape and number of points, and points
#' themselves),
#' 2 - random perspective of data points: toroidal rotation of all data points
#' (micro and macro spatial inter-relations of points are maintained)
#' Other possible n.m. types:
#' 3 - random macro inter-relations of data points: toroidal rotation of
#' data point groups, that lie in the same quadrant of polygon extent (micro
#' inter-relations of points are maintained) # quadrant number would be required
#' as an argument.
#' 4 - random micro inter-relations of data points: jitter of data
#' points. The amount of jittery equal to points' distance to the xth
#' nearest neighbor (density surface (intensity, first moment) as well as macro
#' inter-relations of points are maintained) # xth nearest neighbour would be
#' required as an argument.
#' 5?? - random micro and macro inter-relations of data points: rotation of
#' quadrants and jitter of the points (maintained are point inter-relations that
#' are occur between the given micro and macro scales.)
#' @return list of one or two elements: 1 - vector of estimated performances of
#' the split-line for each of null model simulation. 2 -  list of data frames of
#' simulated data (null models), returned only when n.m.keep is TRUE.
#' @author Liudas Daumantas
#' @noRd
.sp.n.m <- function(data,ribs,n,n.m.keep,type){
    if (type == 1) {
      N <- nrow(data)
    } else {
      N <- 1
    }
    mirror.data <- data
    pseudo.kokybe <- numeric(n.m.N)
    if (n.m.keep){
      sim <- list()
    }
  for (i in 1:n){
    x.shift <- runif(N,min=0,max=dist(range(data$x)))
    y.shift <- runif(N,min=0,max=dist(range(data$y)))
    testx <- data$x+x.shift
    testy <- data$y+y.shift
    tx <- case_when(testx>max(data$x) ~ testx-max(data$x)+min(data$x),
                  TRUE ~ testx)
    ty <- case_when(testy>max(data$y) ~ testy-max(data$y)+min(data$y),
                  TRUE ~ testy)

    mirror.data[,c("x","y")] <- data.frame(x=tx,y=ty)
    I.dat <- .get_data(ribs[[1]], mirror.data)
    II.dat <- .get_data(ribs[[2]], mirror.data)

    pseudo.kokybe[i] <- dif.fun(I.dat, II.dat)

    if (n.m.keep){
    sim <- do.call(c,list(sim1,list(mirror.data)))
    }
  }
    if (n.m.keep){
      return(list(pseudo.kokybe,sim))
    } else {
      return(list(pseudo.kokybe))
    }
}
