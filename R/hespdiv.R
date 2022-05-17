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
#' @param data an R object containing data to be analyzed. Required data
#' structure depends on the method selected (see ...) or created.
#' @param xy.dat a data.frame containing coordinates for objects in \code{data}.
#'  If \code{data} is data frame or matrix that has columns x or y, can be
#' ignored (default NULL).
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
#' @param method String. A name of preset method. Determines which
#' preset \code{generalize.f} and \code{compare.f} functions will be used to
#' perform hierarchical subdivisions of spatial data, if \code{compare.f}
#' function is not provided.
#' @param similarity logical. Does the named metric reflect similarity? Needs
#' only to be provided when custom \code{compare.f} function is used.
#' @param maximize logical. Should the value returned by \code{compare.f} be
#' maximized or minimized? Needs only to be provided when custom
#' \code{compare.f} function is used.
#' @param N.crit Subdivision stopping criteria - number of observations.
#' Minimum number of observations (rows in data) that should
#' be present in areas separated by a split-line in order to establish the
#' split-line. Default is 1.
#' @param S.crit Subdivision stopping criteria - size of plots.
#' Minimum area expressed as a proportion of original study area
#' (provided polygon or estimated as convex hull of \code{xy_dat}) that plots
#' separated by a split-line should have so that the split-line could be
#' established. Default is 0.
#' @param lower.Q.crit integer (default NULL). Only meaningful when c.splits is
#' TRUE. lower.Q.crit determines the minimum performance that the best obtained
#' straight-split line should have so that attempt would be made to generate
#' from it a better non-linear split-line. Recommendation is to leave this value
#' set to default (no performance requirement) and only use different value,
#' when you are quite sure what are the limitations of improvements that
#' non-linear split-lines can make over straight split-lines (e.g. if maximize
#' is TRUE, lower.Q.crit = upper.Q.crit - MAX.c.improv).
#' @param upper.Q.crit Subdivision stopping criteria - upper limit of split-line
#' quality applied to the final split-line. This is a minimum difference as
#' estimated by \code{compare.f} function that separated
#' plots should exhibit, so that a subdivision of a plot using the
#' best split-line would be established. Default is -Inf.
#' @param c.splits Logical (default TRUE).
#' Should non-linear split-lines be estimated?
#' @param C.crit.improv integer. How much non-linear split must be
#' better than straight (in units of provided metric) for it to be selected?
#' When 0 is set (default), then non-linear split-line still won't be selected
#' if it performs the same as straight split line.
#' @param c.X.knots Curve parameter. The number of columns in the net of
#' spline knots. These columns are distributed regularly along to the straight
#' split-line. This parameter controls wiggliness (wave length) of the
#' curvi-linear split-lines. Higher values allow  wigglier curves, thus
#' increasing the fit to the data, but increases the optimization time of
#' curvi-linear split-lines. Default value is 5.
#' @param c.Y.knots Curve parameter. The number of rows in the net of
#' spline knots. These rows are distributed regularly orthogonal to the straight
#' split-line. This parameter controls wiggliness (resolution of tested wave
#' amplitudes) of the curvi-linear split-lines. Higher values allow higher
#' variety of wave amplitudes to be tested, when optimizing the shape of
#' curvi-linear split-lines. Thus higher values increase the fit to the data at
#' the cost of optimization time. Default value is 10.
#' @param c.max.iter.no the maximum number of allowed iterations through the
#' columns of spline knots (default 5). Setting higher values (eg Inf) increases
#' the odds of convergence to the best possible curve. Higher values recommended
#' if c.fast.optim is FALSE, because then only one spline knot can be
#' changed each iteration.
#' @param c.fast.optim Logical (default TRUE). Determines when spline knots
#' are selected: TRUE - when first better curve than before is obtained, FALSE -
#' when the best curve is obtained (after all spline knots are checked).
#' @param c.corr.term The term that defines the correction size of problematic
#' curvi-linear split-lines which intersect the boundary of the polygon.
#' Possible values are between 0 and 1, though small values are recommended
#' (default is 0.05). These values define how much the interval of a generated
#' spline, that crosses the boundary of polygon, should be shifted away from the
#' boundary, inside the polygon, in direction orthogonal to the straight
#' split-line, in  terms of proportion of polygon width where spline intersects
#' the polygon boundary.
#' @param filter.all logical (default TRUE). Defines which data points are
#' filtered by polygons. FALSE - only points strictly inside the polygon are
#' filtered by polygons. TRUE - points located on a split-line are filtered as
#' well in addition to those strictly inside a polygon.
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
#' @param trace.level Sting indicating the trace level. Can be one of the
#' following: "best", "main", "all"
#' @param trace.object String that informs what graphical information to show.
#' Can be either "curves", "straight" or "both".
#' @param pnts.col Color of data points, default is 1. Argument is used when
#' \code{trace} > 0. If set to NULL, data points will not be displayed.
#' @param display logical. Display the resulting polygons at the output?
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
#'   \item \code{root.id} - the ID of \code{hespdiv} iteration, which produced the
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
                  maximize = NULL, method = NULL,
                  compare.f = NULL, N.crit = 0,
                  S.crit = 0, lower.Q.crit, upper.Q.crit,
                  c.splits = TRUE, C.crit.improv = 0, c.X.knots = 5,
                  c.Y.knots = 10, xy.dat = NULL,
                  c.max.iter.no = 5, c.fast.optim = TRUE,
                  c.corr.term = 0.05, filter.all = TRUE,
                  n.m.test = FALSE,
                  n.m.N = 1000, n.m.seed = 1,  n.m.keep = FALSE,
                  study.pol = NULL, trace.level = NULL,
                  trace.object = NULL, pnts.col = 1, display = TRUE){

  if ((is.null(trace.level) & !is.null(trace.object)) |
      (!is.null(trace.level) & is.null(trace.object))){
    stop(paste("Conflicting arguments: trace.level, trace.object"),
         paste("\ntrace.object and trace.level must be both set to 'NULL'"),
         paste( ' or asigned a viable value'))
  }

  args <- sapply(ls(),get,environment())

  if (is.null(xy.dat)){
    if (class(data) %in% c("data.frame", "matrix")){
      if (any(colnames(data) == 'x') & any(colnames(data) == 'y')){
        if (sum(colnames(data) %in% c('x','y')) != 2){
          stop("Data must contain no more than one of columns that are named
               'x' and 'y'")
        }
        xy.dat <- as.data.frame(data[,c('x','y')])
        data <- data[,-which(colnames(data) %in% c('x','y'))]
      } else {
        stop("'xy.dat' was not provided and 'data' is missing 'x', 'y' columns",
             " containing coordinates of observations.")
      }
    } else {
      stop("'xy.dat' was not provided and 'data' is neither a data frame or",
           "a matrix.")
    }
  } else {
    if (!is.data.frame(xy.dat)){
      if (is.matrix(xy.dat)){
        xy.dat <- as.data.frame(xy.dat)
      } else {
        stop("'xy.dat' must be a data.frame.")
      }
    }
    if (ncol(xy.dat) !=2 )
      stop("number of columns in 'xy.dat' must be 2, named 'x' and 'y'.")

    if (any(colnames(xy.dat) != c('x','y')))
      stop("there must be 'x' and 'y' columns in xy.dat positioned in this",
           " order.")

    if (is.matrix(data) | is.data.frame(data)){
      if (nrow(data) != nrow(xy.dat))
        stop("number of rows in 'xy.dat' must identical to the number ",
             "of observations in data.")
    } else {
      if (length(data) != nrow(xy.dat))
        stop("number of rows in 'xy.dat' must identical to the number ",
             "of observations in data.")
    }
  }

  if (!is.null(trace.object)){
    trace.object <- .arg_check("trace.object", trace.object,
                               c("straight", "curve", "both"))
  }
  if (!is.null(trace.level)){
    trace.level <- .arg_check("trace.level", trace.level,c("all","main","best"))
  }

  if (is.null(compare.f)){
    if (is.null(method)){
      stop("Neither 'method' or 'compare.f' argument is specified. ",
           "Please select a preset method or provide a custom method.")
    }
    method.type <- "preset"
    method <- .arg_check("metric", method, names(
      .get_methods()[['biozonation']]))
    if (method == "pielou") {
      variant <- '1'
      metric <- "pielou"
      method <- 'biozonation'
      similarity <- FALSE
      maximize <- TRUE
      if(!is.vector(data))
        stop("Data must be a vector, when using '",method,"' method,",
             " '",metric,"' metric and '", variant,"' variant.")

      generalize.f <- function(plot.dat){
        x <- table(factor(plot.dat))
        p <- x/sum(x)
        -sum(log(p)*p)
      }
      compare.f <- function(eveness1,eveness2) {
        base.eveness <- poly.obj[[testid]]
        (1 - mean(c(eveness1, eveness2)) / base.eveness) * 100 # percent change in eveness
      }

    } else {
      if (method == "sorensen"){
        similarity <- TRUE
        variant <- '1'
        metric <- 'sorensen'
        method <- 'biozonation'
        maximize <- FALSE
        if(!is.vector(data))
          stop("Data must be a vector, when using '",method,"' method,",
               " '",metric,"' metric and '", variant,"' variant.")

        generalize.f <- function(plot.dat){
          unique(plot.dat)
        }
        compare.f <- function(uniq_tax1,uniq_tax2) {
          sum <- length(uniq_tax1) + length(uniq_tax2)
          int_2x <- length(which(duplicated(c(uniq_tax1,uniq_tax2))))*2
          if (length(int_2x)!=0){
            int_2x/sum
          } else {
            0
          }
        }
      } else {
        if (method == "morisita"){
          similarity <- TRUE
          metric <- "morisita"
          maximize <- FALSE
          variant <- '1'
          method <- 'biozonation'
          if(!is.vector(data))
            stop("Data must be a vector, when using '",method,"' method,",
                 " '",metric,"' metric and '", variant,"' variant.")

          generalize.f <- function(plot.dat){
            plot.dat
          }

          compare.f <- function(x,y) {
            all_sp <- unique(c(x,y))
            x_f <- factor(x,levels = all_sp)
            y_f <- factor(y,levels = all_sp)
            (2*sum(table(x_f) * table(y_f)))/
              (length(x) * length(y) *
                 ((sum(table(x_f)*(table(x_f)-1)) /
                     (length(x)* (length(x)-1))) +
                    (sum(table(y_f)*(table(y_f)-1)) /
                       (length(y)* (length(y)-1)))))
          }
        } else {
          if (method == "horn.morisita"){
            similarity <- TRUE
            metric <- "horn.morisita"
            maximize <- FALSE
            variant <- '2'
            method <- 'biozonation'
            if(!is.vector(data))
              stop("Data must be a vector, when using '",method,"' method,",
                   " '",metric,"' metric and '", variant,"' variant.")

            generalize.f <- function(plot.dat){
              plot.dat
            }

            compare.f <- function(x,y) {
              all_sp <- unique(c(x,y))
              x_f <- factor(x,levels = all_sp)
              y_f <- factor(y,levels = all_sp)
              (2*sum(table(x_f) * table(y_f))) /
                (length(x) * length(y) *
                   ( sum(table(x_f)^2) / (length(x)^2)  +
                        sum(table(y_f)^2) / (length(y)^2)))
            }
          }
        }
      }
    }
  } else {
    if (is.null(generalize.f))
      generalize.f <- function(x) x
    method.type <- "custom"
    variant <- NULL
    metric <- NULL
    similarity <- NULL
    if (is.null(maximize)){
      stop("Missing maximize argument value.")
    }
  }

  if (is.data.frame(data) | is.matrix(data)){
    .slicer <- .slicer.table
  } else {
    if (is.list(data)) {
      .slicer <- .slicer.list
    } else {
      if (!is.vector(data))
        warning("Data structure of 'data' is not recognized by 'hespdiv'.",
                " Slicing of 'data' are attempted using '[]' brackets.")
      .slicer <- .slicer.vect
    }
  }


  if (!maximize){
    if (is.null(lower.Q.crit))
      lower.Q.crit <- +Inf
    .comp <- function(x,criteria){ x < criteria}
    c.sign <- "<"
    .minormax <- min
    .which_minormax <- which.min
  } else {
    if (is.null(lower.Q.crit))
      lower.Q.crit <- -Inf
    .comp <- function(x,criteria){ x > criteria}
    c.sign <- ">"
    .minormax <- max
    .which_minormax <- which.max
  }

  if ((c.splits == FALSE & upper.Q.crit != lower.Q.crit) |
      ((upper.Q.crit > lower.Q.crit) & maximize) |
      ((upper.Q.crit < lower.Q.crit) & !maximize) ){
    if ((c.splits == FALSE & upper.Q.crit != lower.Q.crit)){
      warning(paste("Since 'c.splits' is FALSE, 'lower.Q.crit' was set equal to
          'upper.Q.crit'"))
      lower.Q.crit <- upper.Q.crit
    }
    if ((upper.Q.crit > lower.Q.crit) & !maximize){
      warning(paste("upper.Q.crit should be lower or equal to lower.Q.crit,",
                    "when optimization is reached by minimazing the metric."))
    }
    if ((upper.Q.crit < lower.Q.crit) & maximize){
      warning(paste("upper.Q.crit should be higher or equal to lower.Q.crit,",
                    "when optimization is reached by maximising the metric."))
    }
  }

  if (is.null(study.pol)){
    ### is x, y ids ch ever used again ????
    ch <- chull(xy.dat$x, xy.dat$y)
    ids <- c(ch,ch[1])
    x <- xy.dat$x[ids]
    y <- xy.dat$y[ids]
    study.pol <- data.frame(x=x,y=y)
  }
  first.p <- study.pol[1,]
  if (filter.all) {
    # ... is used to ignore split_endpnts, when they are added.
    .get_ids <- function(polygon, xy_dat,...){
      which(sp::point.in.polygon(pol.x = polygon[,1],
                           pol.y = polygon[,2],
                           point.x = xy_dat$x,
                           point.y = xy_dat$y) != 0)
    }
  } else {
    .get_ids <- function(polygon, xy_dat,first.p,split_endpnts){
      if (any((first.p[,1] == split_endpnts[c(1,2),1]) &
              (first.p[,2] == split_endpnts[c(1,2),2])) ){
        del.id <- which(xy_dat$x == first.p[,1] &
                          xy_dat$y == first.p[,2])
        if (length(del.id)>0)
          return(which(sp::point.in.polygon(pol.x = polygon[,1],
                                      pol.y = polygon[,2],
                                      point.x = xy_dat$x[-del.id],
                                      point.y = xy_dat$y[-del.id])
                 != 0 ))

      }
      which(sp::point.in.polygon(pol.x = polygon[,1],
                           pol.y = polygon[,2],
                           point.x = xy_dat$x,
                           point.y = xy_dat$y) != 0)

    }
  }

  rims <- list(study.pol)

  poly.info <- data.frame(root.id = numeric(),
                          n.splits = numeric(),
                          n.obs = numeric(),
                          mean = numeric(), # mean spatial heterogeneity irrespective of split-line position
                          sd = numeric(), # anizotropy of heterogeneity based on straight split-lines
                          str.best = numeric(),
                          str.z.score = numeric(), # level of outstandingness. Are there other competetive candidate splits?,
                          has.split = logical(),
                          is.curve = logical(),
                          crv.best = numeric(),
                          crv.z.score = numeric(),
                          c.improv = numeric())

  poly.obj <- list()
  str.split.quals <- list()

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

  S.cond <- round(abs(pracma::polyarea(study.pol$x,study.pol$y)) * S.crit,2)
  splits <- numeric()

  e <- environment()
  environment(.spatial_div) <- e
  ### obtaining results:
  .spatial_div(data,xy.dat,root.id=0)
  ### results obtained. Formatting them.
  plot.id <- 1:nrow(poly.info)
  names(poly.obj) <- plot.id
  names(rims) <- plot.id
  poly.info <- data.frame(plot.id = plot.id,poly.info)


  # which best is not clear, when no split was above lower.Q.crit
  if (any(poly.info$str.best == lower.Q.crit)){
    plot.id.with.ns <- which(poly.info$n.splits != 0)
    poly.info$str.best[plot.id.with.ns] <- ifelse(poly.info$str.best[
      plot.id.with.ns] == lower.Q.crit,
      lapply(str.split.quals[plot.id.with.ns],.minormax),
      poly.info$str.best[plot.id.with.ns])
    poly.info$str.best[-plot.id.with.ns] <- NA
  }
  if (!c.splits){
    poly.info <- poly.info[,-c(10:13)]
  } else {
    # if no data about curves, then we don't know how curves compare to straight
    poly.info$is.curve <- ifelse(is.na(poly.info$crv.best),NA,
                                 poly.info$is.curve)
    # is no improvement was made by the curve, then we dont know the performance
    # of the best curve since worse values than straight-split performance
    # are not saved.
    poly.info$crv.best <- ifelse(poly.info$c.improv == 0,NA,
                                 poly.info$crv.best)
    poly.info$crv.z.score <- ifelse(poly.info$c.improv == 0,NA,
                                    poly.info$crv.z.score)
  }
  if (any(poly.info$has.split)){
    names(splits) <- 1:sum(poly.info$has.split)
  }

  plot.id <- which(poly.info$has.split)
  e <- environment()
  result <- .format_result(e)

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

  if (display){
    .visualise_splits.end(pnts.col, xy.dat, rims)
  }

  return(print.hespdiv(result))
}
.slicer.table <- function(sample,ind) sample[ind,]
.slicer.list <- function(sample,ind) sample[[ind]]
.slicer.vect <- function(sample,ind) sample[ind]

.format_result <- function(e){
  if (e$method.type == "custom"){
    return(.format_result.general(e))
  } else{
    if (e$metric == "pielou"){
      return(.form.bioz.Pielou(.format_result.general(e)))
    } else {
      return(.format_result.general(e))
    }
  }
}
.structurise <- function(e,split.stats,split.lines){

  structure(list(
    split.lines = split.lines,
    polygons.xy = e$rims,
    poly.stats = e$poly.info,
    poly.obj = structure(e$poly.obj),
    split.stats = structure(split.stats),
    call.info = list(METHOD =
                       list(method = e$method, metric = e$metric,
                            variant = e$variant, similarity = e$similarity,
                            maximize = e$maximize, method.type = e$method.type),
                     Call_ARGS = e$args),
    str.difs = e$str.split.quals
  ),
  class = "hespdiv"
  )
}
.form.bioz.Pielou <- function(result){
  parent.Pe <- unlist(result$poly.obj)
  result$poly.stats$parent.Pe <- parent.Pe
  if (!is.null(result$split.stats) ){
    result$split.stats$parent.Pe <- parent.Pe[result$split.stats$plot.id]
    result$split.stats$delta.Pe <- -result$split.stats$parent.Pe *
      result$split.stats$performance/100
  }
  result
}
.format_result.general <- function(e){
  if (length(e$plot.id) == 0 ){
    split.stats <- NULL
    split.lines <- NULL
  } else {
    split.lines <- e$splits
    if (e$c.splits){
      split.stats <- data.frame(
        plot.id = e$plot.id,
        n.splits = e$poly.info$n.splits[e$plot.id],
        n.obs = e$poly.info$n.obs[e$plot.id],
        mean = e$poly.info$mean[e$plot.id],
        sd = e$poly.info$sd[e$plot.id],
        z.score = apply(e$poly.info[e$plot.id,],1,
                        function(o) if (o[[10]]) {o[[12]]} else {o[[8]]}),
        performance = apply(e$poly.info[e$plot.id,],1,
                     function(o) if (o[[10]]) {o[[11]]} else {o[[7]]}),
        is.curve = e$poly.info$is.curve[e$plot.id]
      )
    } else {
      split.stats <- data.frame(
        plot.id = e$plot.id,
        n.splits = e$poly.info$n.splits[e$plot.id],
        n.obs = e$poly.info$n.obs[e$plot.id],
        mean = e$poly.info$mean.en.p.red[e$plot.id],
        sd = e$poly.info$sd.en.p.red[e$plot.id],
        z.score = e$poly.info[e$plot.id,"str.z.score"],
        performance = e$poly.info[e$plot.id,"str.best"]
      )
    }
  }
  .structurise(e,split.stats,split.lines)
}

.arg_check <- function(name, given,NAMES){
  matched.i <- pmatch(given, NAMES)
  if(is.na(matched.i))
    stop("invalid ",name, " ", paste0("'", given,"'."),
         paste('\nPlease select viable', name,": "),
         paste(NAMES,collapse = ", ",sep = "'"))
  NAMES[matched.i]
}

.get_methods <- function(){
  list(
    "biozonation" = list(
      "pielou" = c("1"),
      "morisita" = c("1"),
      "sorensen" = c("1"),
      "horn.morisita" = c("2")
    )
  )
}
