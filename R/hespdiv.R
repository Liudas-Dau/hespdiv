#' Hierarchically subdivide spatial data
#'
#' This function is an implementation of spatial data analysis method "HespDiv".
#' It performs hierarchical spatial data subdivision by recursively
#' dividing the data using random split-lines, evaluating their comparison
#' values (how well they separate data), and using the best to perform
#' subdivisions.
#
#' @param data An R object that contains the data to be analyzed. The required data
#' structure depends on the selected method (e.g., character vector of taxa names).
#' @param xy.dat  A data.frame containing the coordinates for observations in the
#' \code{data}. This parameter can be ignored if \code{data} is a data frame or
#'  matrix with columns \code{x} or \code{y}.
#' @param n.split.pts integer. The number of split-points - 1. These points are
#' used in creating straight split-lines (see details). The total number of straight
#' split-lines generated can be obtained by \code{sum(1:n.split.pts)}.
#' Increasing the value of \code{n.split.pts} leads to an increase in both the
#' computation time and the fit to the data.
#' @param same.n.split logical. Should the number of split-points
#' (\code{n.split.pts}) remain the same within each lower-ranked polygon?
#' Choosing \code{TRUE} would result in a higher density of straight split-lines
#' within lower-ranked polygons, whereas \code{FALSE} would preserve the same
#' split-line density. Thus, essentially you are choosing between cross-scale
#' and fixed-scale analysis.
#' @param method character. A name (or its abbreviation) of a preset method:
#' 'sorensen', 'pielou', 'morisita' or 'horn.morisita'. Internally determines
#' values for \code{compare.f}, \code{generalize.f} and \code{maximize}.
#' @param generalize.f function. Optional function used in custom methods to
#' prepare input for \code{compare.f} function (see details). It requires a
#' single argument (e.g. "x") to which a spatially filtered subset of the
#' \code{data} could be assigned (see notes for exceptions).
#' @param compare.f function. Only required in custom methods. Employed to
#' quantify the comparison value of a split-line. For this purpose, the
#' \code{compare.f} function requires two arguments where the outputs of the
#' \code{generalize.f} function can be assigned (see notes for exceptions).
#' @param maximize logical. Only required in custom methods. Determines whether
#' the split-line comparison value should be maximized or minimized during the
#' optimization process.
#' @param N.crit number. Minimum required number of observations that should
#' be present in each polygon obtained by a split-line in order for it to meet
#' this criterion.
#' @param N.rel.crit number from 0 to 0.5. Each polygon obtained with a
#' split-line must have at least such proportion of observations to pass this
#' criterion. Equation of the proportion:  (Number of observations in 1st/2nd
#' resulting polygon) / (Number of observations in the polygon being subdivided)
#' @param N.loc.crit number. Minimum required number of different locations that
#' should be present in each polygon obtained by a split-line in order for it to
#' meet this criterion.
#' @param N.loc.rel.crit number from 0 to 0.5. Each polygon obtained with a
#' split-line must have at least such proportion of different locations to pass this
#' criterion. Equation of the proportion:  (Number of different locations in 1st/2nd
#' resulting polygon) / (Number of different locations in the polygon being subdivided)
#' @param S.crit number from 0 to 1. Each polygon obtained with a
#' split-line must have at least such area proportion to pass this
#' criterion. Equation of the proportion:  (Area of 1st/2nd resulting polygon) /
#' (Area of the first polygon). The first polygon is the provided study area
#' or convex hull of observation locations.
#' @param S.rel.crit number from 0 to 0.5. Each polygon obtained with a
#' split-line must have at least such area proportion to pass this
#' criterion. Equation of the proportion:  (Area of 1st/2nd resulting polygon) /
#' (Area of the polygon being subdivided).
#' @param Q.crit number. The threshold for a split-line comparison value to be
#' considered acceptable for a subdivision. When \code{maximize = TRUE}, higher
#' values of the comparison value indicate a better subdivision. Conversely, when
#' \code{maximize = FALSE}, lower values of the comparison value indicate a
#' better subdivision.
#' @param c.splits logical. When set to TRUE, the algorithm will explore
#' nonlinear split-lines in addition to straight split-lines in order to
#' find the optimal subdivision.
#' @param c.Q.crit number.  The threshold for a split-line comparison value to
#' be considered acceptable for generating nonlinear split-lines. It is
#' recommended to use the default value, which does not impose a performance
#' requirement, unless you have a clear understanding of the potential
#' improvements that nonlinear split-lines can achieve over straight split-lines.
#'  If \code{maximize = TRUE}, a suggested value for \code{c.Q.crit} is
#'  \code{Q.crit} minus the maximum potential improvement.
#' @param c.crit.improv number. The threshold for the improvement in a
#' split-line comparison value required for a nonlinear split-line to be
#' selected instead of a straight split-line for subdivision. The default value
#' of 0 means that even if a nonlinear split-line performs equally to a
#' straight split-line, the straight split-line will still be chosen.
#' @param c.X.knots integer. Specifies the number of columns in a network of
#' spline knots used to generate nonlinear split-lines. These knots are evenly
#' distributed along the straight split-line. Adjusting the value of
#' \code{c.X.knots} controls the degree of wiggliness (number of turns) in the
#' resulting nonlinear  split-lines.
#' @param c.Y.knots integer. specifies the number of  rows in a network of
#' spline knots used to generate nonlinear split-lines. These knots are
#' distributed regularly along lines orthogonal to the straight
#' split-line. Adjusting the value of \code{c.Y.knots} controls the number
#' of amplitudes tested in each wiggle of a spline, influencing the shape of
#' nonlinear split-lines.
#' @param c.max.iter.no integer. The maximum number of iterations allowed through
#' the network of spline knots when searching for the optimal shape of a
#' nonlinear split-line. Setting a higher value, such as \code{+Inf}, increases
#' the chances of converging to the best possible curve. It is recommended to
#' use higher values when \code{c.fast.optim = FALSE}, as in this case, only a
#' single spline knot can be selected per iteration.
#' @param c.fast.optim logical. Determines when spline knots are selected. If
#' \code{TRUE}, the algorithm selects the first knot that generates a curve with
#' a better comparison value, ensuring faster convergence. If set to
#' \code{FALSE}, the algorithm completes a full iteration through the network
#' of spline knots, which may result in slower convergence but potentially
#' better overall performance.
#' @param c.corr.term number from 0.01 to 0.2. A correction term for nonlinear
#' split-lines that intersect the boundary of the polygon. Smaller values
#' (default is 0.05) are recommended, as they determine the extent to
#' which the outlying interval of the generated spline, which crosses the
#' polygon boundary, should be shifted away from the boundary and inside the
#' polygon in a direction orthogonal to the straight split-line. This shift is
#' specified as a proportion of the polygon height where the spline intersects
#' the polygon boundary.
#' @param study.pol data frame with two columns, 'x' and 'y'. Rows of this data
#' frame should correspond to vertexes of a polygon that defines the study area
#' encompassing all observation locations (\code{xy.dat}). If not provided
#' (default is NULL), convex hull of \code{xy.dat} will be used as
#' study area polygon.
#' @param use.chull logical. If \code{study.pol} is provided, you can use
#' \code{use.chull = TRUE} (which is default) to use it only for visualization.
#' @param tracing a character vector with two elements. The first element
#' indicates the level of tracing, which can be "best", "main", or "all". The
#' second element specifies the object to be traced, which can be "curves",
#' "straight", or "both". By default, when set to NULL, no tracing will be
#' performed.
#' @param pnts.col character or numeric. Specifies the color of observations
#' in a plot. The argument is used when \code{tracing} is not NULL. If
#' \code{pnts.col = NULL}, observations will not be displayed.
#' @param display logical. Display a simple plot of results at the end of
#' computations?
#' @param pacific.region logical (default is FALSE). When set to TRUE, indicates
#' that the study area is crossed by the 180th meridian, such as being within
#' the Pacific Ocean. In this case, the coordinates of  \code{xy.dat} and
#' \code{study.pol} are transformed to eliminate the artificial abrupt change
#' in x-coordinate values at the 180th meridian.
#' @return The 'hespdiv' class object eith a list with seven elements:
#' \describe{
#'  \item{\code{poly.stats}}{ The data frame containing information about
#'   hespdiv polygons established by the selected split-lines. Columns of the
#'   data frame:}
#'   \itemize{
#'   \item \code{rank} - The rank of a polygon. It corresponds to the rank of the
#'   split-line that produced the polygon and polygon position in the
#'    hierarchical structure of the subdivisions.
#'   \item \code{plot.id} - The ID assigned to the polygon. Corresponds to the
#'   order within hespdiv analysis in which it was processed.
#'   \item \code{root.id} - The ID of the parent polygon whose subdivision
#'   resulted in the current polygon.
#'   \item \code{n.splits} - The count of straight split-lines that were
#'   evaluated for comparison values in an attempt to subdivide the polygon.
#'   This count excludes split-lines that crossed the polygon boundary or did
#'   not meet area, sample size or location number criteria.
#'   \item \code{n.obs} - The number of observations (e.g., data points) inside
#'   the polygon.
#'   \item \code{mean} - The average comparison value of the straight
#'   split-lines used in the attempted subdivision of the polygon. This value
#'   reflects the general spatial heterogeneity of the data
#'   \item \code{sd} - The standard deviation of the comparison values of the
#'   straight split-lines used in the attempted subdivision. It indicates
#'   the extent of anisotropy or variation in the spatial heterogeneity of
#'   the data within the polygon.
#'   \item \code{str.best} - The comparison value of the best straight
#'   split-line produced within the polygon.
#'   \item \code{str.z.score} - The z-score of the comparison value of the best
#'   straight split-line within the polygon. It indicates how outstanding
#'   the best straight split-line is compared to other evaluated random straight
#'   subdivisions. Normally, this value should rarely reach extremely low or
#'   high values (+-2-3) since the comparison values of split-lines are
#'   spatially autocorrelated.
#'   \item \code{has.split} - A Boolean indicator (TRUE or FALSE) that shows
#'   whether a subdivision was established in the polygon.
#'   \item \code{is.curve} - A Boolean indicator (TRUE or FALSE) that shows
#'   whether the established subdivision was obtained using a curve.
#'   If \code{has.split} is FALSE, this column will have an NA value.
#'   \item \code{crv.best} - The same as \code{str.best} but for nonlinear split-lines.
#'   \item \code{crv.z.score} - The same as \code{crv.z.score} but for nonlinear
#'   split-lines.
#'   \item \code{c.improv} - The improvement in comparison values
#'   achieved by using nonlinear split-lines over straight split-lines.
#'   }
#'   \item{\code{split.stats}}{ A data frame containing information about the
#'   established split-lines. Columns of this data frame can be interpreted
#'   the same way as in \code{poly.stats}, expect through the perspective of the
#'   split-line (e.g. \code{rank} - is the rank of the split-line, not the
#'    polygon). Additionally, \code{performance} column is the comparison value
#'    of the split-line, and is the same as \code{str.best} or \code{crv.best} in \code{poly.stats},
#'    depending on the value of is.curve column. The same applies to the z.score
#'    column.}
#'   \item{\code{split.lines}}{ A list containing data frames with the
#'   coordinates of the established split-lines. The order of split-lines
#'   in this list corresponds to their order in the \code{split.stats} data frame..}
#'   \item{\code{polygons.xy}}{ A list containing data frames of the established
#'   polygons resulting from the split-lines. The order of polygons in this
#'   list corresponds to their order in the \code{poly.stats} data frame.}
#'   \item{\code{poly.obj}}{ A list containing the polygons objects (the outputs of the \code{generalize.f}
#'   function for each polygon). The order of elements in this list corresponds
#'   to the row order of \code{poly.stats} and the polygon order in \code{polygons.xy}.}
#'   \item{\code{call.info}}{Information about hespdiv call: method and
#'   arguments used.}
#'   \item{\code{str.difs}}{ A list containing comparison values of
#'   evaluated straight split-lines from each polygon. The elements of this
#'   list correspond to the rows of \code{poly.stats}.
#'   }
#'   }
#' @note Please note that if you use the \code{method} argument, the arguments
#' \code{generalize.f}, \code{compare.f}, and \code{maximize} are determined
#' internally and should not be provided. Therefore, you should only assign
#' values to these arguments when using a custom method, not a predefined one.
#' Additionally, you can ignore the \code{generalize.f} argument even when
#' applying custom methods. If \code{generalize.f} is set to NULL (default), the
#' data remains unchanged, as \code{generalize.f} acts as an identity function.
#' Hence, \code{generalize.f} is only an optional argument that allows to omit
#' the transformation or generalization step in \code{compare.f} function,
#' simplifying it.
#'
#' Both \code{generalize.f} and \code{compare.f} functions inherit parent
#' function environments (e. g. hespdiv, .spatial_div), allowing to use
#' additional variables (such as xy.dat, samp.xy, id1, id2) as arguments from
#' those environments.
#'
#' There is a small possibility for a nonlinear split-line to cross a polygon
#' boundary. If the result contains such a split-line, or if an error was received
#' related to this issue, the recommendation would be to make a small change
#' in some of the arguments (e.g., 'n.split.pts') and re-run hespdiv.
#'
#' @details
#' \subsection{The Algorithm}{
#' \bold{1) Split-point Placement:} The function places a predetermined number of
#' split-points (\code{n.split.pts + 1}) along the perimeter of the study area
#' (\code{study.pol}) or the convex hull of observation locations
#' (\code{xy.dat}) if a study area polygon is not provided. These split-points
#' are evenly spaced, resulting in a distance between points equal to
#' 1/(n.split.pts + 1) fraction (1/16 by default) of a polygon circumference.
#'
#' \bold{2) Straight Split-lines:} Straight split-lines are generated by connecting
#' the split-points. The total number of straight split-lines generated is equal
#' to the value of \code{sum(1:n.split.pts)}.This holds true only in the first iteration when
#' \code{same.n.split} is set to \code{FALSE} or in all iterations when
#' \code{same.n.split} is set to \code{TRUE}. Note that the total number of
#' split-lines generated will not be equal to the number of split-lines
#' evaluated since split-lines that cross polygon boundary or do not pass
#' sample size, area or location number subdivision criteria are not evaluated.
#'
#' \bold{3) Subdivisions:} Each split-line spatially divides the data and study area
#' polygon into two subsets.
#'
#' \bold{4) Criteria:} Both subsets are then checked to see if they meet
#' sample size, area and location number criteria.
#'
#' \bold{5) Obtaining Comparison Values:} Subsets that meet criteria are
#' compared using \code{generalize.f} and \code{compare.f} functions to obtain
#' a comparison value. First, each subset is passed into the
#' \code{generalize.f} function to obtain a generalization value
#' (e.g., Pielou evenness index in 'pielou' method) for it (please refer to the
#''note' section for clarification). Then these values are used in the
#' \code{compare.f} function to compare the subsets, producing a comparison
#' value. This way each split-line that passed all subdivision criteria is
#' assigned a comparison value.
#'
#' \bold{6) Best Straight Split-line Selection:} The best performing straight
#' split-line is determined based on whether the \code{maximize} argument is
#' set to \code{TRUE} or \code{FALSE}. If \code{maximize} is \code{TRUE}, the
#' best split-line is the one with the highest comparison value; if
#' \code{maximize} is \code{FALSE}, the best split-line has the lowest
#' comparison value.
#'
#' The combination of the \code{generalize.f}, \code{compare.f}, and
#' \code{maximize} arguments can be provided, enabling the creation of custom
#' methods, or it can be determined internally based on the chosen preset
#' subdivision methods (see below).
#'
#' \bold{7) Nonlinear Split-lines:} If the best straight split-line meets
#' the quality criteria specified by the \code{c.Q.crit}, it serves as the basis
#' for generating variously shaped curves (nonlinear split-lines). These curves
#' are produced using splines, which are mathematical functions that can create
#' smooth and flexible curves.
#'
#' To generate the curves, a number of knots (control points) for the splines are
#' distributed evenly along a set of lines orthogonal to the best straight
#' split-line. These orthogonal lines are also evenly distributed along the
#' straight split-line itself. The number of knots and lines used can be
#' adjusted through the parameters \code{c.Y.knots} and \code{c.X.knots},
#' respectively.
#'
#' The algorithm then iterates through this network of knots, considering
#' different combinations of knots, to produce curves. By varying the selection
#' and arrangement of knots, different shapes of curves are generated.
#'
#'
#' \bold{8) Nonlinear Split-line Evaluation and Selection:} Curves are then
#' processed in the same manner as straight split-lines in steps 3 to 6.
#'
#' \bold{9) Final Split-line Selection and Establishment of Subdivision:} If the best
#' curve outperforms the best straight split-line by a margin of
#' \code{c.crit.improv}, the best split-line becomes nonlinear. Otherwise, it
#' remains straight. The split-line must also satisfy the criteria established
#' by the \code{Q.crit} argument to be used for subdivision.
#'
#' \bold{10) Recursive Iteration:} The process described above is iteratively applied,
#' resulting in a collection of selected split-lines with their performance
#' values. These split-lines hierarchically subdivide space and data, forming
#' polygons of various shapes.
#' }
#' \subsection{Preset Methods}{
#' Preset methods and their combinations of \code{compare.f},
#' \code{generalize.f} and \code{maximize} arguments:
#' \describe{
#' \itemize{
#'   \item "sorensen" - functions calculate Sorensen similarity index
#'   (Sorensen 1948), maximize = FALSE. Values vary from 0 to 1.
#'   \item "pielou" - functions calculate the mean proportional decrease in
#'   Pielou evenness (Pielou 1966) that occurs after polygon subdivision,
#'   maximize = TRUE. Values vary from 0 to 1.
#'   \item "morisita" - functions calculate Morista overlap index
#'   (Morisita 1959), maximize = FALSE. Values vary from 0 to 1.
#'   \item "horn.morisita" - functions calculate Morista overlap index,
#'   as modified by Horn (1966), maximize = FALSE. Values vary from 0 to 1.
#'   }
#'   }
#' All the preset methods currently available are specifically designed for
#' bioregionalization purposes.These methods necessitate two key inputs: the
#' coordinates of fossil taxa occurrences (\code{xy.dat}) and the names or IDs
#' of the taxa (\code{data}).These names or IDs should be structured as
#' character or numeric vectors, with each element corresponding to a row in the
#' \code{xy.dat} data frame. Each method compares fossil communities on
#' opposite sides of a split-line, aiming to minimize similarity or
#' maximize difference. Their outcome yields biogeographical provinces with a
#' hierarchical structure.
#'}
#' @author Liudas Daumantas
#' @importFrom grDevices chull
#' @importFrom pracma polyarea
#' @references Horn, H. S. (1966). Measurement of" overlap" in comparative ecological studies. The American Naturalist, 100(914), 419-424.
#' @references Morisita, M. (1959). Measuring of interspecific association and similarity between assemblages. Mem Fac Sci Kyushu Univ Ser E Biol, 3, 65-80.
#' @references Pielou, E. C. (1966). The measurement of diversity in different types of biological collections. Journal of theoretical biology, 13, 131-144.
#' @references Sorensen, T. A. (1948). A method of establishing groups of equal amplitude in plant sociology based on similarity of species content and its application to analyses of the vegetation on Danish commons. Biol. Skar., 5, 1-34.
#' @examples
#' ## DO NOT RUN:
#' ## Preparing fossil assemblage data:
#' ## loading hespdiv data package that contains fossil mammal occurrence data from US
#' # library(HDData)
#' # species <- mio_mams$accepted_name # taxa names
#' # sp_coords <- data.frame(x = mio_mams$lng, y = mio_mams$lat)
#' ## Running HespDiv with default arguments
#' # h <- hespdiv(data = species, xy.dat = sp_coords, study.pol = us)
#' # us is just a polygon of contiguous US for visualization, see str(us)
#' # plot_hespdiv(h, n.loc = TRUE)
#'
#' @export

hespdiv<-function(data,
                  xy.dat = NULL,
                  n.split.pts = 15,
                  same.n.split = TRUE,
                  method = 'horn.morisita',
                  generalize.f = NULL,
                  compare.f = NULL,
                  maximize = NULL,
                  N.crit = 1,
                  N.rel.crit = 0.2,
                  N.loc.crit = 1,
                  N.loc.rel.crit = 0.2,
                  S.crit = 0.05,
                  S.rel.crit = 0.2,
                  Q.crit = NULL,
                  c.splits = TRUE,
                  c.Q.crit = NULL,
                  c.crit.improv = 0,
                  c.X.knots = 5,
                  c.Y.knots = 10,
                  c.max.iter.no = +Inf,
                  c.fast.optim = FALSE,
                  c.corr.term = 0.05,
                  study.pol = NULL,
                  use.chull = TRUE,
                  tracing = NULL,
                  pnts.col = 1,
                  display = FALSE,
                  pacific.region = FALSE){

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
    if (any(colnames(data) == 'x') & any(colnames(data) == 'y')){
      stop("x and y columns are present in data, when xy.dat is provided. ",
           "Remove either one from the input.")
    }
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
  args <- sapply(ls(),get,environment())
  n.split.pts <- n.split.pts + 1
  c.Y.knots <- c.Y.knots + 2
  c.X.knots <- c.X.knots + 2

  if (!is.null(tracing)){
    if (is.vector(tracing)){
      if (length(tracing) == 2){
        trace.object <- tracing[1]
        trace.level <- tracing[2]
      } else { stop("tracing must contain exactly two elements")}
    } else { stop("tracing must be a vector")}
  }

  if (exists("trace.level")){
    if (!is.null(trace.object)){
      trace.object <- .arg_check("trace.object", trace.object,
                                 c("straight", "curve", "both"))
    }
    if (!is.null(trace.level)){
      trace.level <- .arg_check("trace.level", trace.level,c("all","main","best"))
    }
  } else {
    trace.object <- NULL
    trace.level <- NULL
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
        stop("Data must be a vector when using '",method,"' method.")

      generalize.f <- function(plot.dat){
        x <- table(plot.dat)
        p <- x/sum(x)
        -sum(log(p)*p)/log(length(p))
      }
      compare.f <- function(x,y) {
        base.eveness <- poly.obj[[testid]]
        (1 - mean(c(x, y)) / base.eveness) # proportional change in mean eveness
      }

    } else {
      if (method == "sorensen"){
        similarity <- TRUE
        variant <- '1'
        metric <- 'sorensen'
        method <- 'biozonation'
        maximize <- FALSE
        if(!is.vector(data))
          stop("Data must be a vector when using '",method,"' method.")

        generalize.f <- function(plot.dat){
          unique(plot.dat)
        }
        compare.f <- function(x,y) {
          sum <- length(x) + length(y)
          int_2x <- length(which(duplicated(c(x,y))))*2
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
            stop("Data must be a vector when using '",method,"' method.")

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
              stop("Data must be a vector when using '",method,"' method.")

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
    args$maximize <- maximize
    args$generalize.f <- generalize.f
    args$compare.f <- compare.f
  } else {
    if (!is.null(method)){
      stop(paste0("Please check the input. One of the method arguments",
                  " ('compare.f' or 'method') must be left NULL."))
    }
    if (is.null(generalize.f)) {
      generalize.f <- function(plot.dat) plot.dat
      args$generalize.f <- generalize.f
      }
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
                " Slicing of 'data' is attempted using '[]' brackets.")
      .slicer <- .slicer.vect
    }
  }


  if (!maximize){
    if (is.null(Q.crit)) Q.crit <- +Inf
    if (is.null(c.Q.crit)) c.Q.crit <- Q.crit
    .comp <- function(x,criteria){ x < criteria}
    c.sign <- "<"
    .minormax <- min
    .which_minormax <- which.min
  } else {
    if (is.null(Q.crit)) Q.crit <- -Inf
    if (is.null(c.Q.crit)) c.Q.crit <- Q.crit

    .comp <- function(x,criteria){ x > criteria}
    c.sign <- ">"
    .minormax <- max
    .which_minormax <- which.max
  }
  args$c.Q.crit <- c.Q.crit
  args$Q.crit <- Q.crit
  if ((c.splits == FALSE & Q.crit != c.Q.crit) |
      ((Q.crit > c.Q.crit) & maximize) |
      ((Q.crit < c.Q.crit) & !maximize) ){
    if (c.splits == FALSE & Q.crit != c.Q.crit & !c.Q.crit %in% c(+Inf,-Inf)){
      warning(paste("Since 'c.splits' is FALSE, 'c.Q.crit' was set equal to
          'Q.crit'"))
      c.Q.crit <- Q.crit
    }
    if ((Q.crit > c.Q.crit) & !maximize){
      warning(paste("Q.crit should be lower or equal to c.Q.crit,",
                    "when optimization is reached by minimazing the metric."))
    }
    if ((Q.crit < c.Q.crit) & maximize){
      warning(paste("Q.crit should be higher or equal to c.Q.crit,",
                    "when optimization is reached by maximising the metric."))
    }
  }
  if (pacific.region){
    if (!is.null(study.pol)) {
      ind <- study.pol[,1] >= -180 & study.pol[,1] <= 0
      study.pol[ind,1] <- study.pol[ind,1] + 360
    }
    ind <- xy.dat$x >= -180 & xy.dat$x <= 0
    xy.dat$x[ind] <- xy.dat$x[ind] + 360
  }

  if (use.chull) {
    if (!is.null(study.pol)) {
      study.pol.viz <- study.pol
    }
    ch <- chull(xy.dat$x, xy.dat$y)
    study.pol <- data.frame(x = xy.dat$x[c(ch,ch[1])],
                            y = xy.dat$y[c(ch,ch[1])])
  } else {
    if (is.null(study.pol)){
      stop("Please provide study.pol arg. or use 'use.chull = TRUE' ")
    }
    study.pol <- .clean_poly(study.pol)
  }
  if (!exists("study.pol.viz")){
    study.pol.viz <- NULL
  }
  if (!same.n.split){
    dst.pts <- .calc.perim(study.pol) / n.split.pts
  } else {
    dst.pts <- NULL
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

  S.cond <- round(abs(pracma::polyarea(study.pol$x,study.pol$y)) * S.crit,2)
  splits <- numeric()

  e <- environment()
  environment(.spatial_div) <- e
  ### obtaining results:
  .spatial_div(data,xy.dat,root.id=0)

  ### results obtained. Formatting them.
  if (length(splits)>0) {
    plot.id <- 1:nrow(poly.info)
    names(poly.obj) <- plot.id
    names(rims) <- plot.id
    poly.info <- data.frame(plot.id = plot.id,poly.info)

    names(str.split.quals) <- 1:length(str.split.quals)


    # which best is not clear, when no split was above c.Q.crit
    # (c.Q.crit is lower performance boundary, Q.crit - upper;
    # if c.Q.crit is beaten with straight split-line, but not Q.crit,
    # then it is possible that non-linear split-line will beat Q.crit;
    # If no split-line beats c.Q.crit, the str.best value in poly.info will be
    # c.Q.crit, but not the best performance of straight split-line)
    if (any(poly.info$str.best == c.Q.crit, na.rm = TRUE)){

      plot.id.with.ns <- which(poly.info$n.splits != 0)
      poly.info$str.best[-plot.id.with.ns] <- NA
      poly.info$str.best[plot.id.with.ns] <- ifelse(poly.info$str.best[
        plot.id.with.ns] == c.Q.crit,
        unlist(lapply(str.split.quals[plot.id.with.ns],.minormax)),
        poly.info$str.best[plot.id.with.ns])
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
    result$poly.stats <- data.frame(rank = .split.rank(
      result$poly.stats),
      result$poly.stats)
    result$split.stats <- data.frame(
      rank = result$poly.stats[result$split.stats$plot.id,"rank"],
      result$split.stats)

    if (display){
      .visualise_splits.end(pnts.col, xy.dat, rims, study.pol.viz)
    }

    return(print.hespdiv(result))
  } else {
    base::message("No split-lines were established.\nTry loosening the split-line selection criteria.")
  }
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
        mean = e$poly.info$mean[e$plot.id],
        sd = e$poly.info$sd[e$plot.id],
        z.score = e$poly.info[e$plot.id,"str.z.score"],
        performance = e$poly.info[e$plot.id,"str.best"]
      )
    }
  }
  .structurise(e,split.stats,split.lines)
}

.arg_check <- function(name, given,NAMES){
  matched.i <- pmatch(tolower(given), tolower(NAMES))
  if(is.na(matched.i))
    stop("invalid argument '",name, "' value: ", paste0('"', given,'".'),
         paste0("\nPlease select viable option: "),
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
# Calculate rank of the split-line
#' @noRd
.split.rank <- function(poly.stats){

  roots <- poly.stats$root.id
  ranks.id <- numeric(length(roots))
  for (split in seq(length(roots))){
    split_rank <- 1
    split_root <- roots[split]
    while (split_root != 0) {
      split_rank <- split_rank + 1
      split_root <- poly.stats$root.id[
        which(poly.stats$plot.id == split_root)]
    }
    ranks.id[split] <- split_rank
  }
  ranks.id
}
#  get ids of observations within polygon
#' @noRd
.get_ids <- function(polygon, xy_dat) {
  which(.point.in.polygon(pol.x = polygon[,1],
                             pol.y = polygon[,2],
                             point.x = xy_dat$x,
                             point.y = xy_dat$y) != 0 )
}
