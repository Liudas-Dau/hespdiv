% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hsa.R
\name{hsa}
\alias{hsa}
\title{HespDiv Sensitivity Analysis}
\usage{
hsa(
  obj,
  n.runs = 100,
  data.paired = TRUE,
  display = FALSE,
  images.path = NULL,
  pnts.col = 1,
  data = NULL,
  xy.dat = NULL,
  same.n.split = NULL,
  n.split.pts = NULL,
  N.crit = NULL,
  N.rel.crit = NULL,
  N.loc.crit = NULL,
  N.loc.rel.crit = NULL,
  S.crit = NULL,
  S.rel.crit = NULL,
  Q.crit = NULL,
  c.splits = NULL,
  c.Q.crit = NULL,
  c.crit.improv = NULL,
  c.X.knots = NULL,
  c.Y.knots = NULL,
  c.max.iter.no = NULL,
  c.fast.optim = NULL,
  c.corr.term = NULL,
  study.pol = NULL,
  use.chull = NULL,
  generalize.f = NULL,
  maximize = NULL,
  method = NULL,
  compare.f = NULL,
  .run.id = NULL
)
}
\arguments{
\item{obj}{A 'hespdiv' class object.}

\item{n.runs}{An integer representing how many alternative hespdiv calls to evaluate.}

\item{data.paired}{A Boolean controlling whether the provided alternative
values of 'data' are paired with 'xy.dat'.}

\item{display}{A Boolean indicating the value of the "display" argument in
each hespdiv call.}

\item{images.path}{A path to an existing directory where PNG images of the
displayed results will be saved. If NULL (default), images won't be saved.}

\item{pnts.col}{the value of the "pnts.col" argument in each hespdiv call.}

\item{data}{A list containing matrices, time-series, lists, data frames,
vectors, or other data structures.}

\item{xy.dat, study.pol}{Lists of data frames with two columns: 'x' and 'y'.}

\item{same.n.split, c.fast.optim, use.chull, c.splits}{A Boolean value (if used, should be different from the
one in the basal hespdiv call).}

\item{n.split.pts, c.max.iter.no, N.crit, N.loc.crit, c.X.knots, c.Y.knots}{Numeric integer vectors.}

\item{N.rel.crit, N.loc.rel.crit, S.crit, S.rel.crit}{Numeric vectors with values
between 0 and 1.}

\item{Q.crit, c.Q.crit, c.crit.improv}{Numeric vectors}

\item{c.corr.term}{A numeric vector with values between 0.01 and 0.2.}

\item{generalize.f, compare.f}{Lists of functions.}

\item{maximize}{A logical vector of the same length as 'compare.f' list.}

\item{method}{A character vector.}

\item{.run.id}{An integer. Runs up to this id will be skipped.}
}
\value{
A 'hsa' class object. It is a list of two items:
\describe{
\item{\bold{Alternatives:}}{ A list containing the produced
alternative hespdiv objects. }
\item{\bold{Basis}}{ The basal hespdiv object whose call was modified to produce
alternative subdivisions.}
}
}
\description{
This function is the main function that performs HespDiv sensitivity analysis.
The function takes a provided hespdiv object as a base and produces a
desired number of alternative versions by randomly sampling new values for
its arguments from a provided selection of values for each argument.
}
\details{
\subsection{Difference Between "hsa" And "hsa_detailed"}{
The major difference between "hsa_detailed" and "hsa" is that the former produces
all possible hespdiv calls from combinations of the provided hespdiv arguments.
Therefore, it samples a much smaller segment of the parameter space but more
densely, requiring much more computation time. Although such behavior may
be desired in some cases, the "hsa" function is generally more suitable for
performing hespdiv sensitivity analysis.

Additionally, alternative values for hespdiv arguments in the "hsa_detailed"
function are provided in lists, whereas in the "hsa" function, they are
provided in vectors or lists (depending on the argument).
}
\subsection{Paired Arguments}{
The argument 'data.paired = TRUE' means that in each produced hespdiv
call, the index that samples lists of 'xy.dat' and 'data' will be the same.
This way, a different number of observations or entirely different data sets
can be used for sensitivity analysis. If 'data.paired = FALSE', then the number
of observations must be kept the same in all provided elements of 'data' and
'xy.dat' lists. This option allows to add noise to the data or coordinates,
or shuffle them entirely in order to test the significance of the detected spatial
structure.

Arguments that define custom methods ('compare.f', 'generalize.f' and
'maximize') are paired by default. Therefore, their lengths must also be
the same.
}
}
\note{
If a particular call produces a warning or error, a list of length 2
will be returned for that call. If a warning was produced, the first element
of the list will hold the created hespdiv object, and the second element
will contain the warning message. In case of an error, the first element will
be a list of arguments used to produce the call, and the second element will
contain the error message.
}
\seealso{
Other {functions for hespdiv sensitivity analysis}: 
\code{\link{change_base}()},
\code{\link{hsa_detailed}()},
\code{\link{hsa_quant}()},
\code{\link{plot_hsa_q}()},
\code{\link{plot_hsa}()}

Other {functions for hespdiv post-prossesing}: 
\code{\link{cross_comp}()},
\code{\link{hsa_detailed}()}
}
\concept{{functions for hespdiv post-prossesing}}
\concept{{functions for hespdiv sensitivity analysis}}
