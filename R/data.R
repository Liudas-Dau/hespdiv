#' Example hespdiv object
#'
#' A small \code{hespdiv} object created from simulated occurrence data with a
#' known boundary at \code{x = 0.45}. The object is intended for examples,
#' tests, and demonstrations of plotting and post-processing functions.
#'
#' @format A \code{hespdiv} object.
#'
#' @details
#' The object was generated using simulated coordinates and taxon names.
#' Occurrences on the left and right sides of the artificial boundary share
#' several common taxa but also include side-specific taxa.
#' The exact code used to simulate the dataset and run \code{hespdiv} is also
#' demonstrated in the examples of \code{hespdiv} function.
#'
#' @examples
#' data(example_hespdiv)
#' plot_hespdiv(example_hespdiv)
#'
#' @name example_hespdiv
#' @docType data
#' @keywords datasets
NULL
