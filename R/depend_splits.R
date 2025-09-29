#' Identify dependent split-lines and polygons
#'
#' @description  For a provided split-line id, function returns a list of
#' all offspring polygons and split-lines, if there are any.
#'
#' @return a list:
#' \enumerate{
#'   \item Dependent split-line IDs, if any
#'   \item Dependent polygon IDs
#' }
#' @param obj hespdiv object
#' @param id id of a split-line
#' @export
depend_splits <- function(obj, id){
  offspring <- .dep_split(obj, id) # obtain dependent split ids and polygon ids

  # add offspring to a total vector of offspring:
  all_depend_splits <- offspring$splits
  all_depend_pols <- offspring$polygons

  # create a vector of split-line still needed to check for offspring
  unchecked <- offspring$splits

  # while there are unchecked offspring split-lines, obtain offspring of each
  # offspring
  while (length(unchecked) > 0){
    new <- numeric()
    for (split.id in unchecked){
      offspring <- .dep_split(obj, split.id)
      new <- c(new, offspring$splits)
      all_depend_splits <- c(all_depend_splits, offspring$splits)
      all_depend_pols <- c(all_depend_pols, offspring$polygons)

    }
    unchecked <- new # reset to new list of offspring split-lines
  }
  return(list(
    splits = all_depend_splits,
    polygons = all_depend_pols
  ))
}

# Identify dependent splits and polygons if any
#' @noRd
.dep_split <- function(obj, id){
  # Step 1: identify produced polygons
  offspring_id <-
    which(obj$poly.stats$root.id == obj$split.stats$plot.id[id])
  # Step 2: does each polygon has a split-line?
  depend_l <- obj$poly.stats[offspring_id, "has.split"]
  # if any TRUE, return split ids + polygon ids; else return only polygon ids
  if (any(depend_l)) {
    list(splits = which(obj$split.stats$plot.id %in% offspring_id[depend_l]),
         polygons = offspring_id)
  } else {
    list(splits = numeric(),
         polygons = offspring_id)
  }
}
