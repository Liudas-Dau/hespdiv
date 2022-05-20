#' @noRd
.generate_cols <- function(n,seed){
  set.seed(seed)
  #sp.n <- length(unique(data[,-which(colnames(data) %in% c('x','y'))]))
  if (n<=74){
    qual_col_pals <- RColorBrewer::brewer.pal.info[
      RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector <- unlist(mapply(
      RColorBrewer::brewer.pal, qual_col_pals$maxcolors,
      rownames(qual_col_pals)))
    if (n <=27) {
    palete <- sample(col_vector[
      -c(3:4,5,15,12,10,8,1,21,23,25,29,30:44,19,13,27,28, 45,22,18,48,24,51,
         58,62,63,64,65,59,70,71,73,74)], n, replace = FALSE)
    } else {
      palete <- sample(col_vector, n, replace = FALSE)
    }
  } else{
    color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    palete <- sample(color,n)
  }
  palete
}
