#' IN PROGRESS...
#' 
#' @param X ...
#' 
#' @importFrom vegan diversity
#' @export
average_pielou <- function(X){
  
  stopifnot(is.matrix(X))
  stopifnot(is.numeric(X))
  stopifnot(all(X>=0))
  
  mean(vegan::diversity(X) / log(ncol(X)))
  
}