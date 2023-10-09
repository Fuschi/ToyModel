#' Centered Log-Ratio Transformation
#' 
#' @description It calculates the centered log-ratio transformation of X 
#' considering the rows as samples. If zero are present in X the function uses
#' the simplest replacement strategy replacing them with a constant value equal 
#' to the 65% of the detection limit of the sample (https://doi.org/10.1023/A:1023866030544).
#'
#' @param X numeric matrix or vector with all elements greater than or equal to 0.
#' 
#' @export
clr <- function(X){
  
  # CHECKS
  stopifnot(exprs={
    is.matrix(X) | is.vector(X)
    is.numeric(X)
    all(X>=0)
    
    })
  
  # ZERO REPLACEMENTS
  if(all(round(X)-X==0)){
    dl <- matrix(.65,nrow=nrow(X),ncol=ncol(X))
  } else {
    dl <- apply(X, 1, function(x)min(x[x>0]))
    dl <- replicate(ncol(X),dl)
    dl <- .65*dl
  }
  
  Y <- X + dl*(X==0)
  
  # CENTERED LOG-RATIO
  ref <- apply(X, 1, function(x) mean(log(x)) )
  return(as.matrix(log(X) - ref))

}