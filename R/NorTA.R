#' @name NorTA
#' @rdname NorTA
#' 
#' @title Normal To Anything
#' 
#' @description statistical technique to produce variables with correlation 
#' structure and type of distribution chosen by the user.
#'
#' @param n number of observations.
#' @param cor correlation structure of resulting data; cor must be a number or a
#' correlation matrix. If it is a number it is understood that it is the 
#' correlation between two variables.
#' @param dist characters indicate the chosen distribution.
#' @param param list containing the parameters of the distribution for each 
#' dimension. All elements in the list must be named as the parameter name of
#' the distribution. If the elements have all lengths 1 the values are repeated
#' for all the dimensions.
#' 
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats pnorm
#' 
#' @export
NorTA <- function(n,cor,dist,param=list()){
  
  #check n
  #-----------------------------------#
  if(!is.numeric(n) | n <=0 | round(n)!=n) stop("n must be positve integer number")
  #check dist
  #-----------------------------------#
  ifelse(is.character(dist), 
         qdist<-paste("q",dist,sep=""),
         stop("dist must be a character string"))
  if(!exists(qdist))stop(paste("quantile function",qdist,
                               "has not been found, maybe you need to include the required package?"))
  #check cor
  #-----------------------------------#
  if(length(cor)==0){
    stop("cor has length 0, what happen?")
  } else if(length(cor)==1){
    if(!is.numeric(cor) | abs(cor)>1) stop("if cor is a number must be in range [-1,1]")
  } else if(length(cor)>1){
    if(!is.numeric(cor) | !isSymmetric(cor) | any(diag(cor)!=1) | any(abs(cor)>1))
      stop("cor must be a symmetric matrix with all elements in range [-1,1] and 
           with all elements values equal to 1")
  }
  # Set parameter
  if(length(cor)==1) cor <- matrix(c(1,cor,cor,1),nrow=2)
  D <- nrow(cor)
  #check param
  #-----------------------------------#
  if(length(param)!=0){
    lengths <- lapply(param, length)
    if(!all(lengths==1 | lengths==D)){
      stop("all elemnts in param must have lengts equal to 1 or equal to the dimension of cor")
    } 
    param <- lapply(param, function(x) if(length(x)==1) rep(x,D) else x )
  }
  
  rnorm <- mvtnorm::rmvnorm(n=n,sigma=cor)
  unif <- stats::pnorm(rnorm)

  ans <- sapply(1:nrow(cor), function(idx){
    sub.param <- c(list(p=unif[,idx]),lapply(param, `[[`, idx))
    do.call(what=paste("q",dist,sep=""),args=sub.param)
  })
  
  return(ans)
}