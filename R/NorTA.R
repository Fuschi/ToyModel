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
#' @param qdist quantile function of the selected distribution, it must the vector
#' of probability as first argument called p. The other parameters names must match with
#' the ones present in param list.
#' @param param named array or matrix/data.frame with parameters of the selected target 
#' distribution indicated in qdist. If params is a matrix then it must contain 
#' the parameters of the final distribution along the columns and rows express 
#' each dimension equal to the dimension of the cor matrix, instead if it is a 
#' vector the parameters are repeated for each dimension. The names of the columns,
#' or of the vector elements, must be equal to the function parameters of qdist.
#' If some parameters are missing then the default values of the qdist function 
#' are chosen. Default value is NULL that indicates to use always the default 
#' parameters of qdist.
#' 
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats pnorm
#' 
#' @export
NorTA <- function(n,cor,qdist,param=NULL){
  
  # START CHECKS
  #----------------------------------------------------------------------------#
  
  # n
  #-----------------------------------#
  if(!is.numeric(n) | n <=0 | round(n)!=n) stop("n must be positve integer number")
  # qdist
  #-----------------------------------#
  if(!is.function(qdist)) stop("qdist must be a function")
  if( !("p"%in%names(formals(qdist))) ) stop("p must be an argument of the quantile function qdist, it's where the probabilities are assigned")
  # cor
  #-----------------------------------#
  if(length(cor)==0){
    stop("cor has length 0, what happen?")
  } else if(length(cor)==1){
    if(!is.numeric(cor) | abs(cor)>1) stop("if cor is a number must be in range [-1,1]")
  } else if(length(cor)>1){
    if(!is.numeric(cor) | !isSymmetric(cor) | any(diag(cor)!=1) | any(abs(cor)>1))
      stop("cor must be a symmetric matrix with all elements in range [-1,1] and all
           elements on the diagonal equal to 1.")
  }
  
  # Set input parameters
  #-----------------------------------#
  if(length(cor)==1) cor <- matrix(c(1,cor,cor,1),nrow=2)
  D <- nrow(cor)
  
  # Check param 
  #-----------------------------------#
  if(!is.null(param)){
    
    if(length(dim(param))>2) stop("param must be a named vector or a matrix/data.frame")
    if(is.null(dim(param))) param <- as.matrix(t(param))
    if(is.null(colnames(param))) stop("The parameter names of the target distribution could not be found in param")
    if( !all(colnames(param) %in% names(formals(qdist))) ) stop("Find at least a parameter in param not present in the arguments of qdist")
    
    if(nrow(param)==1) param <- do.call(rbind, replicate(D, param, simplify=FALSE)) 
    param <- as.matrix(param)
  }
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  
  
  rnorm <- mvtnorm::rmvnorm(n=n,sigma=cor)
  unif <- stats::pnorm(rnorm)

  ans <- matrix(0,nrow=n, ncol=D); colnames(ans) <- colnames(cor)
  for(d in 1:D){
    
    ifelse(is.null(param),
           
           param.d <- list(p=unif[,d]),
           
           param.d <- c(list(p=unif[,d]),
                        lapply(split(param[d,], names(param[d,])), unname)))
    
    ans[,d] <- do.call(what=qdist,args=param.d)
  }

  return(ans)
}