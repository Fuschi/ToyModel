#' @name toy_model
#' @rdname toy_model
#' 
#' @title Meta-genomic Toy Model
#' 
#' @description numerical simulations for meta-genomic data with user-chosen 
#' correlation structure and sample heterogeneity.
#'
#' @param n number of observations.
#' @param cor correlation structure of final data, it can be a matrix or a number.
#' If it's a matrix it must be symmetric, with diagonal elements all equal to 1 and
#' with all values in range \[-1,1\]. The size of cor matrix implies the dimensions 
#' of the resulting simulated data. Instead, if cor is a number, it 
#' indicates the correlation between only two variables.
#' @param M magnification factor, real positive number that modify the 
#' heterogeneity of the samples. In practice the first variable is multiplied by 
#' the factor M.
#' @param qdist quantile function of the selected distribution, it must the vector
#' of probability as first argument called p. The other parameters names must match with
#' the ones present in param list.
#' @param param param named array or matrix/data.frame with parameters of the selected target 
#' distribution indicated in qdist. If params is a matrix then it must contain 
#' the parameters of the final distribution along the columns and rows express 
#' each dimension equal to the dimension of the cor matrix, instead if it is a 
#' vector the parameters are repeated for each dimension. The names of the columns,
#' or of the vector elements, must be equal to the function parameters of qdist.
#' If some parameters are missing then the default values of the qdist function 
#' are chosen. Default value is NULL that indicates to use always the default 
#' parameters of qdist.
#' @param method type of correlation used. Possible choices are "pearson", 
#' "kendall", "spearman" (default pearson). 
#' @param force.positive logical, indicates when to force all generated NorTA 
#' data to positive numbers, adding the minimum to all others (this passage does 
#' not affect correlations).
#'
#' @returns returns an object of class "toy_model" containing:
#' \itemize{
#'  \item cor input correlation matrix.
#'  \item normal generated Gaussian data.
#'  \item cor_normal correlation matrix of normal.
#'  \item NorTA generated data of choosen distribution.
#'  \item cor_NorTA correlation matrix of NorTA.
#'  \item L1 relative abundances of NorTA.
#'  \item cor_L1 correlation matrix of L1.
#'  \item CLR matrix of abundances of clr transformed data from NorTA.
#'  \item cor_CLR correlation matrix of CLR.
#' }
#' 
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats pnorm cor
#' 
#' @export
toy_model <- function(n,cor,M,qdist,param=NULL,method="pearson",
                      force.positive=FALSE){
  
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
  #check M
  #-----------------------------------#
  if(!is.numeric(M) | M<=0) stop("M must be positve number")
  #check method
  #-----------------------------------#
  method  <- match.arg(method,c("pearson", "kendall", "spearman"))
  
  
  if(missing(cor)){
    cor <- matrix(0,nrow=D,ncol=D)
    cor[ceiling(.21*D),ceiling(.61*D)] <- .8
    cor[ceiling(.41*D),ceiling(.81*D)] <- -.8
    cor <- cor+t(cor)
    diag(cor) <- 1
  } else {
    D <- nrow(cor)
  }
  
  rnorm <- mvtnorm::rmvnorm(n=n,sigma=cor)
  unif <- stats::pnorm(rnorm)
  
  NorTA <- matrix(0,nrow=n, ncol=D); colnames(NorTA) <- colnames(cor)
  for(d in 1:D){
    
    ifelse(is.null(param),
           
           param.d <- list(p=unif[,d]),
           
           param.d <- c(list(p=unif[,d]),
                        lapply(split(param[d,], names(param[d,])), unname)))
    
    NorTA[,d] <- do.call(what=qdist,args=param.d)
  }
  
  # if(force.positive && any(NorTA<=0)){
  #   NorTA <- apply(NorTA, 2, function(x){
  #     xpos <- x - min(x)
  #     xpos[xpos==0] <- .65*(min(xpos[xpos>0]))
  #     return(xpos)
  #   })}
  if(force.positive && any(NorTA<=0)){
    NorTA <- NorTA - min(NorTA)
    NorTA[NorTA==0] <- min(NorTA[NorTA>0])
  }
  
  NorTA[,1] <- NorTA[,1]*M
  
  if(any(NorTA<0)) stop("the transformed data NorTA cannot have negative values.")
  
  cor_normal = stats::cor(rnorm,method=method)
  cor_NorTA = stats::cor(NorTA,method=method)
  
  L1 <- NorTA/rowSums(NorTA)
  cor_L1 <- stats::cor(L1,method=method)
  
  .clr <- function(X){
    if(any(X==0)){
      min.x <- min(X[X>0])
      X[X==0] <- .65*min.x
    }
    ref <- apply(X, 1, function(x) mean(log(x)) )
    return(as.matrix(log(X) - ref))
  }
  
  CLR <- .clr(NorTA)
  cor_CLR <- stats::cor(CLR,method=method)
  
  results <- list()
  results$cor <- cor
  results$normal <- rnorm
  results$cor_normal <- cor_normal 
  results$NorTA <- NorTA
  results$cor_NorTA <- cor_NorTA
  results$L1 <- L1
  results$cor_L1 <- cor_L1
  results$CLR <- CLR
  results$cor_CLR <- cor_CLR
  class(results) <- "toy_model"
  
  return(results)
}