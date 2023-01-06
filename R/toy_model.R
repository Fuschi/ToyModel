#' @name toy_model
#' @rdname toy_model
#' 
#' @title Meta-genomic Toy Model
#' 
#' @description numerical simulations for meta-genomic data with user-chosen 
#' correlation structure and sample heterogeneity.
#'
#' @param n number of observations.
#' @param cor correlation structure of resulting data; cor must be a number or a
#' correlation matrix. If it is a number it is understood that it is the 
#' correlation between two variables.
#' @param D alternative choice to cor and indicates dimensionality. Produces a 
#' correlation matrix of dimension D with by default two pairs of strongly 
#' correlated variables, one positively and one negatively with correlations 
#' values equal to +0.8 and -0.8. D cannot take values less than 5.
#' @param M magnification factor, real positive number that modify the 
#' heterogeneity of the samples. In practice the function multiplies by a 
#' variable by the value of M.
#' @param dist characters indicate the chosen distribution.
#' @param param list containing the parameters of the distribution for each 
#' dimension. All elements in the list must be named as the parameter name of
#' the distribution. If the elements have all lengths 1 the values are repeated
#' for all the dimensions.
#' @param method type of correlation used. Possible choices are "pearson", 
#' "kendall", "spearman" (default pearson). 
#' @param seed random seed for reproducibility (default 42).
#' @param force.positive logical indicates when to force all elements of Y
#' to be positive number adding the minimum to all others (this passage does not
#' affect correlations).
#'
#' @returns \code{toy_model} returns and object of class "toy_model" containing:
#' \itemize{
#'  \item cor input correlation matrix.
#'  \item normal generated gaussian data.
#'  \item cor_normal correlation matrix of normal.
#'  \item NorTA generated data of choosen distribution.
#'  \item cor_NorTA correlation matrix of Y.
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
toy_model <- function(n,cor,D,M,dist,param=list(),method="pearson",seed=42,
                      force.positive=FALSE){
  
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
  #check cor and D
  #-----------------------------------#
  if(missing(cor) && missing(D)) stop("cor and D cannot be unassigned togheter")
  if(!missing(cor) && !missing(D)) stop("cor and D cannot be assigned togheter")
  #check cor
  #-----------------------------------#
  if(!missing(cor)){
    if(!is.numeric(cor) | !isSymmetric(cor) | any(diag(cor)!=1) | any(abs(cor)>1))
      stop("cor must be a symmetric matrix with all elements in range [-1,1] and 
           with all elements values equal to 1")
  }
  #check D
  #-----------------------------------#
  if(!missing(D)){
    if(!is.numeric(D) | D<5 | round(D)!=D) stop("D must be positve integer number greater or equal to 5")
  }
  #check param
  #-----------------------------------#
  if(length(param)!=0){
    lengths <- lapply(param, length)
    if(!all(lengths==1 | lengths==D)){
      stop("all elemnts in param must have lengts equal to 1 or equal to the dimension of cor")
    } 
    param <- lapply(param, function(x) if(length(x)==1) rep(x,D) else x )
  }
  #check M
  #-----------------------------------#
  if(!is.numeric(M) | M<=0 | round(M)!=M) stop("M must be positve number")
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
  
  X <- mvtnorm::rmvnorm(n=n,sigma=cor)
  unif <- stats::pnorm(X)
  
  Y <- sapply(1:nrow(cor), function(idx){
    sub.param <- c(list(p=unif[,idx]),lapply(param, `[[`, idx))
    do.call(what=paste("q",dist,sep=""),args=sub.param)
  })
  
  if(force.positive && any(Y<0)) Y <- Y - min(Y)
  Y[,1] <- Y[,1]*M
  
  if(any(Y<0)) stop("the transformed data Y cannot have negative values.")
  
  corXhat = stats::cor(X,method=method)
  corYhat = stats::cor(Y,method=method)
  
  L1 <- Y/rowSums(Y)
  corL1 <- stats::cor(L1,method=method)
  
  clr <- function(X){
    if(any(X==0)) X <- X+1
    ref <- apply(X, 1, function(x) mean(log(x)) )
    return(as.matrix(log(X) - ref))
  }
  
  CLR <- clr(Y)
  corCLR <- stats::cor(CLR,method=method)
  
  results <- list()
  results$cor <- cor
  results$normal <- X
  results$cor_normal <- corXhat 
  results$NorTA <- Y
  results$cor_NorTA <- corYhat
  results$L1 <- L1
  results$cor_L1 <- corL1
  results$CLR <- CLR
  results$cor_CLR <- corCLR
  class(results) <- "toy_model"
  
  return(results)
}


#' @name graph_toy_model
#' @rdname graph_toy_model
#' 
#' @title Construct toy model network
#' 
#' @description construct an igraph object from toy_model object
#' 
#' @param obj toy_model class object.
#' @param what normalization selected.
#' 
#' @importFrom igraph graph_from_adjacency_matrix V<- E<-
#'
#' @export
graph_toy_model <- function(obj,what){
  
  if(class(obj)!="toy_model") stop("obj must belong to toy_model class")
  what=match.arg(what,c("NorTA","L1","CLR"))
  
  x <- obj[[what]]
  cor <- obj[[paste("cor_",what,sep="")]]
  
  
  g <- graph_from_adjacency_matrix(adj=cor,mode="undirected",weighted=T,diag=F)
  
  V(g)$size <- 2 + (( colMeans(x) / max(colMeans(x))) * 4)

  E(g)$width <- abs(E(g)$weight)
  E(g)$color <- ifelse(E(g)$weight>=0,rgb(0,0,1),rgb(1,0,0))
  
  return(g)
}


#' @name plot_toy_model
#' @rdname plot_toy_model
#' 
#' @title Custom plot for toy_model
#' 
#' @description plot ring network with custom graphical properties
#' 
#' @param obj toy_model class object.
#' @param what normalization selected.
#' 
#' @importFrom igraph plot.igraph layout_in_circle
#'
#' @export
plot.toy_model <- function(obj,what,...){
 
  if(class(obj)!="toy_model") stop("obj must belong to toy_model class")
  what=match.arg(what,c("NorTA","L1","CLR"))
  
  g <- graph_toy_model(obj,what)
  
  plot.igraph(g,layout=layout_in_circle(g),vertex.label=NA,...)
}