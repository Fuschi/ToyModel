#' IN PROGRESS...
#' 
#' @param x obiect belong to toy_model S3 class
#' @param what ...
#' @param layout coordinates of vertices of the network, layout_in_circle is set by default.
#' @param vertex.size ...
#' @param edge.width ...
#' @param ... graphical paramaters to igraph plot
#' 
#' @importFrom igraph graph_from_adjacency_matrix V<- E<- V E layout_in_circle
#' @importFrom methods is
#' @importFrom grDevices rgb
#' 
#' @export
plot.toy_model <- function(x, what, 
                           layout=NULL,
                           vertex.size=NULL, edge.width=NULL,
                           ...){
  
  D <- ncol(x$cor)
  
  what <- match.arg(what, c("Normal","NorTA","L1","CLR"))
  stopifnot(is(x,"toy_model"))
  if(!is.null(vertex.size)){
    stopifnot(exprs={
      is.numeric(vertex.size)
      length(vertex.size)==1 | length(vertex.size)==D
    })}
  if(!is.null(edge.width)){
    stopifnot(exprs={
      is.numeric(edge.width)
      length(edge.width)==1 | length(edge.width)==D
    })}
  if(!is.null(layout)){
    stopifnot(exprs={
      is.numeric(layout)
      nrow(layout)==D 
      ncol(layout)==2
    })}
  
  if(what=="Normal"){
    adj <- x$cor_normal
  } else if(what=="NorTA"){
    adj <- x$cor_NorTA
  } else if(what=="L1"){
    adj <- x$cor_L1
  } else if(what=="CLR"){
    adj <- x$cor_CLR
  }
  
  g <- graph_from_adjacency_matrix(adj, 
                                   mode=c("undirected"),
                                   weighted=T, diag=F)
  
  E(g)$color <- ifelse(E(g)$weight>0, 
                       rgb(0,0,1,abs(E(g)$weight)^.5), 
                       rgb(1,0,0,abs(E(g)$weight)^.5))
  if(is.null(edge.width)){
    edge.width <- ifelse(abs(E(g)$weight) >= .05, abs(E(g)$weight)^.5, abs(E(g)$weight))
  }
  if(is.null(vertex.size)){
    vertex.size <- 15*((colMeans(x$NorTA)/ max(colMeans(x$NorTA)))) + 6
  }
  if(is.null(layout)) layout <- layout_in_circle(g)
  
  V(g)$size <- vertex.size
  E(g)$width <- edge.width
  
  plot(g, layout=layout, ...)
}