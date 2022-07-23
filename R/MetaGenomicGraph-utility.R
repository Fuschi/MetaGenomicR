#' Edge Density Threshold
#'
#'@description
#'Get the adjacency matrix considering only the highest absolute value weights
#'up to the desired link density.
#'
#'@param weights correlation matrix
#'@param edge.density threshold value
#'
#'@export
density_link_threshold <- function(weights, edge.density=.05){

  #Checks
  if(nrow(weights)!=ncol(weights)) stop("matrix must be square")
  if(any(rownames(weights)!=colnames(weights))) stop("Weights rownames and colnames must be the same")
  if(edge.density<0 || edge.density>1) stop("edge.density must be in range [0,1]")

  max.edges = ncol(weights) * (ncol(weights)-1) *.5
  weights.list <- unroll(weights,diag=F,symmetric=T)

  weights.list <- cbind(weights.list, 'abs.value'=abs(weights.list$value))
  weights.list <- weights.list[order(weights.list$abs.value, decreasing=TRUE),]
  weights.list <- cbind(weights.list, 'density'= 1:nrow(weights.list) / max.edges )

  if(max(weights.list$density) < edge.density) stop("the edge.density is higher than the links present in the weight matrix ")

  weight.fixed.id <- which(min(abs(edge.density-weights.list$density))==abs(edge.density-weights.list$density))
  weights.list <- weights.list[1:max(weight.fixed.id),]

  adj.weights <- matrix(0, nrow=ncol(weights), ncol=ncol(weights), dimnames=list(colnames(weights),colnames(weights)) )
  for(k in 1:nrow(weights.list)){
    adj.weights[weights.list$row[k],weights.list$col[k]] <- weights.list$value[k]
  }

  return(adj.weights + t(adj.weights))
}

#------------------------------------------------------------------------------#
#' Signed-Weighted Graph Layout
#'
#' @description It elaborates the coordinates for the representation of
#' the vertices of the graph considering only the links with a positive sign.
#'
#'@importFrom igraph layout.fruchterman.reingold subgraph.edges E
#'@export
coord.positive.link <- function(graph, seed=123){

  graph.sub <- subgraph.edges(graph=graph,
                              eids=which(E(graph)$weight>0),
                              delete.vertices=FALSE)

  set.seed(seed)
  layout <- layout.fruchterman.reingold(graph.sub)
  return(layout)
}

#------------------------------------------------------------------------------#
#' Get Communities of a Signed Weighted Graph
#'
#'@description Adaptation in R of the algorithm for detention of weighted communities with sign of a graph developed by Sergio Gomez.
#' https://deim.urv.cat/~sergio.gomez/radatools.php
#'
#'@importFrom igraph graph_from_adjacency_matrix write_graph make_clusters is.weighted is.directed as_adjacency_matrix make_clusters
#'@importFrom stringr str_split
#'
#'@export
community_detection <- function(graph){

  #Check Graph
  if(!is.weighted(graph)) stop("graph must be weighted graph")
  if(is.directed(graph))  stop("graph must be undirected")

  adj <- as_adjacency_matrix(graph, attr="weight", sparse=FALSE)
  #End Checks

  #get path of executable Communities_Detection.exe
  path <- system.file("exec", package="MetaGenomicR", mustWork=TRUE)
  path <- paste(path,"/",sep="")

  #write graph in pajek format as request from executable
  write_graph(graph  = graph_from_adjacency_matrix(adjmatrix=adj, mode="undirected", weighted=TRUE),
              file   = paste(path,"graph.net",sep=""),
              format ="pajek")

  #path for graph and results files
  path.graph <- paste(path,"graph.net",sep="")
  path.result <- paste(path, "res.txt", sep="")

  #save command for the execution
  cmd <- paste(path,"Communities_Detection.exe none WS l 1",sep="")
  cmd <- paste(cmd, path.graph, path.result)

  #make communities detection
  system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)

  #open file of results
  res.file = file(path.result, "r")

  #read and store results
  file.lines <- readLines(res.file)
  info <- file.lines[2]
  modularity <- as.numeric(str_split(file.lines[3]," ")[[1]][3])

  vertex.num <- as.numeric(str_split(file.lines[5]," ")[[1]][4])
  comm.num <- as.numeric(str_split(file.lines[6]," ")[[1]][4])
  comm.vert <- vector(mode="list",length=comm.num)

  for(line in 8:(7+comm.num)){

    i <- line-7
    #comm[[i]]
    tmp <- unlist(str_split(file.lines[line]," "))
    tmp <- as.numeric(tmp[2:length(tmp)])

    comm.vert[[i]] <- tmp
  }

  comm <- vector(mode="integer", length=vertex.num)
  for(c in 1:comm.num){comm[comm.vert[[c]]] <- c}

  #close connection
  close(res.file)

  #remove files
  file.remove(path.graph)
  file.remove(path.result)


  #return the results as communities structure of igraph package
  return(make_clusters(graph=graph,
                       membership=comm,
                       algorithm="signed weights (radatools)",
                       modularity=modularity))
}



# Get Weights
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' Spiec-Easi MB weights extraction
#'@description Extract the weight matrix of the Meinshausen and Bühlmann criterion
#'in the Spiec-Easi algorithm.
#'@importFrom SpiecEasi symBeta getOptBeta
#'
#'@param res.mb Spiec-Easi algorithm results with mb criteria
#'
#'@export
get_mb_weights <- function(res.mb){

  #CHECK ARGUMENTS
  if(class(res.mb)!="pulsar.refit") stop("class must be pulsar.refit")
  if(res.mb$est$method!="mb") stop("must be mb method")
  #END CHECK

  # get weights matrix W
  W <- as.matrix(symBeta(as.matrix(getOptBeta(res.mb)), mode='maxabs'))
  colnames(res.mb$est$data) -> colnames(W) -> rownames(W)
  W[W>0] <- 1
  W[W<0] <- -1

  return(W)
}

#' Spiec-Easi GLASSO weights extraction
#'@description Extract the weight matrix of the GLASSO criterion
#'in the Spiec-Easi algorithm.
#'@importFrom SpiecEasi getOptCov getRefit
#'@importFrom stats cov2cor
#'
#'@export
get_gl_weights <- function(res.gl, Refit=TRUE){

  #CHECK ARGUMENTS
  if(class(res.gl)!="pulsar.refit") stop("class must be pulsar.refit")
  if(res.gl$est$method!="glasso") stop("must be gl method")
  #END CHECK

  # get weights matrix W
  W <- cov2cor(as.matrix(getOptCov(res.gl)))
  if(Refit){W <- W * as.matrix(getRefit(res.gl))}

  colnames(res.gl$est$data) -> colnames(W) -> rownames(W)

  return(W)
}



#'@title communities_color
#'@author Alessandro Fuschi
#'@importFrom qualpalr qualpal
#'@importFrom igraph membership
#'@export
communities_color <- function(comm, n=30){

  assert((round(n)==n && n>0),"n must be an integer positive number")
  assert(class(comm)=="communities","comm must be a communities class object")
  assert(max(membership(comm))<=n,"There are more communities than generated color (increase n)")

  set.seed(1)
  colormap <- qualpal(n=(n+1),colorspace=list(h=c(0,360),s=c(0,1),l=c(.5,1)))
  colormap <- rownames(colormap$RGB)
  names(colormap) <- 0:n

  if("0" %in% names(colormap)){colormap["0"] <- rgb(1,1,1,.8)}
  return(colormap[as.character(membership(comm))])
}






#'@title remove.small.communities
#'@author Alessandro Fuschi
#'@importFrom igraph vcount induced.subgraph V
#'@export
setGeneric("remove.small.communities", function(obj,vertex.num.threshold,keep.discarded) standardGeneric("remove.small.communities"))
setMethod("remove.small.communities", c("MetaGenomicGraph","numeric","logical"), function(obj, vertex.num.threshold,keep.discarded){

  #Checks Arguments
  assert(vertex.num.threshold>0,"vertex.num.threshold must be integer greater than 0")
  #End Checks

  graph <- obj@netw
  comm <- obj@comm

  keep.comm.names <- as.numeric(names(sizes(comm)[sizes(comm)>=vertex.num.threshold]))
  keep.comm.vids <-  which(comm$membership %in% keep.comm.names)

  if(keep.discarded){
    graph.sub <- graph
    comm.sub <- comm
    comm.sub$membership[setdiff(1:comm$vcount, keep.comm.vids)] <- 0
    comm.sub$modularity <- NA
  } else {
    graph.sub <- induced.subgraph(graph, V(graph)[keep.comm.vids])
    comm.sub <- comm
    comm.sub$membership <- comm$membership[keep.comm.vids]
    comm.sub$vcount <- length(comm.sub$membership)
    comm.sub$modularity <- NA
  }


  return(new("MetaGenomicGraph",
             data=obj@data, meta=obj@meta, taxa=obj@taxa,
             netw=graph.sub,comm=comm.sub))
})













