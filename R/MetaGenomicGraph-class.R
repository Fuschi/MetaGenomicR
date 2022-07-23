setOldClass("igraph")
setOldClass("communities")
#' S4 class to manage metagenomic data
#'
#'
#'
#' @slot data table of data
#' @slot meta sample features
#' @slot taxa taxonomic information
#' @slot netw Metagenomic igraph network
#' @slot comm graph communities
#'
#' @name MetaGenomicGraph-class
#' @rdname MetaGenomicGraph-class
#' @exportClass MetaGenomicGraph
MetaGenomicGraph <- setClass(
  Class="MetaGenomicGraph",
  contains="MetaGenomic",
  slot=c(
    netw="igraph",
    comm="communities")
)

#' MeteGenomicGraph Constructor
#' @importFrom igraph graph_from_adjacency_matrix
#'
#' @param MetaGenomic Metagenomic Class Object.
#' @param Corr Adjacency Matrix.
#' @export
MetaGenomicGraph <- function(MetaGenomic, Adjacency){

  assert(!is.null(Adjacency) || !is.null(Adjacency), "Adjacency must contain rows and cols names ")
  assert(all(rownames(Adjacency)==colnames(Adjacency)),"Adjacency rows and cols names must be all equals")
  assert(all(rownames(Adjacency)==taxaID(MetaGenomic)),"Adjacency rows and cols names must be equals to taxaID")

  netw <- graph_from_adjacency_matrix(adjmatrix=Adjacency,
                                      mode="undirected",
                                      weighted=TRUE, diag=FALSE)

  comm <- community_detection(netw)
  isolated.membership <- which(sizes(comm)==1)
  comm$membership[comm$membership %in% isolated.membership] <- 0

  return(new("MetaGenomicGraph",
             data=MetaGenomic@data,
             meta=MetaGenomic@meta,
             taxa=MetaGenomic@taxa,
             netw=netw,
             comm=comm))
}


#'@importFrom igraph vcount
setValidity(Class="MetaGenomicGraph",method=function(object){

  assert(ntaxa(object)==vcount(object@netw), "graph vertex number must equal to taxa number of MetaGenomic")

})


#'@importFrom igraph edge_density ecount sizes
setMethod(f="show",
          signature="MetaGenomicGraph",
          definition=function(object){

            print("MetaGenomicGraph Class Object")
            print(paste("Sample Number:",nsample(object)))
            print(paste("Taxa Number:",ntaxa(object)))
            print(paste("Sample Meta Data:",paste(sample_info(object),collapse="," )))
            print(paste("Taxonomic Ranks:",paste(ranks(object),collapse=",")))
            print(paste("Link Number:",ecount(object@netw)))
            print(paste("Edge Density:",round(edge_density(object@netw),2)))
            print(paste("Signed Communities Number:",max(membership(object@comm))))

            if("0" %in% names(sizes(object@comm))){
              print(paste("Communities Sizes:",paste(sizes(object@comm)[-1],collapse=",")))
              print(paste("Isolated Nodes:", sizes(object@comm)[[1]]))
            } else {
              print(paste("Communities Sizes:",paste(sizes(object@comm)[-1],collapse=",")))
              print("There aren't isolated nodes")
            }
          })

#'@export
setMethod("plot", "MetaGenomicGraph", function(x,...){
  plot(x@netw, layout=coord.positive.link(x@netw),...)
})

#' Method extensions to subset MetaGenomicGraph objects.
#'
#'@importFrom igraph V V<- E E<- subgraph as_adjacency_matrix
#' @export
setMethod("[", "MetaGenomicGraph",
          function(x,i,j,...){

            adj <- as_adjacency_matrix(x@netw, attr="weight", sparse=F)
            adj.sub <- adj[j,j]

            mg.sub <- MetaGenomic(x@data[i,j],x@meta[i,],x@taxa[j,])

            return(MetaGenomicGraph(mg.sub,adj.sub)
            )})


#------------------------------------------------------------------------------#
#' Add graphical decorations
#'
#'@description Add graphical descriptions to network:
#'\itemize{
#'  \item Vertex size proportional to mean clr abundances over samples.
#'  \item Link color could be red (+) or blue (-) respect the weights sign.
#'}
#'
#'@importFrom igraph V V<- E E<- graph_from_adjacency_matrix write_graph make_clusters is.weighted is.directed as_adjacency_matrix make_clusters
#'@importFrom stringr str_split str_replace
#'
#'@export
setGeneric("default_decoration", function(obj) standardGeneric("default_decoration"))
setMethod("default_decoration", "MetaGenomicGraph", function(obj){

  # Vertex size
  V(obj@netw)$size  <- 4 + colMeans(clr(obj)) + abs(min(colMeans(clr(obj))))

  # Edges color and width
  w <- E(obj@netw)$weight
  E(obj@netw)$color <- ifelse(w>0, rgb(0,0,1,.5), rgb(1,0,0,.5))
  E(obj@netw)$width <- abs(w) / max(abs(w))

  return(obj)
})



#------------------------------------------------------------------------------#
#' Arrange the vertex number
#'
#'@description Edit the vertices / taxa of the MetaGenomicGraph object (for example to compare two different networks):
#'@param new_taxonomy data.frame with taxonomic information about the vertices wanted in new graph.
#'  The dataset must have the same form as taxa in MetaGenomic/MetaGenomicGraph object.
#'  Naturally in new_taxonomy all the old vertices must be present.
#'@export
setGeneric(name="arrange_vertices",
           def=function(obj, new_taxonomy) standardGeneric("arrange_vertices"))
setMethod("arrange_vertices",c("MetaGenomicGraph","data.frame"),
          function(obj, new_taxonomy){

            assert(all(taxaID(obj)%in%rownames(new_taxonomy)),"find at least a taxa not present in new_taxonomy")
            assert(all(ranks(obj)==colnames(new_taxonomy)),'new_taxonomy must have the same ranks as obj')


            sample_name <- sample_name(obj)
            ntaxa <- nrow(new_taxonomy)
            taxa_name <- rownames(new_taxonomy)
            nsample <- nsample(obj)

            data <- as.data.frame(matrix(0,nrow=nsample,ncol=ntaxa,
                                         dimnames=list(sample_name,taxa_name)))
            data[,taxaID(obj)] <- obj@data

            adj <- as_adjacency_matrix(obj@netw,sparse=F,attr="weight")
            adj.new <- as.data.frame(matrix(0,nrow=ntaxa,ncol=ntaxa,
                                            dimnames=list(taxa_name,taxa_name)))
            adj.new[taxaID(obj),taxaID(obj)] <- adj
            adj.new <- as.matrix(adj.new)

            new.metagenomic <- MetaGenomic(data,obj@meta,new_taxonomy)
            return(MetaGenomicGraph(new.metagenomic,adj.new))
          })

