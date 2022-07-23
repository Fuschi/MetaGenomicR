#' S4 class to manage metagenomic data
#'
#'
#' @slot data table of data
#' @slot meta sample features
#' @slot taxa taxonomic information
#'
#' @name MetaGenomic-class
#' @rdname MetaGenomic-class
#' @exportClass MetaGenomic
MetaGenomic <- setClass(
  Class="MetaGenomic",

  slot=c(data="data.frame",
         meta="data.frame",
         taxa="data.frame")
)

# Validation
setValidity(Class="MetaGenomic", method=function(object){
  assert(is.data.frame(object@data), "data must be a data.frame")
  assert(is.data.frame(object@meta), "meta must be a data.frame")
  assert(is.data.frame(object@taxa), "taxa must be a data.frame")
  assert(nrow(object@data)==nrow(object@meta), "samples number on data and meta are different")
  assert(ncol(object@data)==nrow(object@taxa), "taxa number on data and taxa are different")
  assert(!is.null(rownames(object@data)) || !is.null(colnames(object@data)), "data must have both row/sample and column/taxa names ")
  assert(!is.null(rownames(object@meta)) || !is.null(colnames(object@meta)), "meta must have both row/sample and column/info names ")
  assert(!is.null(rownames(object@taxa)) || !is.null(colnames(object@taxa)), "data must have both row/taxa and column/ranks names ")
  assert(all(rownames(object@data)==rownames(object@meta)), "at least one of the samples names does not coincide between data and meta")
  assert(all(colnames(object@data)==rownames(object@taxa)), "at least one of the taxa names does not coincide between data and taxa")

  data.available <- object@data[!is.na(object@data)]
  assert( (is.numeric(data.available) && all(data.available>=0)), "all elements in data must be number greater than or equal to 0")
})


#' MetaGenomic Class Constructor
#' @export
MetaGenomic <- function(data, meta, taxa){
  new("MetaGenomic",data=data,meta=meta,taxa=taxa)
}


#' Method extensions to subset MetaGenomic objects.
#' @export
setMethod("[", "MetaGenomic",
  function(x,i,j,...){
    return(new("MetaGenomic",
               data=x@data[i,j,drop=FALSE],
               meta=x@meta[i, ,drop=FALSE],
               taxa=x@taxa[j, ,drop=FALSE]
))})


#' @export
setMethod(Ops, c("MetaGenomic","MetaGenomic"), function(e1,e2) {

  assert(is.empty(intersect(sample_name(e1),sample_name(e2))), "The two obj must have different sample")
  assert(all(ranks(e1)==ranks(e2)), "taxonomic ranks in e1 and e2 must be equal")
  taxa_name_intersect <- intersect(taxa_name(e1), taxa_name(e2))
  assert(all(e1@taxa[taxa_name_intersect,]==e2@taxa[taxa_name_intersect,]), "find taxa shared between the two MetaGenomic objects with different taxonomic description")
  #End Checks

  taxa_name_union <- union(taxa_name(e1), taxa_name(e2))
  ntaxa_union <- length(taxa_name_union)
  sample_name_union <- union(sample_name(e1) ,sample_name(e2))
  nsample_union <- nsample(e1) + nsample(e2)
  sample_info_union <- union(sample_info(e1),sample_info(e2))
  nsample_info <- length(sample_info_union)

  data <- as.data.frame(matrix(data=NA,
                                nrow=nsample_union,ncol=ntaxa_union,
                                dimnames=list(sample_name_union,taxa_name_union)))
  meta  <- as.data.frame(matrix(data=NA,
                                nrow=nsample_union,ncol=nsample_info,
                                dimnames=list(sample_name_union,sample_info_union)))
  taxa  <- as.data.frame(matrix(data=NA,
                                nrow=ntaxa_union,ncol=nrank(e1),
                                dimnames=list(taxa_name_union,ranks(e1))))

  data[sample_name(e1), taxa_name(e1)] <- e1@data
  data[sample_name(e2), taxa_name(e2)] <- e2@data

  meta[sample_name(e1),sample_info(e1)] <- e1@meta
  meta[sample_name(e2),sample_info(e2)] <- e2@meta

  taxa[taxa_name(e1),] <- e1@taxa
  taxa[taxa_name(e2),] <- e2@taxa

  return(new("MetaGenomic",data=data,meta=meta,taxa=taxa))
})


setMethod(f="show",
          signature="MetaGenomic",
          definition=function(object){

            print("MetaGenomic Class Object")
            print(paste("Sample Number:",nsample(object)))
            print(paste("Taxa Number:",ntaxa(object)))
            print(paste("Sample Meta Data:",paste(sample_info(object),collapse="," )))
            print(paste("Taxonomic Ranks:",paste(ranks(object),collapse=",")))

          })


#' Check if object belong to MetaGenomic Class
#'
#' @param object object to be check
#' @export
setGeneric("is.MetaGenomic", function(obj) standardGeneric("is.MetaGenomic"))
setMethod(f="is.MetaGenomic",
          definition=function(object){
            return(class(object)[1]=="MetaGenomic")
          })


#  GETTERS
#------------------------------------------------------------------------------#

#' Get table of counts
#'
#' @param obj MetaGenomic object
#' @export
setGeneric("data", function(obj) standardGeneric("data"))
setMethod("data", "MetaGenomic", function(obj) obj@data)

#' Get sample meta-data
#'
#' @param obj MetaGenomic object
#' @export
setGeneric("meta", function(obj) standardGeneric("meta"))
setMethod("meta", "MetaGenomic", function(obj) obj@meta)

#' Get taxonomic information
#'
#' @param obj MetaGenomic object
#' @export
setGeneric("taxa", function(obj) standardGeneric("taxa"))
setMethod("taxa", "MetaGenomic", function(obj) obj@taxa)

#' Get table of counts as matrix
#'
#' @param obj MetaGenomic object
#' @export
setGeneric("mdata", function(obj) standardGeneric("mdata"))
setMethod("mdata", "MetaGenomic", function(obj) as.matrix(obj@data))

#' Get centered log-ratio transformed data
#'
#' @param obj MetaGenomic object
#' @export
setGeneric("clr", function(obj) standardGeneric("clr"))
setMethod("clr", "MetaGenomic",
          function(obj){
            ref <- apply(obj@data+1, 1, function(x) mean(log(x)) )
            return(as.matrix(log(obj@data+1) - ref))
          })
setMethod("clr", "matrix",
          function(obj){
            ref <- apply(obj@data+1, 1, function(x) mean(log(x)) )
            return(as.matrix(log(obj@data+1) - ref))
          })

#' Get sample number
#'
#' @param obj MetaGenomic object
#' @export
setGeneric("nsample", function(obj) standardGeneric("nsample"))
setMethod("nsample", "MetaGenomic", function(obj) nrow(obj@data))

#' Get taxa number
#' @param obj MetaGenomic object
#' @export
setGeneric("ntaxa", function(obj) standardGeneric("ntaxa"))
setMethod("ntaxa", "MetaGenomic", function(obj) ncol(obj@data))

#' Get sample name
#' @param obj MetaGenomic object
#' @export
setGeneric("sample_name", function(obj) standardGeneric("sample_name"))
setMethod("sample_name", "MetaGenomic", function(obj) rownames(obj@data))

#' Get taxonomic rank levels
#' @param obj MetaGenomic object
#' @export
setGeneric("ranks", function(obj) standardGeneric("ranks"))
setMethod("ranks", "MetaGenomic", function(obj) colnames(obj@taxa))

#' Get number of taxonomic rank levels
#' @param obj MetaGenomic object
#' @export
setGeneric("nrank", function(obj) standardGeneric("nrank"))
setMethod("nrank", "MetaGenomic", function(obj) ncol(obj@taxa))

#' Get sample metadata
#' @param obj MetaGenomic object
#' @export
setGeneric("sample_info", function(obj) standardGeneric("sample_info"))
setMethod("sample_info", "MetaGenomic", function(obj) colnames(obj@meta))

#' Get taxa name
#' @param obj MetaGenomic object
#' @param rank taxonomic level choosen (if not set, the finest taxonomic rank is assumed)
#' @export
setGeneric("taxa_name", function(obj, rank) standardGeneric("taxa_name"))
setMethod("taxa_name", "MetaGenomic", function(obj) obj@taxa[,nrank(obj)])
setMethod("taxa_name", c("MetaGenomic","character"),
          function(obj, rank){
            assert(rank%in%ranks(obj) , paste("rank must be one this possible choises {",toString(ranks(obj)),"}"))
            return(obj@taxa[,rank])
          })


#' Get taxa ID
#' @param obj MetaGenomic object
#' @export
setGeneric("taxaID", function(obj) standardGeneric("taxaID"))
setMethod("taxaID", "MetaGenomic", function(obj) rownames(taxa(obj)))




# FUNCTIONS
#------------------------------------------------------------------------------#


#' Organize data in higher taxonomic level
#' @param obj MetaGenomic object
#' @param rank taxonomic level choosen
#' @export
setGeneric("aggregate_taxa", function(obj, rank) standardGeneric("aggregate_taxa"))
setMethod("aggregate_taxa", c("MetaGenomic","character"),
          function(obj, rank){
            assert(rank%in%ranks(obj) , paste("rank must be one this possible choises {",toString(ranks(obj)),"}"))

            different.taxa <- unique(obj@taxa[,rank])
            data.aggregate <- data.frame(matrix(NA, nrow=nsample(obj), ncol=length(different.taxa),
                                      dimnames=list(sample_name(obj),different.taxa)))

            for(taxa.i in different.taxa){
              idx <- which(taxa.i == obj@taxa[,rank])
              data.aggregate[,taxa.i] <- apply(X=obj@data, MARGIN=1, function(x) sum(x[idx]) )
            }

            taxa.aggregate <- obj@taxa[,1:which(ranks(obj)==rank)]
            taxa.aggregate <- taxa.aggregate[!duplicated(taxa.aggregate), ]
            rownames(taxa.aggregate) <- taxa.aggregate[,rank]
            taxa.aggregate <- taxa.aggregate[colnames(data.aggregate),]

            return(new("MetaGenomic",data=data.aggregate, meta=obj@meta, taxa=taxa.aggregate))
})


#' @title Smart sample selection
#' @description Allows you to choose samples based on the information stored in the meta table
#' @param obj `MetaGenomic` object
#' @param condition in progress
#' @param ... in progress
#' @return `MetaGenomic` object with select samples according to the features choose from meta table.
#' @export
setGeneric("sample_selection",function(obj, condition="AND", ...) standardGeneric("sample_selection"))
setMethod("sample_selection",c("MetaGenomic","character"),
          function(obj, condition="AND",...){

            l <- list(...)
            #checks
            #----------------------------------------#
            assert((condition=="OR" || condition=="AND"),"condition must be 'OR' or 'AND'")

            err <- unique(names(l)[names(l) %ni% sample_info(obj)])
            if(!is.empty(err)){
              err <- paste(err,collapse=",")
              err <- paste("{",err,"} parameters not found in meta info",sep="")
              stop(err)
            }

            err <- list()
            for(n in names(l)){err[[n]] <- unique(l[[n]][l[[n]] %ni% obj@meta[[n]]])}
            if(!all(sapply(err,is.empty))){
              err <- unlist(err)
              err <- paste(err,collapse=",")
              err <- paste("{",err,"} not found in their respective meta variables",sep="")
              stop(err)
            }
            #end checks
            #----------------------------------------#

            idx <- matrix(0, nrow=nsample(obj), ncol=length(l),
                          dimnames=list(sample_name(obj), names(l)))

            for(n in names(l)){
              idx[,n] <- obj@meta[[n]] %in% l[[n]]
            }

            if(condition=="OR"){
              idx <- rowSums(idx)
              idx <- idx>0
            } else if(condition=="AND"){
              idx <- apply(idx,1,prod)
              idx <- idx>0
            }

            return(obj[idx,])
          })


#' @title Smart taxa filtering
#' @description Filter the rarest taxa based on the conditions imposed in the parameters.
#' @param obj MetaGenomic object
#' @param condition in progress
#' @param prevalence in progress
#' @param detection in progress
#' @param median_non_zero in progress
#' @param data in progress
#' @param relative_abundance in progress
#' @return `MetaGenomic` object with select taxa.
#' @export
setGeneric(name="taxa_filtering",
           def=function(obj, condition, prevalence, detection,
                        median_non_zero, min_count, relative_abundance,
                        perc_count_lost, keep_discarded) standardGeneric("taxa_filtering"))
setMethod("taxa_filtering",c("MetaGenomic","ANY","ANY","ANY",
                             "ANY","ANY","ANY",
                             "ANY","ANY"),
          function(obj, condition, prevalence, detection,
                   median_non_zero, min_count, relative_abundance,
                   perc_count_lost,keep_discarded){

  if(missing(condition         )) condition          <-"AND"
  if(missing(prevalence        )) prevalence         <- NULL
  if(missing(detection         )) detection          <- NULL
  if(missing(median_non_zero   )) median_non_zero    <- NULL
  if(missing(min_count         )) min_count          <- NULL
  if(missing(relative_abundance)) relative_abundance <- NULL
  if(missing(perc_count_lost   )) perc_count_lost    <- NULL
  if(missing(keep_discarded    )) keep_discarded     <- FALSE
  #checks
  #----------------------------------------#
  assert( (condition=="OR" || condition=="AND"),"condition must be 'OR' or 'AND'")
  if(!is.null(prevalence))        {assert( (prevalence>0 && prevalence<=1), "prevalence must be a number in range (0,1]")}
  if(!is.null(detection))         {assert( (detection>0 && round(detection)==detection), "detection must be a positive integer")}
  if(!is.null(detection))         {assert( (detection <= nsample(obj)), "detection cannot exceed the number of samples")}
  if(!is.null(median_non_zero))   {assert( median_non_zero>0, "median_non_zero must be a positive number")}
  if(!is.null(min_count))             {assert( (min_count>0 && round(min_count)==min_count), "min_count must be a positive integer")}
  if(!is.null(relative_abundance)){assert( (relative_abundance>0 && relative_abundance<=1), "relative_abundance must be a number in range (0,1]")}
  if(!is.null(perc_count_lost)){assert( (perc_count_lost>0 && perc_count_lost<=1), "perc_count_lost must be a number in range (0,1]")}

  assert(is.logical(keep_discarded), 'keep_discarded must be logical')
  #end checks
  #----------------------------------------#

  # filtering type allowed
  filter.type <- c("prevalence","detection","min_count","median_non_zero","relative_abundance")

  # final data.frame of indices
  if(condition=="AND"){
    idx <- as.data.frame(matrix(1,nrow=length(filter.type),ncol=ntaxa(obj),
                                dimnames=list(filter.type,taxa_name(obj)) ))
  } else {
    idx <- as.data.frame(matrix(0,nrow=length(filter.type),ncol=ntaxa(obj),
                                dimnames=list(filter.type,taxa_name(obj)) ))
  }


  # start filter
  if(!is.null(prevalence)){
    idx["prevalence",] <- colSums(obj@data>0)/nsample(obj)
    idx["prevalence",] <- idx["prevalence",]>=prevalence
  }

  if(!is.null(detection)){
    idx["detection",] <- colSums(obj@data>0)
    idx["detection",] <- idx["detection",]>=detection
  }

  if(!is.null(median_non_zero)){
    idx.median <- apply(obj@data, 2, function(x){median(x[x>0])})
    idx.median[is.na(idx.median)] <- 0
    idx["median_non_zero",] <- idx.median>=median_non_zero
  }

  if(!is.null(min_count)){
    idx["min_count",] <- apply(obj@data, 2, function(x){any(x>=min_count)})
  }

  if(!is.null(relative_abundance)){
    idx["relative_abundance",] <- apply(obj@data/rowSums(obj@data), 2, function(x){any(x>=relative_abundance)})
  }



  # Final subset
  if(condition=="OR"){
    idx <- colSums(idx)
    idx <- idx>0
  } else if(condition=="AND"){
    idx <- apply(idx,2,prod)
    idx <- idx>0
  }

  # removing samples where too many counts have been lost
  if(!is.null(perc_count_lost)){
    lost.counts.per.sample <- 1 - (rowSums(obj@data[,idx])/rowSums(obj@data))
    idx.sample <- lost.counts.per.sample <= perc_count_lost
  } else {
    idx.sample <- 1:nsample(obj)
  }


  result <- obj[idx.sample ,idx]
  if(keep_discarded){
    result@data <- cbind(result@data, "discarded"=rowSums(result@data[idx.sample ,-idx]))
    result@taxa <- rbind(result@taxa, "discarded"=rep("merge",nrank(result)))
  }

  return(result)

})


#' @title Enlarge taxa in MetaGenomic
#' @description Replaces a larger taxonomy data.frame than the present one. All taxa present in obj must also be present in the new taxa_new object. It is used to compare different MetaGenomic object.
#' @param obj `MetaGenomic` object
#' @param new_taxa data.frame with larger taxonomic information
#' @return `MetaGenomic` object with enlarged taxonomic information. The new taxa in the data.frame data will be represented by columns of zeros.
#' @export
setGeneric("arrange_taxa",function(obj, new_taxa) standardGeneric("arrange_taxa"))
setMethod("arrange_taxa",c("MetaGenomic","data.frame"),
          function(obj, new_taxa){
            assert(all(unlist(obj@taxa)%in%unlist(new_taxa)), "at least a taxa obj is not present in new_taxa")
            assert(identical(obj@taxa,new_taxa[taxaID(obj),]),"Found different taxonomic information between the obj taxa and the new taxa")

            new_data <- matrix(0,nrow=nsample(obj),ncol=nrow(new_taxa),
                               dimnames=list(sample_name(obj),rownames(new_taxa)) )
            new_data <- as.data.frame(new_data)

            new_data[,taxaID(obj)] <- obj@data[,taxaID(obj)]

            return(MetaGenomic(data=new_data,meta=obj@meta,taxa=new_taxa))
          })





