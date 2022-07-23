#'@title load.HMP2
#'@description Create MetaGenomic object with HMP2 data
#'@export
load.HMP2 <- function(){
  return(new("MetaGenomic",data=otu.HMP2.16S,meta=meta.HMP2,taxa=taxa.otu.HMP2.16S))
}

#'@title load.2DG.16S
#'@description Create MetaGenomic object with 16S data of 2DG mice samples
#'@export
load.2DG.16S <- function(){
  return(new("MetaGenomic",data=otu.2DG.16S,meta=meta.2DG,taxa=taxa.otu.2DG.16S))
}

#'@title load.2DG.WGS
#'@description Create MetaGenomic object with shotgun data of 2DG mice samples
#'@export
load.2DG.WGS <- function(){
  return(new("MetaGenomic",data=count.2DG.WGS,meta=meta.2DG,taxa=taxa.2DG.WGS))
}

#'@title load.POMP.DO
#'@description Create MetaGenomic object with 16S data of POMP DO mice samples
#'@export
load.POMP.DO <- function(){
  return(new("MetaGenomic",data=otu_POMP_DO,meta=meta_POMP_DO,taxa=taxa_POMP_DO))
}
