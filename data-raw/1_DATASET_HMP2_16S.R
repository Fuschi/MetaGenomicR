## code to prepare `DATASET` dataset goes here




#HMP2
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#Originals
#------------------------------------------------------------------------------#
path <- system.file("extdata", "otu_HMP2_16S.csv" ,package="MetaGenomicR", mustWork=TRUE)
otu.HMP2.16S <- read.csv(file=path, row.names=1)
usethis::use_data(otu.HMP2.16S, overwrite = TRUE)

path <- system.file("extdata", "taxa_otu_HMP2_16S.csv" ,package="MetaGenomicR", mustWork=TRUE)
taxa.otu.HMP2.16S <- read.csv(file=path, row.names=1)
usethis::use_data(taxa.otu.HMP2.16S, overwrite = TRUE)

path <- system.file("extdata", "meta_HMP2.csv" ,package="MetaGenomicR", mustWork=TRUE)
meta.HMP2 <- read.csv(file=path, row.names=1)
usethis::use_data(meta.HMP2, overwrite = TRUE)

#check HMP2 data
if(any(rownames(otu.HMP2.16S)!=rownames(meta.HMP2)))stop("daiiii")
if(any(colnames(otu.HMP2.16S)!=rownames(taxa.otu.HMP2.16S)))stop("daiiii2")
if(any(rownames(otu.HMP2.16S)!=unique(rownames(otu.HMP2.16S))))stop("daiiii3")
if(any(colnames(otu.HMP2.16S)!=unique(colnames(otu.HMP2.16S))))stop("daiiii4")
if(any(rownames(taxa.otu.HMP2.16S)!=unique(rownames(taxa.otu.HMP2.16S))))stop("daiiii5")
if(any(colnames(taxa.otu.HMP2.16S)!=unique(colnames(taxa.otu.HMP2.16S))))stop("daiiii6")
if(any(rownames(meta.HMP2)!=unique(rownames(meta.HMP2))))stop("daiiii7")
if(any(colnames(meta.HMP2)!=unique(colnames(meta.HMP2))))stop("daiiii8")
rm(path)

# #Different taxanomic levels
# #------------------------------------------------------------------------------#
# #------------------------------------------------------------------------------#
# #------------------------------------------------------------------------------#
# sample.id <- row.names(otu.HMP2.16S)
#
# #GENUS
# #------------------------------------------------------------------------------#
# taxa.level <- "genus"
# different.taxa <- unique(get(taxa.level,taxa.otu.HMP2.16S))
# genus.HMP2.16S <- data.frame(matrix(NA, nrow=nrow(otu.HMP2.16S), ncol=length(different.taxa),
# 							  dimnames=list(sample.id,different.taxa)))
#
# for(i in 1:length(different.taxa)){
# 	idx <- which(different.taxa[i] == get(taxa.level,taxa.otu.HMP2.16S))
# 	genus.HMP2.16S[,i] <- apply(X=otu.HMP2.16S, MARGIN=1, function(x) sum(x[idx]) )
# }
#
# taxa.genus.HMP2.16S <- unique(subset(taxa.otu.HMP2.16S, select=c(-otu)))
# rownames(taxa.genus.HMP2.16S) <- taxa.genus.HMP2.16S$genus
#
# if(any(rownames(genus.HMP2.16S)!=rownames(otu.HMP2.16S))) stop("a beh")
# if(!any(colnames(genus.HMP2.16S) %in% get(taxa.level,taxa.otu.HMP2.16S))) stop("a boh")
# if(any(colnames(genus.HMP2.16S)!=rownames(taxa.genus.HMP2.16S))) stop("a buh")
#
#
# usethis::use_data(genus.HMP2.16S     , overwrite = TRUE)
# usethis::use_data(taxa.genus.HMP2.16S, overwrite = TRUE)
# rm(genus.HMP2.16S,taxa.genus.HMP2.16S,i,idx,different.taxa,taxa.level)
#
#
# #FAMILY
# #------------------------------------------------------------------------------#
# taxa.level <- "family"
# different.taxa <- unique(get(taxa.level,taxa.otu.HMP2.16S))
# family.HMP2.16S <- data.frame(matrix(NA, nrow=nrow(otu.HMP2.16S), ncol=length(different.taxa),
# 									dimnames=list(sample.id,different.taxa)))
#
# for(i in 1:length(different.taxa)){
# 	idx <- which(different.taxa[i] == get(taxa.level,taxa.otu.HMP2.16S))
# 	family.HMP2.16S[,i] <- apply(X=otu.HMP2.16S, MARGIN=1, function(x) sum(x[idx]) )
# }
#
# taxa.family.HMP2.16S <- unique(subset(taxa.otu.HMP2.16S, select=c(-otu,-genus)))
# rownames(taxa.family.HMP2.16S) <- taxa.family.HMP2.16S$family
#
# if(any(rownames(family.HMP2.16S)!=rownames(otu.HMP2.16S))) stop("a beh")
# if(!any(colnames(family.HMP2.16S) %in% get(taxa.level,taxa.otu.HMP2.16S))) stop("a boh")
# if(any(colnames(family.HMP2.16S)!=rownames(taxa.family.HMP2.16S))) stop("a buh")
#
#
# usethis::use_data(family.HMP2.16S     , overwrite = TRUE)
# usethis::use_data(taxa.family.HMP2.16S, overwrite = TRUE)
# rm(family.HMP2.16S,taxa.family.HMP2.16S,i,idx,different.taxa,taxa.level)
#
#
# #ORDER
# #------------------------------------------------------------------------------#
# taxa.level <- "order"
# different.taxa <- unique(get(taxa.level,taxa.otu.HMP2.16S))
# order.HMP2.16S <- data.frame(matrix(NA, nrow=nrow(otu.HMP2.16S), ncol=length(different.taxa),
# 									 dimnames=list(sample.id,different.taxa)))
#
# for(i in 1:length(different.taxa)){
# 	idx <- which(different.taxa[i] == get(taxa.level,taxa.otu.HMP2.16S))
# 	order.HMP2.16S[,i] <- apply(X=otu.HMP2.16S, MARGIN=1, function(x) sum(x[idx]) )
# }
#
# taxa.order.HMP2.16S <- unique(subset(taxa.otu.HMP2.16S, select=c(-otu,-genus,-family)))
# rownames(taxa.order.HMP2.16S) <- taxa.order.HMP2.16S$order
#
# if(any(rownames(order.HMP2.16S)!=rownames(otu.HMP2.16S))) stop("a beh")
# if(!any(colnames(order.HMP2.16S) %in% get(taxa.level,taxa.otu.HMP2.16S))) stop("a boh")
# if(any(colnames(order.HMP2.16S)!=rownames(taxa.order.HMP2.16S))) stop("a buh")
#
#
# usethis::use_data(order.HMP2.16S     , overwrite = TRUE)
# usethis::use_data(taxa.order.HMP2.16S, overwrite = TRUE)
# rm(order.HMP2.16S,taxa.order.HMP2.16S,i,idx,different.taxa,taxa.level)
#
#
# #CLASS
# #------------------------------------------------------------------------------#
# taxa.level <- "class"
# different.taxa <- unique(get(taxa.level,taxa.otu.HMP2.16S))
# class.HMP2.16S <- data.frame(matrix(NA, nrow=nrow(otu.HMP2.16S), ncol=length(different.taxa),
# 									dimnames=list(sample.id,different.taxa)))
#
# for(i in 1:length(different.taxa)){
# 	idx <- which(different.taxa[i] == get(taxa.level,taxa.otu.HMP2.16S))
# 	class.HMP2.16S[,i] <- apply(X=otu.HMP2.16S, MARGIN=1, function(x) sum(x[idx]) )
# }
#
# taxa.class.HMP2.16S <- unique(subset(taxa.otu.HMP2.16S, select=c(-otu,-genus,-family,-order)))
# rownames(taxa.class.HMP2.16S) <- taxa.class.HMP2.16S$class
#
# if(any(rownames(class.HMP2.16S)!=rownames(otu.HMP2.16S))) stop("a beh")
# if(!any(colnames(class.HMP2.16S) %in% get(taxa.level,taxa.otu.HMP2.16S))) stop("a boh")
# if(any(colnames(class.HMP2.16S)!=rownames(taxa.class.HMP2.16S))) stop("a buh")
#
#
# usethis::use_data(class.HMP2.16S     , overwrite = TRUE)
# usethis::use_data(taxa.class.HMP2.16S, overwrite = TRUE)
# rm(class.HMP2.16S,taxa.class.HMP2.16S,i,idx,different.taxa,taxa.level)
#
#
# #PHYLUM
# #------------------------------------------------------------------------------#
# taxa.level <- "phylum"
# different.taxa <- unique(get(taxa.level,taxa.otu.HMP2.16S))
# phylum.HMP2.16S <- data.frame(matrix(NA, nrow=nrow(otu.HMP2.16S), ncol=length(different.taxa),
# 									dimnames=list(sample.id,different.taxa)))
#
# for(i in 1:length(different.taxa)){
# 	idx <- which(different.taxa[i] == get(taxa.level,taxa.otu.HMP2.16S))
# 	phylum.HMP2.16S[,i] <- apply(X=otu.HMP2.16S, MARGIN=1, function(x) sum(x[idx]) )
# }
#
# taxa.phylum.HMP2.16S <- unique(subset(taxa.otu.HMP2.16S, select=c(-otu,-genus,-family,-order,-class)))
# rownames(taxa.phylum.HMP2.16S) <- taxa.phylum.HMP2.16S$phylum
#
# if(any(rownames(phylum.HMP2.16S)!=rownames(otu.HMP2.16S))) stop("a beh")
# if(!any(colnames(phylum.HMP2.16S) %in% get(taxa.level,taxa.otu.HMP2.16S))) stop("a boh")
# if(any(colnames(phylum.HMP2.16S)!=rownames(taxa.phylum.HMP2.16S))) stop("a buh")
#
#
# usethis::use_data(phylum.HMP2.16S     , overwrite = TRUE)
# usethis::use_data(taxa.phylum.HMP2.16S, overwrite = TRUE)
# rm(phylum.HMP2.16S,taxa.phylum.HMP2.16S,i,idx,different.taxa,taxa.level)
#
#
rm(list=ls())
