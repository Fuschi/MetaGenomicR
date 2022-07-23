## code to prepare `DATASET` dataset goes here

#16S
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
path <- system.file("extdata", "otu_2DG_16S.tsv" ,package="MetaGenomicR", mustWork=TRUE)
otu.2DG.16S <- read.table(file=path,header=TRUE,row.names=1,sep="\t")
otu.2DG.16S <- as.data.frame(t(otu.2DG.16S))

path <- system.file("extdata", "taxa_otu_2DG_16S.tsv" ,package="MetaGenomicR", mustWork=TRUE)
taxa.otu.2DG.16S <- read.table(file=path,header=TRUE,row.names=1,sep="\t")
taxa.otu.2DG.16S <- taxa.otu.2DG.16S[colnames(otu.2DG.16S),]
taxa.otu.2DG.16S <- cbind(taxa.otu.2DG.16S,"otu"=rownames(taxa.otu.2DG.16S))
require(stringr)
taxa.otu.2DG.16S <- as.data.frame(apply(taxa.otu.2DG.16S, c(1,2), function(x) str_replace(x," ","_")))
taxa.otu.2DG.16S <- as.data.frame(apply(taxa.otu.2DG.16S, c(1,2), function(x) str_replace(x," ","_")))

path <- system.file("extdata", "meta_2DG.txt" ,package="MetaGenomicR", mustWork=TRUE)
meta.2DG <- read.table(file=path, header=TRUE,sep="\t")
meta.2DG <- cbind(meta.2DG, "long.SampleID.16S"=rownames(otu.2DG.16S))
rownames(meta.2DG) <- meta.2DG$Unique.Collaborator.ID

#correct rownames in otu.2DG.16S
rownames(otu.2DG.16S) <- substr(rownames(otu.2DG.16S),29,35)

if(any(rownames(otu.2DG.16S)!=rownames(meta.2DG))) stop("a bah")
if(any(colnames(otu.2DG.16S)!=rownames(taxa.otu.2DG.16S))) stop("a beh")
rm(path)

usethis::use_data(otu.2DG.16S     , overwrite = TRUE)
usethis::use_data(taxa.otu.2DG.16S, overwrite = TRUE)
usethis::use_data(meta.2DG,         overwrite = TRUE)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#



#SHOTGUN
#--------------------------#
#-------------------------------------------------#
#--------------------------#
samples.removed.2DG.WGS <- list()

#ARCHEA
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#GENUS
#------------------------------#
path <- system.file("extdata","pathseq.archaea.genus.readcount.txt",package="MetaGenomicR",mustWork=TRUE)
res.archaea.genus.count.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

path <- system.file("extdata","pathseq.archaea.genus.score.txt",package="MetaGenomicR",mustWork=TRUE)
res.archaea.genus.score.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

# verify colnames and rownames of counts and score are equal
if(any(colnames(res.archaea.genus.count.2DG.WGS)!=colnames(res.archaea.genus.score.2DG.WGS)))stop("so sad 1")
if(any(rownames(res.archaea.genus.count.2DG.WGS)!=rownames(res.archaea.genus.score.2DG.WGS)))stop("so sad 2")

# take only 20 sample and remove one without meta
# transpose data to obtain sample as row
archaea.genus.count.2DG.WGS <- as.data.frame(t(res.archaea.genus.count.2DG.WGS[,5:24]))
archaea.genus.score.2DG.WGS <- as.data.frame(t(res.archaea.genus.score.2DG.WGS[,5:24]))

# take sample name verify are in the same order as meta data
long.sampleID.archaea.genus <- rownames(archaea.genus.count.2DG.WGS)
for(i in 1:length(long.sampleID.archaea.genus)){
  if(!grepl(meta.2DG$Unique.Collaborator.ID[i],long.sampleID.archaea.genus[i],fixed=TRUE)) stop(print(paste("sad 2.1 -- i",i)))
}

# add long names to meta and substitute the short name in rownames of data
meta.2DG <- cbind(meta.2DG,"long.sampleID.archaea.genus"=long.sampleID.archaea.genus)
rownames(archaea.genus.count.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID
rownames(archaea.genus.score.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID

# verify taxonomic information are equal for score and genus
if(any(res.archaea.genus.count.2DG.WGS[,1:4]!=res.archaea.genus.score.2DG.WGS[,1:4])) stop("so sad 3")

# take taxonomic informations
taxa.archaea.genus.2DG.WGS.raw <- res.archaea.genus.count.2DG.WGS[,1:4]

# verify taxonomix information are equal to data columns
if(any(rownames(taxa.archaea.genus.2DG.WGS.raw)!=colnames(archaea.genus.count.2DG.WGS))) stop("so sad 4")
if(any(rownames(taxa.archaea.genus.2DG.WGS.raw)!=colnames(archaea.genus.score.2DG.WGS))) stop("so sad 5")

#collect one sample removed
samples.removed.2DG.WGS$archaea.genus.count.2DG.WGS <- res.archaea.genus.count.2DG.WGS[,25]
samples.removed.2DG.WGS$archaea.genus.score.2DG.WGS <- res.archaea.genus.count.2DG.WGS[,25]

idx <- which(colSums(archaea.genus.count.2DG.WGS)!=0)
archaea.genus.count.2DG.WGS <- archaea.genus.count.2DG.WGS[,idx]
archaea.genus.score.2DG.WGS <- archaea.genus.score.2DG.WGS[,idx]
taxa.archaea.genus.2DG.WGS.raw <- taxa.archaea.genus.2DG.WGS.raw[idx,]
usethis::use_data(archaea.genus.count.2DG.WGS , overwrite = TRUE)
usethis::use_data(archaea.genus.score.2DG.WGS, overwrite = TRUE)
usethis::use_data(taxa.archaea.genus.2DG.WGS.raw , overwrite = TRUE)
usethis::use_data(samples.removed.2DG.WGS, overwrite=TRUE )
#------------------------------#


#SPECIES
#------------------------------#
path <- system.file("extdata","pathseq.archaea.species.readcount.txt",package="MetaGenomicR",mustWork=TRUE)
res.archaea.species.count.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

path <- system.file("extdata","pathseq.archaea.species.score.txt",package="MetaGenomicR",mustWork=TRUE)
res.archaea.species.score.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

# verify colnames and rownames of counts and score are equal
if(any(colnames(res.archaea.species.count.2DG.WGS)!=colnames(res.archaea.species.score.2DG.WGS)))stop("so sad 1")
if(any(rownames(res.archaea.species.count.2DG.WGS)!=rownames(res.archaea.species.score.2DG.WGS)))stop("so sad 2")

# take only 20 sample and remove one without meta
# transpose data to obtain sample as row
archaea.species.count.2DG.WGS <- as.data.frame(t(res.archaea.species.count.2DG.WGS[,5:24]))
archaea.species.score.2DG.WGS <- as.data.frame(t(res.archaea.species.score.2DG.WGS[,5:24]))

# take sample name verify are in the same order as meta data
long.sampleID.archaea.species <- rownames(archaea.species.count.2DG.WGS)
for(i in 1:length(long.sampleID.archaea.species)){
  if(!grepl(meta.2DG$Unique.Collaborator.ID[i],long.sampleID.archaea.species[i],fixed=TRUE)) stop(print(paste("sad 2.1 -- i",i)))
}

# add long names to meta and substitute the short name in rownames of data
if(any(long.sampleID.archaea.genus!=long.sampleID.archaea.species))stop("so sad")
rownames(archaea.species.count.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID
rownames(archaea.species.score.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID

# verify taxonomic information are equal for score and genus
if(any(res.archaea.species.count.2DG.WGS[,1:4]!=res.archaea.species.score.2DG.WGS[,1:4])) stop("so sad 3")

# take taxonomic informations
taxa.archaea.species.2DG.WGS.raw <- res.archaea.species.count.2DG.WGS[,1:4]

# verify taxonomic information are equal to data columns
if(any(rownames(taxa.archaea.species.2DG.WGS.raw)!=colnames(archaea.species.count.2DG.WGS))) stop("so sad 4")
if(any(rownames(taxa.archaea.species.2DG.WGS.raw)!=colnames(archaea.species.score.2DG.WGS))) stop("so sad 5")

#collect one sample removed
samples.removed.2DG.WGS$archaea.species.count.2DG.WGS <- res.archaea.species.count.2DG.WGS[,25]
samples.removed.2DG.WGS$archaea.species.score.2DG.WGS <- res.archaea.species.count.2DG.WGS[,25]


idx <- which(colSums(archaea.species.count.2DG.WGS)!=0)
archaea.species.count.2DG.WGS <- archaea.species.count.2DG.WGS[,idx]
archaea.species.score.2DG.WGS <- archaea.species.score.2DG.WGS[,idx]
taxa.archaea.species.2DG.WGS.raw <- taxa.archaea.species.2DG.WGS.raw[idx,]
all(colnames(archaea.species.count.2DG.WGS)==rownames(taxa.archaea.species.2DG.WGS.raw))
usethis::use_data(archaea.species.count.2DG.WGS , overwrite = TRUE)
usethis::use_data(archaea.species.score.2DG.WGS, overwrite = TRUE)
usethis::use_data(taxa.archaea.species.2DG.WGS.raw , overwrite = TRUE)
usethis::use_data(samples.removed.2DG.WGS, overwrite=TRUE )
#------------------------------#





#FUNGI
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#GENUS
#------------------------------#
path <- system.file("extdata","pathseq.fungi.genus.readcount.txt",package="MetaGenomicR",mustWork=TRUE)
res.fungi.genus.count.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

path <- system.file("extdata","pathseq.fungi.genus.score.txt",package="MetaGenomicR",mustWork=TRUE)
res.fungi.genus.score.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

# verify colnames and rownames of counts and score are equal
if(any(colnames(res.fungi.genus.count.2DG.WGS)!=colnames(res.fungi.genus.score.2DG.WGS)))stop("so sad 1")
if(any(rownames(res.fungi.genus.count.2DG.WGS)!=rownames(res.fungi.genus.score.2DG.WGS)))stop("so sad 2")

# take only 20 sample and remove one without meta
# transpose data to obtain sample as row
fungi.genus.count.2DG.WGS <- as.data.frame(t(res.fungi.genus.count.2DG.WGS[,5:24]))
fungi.genus.score.2DG.WGS <- as.data.frame(t(res.fungi.genus.score.2DG.WGS[,5:24]))

# take sample name verify are in the same order as meta data
long.sampleID.fungi.genus <- rownames(fungi.genus.count.2DG.WGS)
for(i in 1:length(long.sampleID.fungi.genus)){
  if(!grepl(meta.2DG$Unique.Collaborator.ID[i],long.sampleID.fungi.genus[i],fixed=TRUE)) stop(print(paste("sad 2.1 -- i",i)))
}

# add long names to meta and substitute the short name in rownames of data
meta.2DG <- cbind(meta.2DG,"long.sampleID.fungi.genus"=long.sampleID.fungi.genus)
rownames(fungi.genus.count.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID
rownames(fungi.genus.score.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID

# verify taxonomic information are equal for score and genus
if(any(res.fungi.genus.count.2DG.WGS[,1:4]!=res.fungi.genus.score.2DG.WGS[,1:4])) stop("so sad 3")

# take taxonomic informations
taxa.fungi.genus.2DG.WGS.raw <- res.fungi.genus.count.2DG.WGS[,1:4]

# verify taxonomix information are equal to data columns
if(any(rownames(taxa.fungi.genus.2DG.WGS.raw)!=colnames(fungi.genus.count.2DG.WGS))) stop("so sad 4")
if(any(rownames(taxa.fungi.genus.2DG.WGS.raw)!=colnames(fungi.genus.score.2DG.WGS))) stop("so sad 5")

#collect one sample removed
samples.removed.2DG.WGS$fungi.genus.count.2DG.WGS <- res.fungi.genus.count.2DG.WGS[,25]
samples.removed.2DG.WGS$fungi.genus.score.2DG.WGS <- res.fungi.genus.count.2DG.WGS[,25]


idx <- which(colSums(fungi.genus.count.2DG.WGS)!=0)
fungi.genus.count.2DG.WGS <- fungi.genus.count.2DG.WGS[,idx]
fungi.genus.score.2DG.WGS <- fungi.genus.score.2DG.WGS[,idx]
taxa.fungi.genus.2DG.WGS.raw <- taxa.fungi.genus.2DG.WGS.raw[idx,]
all(colnames(archaea.genus.count.2DG.WGS)==rownames(taxa.archaea.genus.2DG.WGS.raw))
usethis::use_data(fungi.genus.count.2DG.WGS , overwrite = TRUE)
usethis::use_data(fungi.genus.score.2DG.WGS, overwrite = TRUE)
usethis::use_data(taxa.fungi.genus.2DG.WGS.raw , overwrite = TRUE)
usethis::use_data(samples.removed.2DG.WGS, overwrite=TRUE )
#------------------------------#


#SPECIES
#------------------------------#
path <- system.file("extdata","pathseq.fungi.species.readcount.txt",package="MetaGenomicR",mustWork=TRUE)
res.fungi.species.count.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

path <- system.file("extdata","pathseq.fungi.species.score.txt",package="MetaGenomicR",mustWork=TRUE)
res.fungi.species.score.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

# verify colnames and rownames of counts and score are equal
if(any(colnames(res.fungi.species.count.2DG.WGS)!=colnames(res.fungi.species.score.2DG.WGS)))stop("so sad 1")
if(any(rownames(res.fungi.species.count.2DG.WGS)!=rownames(res.fungi.species.score.2DG.WGS)))stop("so sad 2")

# take only 20 sample and remove one without meta
# transpose data to obtain sample as row
fungi.species.count.2DG.WGS <- as.data.frame(t(res.fungi.species.count.2DG.WGS[,5:24]))
fungi.species.score.2DG.WGS <- as.data.frame(t(res.fungi.species.score.2DG.WGS[,5:24]))

# take sample name verify are in the same order as meta data
long.sampleID.fungi.species <- rownames(fungi.species.count.2DG.WGS)
for(i in 1:length(long.sampleID.fungi.species)){
  if(!grepl(meta.2DG$Unique.Collaborator.ID[i],long.sampleID.fungi.species[i],fixed=TRUE)) stop(print(paste("sad 2.1 -- i",i)))
}

# add long names to meta and substitute the short name in rownames of data
if(any(long.sampleID.fungi.genus!=long.sampleID.fungi.species))stop("so sad")
rownames(fungi.species.count.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID
rownames(fungi.species.score.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID

# verify taxonomic information are equal for score and genus
if(any(res.fungi.species.count.2DG.WGS[,1:4]!=res.fungi.species.score.2DG.WGS[,1:4])) stop("so sad 3")

# take taxonomic informations
taxa.fungi.species.2DG.WGS.raw <- res.fungi.species.count.2DG.WGS[,1:4]

# verify taxonomic information are equal to data columns
if(any(rownames(taxa.fungi.species.2DG.WGS.raw)!=colnames(fungi.species.count.2DG.WGS))) stop("so sad 4")
if(any(rownames(taxa.fungi.species.2DG.WGS.raw)!=colnames(fungi.species.score.2DG.WGS))) stop("so sad 5")

#collect one sample removed
samples.removed.2DG.WGS$fungi.species.count.2DG.WGS <- res.fungi.species.count.2DG.WGS[,25]
samples.removed.2DG.WGS$fungi.species.score.2DG.WGS <- res.fungi.species.count.2DG.WGS[,25]


idx <- which(colSums(fungi.species.count.2DG.WGS)!=0)
fungi.species.count.2DG.WGS <- fungi.species.count.2DG.WGS[,idx]
fungi.species.score.2DG.WGS <- fungi.species.score.2DG.WGS[,idx]
taxa.fungi.species.2DG.WGS.raw <- taxa.fungi.species.2DG.WGS.raw[idx,]
all(colnames(archaea.species.count.2DG.WGS)==rownames(taxa.archaea.species.2DG.WGS.raw))
usethis::use_data(fungi.species.count.2DG.WGS , overwrite = TRUE)
usethis::use_data(fungi.species.score.2DG.WGS, overwrite = TRUE)
usethis::use_data(taxa.fungi.species.2DG.WGS.raw , overwrite = TRUE)
usethis::use_data(samples.removed.2DG.WGS, overwrite=TRUE )
#------------------------------#







#VIRUSES
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#GENUS
#------------------------------#
path <- system.file("extdata","pathseq.viruses.genus.readcount.txt",package="MetaGenomicR",mustWork=TRUE)
res.viruses.genus.count.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

path <- system.file("extdata","pathseq.viruses.genus.score.txt",package="MetaGenomicR",mustWork=TRUE)
res.viruses.genus.score.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

# verify colnames and rownames of counts and score are equal
if(any(colnames(res.viruses.genus.count.2DG.WGS)!=colnames(res.viruses.genus.score.2DG.WGS)))stop("so sad 1")
if(any(rownames(res.viruses.genus.count.2DG.WGS)!=rownames(res.viruses.genus.score.2DG.WGS)))stop("so sad 2")

# take only 20 sample and remove one without meta
# transpose data to obtain sample as row
viruses.genus.count.2DG.WGS <- as.data.frame(t(res.viruses.genus.count.2DG.WGS[,5:24]))
viruses.genus.score.2DG.WGS <- as.data.frame(t(res.viruses.genus.score.2DG.WGS[,5:24]))

# take sample name verify are in the same order as meta data
long.sampleID.viruses.genus <- rownames(viruses.genus.count.2DG.WGS)
for(i in 1:length(long.sampleID.viruses.genus)){
  if(!grepl(meta.2DG$Unique.Collaborator.ID[i],long.sampleID.viruses.genus[i],fixed=TRUE)) stop(print(paste("sad 2.1 -- i",i)))
}

# add long names to meta and substitute the short name in rownames of data
meta.2DG <- cbind(meta.2DG,"long.sampleID.viruses.genus"=long.sampleID.viruses.genus)
rownames(viruses.genus.count.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID
rownames(viruses.genus.score.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID

# verify taxonomic information are equal for score and genus
if(any(res.viruses.genus.count.2DG.WGS[,1:4]!=res.viruses.genus.score.2DG.WGS[,1:4])) stop("so sad 3")

# take taxonomic informations
taxa.viruses.genus.2DG.WGS.raw <- res.viruses.genus.count.2DG.WGS[,1:4]

# verify taxonomix information are equal to data columns
if(any(rownames(taxa.viruses.genus.2DG.WGS.raw)!=colnames(viruses.genus.count.2DG.WGS))) stop("so sad 4")
if(any(rownames(taxa.viruses.genus.2DG.WGS.raw)!=colnames(viruses.genus.score.2DG.WGS))) stop("so sad 5")

#collect one sample removed
samples.removed.2DG.WGS$viruses.genus.count.2DG.WGS <- res.viruses.genus.count.2DG.WGS[,25]
samples.removed.2DG.WGS$viruses.genus.score.2DG.WGS <- res.viruses.genus.count.2DG.WGS[,25]


idx <- which(colSums(viruses.genus.count.2DG.WGS)!=0)
viruses.genus.count.2DG.WGS <- viruses.genus.count.2DG.WGS[,idx]
viruses.genus.score.2DG.WGS <- viruses.genus.score.2DG.WGS[,idx]
taxa.viruses.genus.2DG.WGS.raw <- taxa.viruses.genus.2DG.WGS.raw[idx,]
all(colnames(archaea.genus.count.2DG.WGS)==rownames(taxa.archaea.genus.2DG.WGS.raw))
usethis::use_data(viruses.genus.count.2DG.WGS , overwrite = TRUE)
usethis::use_data(viruses.genus.score.2DG.WGS, overwrite = TRUE)
usethis::use_data(taxa.viruses.genus.2DG.WGS.raw , overwrite = TRUE)
usethis::use_data(samples.removed.2DG.WGS, overwrite=TRUE )
#------------------------------#


#SPECIES
#------------------------------#
path <- system.file("extdata","pathseq.viruses.species.readcount.txt",package="MetaGenomicR",mustWork=TRUE)
res.viruses.species.count.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

path <- system.file("extdata","pathseq.viruses.species.score.txt",package="MetaGenomicR",mustWork=TRUE)
res.viruses.species.score.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

# verify colnames and rownames of counts and score are equal
if(any(colnames(res.viruses.species.count.2DG.WGS)!=colnames(res.viruses.species.score.2DG.WGS)))stop("so sad 1")
if(any(rownames(res.viruses.species.count.2DG.WGS)!=rownames(res.viruses.species.score.2DG.WGS)))stop("so sad 2")

# take only 20 sample and remove one without meta
# transpose data to obtain sample as row
viruses.species.count.2DG.WGS <- as.data.frame(t(res.viruses.species.count.2DG.WGS[,5:24]))
viruses.species.score.2DG.WGS <- as.data.frame(t(res.viruses.species.score.2DG.WGS[,5:24]))

# take sample name verify are in the same order as meta data
long.sampleID.viruses.species <- rownames(viruses.species.count.2DG.WGS)
for(i in 1:length(long.sampleID.viruses.species)){
  if(!grepl(meta.2DG$Unique.Collaborator.ID[i],long.sampleID.viruses.species[i],fixed=TRUE)) stop(print(paste("sad 2.1 -- i",i)))
}

# add long names to meta and substitute the short name in rownames of data
if(any(long.sampleID.viruses.genus!=long.sampleID.viruses.species))stop("so sad")
rownames(viruses.species.count.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID
rownames(viruses.species.score.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID

# verify taxonomic information are equal for score and genus
if(any(res.viruses.species.count.2DG.WGS[,1:4]!=res.viruses.species.score.2DG.WGS[,1:4])) stop("so sad 3")

# take taxonomic informations
taxa.viruses.species.2DG.WGS.raw <- res.viruses.species.count.2DG.WGS[,1:4]

# verify taxonomic information are equal to data columns
if(any(rownames(taxa.viruses.species.2DG.WGS.raw)!=colnames(viruses.species.count.2DG.WGS))) stop("so sad 4")
if(any(rownames(taxa.viruses.species.2DG.WGS.raw)!=colnames(viruses.species.score.2DG.WGS))) stop("so sad 5")

#collect one sample removed
samples.removed.2DG.WGS$viruses.species.count.2DG.WGS <- res.viruses.species.count.2DG.WGS[,25]
samples.removed.2DG.WGS$viruses.species.score.2DG.WGS <- res.viruses.species.count.2DG.WGS[,25]


idx <- which(colSums(viruses.species.count.2DG.WGS)!=0)
viruses.species.count.2DG.WGS <- viruses.species.count.2DG.WGS[,idx]
viruses.species.score.2DG.WGS <- viruses.species.score.2DG.WGS[,idx]
taxa.viruses.species.2DG.WGS.raw <- taxa.viruses.species.2DG.WGS.raw[idx,]
all(colnames(archaea.species.count.2DG.WGS)==rownames(taxa.archaea.species.2DG.WGS.raw))
usethis::use_data(viruses.species.count.2DG.WGS , overwrite = TRUE)
usethis::use_data(viruses.species.score.2DG.WGS, overwrite = TRUE)
usethis::use_data(taxa.viruses.species.2DG.WGS.raw , overwrite = TRUE)
usethis::use_data(samples.removed.2DG.WGS, overwrite=TRUE )
#------------------------------#













#BACTERIA
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#GENUS
#------------------------------#
path <- system.file("extdata","pathseq.bacteria.genus.readcount.txt",package="MetaGenomicR",mustWork=TRUE)
res.bacteria.genus.count.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

#-------------------------------------------#
#STRANGE SITUATION HERE, the seventh row of bacteria.genus.readcount has 22 different sample????
#I remove the last one
res.bacteria.genus.count.2DG.WGS$X <- NULL
#-------------------------------------------#

path <- system.file("extdata","pathseq.bacteria.genus.score.txt",package="MetaGenomicR",mustWork=TRUE)
res.bacteria.genus.score.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

# verify colnames and rownames of counts and score are equal
if(any(colnames(res.bacteria.genus.count.2DG.WGS)!=colnames(res.bacteria.genus.score.2DG.WGS)))stop("so sad 1")
if(any(rownames(res.bacteria.genus.count.2DG.WGS)!=rownames(res.bacteria.genus.score.2DG.WGS)))stop("so sad 2")

# take only 20 sample and remove one without meta
# transpose data to obtain sample as row
bacteria.genus.count.2DG.WGS <- as.data.frame(t(res.bacteria.genus.count.2DG.WGS[,5:24]))
bacteria.genus.score.2DG.WGS <- as.data.frame(t(res.bacteria.genus.score.2DG.WGS[,5:24]))

# take sample name verify are in the same order as meta data
long.sampleID.bacteria.genus <- rownames(bacteria.genus.count.2DG.WGS)
for(i in 1:length(long.sampleID.bacteria.genus)){
  if(!grepl(meta.2DG$Unique.Collaborator.ID[i],long.sampleID.bacteria.genus[i],fixed=TRUE)) stop(print(paste("sad 2.1 -- i",i)))
}

# add long names to meta and substitute the short name in rownames of data
meta.2DG <- cbind(meta.2DG,"long.sampleID.bacteria.genus"=long.sampleID.bacteria.genus)
rownames(bacteria.genus.count.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID
rownames(bacteria.genus.score.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID

# verify taxonomic information are equal for score and genus
if(any(res.bacteria.genus.count.2DG.WGS[,1:4]!=res.bacteria.genus.score.2DG.WGS[,1:4])) stop("so sad 3")

# take taxonomic informations
taxa.bacteria.genus.2DG.WGS.raw <- res.bacteria.genus.count.2DG.WGS[,1:4]

# verify taxonomix information are equal to data columns
if(any(rownames(taxa.bacteria.genus.2DG.WGS.raw)!=colnames(bacteria.genus.count.2DG.WGS))) stop("so sad 4")
if(any(rownames(taxa.bacteria.genus.2DG.WGS.raw)!=colnames(bacteria.genus.score.2DG.WGS))) stop("so sad 5")

#collect one sample removed
samples.removed.2DG.WGS$bacteria.genus.count.2DG.WGS <- res.bacteria.genus.count.2DG.WGS[,25]
samples.removed.2DG.WGS$bacteria.genus.score.2DG.WGS <- res.bacteria.genus.count.2DG.WGS[,25]


idx <- which(colSums(bacteria.genus.count.2DG.WGS)!=0)
bacteria.genus.count.2DG.WGS <- bacteria.genus.count.2DG.WGS[,idx]
bacteria.genus.score.2DG.WGS <- bacteria.genus.score.2DG.WGS[,idx]
taxa.bacteria.genus.2DG.WGS.raw <- taxa.bacteria.genus.2DG.WGS.raw[idx,]
all(colnames(archaea.genus.count.2DG.WGS)==rownames(taxa.archaea.genus.2DG.WGS.raw))
usethis::use_data(bacteria.genus.count.2DG.WGS , overwrite = TRUE)
usethis::use_data(bacteria.genus.score.2DG.WGS, overwrite = TRUE)
usethis::use_data(taxa.bacteria.genus.2DG.WGS.raw , overwrite = TRUE)
usethis::use_data(samples.removed.2DG.WGS, overwrite=TRUE )
#------------------------------#


#SPECIES
#------------------------------#
path <- system.file("extdata","pathseq.bacteria.species.readcount.txt",package="MetaGenomicR",mustWork=TRUE)
res.bacteria.species.count.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

path <- system.file("extdata","pathseq.bacteria.species.score.txt",package="MetaGenomicR",mustWork=TRUE)
res.bacteria.species.score.2DG.WGS <- read.table(path,header=TRUE,row.names=1,sep="\t")

# verify colnames and rownames of counts and score are equal
if(any(colnames(res.bacteria.species.count.2DG.WGS)!=colnames(res.bacteria.species.score.2DG.WGS)))stop("so sad 1")
if(any(rownames(res.bacteria.species.count.2DG.WGS)!=rownames(res.bacteria.species.score.2DG.WGS)))stop("so sad 2")

# take only 20 sample and remove one without meta
# transpose data to obtain sample as row
bacteria.species.count.2DG.WGS <- as.data.frame(t(res.bacteria.species.count.2DG.WGS[,5:24]))
bacteria.species.score.2DG.WGS <- as.data.frame(t(res.bacteria.species.score.2DG.WGS[,5:24]))

# take sample name verify are in the same order as meta data
long.sampleID.bacteria.species <- rownames(bacteria.species.count.2DG.WGS)
for(i in 1:length(long.sampleID.bacteria.species)){
  if(!grepl(meta.2DG$Unique.Collaborator.ID[i],long.sampleID.bacteria.species[i],fixed=TRUE)) stop(print(paste("sad 2.1 -- i",i)))
}

# add long names to meta and substitute the short name in rownames of data
if(any(long.sampleID.bacteria.genus!=long.sampleID.bacteria.species))stop("so sad")
rownames(bacteria.species.count.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID
rownames(bacteria.species.score.2DG.WGS) <- meta.2DG$Unique.Collaborator.ID

# verify taxonomic information are equal for score and genus
if(any(res.bacteria.species.count.2DG.WGS[,1:4]!=res.bacteria.species.score.2DG.WGS[,1:4])) stop("so sad 3")

# take taxonomic informations
taxa.bacteria.species.2DG.WGS.raw <- res.bacteria.species.count.2DG.WGS[,1:4]

# verify taxonomic information are equal to data columns
if(any(rownames(taxa.bacteria.species.2DG.WGS.raw)!=colnames(bacteria.species.count.2DG.WGS))) stop("so sad 4")
if(any(rownames(taxa.bacteria.species.2DG.WGS.raw)!=colnames(bacteria.species.score.2DG.WGS))) stop("so sad 5")

#collect one sample removed
samples.removed.2DG.WGS$bacteria.species.count.2DG.WGS <- res.bacteria.species.count.2DG.WGS[,25]
samples.removed.2DG.WGS$bacteria.species.score.2DG.WGS <- res.bacteria.species.count.2DG.WGS[,25]


idx <- which(colSums(bacteria.species.count.2DG.WGS)!=0)
bacteria.species.count.2DG.WGS <- bacteria.species.count.2DG.WGS[,idx]
bacteria.species.score.2DG.WGS <- bacteria.species.score.2DG.WGS[,idx]
taxa.bacteria.species.2DG.WGS.raw <- taxa.bacteria.species.2DG.WGS.raw[idx,]
all(colnames(archaea.species.count.2DG.WGS)==rownames(taxa.archaea.species.2DG.WGS.raw))
usethis::use_data(bacteria.species.count.2DG.WGS , overwrite = TRUE)
usethis::use_data(bacteria.species.score.2DG.WGS, overwrite = TRUE)
usethis::use_data(taxa.bacteria.species.2DG.WGS.raw , overwrite = TRUE)
usethis::use_data(samples.removed.2DG.WGS, overwrite=TRUE )
#------------------------------#




#Verify long name in shotgun data are the same
#------------------------------#
if(any(meta.2DG$long.sampleID.archaea.genus!=meta.2DG$long.sampleID.fungi.genus)) stop("so sad 1")
if(any(meta.2DG$long.sampleID.archaea.genus!=meta.2DG$long.sampleID.viruses.genus)) stop("so sad 2")
if(any(meta.2DG$long.sampleID.archaea.genus!=meta.2DG$long.sampleID.bacteria.genus)) stop("so sad 3")

# are all the same and resume in only a variable/column
meta.2DG$long.SampleID.WGS <- meta.2DG$long.sampleID.archaea.genus
meta.2DG$long.sampleID.archaea.genus  <- NULL
meta.2DG$long.sampleID.fungi.genus    <- NULL
meta.2DG$long.sampleID.viruses.genus  <- NULL
meta.2DG$long.sampleID.bacteria.genus <- NULL
rownames(meta.2DG) <- meta.2DG$Unique.Collaborator.ID
usethis::use_data(meta.2DG, overwrite=TRUE )













#AT THE END COLLECT 16S 2DG IN HIGHERE TAXONOMIC LEVELS
#------------------------------------------------------------------------------#
sample.id <- row.names(otu.2DG.16S)

# #GENUS
# #------------------------------------------------------------------------------#
# taxa.level <- "genus"
# different.taxa <- unique(get(taxa.level,taxa.otu.2DG.16S))
# genus.2DG.16S <- data.frame(matrix(NA, nrow=nrow(otu.2DG.16S), ncol=length(different.taxa),
#                                     dimnames=list(sample.id,different.taxa)))
#
# for(i in 1:length(different.taxa)){
#   idx <- which(different.taxa[i] == get(taxa.level,taxa.otu.2DG.16S))
#   genus.2DG.16S[,i] <- apply(X=otu.2DG.16S, MARGIN=1, function(x) sum(x[idx]) )
# }
#
# taxa.genus.2DG.16S <- unique(subset(taxa.otu.2DG.16S, select=c(-otu)))
# rownames(taxa.genus.2DG.16S) <- taxa.genus.2DG.16S$genus
#
# if(any(rownames(genus.2DG.16S)!=rownames(otu.2DG.16S))) stop("a beh")
# if(!any(colnames(genus.2DG.16S) %in% get(taxa.level,taxa.otu.2DG.16S))) stop("a boh")
# if(any(colnames(genus.2DG.16S)!=rownames(taxa.genus.2DG.16S))) stop("a buh")
#
#
# usethis::use_data(genus.2DG.16S     , overwrite = TRUE)
# usethis::use_data(taxa.genus.2DG.16S, overwrite = TRUE)
# rm(genus.2DG.16S,taxa.genus.2DG.16S,i,idx,different.taxa,taxa.level)
#
#
# #FAMILY
# #------------------------------------------------------------------------------#
# taxa.level <- "family"
# different.taxa <- unique(get(taxa.level,taxa.otu.2DG.16S))
# family.2DG.16S <- data.frame(matrix(NA, nrow=nrow(otu.2DG.16S), ncol=length(different.taxa),
#                                    dimnames=list(sample.id,different.taxa)))
#
# for(i in 1:length(different.taxa)){
#   idx <- which(different.taxa[i] == get(taxa.level,taxa.otu.2DG.16S))
#   family.2DG.16S[,i] <- apply(X=otu.2DG.16S, MARGIN=1, function(x) sum(x[idx]) )
# }
#
# taxa.family.2DG.16S <- unique(subset(taxa.otu.2DG.16S, select=c(-otu,-genus)))
# rownames(taxa.family.2DG.16S) <- taxa.family.2DG.16S$family
#
# if(any(rownames(family.2DG.16S)!=rownames(otu.2DG.16S))) stop("a beh")
# if(!any(colnames(family.2DG.16S) %in% get(taxa.level,taxa.otu.2DG.16S))) stop("a boh")
# if(any(colnames(family.2DG.16S)!=rownames(taxa.family.2DG.16S))) stop("a buh")
#
#
# usethis::use_data(family.2DG.16S     , overwrite = TRUE)
# usethis::use_data(taxa.family.2DG.16S, overwrite = TRUE)
# rm(family.2DG.16S,taxa.family.2DG.16S,i,idx,different.taxa,taxa.level)
#
#
# #ORDER
# #------------------------------------------------------------------------------#
# taxa.level <- "order"
# different.taxa <- unique(get(taxa.level,taxa.otu.2DG.16S))
# order.2DG.16S <- data.frame(matrix(NA, nrow=nrow(otu.2DG.16S), ncol=length(different.taxa),
#                                     dimnames=list(sample.id,different.taxa)))
#
# for(i in 1:length(different.taxa)){
#   idx <- which(different.taxa[i] == get(taxa.level,taxa.otu.2DG.16S))
#   order.2DG.16S[,i] <- apply(X=otu.2DG.16S, MARGIN=1, function(x) sum(x[idx]) )
# }
#
# taxa.order.2DG.16S <- unique(subset(taxa.otu.2DG.16S, select=c(-otu,-genus,-family)))
# rownames(taxa.order.2DG.16S) <- taxa.order.2DG.16S$order
#
# if(any(rownames(order.2DG.16S)!=rownames(otu.2DG.16S))) stop("a beh")
# if(!any(colnames(order.2DG.16S) %in% get(taxa.level,taxa.otu.2DG.16S))) stop("a boh")
# if(any(colnames(order.2DG.16S)!=rownames(taxa.order.2DG.16S))) stop("a buh")
#
#
# usethis::use_data(order.2DG.16S     , overwrite = TRUE)
# usethis::use_data(taxa.order.2DG.16S, overwrite = TRUE)
# rm(order.2DG.16S,taxa.order.2DG.16S,i,idx,different.taxa,taxa.level)
#
#
# #CLASS
# #------------------------------------------------------------------------------#
# taxa.level <- "class"
# different.taxa <- unique(get(taxa.level,taxa.otu.2DG.16S))
# class.2DG.16S <- data.frame(matrix(NA, nrow=nrow(otu.2DG.16S), ncol=length(different.taxa),
#                                     dimnames=list(sample.id,different.taxa)))
#
# for(i in 1:length(different.taxa)){
#   idx <- which(different.taxa[i] == get(taxa.level,taxa.otu.2DG.16S))
#   class.2DG.16S[,i] <- apply(X=otu.2DG.16S, MARGIN=1, function(x) sum(x[idx]) )
# }
#
# taxa.class.2DG.16S <- unique(subset(taxa.otu.2DG.16S, select=c(-otu,-genus,-family,-order)))
# rownames(taxa.class.2DG.16S) <- taxa.class.2DG.16S$class
#
# if(any(rownames(class.2DG.16S)!=rownames(otu.2DG.16S))) stop("a beh")
# if(!any(colnames(class.2DG.16S) %in% get(taxa.level,taxa.otu.2DG.16S))) stop("a boh")
# if(any(colnames(class.2DG.16S)!=rownames(taxa.class.2DG.16S))) stop("a buh")
#
#
# usethis::use_data(class.2DG.16S     , overwrite = TRUE)
# usethis::use_data(taxa.class.2DG.16S, overwrite = TRUE)
# rm(class.2DG.16S,taxa.class.2DG.16S,i,idx,different.taxa,taxa.level)
#
#
# #PHYLUM
# #------------------------------------------------------------------------------#
# taxa.level <- "phylum"
# different.taxa <- unique(get(taxa.level,taxa.otu.2DG.16S))
# phylum.2DG.16S <- data.frame(matrix(NA, nrow=nrow(otu.2DG.16S), ncol=length(different.taxa),
#                                      dimnames=list(sample.id,different.taxa)))
#
# for(i in 1:length(different.taxa)){
#   idx <- which(different.taxa[i] == get(taxa.level,taxa.otu.2DG.16S))
#   phylum.2DG.16S[,i] <- apply(X=otu.2DG.16S, MARGIN=1, function(x) sum(x[idx]) )
# }
#
# taxa.phylum.2DG.16S <- unique(subset(taxa.otu.2DG.16S, select=c(-otu,-genus,-family,-order,-class)))
# rownames(taxa.phylum.2DG.16S) <- taxa.phylum.2DG.16S$phylum
#
# if(any(rownames(phylum.2DG.16S)!=rownames(otu.2DG.16S))) stop("a beh")
# if(!any(colnames(phylum.2DG.16S) %in% get(taxa.level,taxa.otu.2DG.16S))) stop("a boh")
# if(any(colnames(phylum.2DG.16S)!=rownames(taxa.phylum.2DG.16S))) stop("a buh")
#
#
# usethis::use_data(phylum.2DG.16S     , overwrite = TRUE)
# usethis::use_data(taxa.phylum.2DG.16S, overwrite = TRUE)
# rm(phylum.2DG.16S,taxa.phylum.2DG.16S,i,idx,different.taxa,taxa.level)
#
#
#
rm(list=ls())

