
bacteria.taxa <- taxa.bacteria.species.2DG.WGS
viruses.taxa <- taxa.viruses.species.2DG.WGS
archaea.taxa <- taxa.archaea.species.2DG.WGS
fungi.taxa <- taxa.fungi.species.2DG.WGS

taxa.2DG.WGS <- rbind(bacteria.taxa, viruses.taxa, archaea.taxa, fungi.taxa)

bacteria.count <- bacteria.species.count.2DG.WGS
viruses.count <- viruses.species.count.2DG.WGS
archaea.count <- archaea.species.count.2DG.WGS
fungi.count <- fungi.species.count.2DG.WGS

count.2DG.WGS <- cbind(bacteria.count, viruses.count, archaea.count, fungi.count)


bacteria.score <- bacteria.species.score.2DG.WGS
viruses.score <- viruses.species.score.2DG.WGS
archaea.score <- archaea.species.score.2DG.WGS
fungi.score <- fungi.species.score.2DG.WGS

score.2DG.WGS <- cbind(bacteria.score, viruses.score, archaea.score, fungi.score)


all(rownames(count.2DG.WGS)==rownames(score.2DG.WGS))
all(colnames(count.2DG.WGS)==colnames(score.2DG.WGS))
all(colnames(count.2DG.WGS)==rownames(taxa.2DG.WGS))
all(rownames(count.2DG.WGS)==rownames(meta.2DG))

usethis::use_data(count.2DG.WGS , overwrite = TRUE)
usethis::use_data(score.2DG.WGS , overwrite = TRUE)
usethis::use_data(taxa.2DG.WGS , overwrite = TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data.dir <- "../data/"

delfiles <- dir(path=data.dir ,pattern="archaea*")
file.remove(file.path(data.dir, delfiles))

delfiles <- dir(path=data.dir ,pattern="bacteria*")
file.remove(file.path(data.dir, delfiles))

delfiles <- dir(path=data.dir ,pattern="fungi*")
file.remove(file.path(data.dir, delfiles))


delfiles <- dir(path=data.dir ,pattern="viruses*")
file.remove(file.path(data.dir, delfiles))


delfiles <- file.remove("../data/samples.removed.2DG.WGS.rda")

