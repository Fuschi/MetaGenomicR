library(stringr)

path <- system.file("extdata", "ben1_16s_otu_table.tsv.tsv" ,package="MetaGenomicR", mustWork=TRUE)
otu_POMP_DO <- read.table(path, header=T, row.names=1, sep="\t")
otu_POMP_DO <- as.data.frame(t(otu_POMP_DO))

path <- system.file("extdata", "ben1_16s_otu_taxa.tsv.tsv" ,package="MetaGenomicR", mustWork=TRUE)
taxa_POMP_DO <- read.table(path, header=T, row.names=1, sep="\t")
taxa_POMP_DO <- as.data.frame(apply(taxa_POMP_DO,c(1,2),function(x)str_replace_all(x,'"','')))
taxa_POMP_DO <- as.data.frame(apply(taxa_POMP_DO,c(1,2),function(x)str_replace_all(x,' ','_')))
taxa_POMP_DO <- as.data.frame(apply(taxa_POMP_DO,c(1,2),function(x)str_replace_all(x,'/','_')))
taxa_POMP_DO$otu <- rownames(taxa_POMP_DO)
taxa_POMP_DO <- taxa_POMP_DO[colnames(otu_POMP_DO),]

sample_name <- rownames(otu_POMP_DO)
meta_POMP_DO <- str_split(sample_name,"[.]", simplify=T)
meta_POMP_DO <- as.data.frame(meta_POMP_DO[,2:4])
colnames(meta_POMP_DO) <- c("diet","age","R")
rownames(meta_POMP_DO) <- rownames(otu_POMP_DO)

stopifnot(is.MetaGenomic(new('MetaGenomic',data=otu_POMP_DO,meta=meta_POMP_DO,taxa=taxa_POMP_DO)))
usethis::use_data(otu_POMP_DO, overwrite = TRUE)
usethis::use_data(meta_POMP_DO, overwrite = TRUE)
usethis::use_data(taxa_POMP_DO, overwrite = TRUE)
