library(MetaGenomicR)
library(stringr)

data <- bacteria.species.count.2DG.WGS
taxa <- taxa.bacteria.species.2DG.WGS.raw[colnames(data),]

idx <- which(colSums(data)!=0)

data <- data[,idx]
taxa <- taxa.bacteria.species.2DG.WGS.raw[idx,]

nlvls <- sapply(taxa$taxonomy, function(x) str_split(x,'\\|'))
nlvls <- sapply(nlvls, length)
hist(nlvls, main="0")
#-----------------#

#-----------------#
all(str_detect(taxa$taxonomy,'root\\|cellular_organisms\\|'))
taxa$taxonomy <- str_remove(taxa$taxonomy,'root\\|cellular_organisms\\|')
#-----------------#

#-----------------#
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Clostridiales_Family_XIII._Incertae_Sedis')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Roseiflexineae')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Rickettsieae')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Bacillales_Family_XII._Incertae_Sedis')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|delta/epsilon_subdivisions')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Cystobacterineae')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Bacillales_Family_XI._Incertae_Sedis')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Agrobacterium_tumefaciens_complex')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|unclassified_Ruminococcaceae_\\(miscellaneous\\)')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Mycobacterium_avium_complex_\\(MAC\\)')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Clostridiales_Family_XVI._Incertae_Sedis')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|unclassified_Erysipelotrichaceae_\\(miscellaneous\\)')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Sphaerobacterineae')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Sphaerobacteridae')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Burkholderiales_Genera_incertae_sedis')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Clostridiales_Family_XVII._Incertae_Sedis')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Bacillus_altitudinis_complex')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Sorangiineae')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Mycobacterium_simiae_complex')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Citrobacter_freundii_complex')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Enterobacter_cloacae_complex')
taxa$taxonomy <- str_remove(taxa$taxonomy,'\\|Nannocystineae')
#-----------------#

#-----------------#
res <- list()
for(i in 1:length(taxa$taxonomy)){
	tmp <- unlist(str_split(taxa$taxonomy[i],"\\|"))
	tmp <- tmp[!str_detect(tmp,"group")]
	res[[i]] <- tmp
}


#-----------------#
taxonomy <- as.data.frame(matrix(NA,nrow=nrow(taxa), ncol=12))
rownames(taxonomy) <- rownames(taxa)
colnames(taxonomy)[1:7] <- c("domain","phylum","class","order","family","genus","species")
for(i in 1:length(taxa$taxonomy)){taxonomy[i,] <- res[[i]][1:12]}
#-----------------#


#Unclassified
#----------------------#
idx <- which(str_detect(taxonomy[,8],"sp."))
taxonomy[idx, 7] <- taxonomy[idx, 8]
taxonomy[idx, 8] <- NA

#Less than 7 step
#----------------------#
idx <- which(is.na(taxonomy[,5]))
View(taxonomy[idx,])
taxonomy[idx,7] <- taxonomy[idx,4]
taxonomy[idx,4:6] <- taxonomy[idx,3]



idx <- which(is.na(taxonomy[,6]))
View(taxonomy[idx,])

taxonomy["1457365",6:7] <- taxonomy["1457365",3:4]
taxonomy["1457365",3] -> taxonomy["1457365",4:5]

taxonomy["1541959",6] <- "Candidatus_Izimaplasma"
taxonomy["1541959",7] <- "Candidatus_Izimaplasma_sp._HR1"
taxonomy["1541959",3:5] <- "unclassified_Candidatus_Izimaplasma"

taxonomy["1637999",7] <- "Verrucomicrobia_bacterium_IMCC26134"
taxonomy["1637999",3:6] <- "unclassified_Verrucomicrobia"

taxonomy["1839801",2:7] <- c("Chloroflexi","Dehalococcoidia","Dehalococcoidales",
							 "Dehalococcoidia","Dehalogenimonas","Dehalogenimonas_formicexedens")

taxonomy["1871025",6:7] <- taxonomy["1871025",4:5]
taxonomy["1871025",3] -> taxonomy["1871025",4:5]

taxonomy["243899",7] <- taxonomy["243899",5]
taxonomy["243899",3:6] <- "unclassified_Firmicutes_sensu_stricto"

taxonomy["2488809",7] <- taxonomy["2488809",5]
taxonomy["2488809",3:6] <- "unclassified_Verrucomicrobia"

taxonomy["552810",2:7] <- c("Chloroflexi","Dehalococcoidia","Dehalococcoidales",
							"Dehalococcoidia","Dehalogenimonas","Dehalogenimonas_lykanthroporepellens")


idx <- which(is.na(taxonomy[,7]))
View(taxonomy[idx,])

taxonomy["1248727",1:7] <- c("Bacteria","Proteobacteria",
							 "Gammaproteobacteria",
							 "unclassified_Gammaproteobacteria",
							 "unclassified_Gammaproteobacteria",
							 "unclassified_Gammaproteobacteria",
							 "endosymbiont_of_unidentified_scaly_snail_isolate_Monju")

taxonomy["143393",1:7]  <- c("Bacteria","Firmicutes",
							 "Clostridia",
							 "Eubacteriales",
							 "Eubacteriales incertae sedis",
							 "Eubacteriales incertae sedis",
							 "[Eubacterium]_sulci")

taxonomy["146919",1:7]  <- c("Bacteria"," Bacteroidetes",
							 "Bacteroidia",
							 "Bacteroidetes_Order_II._Incertae_sedis",
							 "Rhodothermaceae",
							 "Salinibacter",
							 "Salinibacter_ruber")

taxonomy["1852374",1:7]  <- c("Bacteria"," Firmicutes",
							  "Tissierellia",
							  "Tissierellia_incertae_sedis",
							  "Tissierellia_incertae_sedis",
							  "Ezakiella",
							  "Ezakiella_massiliensis")


taxonomy["1855912",1:7]  <- c("Bacteria"," Acidobacteria",
							  "Vicinamibacteria",
							  "Bacteroidetes Order_II._Incertae_sedis",
							  "Vicinamibacteraceae",
							  "Luteitalea",
							  "Luteitalea_pratensis")

taxonomy["1868589",1:7]  <- c("Bacteria","Proteobacteria",
							  "Alphaproteobacteria",
							  "Hyphomicrobiales",
							  "Phreatobacteraceae",
							  "Phreatobacter",
							  "Phreatobacter_cathodiphilus")

taxonomy["1971488",1:7]  <- c("Bacteria","Proteobacteria",
							  "Gammaproteobacteria",
							  "Gammaproteobacteria_incertae_sedis",
							  "Gammaproteobacteria_incertae_sedis",
							  "Candidatus_Nardonella",
							  "endosymbiont_of_Pachyrhynchus_infernalis")

taxonomy["2006",1:7]  <- c("Bacteria","Actinobacteria",
						   "Actinomycetia",
						   "Streptosporangiales",
						   "Streptosporangiaceae",
						   "Thermobispora",
						   "Thermobispora_bispora")

taxonomy["2026787",1:7]  <- c("Bacteria","Bacteroidetes",
							  "Bacteroidetes_Order_II._Incertae_sedis",
							  "Bacteroidetes_Order_II._Incertae_sedis",
							  "Rhodothermaceae",
							  "unclassified_Rhodothermaceae",
							  "Rhodothermaceae_bacterium")

taxonomy["2066483",1:7]  <- c("Bacteria"," unclassified_Candidatus_Dependentiae",
							  "unclassified_Candidatus_Dependentiae",
							  "unclassified_Candidatus_Dependentiae",
							  "unclassified_Candidatus_Dependentiae",
							  "unclassified_Candidatus_Dependentiae",
							  "Candidatus_Dependentiae_bacterium_(ex_Spumella_elongata_CCAP_955/1)")

taxonomy["2292766",1:7]  <- c("Bacteria","Proteobacteria",
							  "Deltaproteobacteria",
							  "Bradymonadales",
							  "Bradymonadaceae",
							  "Persicimonas",
							  "Bradymonadales_bacterium_YN101")

taxonomy["2600177",1:7]  <- c("Bacteria","Proteobacteria",
							  "Deltaproteobacteria",
							  "Bradymonadales",
							  "unclassified_Bradymonadales",
							  "unclassified_Bradymonadales",
							  "Bradymonadales_bacterium_V1718")

taxonomy["29549",1:7]  <- c("Bacteria","Bacteroidetes",
							"Bacteroidetes_Order_II._Incertae_sedis",
							"Bacteroidetes_Order_II._Incertae_sedis",
							"Rhodothermaceae",
							"Rhodothermus",
							"Rhodothermus_marinus")

taxonomy["446679",1:7]  <- c("Bacteria","Cyanobacteria",
							 "Cyanophyceae",
							 "Nostocales",
							 "Nostocaceae",
							 "Nostoc",
							 "Nostoc_sphaeroides")

taxonomy["59930",1:7]  <- c("Bacteria","Cyanobacteria",
							"Cyanophyceae",
							"Synechococcales",
							"Synechococcaceae",
							"Cyanobium",
							"Cyanobium_gracile")

taxonomy["947516",1:7]  <- c("Bacteria","Proteobacteria",
							 "Gammaproteobacteria",
							 "unclassified_Gammaproteobacteria",
							 "unclassified_Gammaproteobacteria",
							 "unclassified_Gammaproteobacteria",
							 "gamma_proteobacterium_SS-5")


View(taxonomy[idx,])

!any(is.na(taxonomy[,1]))
!any(is.na(taxonomy[,2]))
!any(is.na(taxonomy[,3]))
!any(is.na(taxonomy[,4]))
!any(is.na(taxonomy[,5]))
!any(is.na(taxonomy[,6]))
!any(is.na(taxonomy[,7]))
all(is.na(taxonomy[,8]))
all(is.na(taxonomy[,9]))
all(is.na(taxonomy[,10]))
all(is.na(taxonomy[,11]))
all(is.na(taxonomy[,12]))



taxonomy <- taxonomy[,1:7]
taxa.bacteria.species.2DG.WGS <- taxonomy
View(taxonomy)


all(colnames(data)==rownames(taxonomy))
usethis::use_data(taxa.bacteria.species.2DG.WGS , overwrite = TRUE)












#Viruses
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
taxa <- taxa.viruses.species.2DG.WGS.raw

#-----------------#
all(str_detect(taxa$taxonomy,'root\\|'))
taxa$taxonomy <- str_remove(taxa$taxonomy,'root\\|')
#-----------------#



#-----------------#
res <- list()
for(i in 1:length(taxa$taxonomy)){
	tmp <- unlist(str_split(taxa$taxonomy[i],"\\|"))
	tmp <- tmp[!str_detect(tmp,"group")]
	res[[i]] <- tmp
}


#-----------------#
taxonomy <- as.data.frame(matrix(NA,nrow=nrow(taxa), ncol=7))
rownames(taxonomy) <- rownames(taxa)
colnames(taxonomy)[1:7] <- c("domain","phylum","class","order","family","genus","species")

taxonomy[1,] <- c("Viruses","Artverviricota","Revtraviricetes","Ortervirales",
				  "Retroviridae","Betaretrovirus","Mouse_mammary_tumor_virus")
taxonomy[2,] <- c("Viruses","Artverviricota","Revtraviricetes","Ortervirales",
				  "Retroviridae", "Gammaretrovirus","Murine_leukemia_virus")
taxonomy[3,] <- c("Viruses","Artverviricota","Revtraviricetes","Ortervirales",
				  "Retroviridae", "Gammaretrovirus","Abelson_murine_leukemia_virus")
taxonomy[4,] <- c("Viruses","Artverviricota","Revtraviricetes","Ortervirales",
				  "Retroviridae", "Gammaretrovirus","Moloney_murine_sarcoma_virus")
taxonomy[5,] <- c("Viruses","Artverviricota","Revtraviricetes","Ortervirales",
				  "Retroviridae", "Gammaretrovirus","Spleen_focus-forming_virus")
taxonomy[6,] <- c("Viruses","Uroviricota","Caudoviricetes","Caudovirales",
				  "Myoviridae", "unclassified_Myoviridae","Lactobacillus_prophage_Lj771")
taxonomy[7,] <- c("Viruses","Uroviricota","Caudoviricetes","Caudovirales",
				  "Myoviridae", "unclassified_Myoviridae","Lactobacillus_phage_phi_jlb1")
taxonomy[8,] <- c("Viruses","Artverviricota","Revtraviricetes","Ortervirales",
				  "Retroviridae", "Gammaretrovirus","Finkel-Biskis-Jinkins_murine_sarcoma_virus")
taxonomy[9,] <- c("Viruses","Artverviricota","Revtraviricetes","Ortervirales",
				  "Retroviridae", "Gammaretrovirus","Mus_musculus_mobilized_endogenous_polytropic_provirus")
taxonomy[10,]<- c("Viruses","Artverviricota","Revtraviricetes","Ortervirales",
				  "Retroviridae", "Gammaretrovirus","Murine_leukemia-related_retroviruses")

taxa.viruses.species.2DG.WGS <- taxonomy
usethis::use_data(taxa.viruses.species.2DG.WGS , overwrite = TRUE)







#Viruses
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
taxa <- taxa.archaea.species.2DG.WGS.raw

#-----------------#
all(str_detect(taxa$taxonomy,'root\\|'))
taxa$taxonomy <- str_remove(taxa$taxonomy,'root\\|')
#-----------------#



#-----------------#
res <- list()
for(i in 1:length(taxa$taxonomy)){
	tmp <- unlist(str_split(taxa$taxonomy[i],"\\|"))
	tmp <- tmp[!str_detect(tmp,"group")]
	res[[i]] <- tmp
}


#-----------------#
taxonomy <- as.data.frame(matrix(NA,nrow=nrow(taxa), ncol=7))
rownames(taxonomy) <- rownames(taxa)
colnames(taxonomy)[1:7] <- c("domain","phylum","class","order","family","genus","species")

taxonomy[1 ,] <- res[[1 ]][2:8]
taxonomy[2 ,] <- res[[2 ]][2:8]
taxonomy[3 ,] <- res[[3 ]][2:8]
taxonomy[4 ,] <- res[[4 ]][2:8]
taxonomy[5 ,] <- res[[5 ]][2:8]
taxonomy[6 ,] <- res[[6 ]][c(2,3,4,5,6,7,9)]
taxonomy[7 ,] <- res[[7 ]][2:8]
taxonomy[8 ,] <- res[[8 ]][2:8]
taxonomy[9 ,] <- res[[9 ]][2:8]
taxonomy[10,] <- res[[10]][2:8]
taxonomy[11,] <- res[[11]][2:8]
taxonomy[12,] <- res[[12]][c(2,3,4,5,6,7,9)]
taxonomy[13,] <- res[[13]][c(2,3,4,5,6,7,9)]
taxonomy[14,] <- res[[14]][2:8]
taxonomy[15,] <- res[[15]][2:8]
taxonomy[16,] <- res[[16]][2:8]
taxonomy[17,] <- res[[17]][2:8]
taxonomy[18,] <- res[[18]][2:8]
taxonomy[19,] <- res[[19]][2:8]
taxonomy[20,] <- res[[20]][2:8]


taxa.archaea.species.2DG.WGS <- taxonomy
usethis::use_data(taxa.archaea.species.2DG.WGS , overwrite = TRUE)












#Fungi
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
taxa <- taxa.fungi.species.2DG.WGS.raw

#-----------------#
all(str_detect(taxa$taxonomy,'root\\|'))
taxa$taxonomy <- str_remove(taxa$taxonomy,'root\\|')
#-----------------#
all(str_detect(taxa$taxonomy,'cellular_organisms\\|'))
taxa$taxonomy <- str_remove(taxa$taxonomy,'cellular_organisms\\|')
#-----------------#
all(str_detect(taxa$taxonomy,'Eukaryota\\|'))
taxa$taxonomy <- str_remove(taxa$taxonomy,'Eukaryota\\|')
#-----------------#
all(str_detect(taxa$taxonomy,'Opisthokonta\\|'))
taxa$taxonomy <- str_remove(taxa$taxonomy,'Opisthokonta\\|')
#-----------------#
all(str_detect(taxa$taxonomy,'Dikarya\\|'))
taxa$taxonomy <- str_remove(taxa$taxonomy,'Dikarya\\|')
taxa$taxonomy[25] <- str_remove(taxa$taxonomy[25],'Fungi_incertae_sedis\\|')
#-----------------#
taxa$taxonomy <- str_remove(taxa$taxonomy,'Saccharomycotina\\|')
taxa$taxonomy <- str_remove(taxa$taxonomy,'saccharomyceta\\|')
taxa$taxonomy <- str_remove(taxa$taxonomy,"Pleosporineae\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Pleosporomycetidae\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"dothideomyceta\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"leotiomyceta\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"leotiomyceta\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Pezizomycotina\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Eurotiomycetidae\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Alternaria_sect._Alternaria\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Wallemiomycotina\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Chaetothyriomycetidae\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Dothideomycetes_incertae_sedis\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Sordariomycetidae\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"sordariomyceta\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Pucciniomycotina\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Xylariomycetidae\\|")
taxa$taxonomy[19] <- str_remove(taxa$taxonomy[19],"Aspergillus\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Nakaseomyces\\/Candida_clade\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Fusarium_sambucinum_species_complex\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Hypocreomycetidae\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Mucoromycotina\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Dothideomycetidae\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Alternaria_alternata_complex\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Fumigati\\|")
taxa$taxonomy <- str_remove(taxa$taxonomy,"Ustilaginomycotina\\|")
#-----------------#


res <- list()
for(i in 1:length(taxa$taxonomy)){
	tmp <- unlist(str_split(taxa$taxonomy[i],"\\|"))
	tmp <- tmp[!str_detect(tmp,"group")]
	res[[i]] <- tmp
}
rm(i,tmp)
which(sapply(res,length)>7)




#-----------------#
taxonomy <- as.data.frame(matrix(NA,nrow=nrow(taxa), ncol=7))
rownames(taxonomy) <- rownames(taxa)
colnames(taxonomy) <- c("domain","phylum","class","order","family","genus","species")

for(i in 1:nrow(taxa)){
	taxonomy[i,] <- res[[i]]
}






taxa.fungi.species.2DG.WGS <- taxonomy
usethis::use_data(taxa.fungi.species.2DG.WGS , overwrite = TRUE)



