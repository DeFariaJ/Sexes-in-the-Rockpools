library(stringr)
#importing data
SCOs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/OrthologsIDS_test.txt", header = FALSE)
#SCOs as rows
rownames(SCOs) <- SCOs[,1]
SCOs[,1] <- NULL

Dictyo_SBG_OGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Dictyo_SB_SCOs", header = FALSE)
Myrio_SBG_OGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Myrio_SB_SCOs", header = FALSE)
Desma_SBG_OGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Desma_SB_SCOs2", header = FALSE)
Saccho_SBG_OGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Saccho_SB_SCOs2", header = FALSE)

Dictyo_vector <- as.vector(Dictyo_SBG_OGs$V1)
Myrio_vector <- as.vector(Myrio_SBG_OGs$V1)
Desma_vector <- as.vector(Desma_SBG_OGs$V1)
Saccho_vector <- as.vector(Saccho_SBG_OGs$V1)

#dictyo upset
Dictyo_upset <- as.vector(SCOs$V2%in%(Dictyo_vector))
SCOs$Dictyota <- Dictyo_upset
SCOs$Dictyota <- as.numeric(SCOs$Dictyota)
sum(SCOs$Dictyo)

#myrio_upset
Myrio_upset <- as.vector(SCOs$V3%in%(Myrio_vector))
SCOs$Myriotrichia <- Myrio_upset
SCOs$Myriotrichia <- as.numeric(SCOs$Myrio)
sum(SCOs$Myriotrichia)

#desma upset
Desma_upset <- as.vector(SCOs$V4%in%(Desma_vector))
SCOs$Desmarestia <- Desma_upset
SCOs$Desmarestia <- as.numeric(SCOs$Desmarestia)
sum(SCOs$Desmarestia)

#saccho upset
Saccho_upset <- as.vector(SCOs$V5%in%(Saccho_vector))
SCOs$Sacchoriza <- Saccho_upset
SCOs$Sacchoriza <- as.numeric(SCOs$Sacchoriza)
sum(SCOs$Saccho)

SCOs_binary <- SCOs[, c(5:8)]

upset(SCOs_binary,nsets=4)


#stuff
SCOs_4species <- SCOs[(SCOs$V2 %in% Dictyo_vector)&(SCOs$V3 %in% Myrio_vector)&(SCOs$V4 %in% Desma_vector)&(SCOs$V5 %in% Saccho_vector),] 
SCOs_3species_NoDictyo <- SCOs[(SCOs$V3 %in% Myrio_vector)&(SCOs$V4 %in% Desma_vector)&(SCOs$V5 %in% Saccho_vector),]
SCOs_3species_NoMyrio <- SCOs[(SCOs$V2 %in% Dictyo_vector)&(SCOs$V4 %in% Desma_vector)&(SCOs$V5 %in% Saccho_vector),] 
SCOs_3species_NoDesma <- SCOs[(SCOs$V2 %in% Dictyo_vector)&(SCOs$V3 %in% Myrio_vector)&(SCOs$V5 %in% Saccho_vector),] 
SCOs_3species_NoSaccho <- SCOs[(SCOs$V2 %in% Dictyo_vector)&(SCOs$V3 %in% Myrio_vector)&(SCOs$V4 %in% Desma_vector),] 

#male or female? **Dictyota**
conserved_Female <- fksake_Female[row.names(fksake_Female) %in% (SCOs_3species_NoDictyo$V2),]
conserved_Male <- fksake_Male[row.names(fksake_Male) %in% (SCOs_3species_NoDictyo$V2),]

#male or female? **Myriotrichia**
conserved_SBGs <- filter_df3[row.names(filter_df3)%in%(SCOs_3species_NoDictyo$V3),]
conserved_Female <- fksake_Female[row.names(fksake_Female)%in%(row.names(conserved_SBGs)),]
conserved_Male <- fksake_Male[row.names(fksake_Male)%in%(row.names(conserved_SBGs)),]

#male or female? **Desma**
desma_sorted <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Desma_4speciesSorted", header = FALSE)
#SCOs as rows
conserved_SBGs <- filter_df3[row.names(filter_df3)%in%(SCOs_3species_NoDictyo$V4),]
conserved_Female <- fksake_Female[row.names(fksake_Female)%in%(row.names(conserved_SBGs)),]
conserved_Male <- fksake_Male[row.names(fksake_Male)%in%(row.names(conserved_SBGs)),]

SCOs_3species_NoDictyo$V4 <- str_replace_all(SCOs_3species_NoDictyo$V4, "prot", "mRNA")


#male or female? **Saccho**
saccho_sorted <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Saccho_4speciesSorted", header = FALSE)
SCOs_3species_NoDictyo$V5 <- str_replace_all(SCOs_3species_NoDictyo$V5, "prot", "mRNA")
#SCOs as rows
conserved_SBGs <- filter_df3[row.names(filter_df3)%in%(SCOs_3species_NoDictyo$V5),]
conserved_Female <- fksake_Female[row.names(fksake_Female)%in%(row.names(conserved_SBGs)),]
conserved_Male <- fksake_Male[row.names(fksake_Male)%in%(row.names(conserved_SBGs)),]

###########################################################################################################################################
#FUNCTION

#functio of these genes Dictio
Dictyo_conserved_F <- ego_prep[(ego_prep$V1 %in% SCOs_4species$V2),]
Dictyo_conserved_FGO <- as.vector(Dictyo_conserved_F[3,])
(allResM$GO.ID)%in%(Dictyo_conserved_FGO)

#functio of these genes Myrio
Myrio_conserved_F <- ego_prep[(ego_prep$V1 %in% SCOs_4species$V3),]
Myrio_conserved_FGO <- as.vector(Myrio_conserved_F[3,])
(allResM$GO.ID)%in%(Myrio_conserved_FGO)


#functio of these genes Desma
SCOs_4species$V4 <- str_replace_all(SCOs_4species$V4, "prot", "mRNA")
Desma_conserved_F <- ego_prep[(ego_prep$V1 %in% SCOs_4species$V4),]
Desma_conserved_FGO <- as.vector(Desma_conserved_F[3,])
(allResM$GO.ID)%in%(Desma_conserved_FGO)

#functio of these genes Saccho
SCOs_4species$V5 <- str_replace_all(SCOs_4species$V5, "prot", "mRNA")
Saccho_conserved_F <- ego_prep[(ego_prep$V1 %in% SCOs_4species$V5),]
Desma_conserved_FGO <- as.vector(Desma_conserved_F[3,])
(allResM$GO.ID)%in%(Desma_conserved_FGO)
