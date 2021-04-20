#This script does all data wrangling needed to have TPM table of all DGE genes
#Gets data ready to build heatmaps
#Buils heatmap of all DEG genes across the 6 samples
#Build heatmap of the average TPMs of both males and females

library("gplots")
library("pheatmap")
library("RColorBrewer")

#importing data
fc_MyrioTPM <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Myrio_featureMatrixTPM.txt", header = TRUE)

#genes as rows
rownames(fc_MyrioTPM) <- fc_MyrioTPM[,1]
fc_MyrioTPM[,1] <- NULL

#columns names a bit more handy
names(fc_MyrioTPM)[2] <- "Female_1"
names(fc_MyrioTPM)[3] <- "Female_2"
names(fc_MyrioTPM)[4] <- "Female_3"
names(fc_MyrioTPM)[5] <- "Male_1"
names(fc_MyrioTPM)[6] <- "Male_2"
names(fc_MyrioTPM)[7] <- "Male_3"

# Step 1 - Normalize for gene length
#Divide each gene count by gene legth (I am adding results as new variables in the sabe df(not very sof))
fc_MyrioTPM$F1<- fc_MyrioTPM[,"Female_1"]/fc_MyrioTPM[,"Length"]
fc_MyrioTPM$F2<- fc_MyrioTPM[,"Female_2"]/fc_MyrioTPM[,"Length"]
fc_MyrioTPM$F3<- fc_MyrioTPM[,"Female_3"]/fc_MyrioTPM[,"Length"]
fc_MyrioTPM$M1<- fc_MyrioTPM[,"Male_1"]/fc_MyrioTPM[,"Length"]
fc_MyrioTPM$M2<- fc_MyrioTPM[,"Male_2"]/fc_MyrioTPM[,"Length"]
fc_MyrioTPM$M3<- fc_MyrioTPM[,"Male_3"]/fc_MyrioTPM[,"Length"]

#put results in a new df 
fc_MyrioTPM_1 <- fc_MyrioTPM[, c(1, 8:13)]

#Step 2 - Normalize for sequence depth
# 2.1 - Add up read counts for each sample and divide by 1 million
F1_sc <- sum(fc_MyrioTPM_1$F1/1000000)
F2_sc <- sum(fc_MyrioTPM_1$F2/1000000)
F3_sc <- sum(fc_MyrioTPM_1$F3/1000000)
M1_sc <- sum(fc_MyrioTPM_1$M1/1000000)
M2_sc <- sum(fc_MyrioTPM_1$M2/1000000)
M3_sc <- sum(fc_MyrioTPM_1$M3/1000000)


# 2.2 - normalized gene counts / new scalling factor
fc_MyrioTPM_1$Female_1 <- fc_MyrioTPM_1[,"F1"]/F1_sc
fc_MyrioTPM_1$Female_2 <- fc_MyrioTPM_1[,"F2"]/F2_sc
fc_MyrioTPM_1$Female_3 <- fc_MyrioTPM_1[,"F3"]/F3_sc
fc_MyrioTPM_1$Male_1 <- fc_MyrioTPM_1[,"M1"]/M1_sc
fc_MyrioTPM_1$Male_2 <- fc_MyrioTPM_1[,"M2"]/M2_sc
fc_MyrioTPM_1$Male_3 <- fc_MyrioTPM_1[,"M2"]/M3_sc
#put final results in a new df
fc_MyrioTPM_2 <- fc_MyrioTPM_1[, c(1, 8:13)]


###########################################
### save fc_DictyoTPM_2 like it is here
###########################################

#this creates a new df with genes (rows) which must be present in two dif DFs 
#fc_DictyoTPM2 and filter_DF3
(rownames(fc_MyrioTPM_2) %in% rownames(filter_df3)) #logical matrix
fc_MyrioTPM_DEGs <- fc_MyrioTPM_2[(rownames(fc_MyrioTPM_2) %in% rownames(filter_df3)),] #apply logical filtering

#get rid of length (not needed for heatMap)
fc_MyrioTPM_DEGs <- subset (fc_MyrioTPM_DEGs, select = -Length)

#log2 +1
fc_MyrioTPM_DEGsLog <- log2(fc_MyrioTPM_DEGs+1)
#get matrix for heatMap
Myrio_heatMatrix <- as.matrix(fc_MyrioTPM_DEGsLog)


#heatMap 1 (2)
heatmap.2(x=Myrio_heatMatrix, col = brewer.pal(11,"Spectral"), keysize = 1.5, margins = c(7,3), trace = "none",labRow=FALSE)

#data wrangling for *binary heatmap*
fc_MyrioTPM_DEGsLog$Females <- rowMeans(fc_MyrioTPM_DEGsLog[,c("Female_1", "Female_2","Female_3")])
fc_MyrioTPM_DEGsLog$Males <- rowMeans(fc_MyrioTPM_DEGsLog[,c("Male_1", "Male_2","Male_3")])
#new df
binary_heatMap <- fc_MyrioTPM_DEGsLog[,c(7,8)]
# make matrix for heatMap
binary_heatMapMatrix <- as.matrix(binary_heatMap)

heatmap.2(x=binary_heatMapMatrix, col = brewer.pal(11,"Spectral"),trace = "none", labRow=FALSE, cexCol = 1 )


logicalTPM <- binary_heatMap$Females<1 & binary_heatMap$Males<1
length(logicalTPM[logicalTPM== TRUE])


################
##Getting rid of genes where TPM <1 in both males and females
#create not-in function
`%nin%` = Negate(`%in%`)
#new df based on condition
new_TPM <- subset(binary_heatMap, binary_heatMap$Females <1 & binary_heatMap$Males <1)
#get logical
(rownames(binary_heatMap) %nin% rownames(new_TPM)) 
#applying logical filtering
binary_heatMap2 <- binary_heatMap[(rownames(binary_heatMap) %nin% rownames(new_TPM)),] 
#as matrix 
binary_heatMapMatrix2 <- as.matrix(binary_heatMap2)

pheatmap(binary_heatMapMatrix2, show_rownames = F)
clr <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(100)  
pheatmap(binary_heatMapMatrix2, show_rownames = F, color = clr, cutree_cols = 2)

####################
#FC >4
####################

#logical matrix, genes which are both present in Female_up and filter_DF1
logical_vector <-(rownames(Male_Up) %in% (fold_4_genes))
#count number of TRUES which (I think, hope, presume) = number of genes FBG FC >2
sum(logical_vector, na.rm = TRUE) # best way to count TRUE values

####################
#FC >6
####################
#logical matrix, genes which are both present in Female_up and filter_DF1
logical_vector <-(rownames(Male_Up) %in% (fold_6_genes))
#count number of TRUES which (I think, hope, presume) = number of genes FBG FC >2
sum(logical_vector, na.rm = TRUE) # best way to count TRUE values

#Sanity check
sanity_check <- row.names(filter_df3[filter_df3$padj<.05 & abs(filter_df3$log2FoldChange) <=2,])
