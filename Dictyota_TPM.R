#This script does all data wrangling needed to have TPM table of all DGE genes
#Gets data ready to build heatmaps
#Buils heatmap of all DEG genes across the 6 samples
#Build heatmap of the average TPMs of both males and females

library("gplots")
library("pheatmap")
library("RColorBrewer")
#importing data
fc_DictyoTPM <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Dictyo_featureMatrixTPM.txt", header = TRUE)

#genes as rows
rownames(fc_DictyoTPM) <- fc_DictyoTPM[,1]
fc_DictyoTPM[,1] <- NULL

#columns names a bit more handy
names(fc_DictyoTPM)[2] <- "Female_1"
names(fc_DictyoTPM)[3] <- "Female_2"
names(fc_DictyoTPM)[4] <- "Female_3"
names(fc_DictyoTPM)[5] <- "Male_1"
names(fc_DictyoTPM)[6] <- "Male_2"
names(fc_DictyoTPM)[7] <- "Male_3"

# Step 1 - Normalize for gene length
#Divide each gene count by gene legth (I am adding results as new variables in the sabe df(not very sof))
fc_DictyoTPM$F1<- fc_DictyoTPM[,"Female_1"]/fc_DictyoTPM[,"Length"]
fc_DictyoTPM$F2<- fc_DictyoTPM[,"Female_2"]/fc_DictyoTPM[,"Length"]
fc_DictyoTPM$F3<- fc_DictyoTPM[,"Female_3"]/fc_DictyoTPM[,"Length"]
fc_DictyoTPM$M1<- fc_DictyoTPM[,"Male_1"]/fc_DictyoTPM[,"Length"]
fc_DictyoTPM$M2<- fc_DictyoTPM[,"Male_2"]/fc_DictyoTPM[,"Length"]
fc_DictyoTPM$M3<- fc_DictyoTPM[,"Male_3"]/fc_DictyoTPM[,"Length"]
#put results in a new df 
fc_DictyoTPM_1 <- fc_DictyoTPM[, c(1, 8:13)]

#Step 2 - Normalize for sequence depth
# 2.1 - Add up read counts for each sample and divide by 1 million
F1_sc <- sum(fc_DictyoTPM_1$F1/1000000)
F2_sc <- sum(fc_DictyoTPM_1$F2/1000000)
F3_sc <- sum(fc_DictyoTPM_1$F3/1000000)
M1_sc <- sum(fc_DictyoTPM_1$M1/1000000)
M2_sc <- sum(fc_DictyoTPM_1$M2/1000000)
M3_sc <- sum(fc_DictyoTPM_1$M3/1000000)

# 2.2 - normalized gene counts / new scalling factor
fc_DictyoTPM_1$Female_1 <- fc_DictyoTPM_1[,"F1"]/F1_sc
fc_DictyoTPM_1$Female_2 <- fc_DictyoTPM_1[,"F2"]/F2_sc
fc_DictyoTPM_1$Female_3 <- fc_DictyoTPM_1[,"F3"]/F3_sc
fc_DictyoTPM_1$Male_1 <- fc_DictyoTPM_1[,"M1"]/M1_sc
fc_DictyoTPM_1$Male_2 <- fc_DictyoTPM_1[,"M2"]/M2_sc
fc_DictyoTPM_1$Male_3 <- fc_DictyoTPM_1[,"M2"]/M3_sc
#put final results in a new df
fc_DictyoTPM_2 <- fc_DictyoTPM_1[, c(1, 8:13)]

###########################################
### save fc_DictyoTPM_2 like it is here
###########################################


#this creates a new df with genes (rows) which must be present in two dif DFs 
#fc_DictyoTPM2 and filter_DF3
(rownames(fc_DictyoTPM_2) %in% rownames(filter_df3)) #logical matrix
fc_DictyoTPM_DEGs <- fc_DictyoTPM_2[(rownames(fc_DictyoTPM_2) %in% rownames(filter_df3)),] #apply logical filtering

#get rid of length (not needed for heatMap)
fc_DictyoTPM_DEGs <- subset (fc_DictyoTPM_DEGs, select = -Length)

#log2 +1
fc_DictyoTPM_DEGsLog <- log2(fc_DictyoTPM_DEGs+1)
#get matrix for heatMap
dictyo_heatMatrix <- as.matrix(fc_DictyoTPM_DEGsLog)

#heatMap 1 (2)
heatmap.2(x=dictyo_heatMatrix, col = brewer.pal(11,"Spectral"), keysize = 1.5, margins = c(7,3), trace = "none", labRow=FALSE)

#data wrangling for *binary heatmap*
fc_DictyoTPM_DEGsLog$Females <- rowMeans(fc_DictyoTPM_DEGsLog[,c("Female_1", "Female_2","Female_3")])
fc_DictyoTPM_DEGsLog$Males <- rowMeans(fc_DictyoTPM_DEGsLog[,c("Male_1", "Male_2","Male_3")])
#new df
binary_heatMap <- fc_DictyoTPM_DEGsLog[,c(7,8)]
# make matrix for heatMap
binary_heatMapMatrix <- as.matrix(binary_heatMap)

heatmap.2(x=binary_heatMapMatrix, col = brewer.pal(12,"Spectral"),trace = "none", labRow=FALSE, cexCol = 1 )


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



######### final heat map Male vs Female averages #############
heatmap.2(x=binary_heatMapMatrix2, col = brewer.pal(9,"RdBu"),trace = "none", labRow=FALSE, cexCol = 1, breaks = c(0,1,1.5,2,2.5,3,3.5,4,4.5,5.5))

#pheatmap
#pheatmap(binary_heatMapMatrix2, show_rownames = F)
clr <-colorRampPalette(brewer.pal(7, "Greys"))(100)  
pheatmap(binary_heatMapMatrix2, show_rownames = F, color = clr, cutree_cols = 2)

#xxx colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(100) #takes middle 7 and extrapulate it into a much larger pallet
#assign it to a variable and pass it to the pheatmapo function
