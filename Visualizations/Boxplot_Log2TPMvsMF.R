library(tidyverse)
library(ggplot2)
library(reshape2)
#importing data
fc_Dictyo <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/DictyoOutput.Rmatrix.tab", header = TRUE)
Dictyo_FBGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/FBGs_Dictyota.tsv", header = TRUE)
Dictyo_MBGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/MBGs_Dictyota.tsv", header = TRUE)
Myrio_FBGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/FBGs_Myrio.tsv", header = TRUE)
Myrio_MBGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/MBGs_Myrio.tsv", header = TRUE)
Desma_FBGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/FBGs_Desma.tsv", header = TRUE)
Desma_MBGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/MBGs_Desma.tsv", header = TRUE)
Saccho_FBGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/FBGs_Saccho.tsv", header = TRUE)
Saccho_MBGs <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/MBGs_Saccho.tsv", header = TRUE)


#genes as rows
rownames(fc_Dictyo) <- fc_Dictyo[,1]
fc_Dictyo[,1] <- NULL

#columns names a bit more handy
names(fc_Dictyo)[1] <- "Female_1"
names(fc_Dictyo)[2] <- "Female_2"
names(fc_Dictyo)[3] <- "Female_3"
names(fc_Dictyo)[4] <- "Male_1"
names(fc_Dictyo)[5] <- "Male_2"
names(fc_Dictyo)[6] <- "Male_3"

#create metadata (data about columns in count data)
metaData_df <-  data.frame(first_column  = c("Female_1", "Female_2", "Female_3", "Male_1", "Male_2", "Male_3"),
                           Sex = c("Female", "Female", "Female", "Male", "Male", "Male")
)
rownames(metaData_df) <- metaData_df$first_column
metaData_df[,1] <- NULL

#getting main DESeq2 object
deseq_dataset = DESeqDataSetFromMatrix(countData = fc_Dictyo, colData = metaData_df, design = ~ Sex)

#prep steps for deseq2
deseq_dataset = estimateSizeFactors(deseq_dataset)
vst_data = varianceStabilizingTransformation(deseq_dataset) #transforms data to make it normally dist

###########**************************************#################
###########
#estimate dispersions (step2) 
deseq_dataset <- estimateDispersions(deseq_dataset) #applies bayesian shrinkage etc
#visualize results
plotDispEsts(deseq_dataset)

#black dots are the gene wise estimates of dispersion 
# red line is the fitted relationship between the mean of the counts and the mean of the dispersion
# blue dots are the result of the bayesian shrinkage 
# blue circles around black dots are outliers. DESeq2 does not shrink them

# 3) apply statistics
deseq_dataset <- nbinomWaldTest(deseq_dataset) #stats being applied to ~20000 genes takes a bit

#DESeq2 is done. Next time we can take the shortcut which does these 3 steps for us and put us here straight away (DESeq function)

#check results
result_table = results(deseq_dataset)
summary(result_table)
#putting results table into a nice data.frame
results_df <- as.data.frame(result_table)

#######FILTERING################

#testing filtering with Aga's script
b = subset(results_df, padj<.05 & abs(log2FoldChange)>1)#select the desired p value and FC


## Filter rows (genes) with all columns == NA
complete.cases(results_df) #complete.cases returns TRUE for columns data have all data (no NAs)
filter_df1 <- results_df[complete.cases(results_df),]

#pvclust was here

# adj p-value <0.05
filter_df1$padj <0.05
filter_df2 <- filter_df1[filter_df1$padj <0.05,]

#p-value <0.05 is not enough, we also want to have a large effect size
#log2FoldChange >= 1 <= -1, we will filter it on magnitude and not sign, regardless of direction
abs(filter_df2$log2FoldChange) > 1
filter_df3 <- filter_df2[abs(filter_df2$log2FoldChange) > 1,]


###                                                                 ###
####################### VISULALISATION ################################

#MA plot
plotMA(result_table, alpha = 0.05) #takes as input original table produced by DESeq2 (not data.frame) #CAN CHANGE TO 0.05
# the aim here is to produce a comparison of the mean normalized counts (baseMean) and fold change


######## VOLCANO PLOT ############
#for all the genes we have analysied, we need the fold change and the adj p-value

#creating a new column on filter_df1 (up, neutra, down)
filter_df1$DiffExpressed <-"Unbiased"
filter_df1$DiffExpressed[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >1] <- "FC >2"
filter_df1$DiffExpressed[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >2] <- "FC >4"
filter_df1$DiffExpressed[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >4] <- "FC >16"                                                


#Plot plot
ggplot(filter_df1, aes(x=log2FoldChange, y= -log10(padj))) +
  geom_point(aes(colour = DiffExpressed)) +
  geom_vline(xintercept = 1, linetype = 3) + 
  geom_vline(xintercept = -1, linetype = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 3) +
  ylim(-5, 110)

######## ************* TPM *************** ##############

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

#creating new variable LengthKB
fc_DictyoTPM$LengthKB <- fc_DictyoTPM$Length/1000

# Step 1 - Normalize for gene length
#Divide each gene count by gene legth (I am adding results as new variables in the sabe df(not very sof))
fc_DictyoTPM$F1<- fc_DictyoTPM[,"Female_1"]/fc_DictyoTPM[,"LengthKB"]
fc_DictyoTPM$F2<- fc_DictyoTPM[,"Female_2"]/fc_DictyoTPM[,"LengthKB"]
fc_DictyoTPM$F3<- fc_DictyoTPM[,"Female_3"]/fc_DictyoTPM[,"LengthKB"]
fc_DictyoTPM$M1<- fc_DictyoTPM[,"Male_1"]/fc_DictyoTPM[,"LengthKB"]
fc_DictyoTPM$M2<- fc_DictyoTPM[,"Male_2"]/fc_DictyoTPM[,"LengthKB"]
fc_DictyoTPM$M3<- fc_DictyoTPM[,"Male_3"]/fc_DictyoTPM[,"LengthKB"]
#put results in a new df 
fc_DictyoTPM_1 <- fc_DictyoTPM[, c(1, 9:14)]

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
#these are the TPM of the GEGs

#get rid of length (not needed for heatMap)
fc_DictyoTPM_DEGs <- subset (fc_DictyoTPM_DEGs, select = -Length)

#log2 +1
fc_DictyoTPM_DEGsLog <- log2(fc_DictyoTPM_DEGs+1)
boxplot(fc_DictyoTPM_DEGsLog)
#get matrix for heatMap
dictyo_heatMatrix <- as.matrix(fc_DictyoTPM_DEGsLog)

#heatMap 1 (2)
heatmap.2(x=dictyo_heatMatrix, col = brewer.pal(11,"Spectral"), keysize = 1.5, margins = c(7,3), trace = "none", labRow=FALSE)

clr85 <- colorRampPalette(brewer.pal(9,"OrRd"))(100)  
pheatmap(dictyo_heatMatrix, show_rownames = F, color = clr85, cutree_cols = 2)

#data wrangling for *binary heatmap*
fc_DictyoTPM_DEGsLog$Females <- rowMeans(fc_DictyoTPM_DEGsLog[,c("Female_1", "Female_2","Female_3")])
fc_DictyoTPM_DEGsLog$Males <- rowMeans(fc_DictyoTPM_DEGsLog[,c("Male_1", "Male_2","Male_3")])

#new df Dictyo *******************************************
binary_heatMap_Dictyo <- fc_DictyoTPM_DEGsLog[,c(7,8)] #****************************
colnames(binary_heatMap_Dictyo) <- c("Female","Male") #*******************
binary_heatMap_Dictyo$Condition <- "condition"
binary_heatMap_Dictyo$Species <- "Dictyota"
binary_heatMap_Dictyo$Condition[rownames(binary_heatMap_Dictyo) %in% rownames(Dictyo_FBGs)] <- "FB"
binary_heatMap_Dictyo$Condition[rownames(binary_heatMap_Dictyo) %in% rownames(Dictyo_MBGs)] <- "MB"


#new df    Myrio   ***************************************************************
binary_heatMap_Myrio <- fc_MyrioTPM_DEGsLog[,c(7,8)]#***********************
colnames(binary_heatMap_Myrio) <- c("Female","Male")#****************************
binary_heatMap_Myrio$Condition <- "condition"
binary_heatMap_Myrio$Species <- "Myriotrichia"
binary_heatMap_Myrio$Condition[rownames(binary_heatMap_Myrio) %in% rownames(Myrio_FBGs)] <- "FB"
binary_heatMap_Myrio$Condition[rownames(binary_heatMap_Myrio) %in% rownames(Myrio_MBGs)] <- "MB"

merged_dfs <- bind_rows(binary_heatMap_Dictyo,binary_heatMap_Myrio)

#new df    Desma  ***************************************************************
#new df
binary_heatMap_Desma <- fc_DesmaTPM_DEGsLog[,c(7,8)]
colnames(binary_heatMap_Desma) <- c("Female","Male")#****************************
binary_heatMap_Desma$Condition <- "condition"
binary_heatMap_Desma$Species <- "Desmarestia"
binary_heatMap_Desma$Condition[rownames(binary_heatMap_Desma) %in% rownames(Desma_FBGs)] <- "FB"
binary_heatMap_Desma$Condition[rownames(binary_heatMap_Desma) %in% rownames(Desma_MBGs)] <- "MB"


merged_dfs_2 <- bind_rows(merged_dfs,binary_heatMap_Desma)

#new df   Saccho  ***************************************************************
#new df
binary_heatMap_Saccho <- fc_SacchoTPM_DEGsLog[,c(5,6)]
colnames(binary_heatMap_Saccho) <- c("Female","Male")#****************************
binary_heatMap_Saccho$Condition <- "condition"
binary_heatMap_Saccho$Species <- "Sacchoriza"
binary_heatMap_Saccho$Condition[rownames(binary_heatMap_Saccho) %in% rownames(Saccho_FBGs)] <- "FB"
binary_heatMap_Saccho$Condition[rownames(binary_heatMap_Saccho) %in% rownames(Saccho_MBGs)] <- "MB"

merged_dfs_all <- bind_rows(merged_dfs_2,binary_heatMap_Saccho)


merged_dfs_all_melted <- melt(merged_dfs_all)
colnames(merged_dfs_all_melted)[3] <- "Sex"
colnames(merged_dfs_all_melted)[4] <- "Log2TPM"

ggplot(merged_dfs_all_melted, aes(x=Species,y=Log2TPM, fill=Sex)) +
  geom_boxplot() +
  theme(axis.title.x = element_text(face="bold", size=16),
          axis.title.y = element_text(face="bold", size=16),
        axis.text=element_text(size=14),
        legend.title = element_text(size=16, face="bold"),
        legend.text = element_text(size = 14))



#Plot plot
ggplot(merged_dfs_all, aes() +
         geom_boxplot())

ggplot(merged_dfs_all, aes(x=Species,y=Female, fill=Species)) +
  geom_boxplot() 
