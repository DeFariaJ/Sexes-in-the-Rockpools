#invitations to the party
library(DESeq2)
library(ggplot2)
library(pvclust)
library("gplots")
library("pheatmap")
library("RColorBrewer")
library(gplots)
library(topGO)
# 1
#importing data
fc_Desma <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/output2X_R.tab", header = TRUE)

#genes as rows
rownames(fc_Desma) <- fc_Desma[,1]
fc_Desma[,1] <- NULL

#columns names a bit more handy
names(fc_Desma)[1] <- "Female_1"
names(fc_Desma)[2] <- "Female_2"
names(fc_Desma)[3] <- "Female_3"
names(fc_Desma)[4] <- "Male_1"
names(fc_Desma)[5] <- "Male_2"
names(fc_Desma)[6] <- "Male_3"

#create metadata (data about columns in count data)
metaData_df <-  data.frame(first_column  = c("Female_1", "Female_2", "Female_3", "Male_1", "Male_2", "Male_3"),
                           Sex = c("Female", "Female", "Female", "Male", "Male", "Male")
)
rownames(metaData_df) <- metaData_df$first_column
metaData_df[,1] <- NULL

#################### ******** Deseq2 object ************* ########################

#getting main DESeq2 object
deseq_dataset = DESeqDataSetFromMatrix(countData = fc_Desma, colData = metaData_df, design = ~ Sex )

#prep steps for deseq2
deseq_dataset <- estimateSizeFactors(deseq_dataset) #normalization (geometric mean)

vst_data = varianceStabilizingTransformation(deseq_dataset) #transforms data to make it normally dist
#this is to perform PCA, it assumes the data is normally distributed

#boxplots before and after transformation 
boxplot(counts(deseq_dataset, normalized = TRUE)) #before transformation
boxplot(assay(vst_data)) #after transformation

#PCA
plotPCA(vst_data, intgroup = "Sex")

#vst_data is more normally distributed
d_hca <- assay(vst_data)
d_hca <- t(d_hca)
d_hca <- dist(d_hca)
hca <- hclust(d_hca) 
plot(hca)

#Dendogram with p-values
pValue_Dendo <- pvclust(assay(vst_data), method.dist="cor", method.hclust="average", nboot=1000, parallel=TRUE)         
plot(pValue_Dendo)
print(pValue_Dendo) #to check p-values

################ ********** DGE ********* #####################

##### Under the wood DESeq2 does 3 things ##########

# 1) estimate size factors (normalization) line30 above

# 2) estimate dispersions (key for negative binomial distribution)

# 3) apply statistics (Wald Test)
####################################################

#estimate dispersions (step2) 
deseq_dataset <- estimateDispersions(deseq_dataset) #applies empirical bayesian shrinkage 
# bayesian shrinkage improves rubusteness of observations in experiments with not much replicates
#visualize results
plotDispEsts(deseq_dataset)
# black points are raw estimates of dispersion
#estimated relationship between the mean counts and the average dispersion
# blue points are the final estimate of dispersion
# black dots inside blue circles are genes considered outliers in terms of dispersion

# 3) apply statistics
deseq_dataset <- nbinomWaldTest(deseq_dataset) #stats being applied to ~20000 genes takes a bit

#check results
result_table = results(deseq_dataset)
summary(result_table)
#putting results table into a nice data.frame
results_df <- as.data.frame(result_table)

## Filter rows (genes) with all columns == NA (might be outliers or did not pass filter tests)
complete.cases(results_df) #complete.cases returns TRUE for columns data have all data (no NAs)
filter_df1 <- results_df[complete.cases(results_df),]

# adj p-value <0.05
filter_df1$padj <0.05
filter_df2 <- filter_df1[filter_df1$padj <0.05,] 

#p-value <0.05 is not enough, we also want to have a large effect size
#log2FoldChange >= 1 <= -1, we will filter it on magnitude and not sign, regardless of direction
abs(filter_df2$log2FoldChange) > 1 # 2 fold regulation in either direction 
filter_df3 <- filter_df2[abs(filter_df2$log2FoldChange) > 1,]

#MA plot
plotMA(result_table, alpha = 0.05) #takes as input original table produced by DESeq2 (not data.frame) #CAN CHANGE TO 0.05
# the aim here is to produce a comparison of the mean normalized counts (baseMean) and fold change

######## VOLCANO PLOT ############
#for all the genes we have analysied, we need the fold change and the adj p-value

#creating a new column on filter_df1 (up, neutra, down)
filter_df1$DiffExpressed <-"Unbiased"
filter_df1$DiffExpressed[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >1] <- "FC >2 "
filter_df1$DiffExpressed[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >2] <- "FC >4"
filter_df1$DiffExpressed[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >4] <- "FC >16"  

#Plot plot
ggplot(filter_df1, aes(x=log2FoldChange, y= -log10(padj))) +
  geom_point(aes(colour = DiffExpressed)) +
  geom_vline(xintercept = 1, linetype = 3) + 
  geom_vline(xintercept = -1, linetype = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 3) 



################################################################################################################
#    ********************************************************** TPMs *************************************
################################################################################################################

#importing data
fc_DesmaTPM <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/output2X_TPM.tab", header = TRUE)
#same as output from feature counts used for DGE analysis plus *gene length*


#genes as rows
rownames(fc_DesmaTPM) <- fc_DesmaTPM[,1]
fc_DesmaTPM[,1] <- NULL

#columns names a bit more handy
names(fc_DesmaTPM)[2] <- "Female_1"
names(fc_DesmaTPM)[3] <- "Female_2"
names(fc_DesmaTPM)[4] <- "Female_3"
names(fc_DesmaTPM)[5] <- "Male_1"
names(fc_DesmaTPM)[6] <- "Male_2"
names(fc_DesmaTPM)[7] <- "Male_3"

fc_DesmaTPM$LengthKB <- fc_DesmaTPM$Length/1000

# Step 1 - Normalize for gene length #########!!!!!!!

#fc_DictyoTPM$F1<- fc_DictyoTPM[,"Female_1"]/fc_DictyoTPM[,"Length"]
#Divide each gene count by gene length (I am adding results as new variables in the same df(not very sofisticated))
fc_DesmaTPM$F1<- fc_DesmaTPM[,"Female_1"]/fc_DesmaTPM[,"LengthKB"]
fc_DesmaTPM$F2<- fc_DesmaTPM[,"Female_2"]/fc_DesmaTPM[,"LengthKB"]
fc_DesmaTPM$F3<- fc_DesmaTPM[,"Female_3"]/fc_DesmaTPM[,"LengthKB"]
fc_DesmaTPM$M1<- fc_DesmaTPM[,"Male_1"]/fc_DesmaTPM[,"LengthKB"]
fc_DesmaTPM$M2<- fc_DesmaTPM[,"Male_2"]/fc_DesmaTPM[,"LengthKB"]
fc_DesmaTPM$M3<- fc_DesmaTPM[,"Male_3"]/fc_DesmaTPM[,"LengthKB"]

#put results in a new df 
fc_DesmaTPM_1 <- fc_DesmaTPM[, c(1, 9:14)]

#Step 2 - Normalize for sequence depth
# 2.1 - Add up read counts for each sample and divide by 1 million
F1_sc <- sum(fc_DesmaTPM_1$F1/1000000)
F2_sc <- sum(fc_DesmaTPM_1$F2/1000000)
F3_sc <- sum(fc_DesmaTPM_1$F3/1000000)
M1_sc <- sum(fc_DesmaTPM_1$M1/1000000)
M2_sc <- sum(fc_DesmaTPM_1$M2/1000000)
M3_sc <- sum(fc_DesmaTPM_1$M3/1000000)


# 2.2 - normalized gene counts / new scalling factor
fc_DesmaTPM_1$Female_1 <- fc_DesmaTPM_1[,"F1"]/F1_sc
fc_DesmaTPM_1$Female_2 <- fc_DesmaTPM_1[,"F2"]/F2_sc
fc_DesmaTPM_1$Female_3 <- fc_DesmaTPM_1[,"F3"]/F3_sc
fc_DesmaTPM_1$Male_1 <- fc_DesmaTPM_1[,"M1"]/M1_sc
fc_DesmaTPM_1$Male_2 <- fc_DesmaTPM_1[,"M2"]/M2_sc
fc_DesmaTPM_1$Male_3 <- fc_DesmaTPM_1[,"M2"]/M3_sc
#put final results in a new df
fc_DesmaTPM_2 <- fc_DesmaTPM_1[, c(1, 8:13)] #******************** TPM table ****************
###########################################
### save fc_MyrioTPM_2 like it is here
###########################################

#this creates a new df with genes (rows) which must be present in two dif DFs 
#fc_DictyoTPM2 and filter_DF3
(rownames(fc_DesmaTPM_2) %in% rownames(filter_df3)) #logical matrix
fc_DesmaTPM_DEGs <- fc_DesmaTPM_2[(rownames(fc_DesmaTPM_2) %in% rownames(filter_df3)),] #apply logical filtering

#get rid of length (not needed for heatMap)
fc_DesmaTPM_DEGs <- subset (fc_DesmaTPM_DEGs, select = -Length)

#log2 +1
fc_DesmaTPM_DEGsLog <- log2(fc_DesmaTPM_DEGs+1)
#get matrix for heatMap
Desma_heatMatrix <- as.matrix(fc_DesmaTPM_DEGsLog)

#heatMap 1 (2)
heatmap.2(x=Desma_heatMatrix, col = brewer.pal(11,"Spectral"), keysize = 1.5, margins = c(7,3), trace = "none",labRow=FALSE)

#test test test
#pheatmap(binary_heatMapMatrix2, show_rownames = F)
clr85 <- colorRampPalette(brewer.pal(9,"OrRd"))(100)  
pheatmap(Desma_heatMatrix, show_rownames = F, color = clr85, cutree_cols = 2)

#data wrangling for *binary heatmap*
fc_DesmaTPM_DEGsLog$Females <- rowMeans(fc_DesmaTPM_DEGsLog[,c("Female_1", "Female_2","Female_3")])
fc_DesmaTPM_DEGsLog$Males <- rowMeans(fc_DesmaTPM_DEGsLog[,c("Male_1", "Male_2","Male_3")])
#new df
binary_heatMap <- fc_DesmaTPM_DEGsLog[,c(7,8)]
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

#pheatmap(binary_heatMapMatrix2, show_rownames = F)
clr <- colorRampPalette(brewer.pal(9,"OrRd"))(100)  
pheatmap(binary_heatMapMatrix2, show_rownames = F, color = clr, cutree_cols = 2)


#****

#################################################################################################
#################################################################################################
#trying to further analysis SBGE number in males vs females
#talk to Aga, not sure if this is correct but here we go

#create a female_up df which is a subset of binary heat_map females > males
Female_Up <- subset(binary_heatMap2[binary_heatMap2$Females>binary_heatMap2$Males,])
#same for males
Male_Up <- subset(binary_heatMap2[binary_heatMap2$Females<binary_heatMap2$Males,])

#######################################FEMALE BIAS GENES######################################
####################
#FC >2
####################

#extract row names (genes) from filter_DF1 with FC > 2
fold_2_genes <- row.names(filter_df1[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >1,])
#logical matrix, genes which are both present in Female_up and filter_DF1
logical_vector <-(rownames(Female_Up) %in% (fold_2_genes))
#count number of TRUES which (I think, hope, presume) = number of genes FBG FC >2
sum(logical_vector, na.rm = TRUE) # best way to count TRUE values

######################
#FC >4
#####################

#extract row names (genes) from filter_DF1 with FC > 4
fold_4_genes <- row.names(filter_df1[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >2,])
#logical matrix, genes which are both present in Female_up and filter_DF1
logical_vector <-(rownames(Female_Up) %in% (fold_4_genes))
#count number of TRUES which (I think, hope, presume) = number of genes FBG FC >4
sum(logical_vector, na.rm = TRUE) # best way to count TRUE values

######################
#FC >6
#####################

#extract row names (genes) from filter_DF1 with FC > 4
fold_6_genes <- row.names(filter_df1[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >4,])
#logical matrix, genes which are both present in Female_up and filter_DF1
logical_vector <-(rownames(Female_Up) %in% (fold_6_genes))
#count number of TRUES which (I think, hope, presume) = number of genes FBG FC >4
sum(logical_vector, na.rm = TRUE) # best way to count TRUE values

#######################################MALE BIAS GENES######################################

####################
#FC >2
####################

#logical matrix, genes which are both present in Female_up and filter_DF1
logical_vector <-(rownames(Male_Up) %in% (fold_2_genes))
#count number of TRUES which (I think, hope, presume) = number of genes FBG FC >2
sum(logical_vector, na.rm = TRUE) # best way to count TRUE values

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
sanity_check <- row.names(filter_df1[filter_df1$padj>.05 & abs(filter_df1$log2FoldChange) <=1,])


#3)
###########################################################################################################
# FUNCTIONAL ANALYSIS
###########################################################################################################

#importing data
ego_prep <- read.table(file = "/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/xxx_Desma.annotations", sep="\t")


##Females
#get vector of DEGs IDs
gene_vector_Females <- factor(rownames(Female_Up))
universal_genes <- rownames(filter_df1)
#logical universal genes in female vector
geneList <- factor(as.integer(universal_genes%in%gene_vector_Females))
names(geneList) <- universal_genes

geneID2GO <- readMappings(file = "/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/xxx_Desma.annotations")

#object of class topGOdata
#This an object so we can extract information from it (use its methods)
topObject_DesmaF <- new("topGOdata",
                        description = "Functional analysis Desmarestia Females", ontology = "BP",
                        allGenes = geneList, geneSel = gene_vector_Females,
                        annot = annFUN.gene2GO,
                        gene2GO = geneID2GO
)

#stats object
teststat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

#*********Main function to run enichment analysys**********#
resultFisher <- getSigGroups(topObject_DesmaF, test.stat = teststat)

#restuls in df
allRes <- GenTable(topObject_DesmaF, classic = resultFisher, topNodes = 50)
allResTop <- GenTable(topObject_DesmaF, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = 10)

#basic info
geneData(resultFisher)

########Males##########

#get vector of DEGs IDs
gene_vector_Males <- factor(rownames(Male_Up))

#logical universal genes in male vector
geneList <- factor(as.integer(universal_genes%in%gene_vector_Males))
universal_genes <- rownames(filter_df1)
names(geneList) <- universal_genes

#object of class topGOdata
#This an object so we can extract information from it (use its methods)
topObject_DesmaM <- new("topGOdata",
                        description = "Functional analysis Desmarestia Males", ontology = "BP",
                        allGenes = geneList, geneSel = gene_vector_Males,
                        annot = annFUN.gene2GO,
                        gene2GO = geneID2GO
)


#stats object
teststat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

#*********Main function to run enichment analysys**********#
resultFisherM <- getSigGroups(topObject_DesmaM, test.stat = teststat)

#restuls in df Males
allResM <- GenTable(topObject_DesmaM, classic = resultFisherM, topNodes = 50)
allResTopM <- GenTable(topObject_DesmaM, classicFisher = resultFisherM, ranksOf = "classicFisher", topNodes = 10)
