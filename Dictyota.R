#This script takes the output from featureCounts out of the command line and:
# 1) brings it to R
# 2) Creates handy data frames
# 3) Does the data wrangling needed to create an DESeq2 object
# 4) Performs exploratory data analysis (PCA and HCA)
# 5) Performs DGE analysis 
# 6) Puts the results in nice handy tables
# 7) Performs visualisation of the results (MA and Volcano plots)
# 8) Performs TPM calculations
#9) Performs functional analysis

#packages
library(DESeq2)
library(ggplot2)
library(pvclust)
library(clusterProfiler)
library(tidyr)
library(topGO)
#importing data
fc_Dictyo <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/DictyoOutput.Rmatrix.tab", header = TRUE)

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

##### Under the wood DESeq2 does 3 things ##########

# 1) estimate size factors (normalization) line30 above

# 2) estimate dispersions (key for negative binomial distribution)

# 3) apply statistics (Wald Test)
####################################################


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
         ylim(-5, 150)

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

clr85 <- colorRampPalette(brewer.pal(9,"OrRd"))(100)  
pheatmap(dictyo_heatMatrix, show_rownames = F, color = clr85, cutree_cols = 2)

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
clr <-colorRampPalette(brewer.pal(9, "OrRd"))(100)  
pheatmap(binary_heatMapMatrix2, show_rownames = F, color = clr, cutree_cols = 2)



################################################################################################################
################################################################################################################

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

######################
#FC >10
######################
#extract row names (genes) from filter_DF1 with FC > 4
fold_10_genes <- row.names(filter_df1[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >3.322,])
#logical matrix, genes which are both present in Female_up and filter_DF1
logical_vector <-(rownames(Female_Up) %in% (fold_10_genes))
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

######################
#FC >10
######################
#logical matrix, genes which are both present in Female_up and filter_DF1
logical_vector <-(rownames(Male_Up) %in% (fold_10_genes))
#count number of TRUES which (I think, hope, presume) = number of genes FBG FC >4
sum(logical_vector, na.rm = TRUE) # best way to count TRUE values


#Sanity check
sanity_check <- row.names(filter_df1[filter_df1$padj>.05 & abs(filter_df1$log2FoldChange) <=1,])
#### correct!!!! I think...

###########################################################################################################
# FUNCTIONAL ANALYSIS
###########################################################################################################

#importing data
ego_prep <- read.table(file = "/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/new.annotations", sep="\t")
ego_tidy <- tidyr::separate_rows(data=ego_prep,"V2",sep = ",")

##Females
#get vector of DEGs IDs
gene_vector_Females <- factor(rownames(Female_Up))
universal_genes <- rownames(filter_df1)
#logical universal genes in female vector
geneList <- factor(as.integer(universal_genes%in%gene_vector_Females))
names(geneList) <- universal_genes

geneID2GO <- readMappings(file = "/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/new.annotations")

#object of class topGOdata
#This an object so we can extract information from it (use its methods)
topObject_DictyoF <- new("topGOdata",
                description = "Functional analysis Dictyota Females", ontology = "BP",
                allGenes = geneList, geneSel = gene_vector_Females,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO
              )
#stats object
teststat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

#*********Main function to run enichment analysys**********#
resultFisher <- getSigGroups(topObject_DictyoF, test.stat = teststat)

#restuls in df
allRes <- GenTable(topObject_DictyoF, classic = resultFisher, topNodes = 50)
allResTop <- GenTable(topObject_DictyoF, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = 10)

#basic info
geneData(resultFisher)

########Males##########

#get vector of DEGs IDs
gene_vector_Males <- factor(rownames(Male_Up))

#logical universal genes in female vector
geneList <- factor(as.integer(universal_genes%in%gene_vector_Males))
universal_genes <- rownames(filter_df1)
names(geneList) <- universal_genes

#object of class topGOdata
#This an object so we can extract information from it (use its methods)
topObject_DictyoM <- new("topGOdata",
                         description = "Functional analysis Dictyota Males", ontology = "BP",
                         allGenes = geneList, geneSel = gene_vector_Males,
                         annot = annFUN.gene2GO,
                         gene2GO = geneID2GO
                         )

#stats object
teststat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

#*********Main function to run enichment analysys**********#
resultFisherM <- getSigGroups(topObject_DictyoM, test.stat = teststat)

#restuls in df Males
allResM <- GenTable(topObject_DictyoM, classic = resultFisherM, topNodes = 50)
allResTopM <- GenTable(topObject_DictyoM, classicFisher = resultFisherM, ranksOf = "classicFisher", topNodes = 10)
