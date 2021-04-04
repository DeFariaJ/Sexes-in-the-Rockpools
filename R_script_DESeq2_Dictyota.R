#packages
library(DESeq2)
library(ggplot2)
library(pvclust)

#importing data
fc_Dictyo <- read.table("DictyoOutput.Rmatrix.tab", header = TRUE)

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
deseq_dataset = DESeqDataSetFromMatrix(countData = fc_Dictyo, colData = metaData_df, design = ~ Sex )

#prep steps for deseq2
deseq_dataset = estimateSizeFactors(deseq_dataset)
vst_data = varianceStabilizingTransformation(deseq_dataset) #transforms data to make it normally dist

#boxplots before and after transformation 
boxplot(counts(deseq_dataset, normalized = TRUE)) #before transformation
boxplot(assay(vst_data)) #after transformation

#PCA
plotPCA(vst_data, intgroup = "Sex")

#HCA
vst_clust = dist(t(assay(vst_data)))
plot(pvclust(vst_clust))


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
filter_df1$DiffExpressed <-"Neutral"
filter_df1$DiffExpressed[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >1] <- "DEGs"
filter_df1$DiffExpressed[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >2] <- "DEGs+"
filter_df1$DiffExpressed[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >4] <- "DEGs++"                                                
                                                 

#Plot plot
ggplot(filter_df1, aes(x=log2FoldChange, y= -log10(padj))) +
         geom_point(aes(colour = DiffExpressed)) +
         geom_vline(xintercept = 1, linetype = 3) + 
         geom_vline(xintercept = -1, linetype = 3) +
         geom_hline(yintercept = -log10(0.05), linetype = 3) +
         ylim(-5, 125)
