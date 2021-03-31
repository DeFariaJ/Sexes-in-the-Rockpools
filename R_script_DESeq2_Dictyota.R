#packages
library(DESeq2)

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
plot(hclust(vst_clust))

##### Under the wood DESeq2 does 3 things ##########

# 1) estimate size factors (normalization) line30 above

# 2) estimate dispersions (key for negative binomial distribution)

# 3) apply statistics (Wald Test)
####################################################


#estimate dispersions (step2) 
deseq_dataset <- estimateDispersions(deseq_dataset) #applies baysian shrinkage etc
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

## this would throw away rows (genes) with baseline NA or other NA
#complete.cases(results_df) #complete.cases returns TRUE for columns data have all data (no NAs)
#filter_df <- results_df[complete.cases(results_df),]


# so we go directly to p-value and foldChange filtering
# adj p-value <0.05
results_df$padj <0.05
filter_df1 <- results_df[results_df$padj <0.05,]

#p-value <0.05 is not enough, we also want to have a large effect size
#log2FoldChange >= 1 <= -1, we will filter it on magnitude and not sign, regardless of direction
abs(filter_df1$log2FoldChange) > 1
filter_df2 <- filter_df1[abs(filter_df1$log2FoldChange) > 1,]
