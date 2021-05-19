#libraries
library(ggplot2)

#imports
F_to_M<- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Agamosdepth_FEMALEtoMALE.regions.bed")
contig_size <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/contig_Size_Male")
male_gff <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Myriotrichia-clavaeformis_MALE.gff")
scaffolds_genes_clean <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/scaffold_genes_Male_clean.gff")

#################################### FEMALE to MALE ###################################
#Ggetting rid of outliers
F_to_M_2 <- subset(F_to_M, F_to_M$V4 < 139) 
#find intervals with zero coverage
zero_Coverage <- subset(F_to_M, F_to_M$V4 < 1) # 6855 intervals with 0 coverage (16%)

#isolating singles and put them into a dataframe 
singles <- names(which(table(zero_Coverage$V1) == 1)) # get frequency table
unique_zero_Coverage <- zero_Coverage[zero_Coverage$V1 %in% singles, ] # use frequency table to subset

#column names 
names(unique_zero_Coverage)[1] <- "Interval"
names(unique_zero_Coverage)[2] <- "Start"
names(unique_zero_Coverage)[3] <- "End"
names(unique_zero_Coverage)[4] <- "Coverage_depth"

#adding interval size columns (end - start)
unique_zero_Coverage$Size_interval <- (unique_zero_Coverage$End - unique_zero_Coverage$Start)

#genes as rows
rownames(unique_zero_Coverage) <- unique_zero_Coverage[,1]
unique_zero_Coverage[,1] <- NULL
#genes as rows
rownames(contig_size) <- contig_size[,1]
contig_size[,1] <- NULL

#this creates a new df with genes (rows) which must be present in two dif DFs 
(rownames(unique_zero_Coverage) %in% rownames(contig_size)) #logical matrix
fused_df<- unique_zero_Coverage[(rownames(unique_zero_Coverage) %in% rownames(contig_size)),] #apply logical filtering

ppp <- subset(contig_size, rownames(contig_size) %in% rownames(fused_df)) 

#merging total length to fused dataframe
fused_df$total_length <- ppp$V2
#getting fractions
fused_df$fraction <- fused_df$Size_interval/fused_df$total_length

# here I am getting rid of intervals which account for < 50% of the full scaffold
fused_df_clean <- fused_df[(fused_df$fraction > 0.50),]

#lets put the 46 isolated 0 intervals into a vector
zero_islands <- row.names(fused_df[(fused_df$fraction < 0.50),])

#remove duplicates from zero coverage (we just need the scaffolds)
scaffolds <- zero_Coverage[!duplicated(zero_Coverage$V1), ]
#scaffolds as rows
rownames(scaffolds) <- scaffolds[,1]
scaffolds[,1] <- NULL

#powerful trick
`%notin%` <- Negate(`%in%`)

# this df has all scaffolds == 0
# duplicate intervals were deleted (keeping one of them)
# large unique intervals kept
scaffolds_unique <- subset(scaffolds, row.names(scaffolds) %notin% (zero_islands))

#sanity check
unique_vector <- row.names(scaffolds_unique)
xxx <- subset(contig_size, row.names(contig_size) %in% (unique_vector))
sum(xxx$V2)
(7700389/157455168)*100 #looks like, in terms of base pairs, the candidate saffolds add up to 4.89% of the genome

############################################################################
# scaffolds vf gff file
###########################3 trouble 

# Remove duplicates (based on Sepal.Width columns ????)
male_gff_unique <- male_gff[!duplicated(male_gff$V1), ]

#scaffolds as rows
rownames(male_gff_unique) <- male_gff_unique[,1]
male_gff_unique[,1] <- NULL

#number of dif entries in column 3 (type)
length(unique(male_gff[["V3"]]))

# get scaffolds which may be part of genes
scaffolds_genes <- subset(male_gff_unique, rownames(male_gff_unique) %in% rownames(fused_df)) 

# re-organize columns and rows 
scaffolds_genes$V1 <- rownames(scaffolds_genes)
row.names(scaffolds_genes) <- NULL
scaffolds_genes <- scaffolds_genes[c("V1","V2","V3","V4","V5","V6","V7","V8","V9")]

#crossing info  with TPM list
g <- subset(scaffolds_genes, row.names(scaffolds_genes) %in% (row.names(fc_MyrioTPM_DEGs)))

vectorXXX <- as.vector(row.names(fc_MyrioTPM_DEGs))
grep(vectorXXX, scaffolds_genes$V9)

####### clean scaffold genes name (Aga) ###########
#genes as rows
rownames(scaffolds_genes_clean) <- scaffolds_genes_clean[,9]
scaffolds_genes_clean[,9] <- NULL

#getting there
XXX <- subset(fc_MyrioTPM_DEGs, row.names(fc_MyrioTPM_DEGs) %in% (row.names(scaffolds_genes_clean)))


# compare again against full list of TPMs (done)
XXX_2 <- subset(fc_MyrioTPM_2, row.names(fc_MyrioTPM_2) %in% (row.names(scaffolds_genes_clean)))

# from there create table
df <- as.data.frame(scaffolds_genes$V1)
names(df)[1] <- "contig"
df$gene <- row.names(scaffolds_genes_clean)

#getting means
df$TPM_Female <- rowMeans(XXX_2[,c("Female_1", "Female_2","Female_3")])
df$TPM_Male <- rowMeans(XXX_2[,c("Male_1", "Male_2","Male_3")])

DEG_vector <- row.names(XXX_2) %in% (row.names(XXX))
df$DEG <- DEG_vector
