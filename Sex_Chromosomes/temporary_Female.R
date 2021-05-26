#imports
M_to_F<- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Agamosdepth_MALEtoFEMALE.regions.bed")
contig_size_F <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/contig_Size_Female")
female_gff <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/Myriotrichia-clavaeformis_FEMALE.gff")
scaffolds_genes_clean_F <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/scaffold_genes_Female_clean.gff")

#################################### MALE to FEMALE ###################################
#Ggetting rid of outliers
M_to_F_2 <- subset(M_to_F, M_to_F$V4 < 165) 
#find intervals with zero coverage
zero_Coverage_F <- subset(M_to_F, M_to_F$V4 < 1) # 1073 intervals with 0 coverage (XX%)

#isolating singles and put them into a dataframe 
singles_F <- names(which(table(zero_Coverage_F$V1) == 1)) # get frequency table
unique_zero_Coverage_F <- zero_Coverage_F[zero_Coverage_F$V1 %in% singles_F, ] # use frequency table to subset

#column names 
names(unique_zero_Coverage_F)[1] <- "Interval"
names(unique_zero_Coverage_F)[2] <- "Start"
names(unique_zero_Coverage_F)[3] <- "End"
names(unique_zero_Coverage_F)[4] <- "Coverage_depth"

#adding interval size columns (end - start)
unique_zero_Coverage_F$Size_interval <- (unique_zero_Coverage_F$End - unique_zero_Coverage_F$Start)

#genes as rows
rownames(unique_zero_Coverage_F) <- unique_zero_Coverage_F[,1]
unique_zero_Coverage_F[,1] <- NULL
#genes as rows
rownames(contig_size_F) <- contig_size_F[,1]
contig_size_F[,1] <- NULL

#this creates a new df with genes (rows) which must be present in two dif DFs 
(rownames(unique_zero_Coverage_F) %in% rownames(contig_size_F)) #logical matrix
fused_df_F<- unique_zero_Coverage_F[(rownames(unique_zero_Coverage_F) %in% rownames(contig_size_F)),] #apply logical filtering

ppp_F <- subset(contig_size_F, rownames(contig_size_F) %in% rownames(fused_df_F))

#merging total length to fused dataframe
fused_df_F$total_length <- ppp_F$V2
#getting fractions
fused_df_F$fraction <- fused_df_F$Size_interval/fused_df_F$total_length

# here I am getting rid of intervals which account for < 50% of the full scaffold
fused_df_clean_F <- fused_df_F[(fused_df_F$fraction > 0.50),]

#lets put the 36 isolated 0 intervals into a vector
zero_islands_F <- row.names(fused_df_F[(fused_df_F$fraction < 0.50),])

#remove duplicates from zero coverage (we just need the scaffolds)
scaffolds_F <- zero_Coverage_F[!duplicated(zero_Coverage_F$V1), ]
#scaffolds as rows
rownames(scaffolds_F) <- scaffolds_F[,1]
scaffolds_F[,1] <- NULL

#powerful trick
`%notin%` <- Negate(`%in%`)

# this df has all scaffolds == 0
# duplicate intervals were deleted (keeping one of them)
# large unique intervals kept
scaffolds_unique_F <- subset(scaffolds_F, row.names(scaffolds_F) %notin% (zero_islands_F))

#sanity check
unique_vector_F <- row.names(scaffolds_unique_F)
xxx_F <- subset(contig_size_F, row.names(contig_size_F) %in% (unique_vector_F))
sum(xxx_F$V2) #3981152
(3981152/157455168)*100 #looks like, in terms of base pairs, the candidate saffolds add up to 2.53% of the genome

############################################################################
# scaffolds vf gff file FEMALES
########################################## 3 trouble

## Remove duplicates (based on Sepal.Width columns ????)
female_gff_unique <- female_gff[!duplicated(female_gff$V1), ]

#scaffolds as rows
rownames(female_gff_unique) <- female_gff_unique[,1]
female_gff_unique[,1] <- NULL

#number of dif entries in column 3 (type) #just checking
length(unique(female_gff[["V3"]]))

# get scaffolds which may be part of genes
scaffolds_genes_F <- subset(female_gff_unique, rownames(female_gff_unique) %in% rownames(fused_df_F)) 

# re-organize columns and rows 
scaffolds_genes_F$V1 <- rownames(scaffolds_genes_F)
row.names(scaffolds_genes_F) <- NULL
scaffolds_genes_F <- scaffolds_genes_F[c("V1","V2","V3","V4","V5","V6","V7","V8","V9")]


#crossing info  with TPM list
g_f <- subset(scaffolds_genes_F, row.names(scaffolds_genes_F) %in% (row.names(fc_MyrioTPM_DEGs)))

vectorXXX_f <- as.vector(row.names(fc_MyrioTPM_DEGs))
grep(vectorXXX_f, scaffolds_genes_F$V9)


####### clean scaffold genes name (Aga) ###########
#genes as rows
rownames(scaffolds_genes_clean_F) <- scaffolds_genes_clean_F[,9]
scaffolds_genes_clean_F[,9] <- NULL

#getting there
XXX_F <- subset(fc_MyrioTPM_DEGs_F, row.names(fc_MyrioTPM_DEGs_F) %in% (row.names(scaffolds_genes_clean_F)))

# compare again against full list of TPMs (done)
XXX_2_F <- subset(fc_MyrioTPM_2_F, row.names(fc_MyrioTPM_2_F) %in% (row.names(scaffolds_genes_clean_F)))

# from there create table
df_F <- as.data.frame(scaffolds_genes_F$V1)
names(df_F)[1] <- "contig"
df_F$gene <- row.names(scaffolds_genes_clean_F)

#getting means
df_F$TPM_Female_F <- rowMeans(XXX_2_F[,c("Female_1", "Female_2","Female_3")])
df_F$TPM_Male_F <- rowMeans(XXX_2_F[,c("Male_1", "Male_2","Male_3")])

DEG_vector_F <- row.names(XXX_2_F) %in% (row.names(XXX_F))
df_F$DEG <- DEG_vector_F
