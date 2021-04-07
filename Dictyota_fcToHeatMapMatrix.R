#importing data
fc_DictyoTPM <- read.table("Dictyo_featureMatrixTPM.txt", header = TRUE)

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
fc_DictyoTPM$F3<- fc_DictyoTPM[,"Female_2"]/fc_DictyoTPM[,"Length"]
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
