################################################################################################################
#trying to further analysis SBGE number in males vs females
#not sure if this is correct but here we go


#create a female_up df which is a subset of binary heat_map females > males
Female_Up <- subset(binary_heatMap2[binary_heatMap2$Females>binary_heatMap2$Males,])
#same for males
Male_Up <- subset(binary_heatMap2[binary_heatMap2$Females<binary_heatMap2$Males,])

#######################################FEMALE BIAS GENES######################################
####################
#FC >2
####################

#extract row names (genes) from filter_DF1 with FC > 2
fold_2_genes <- row.names(filter_df1[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >2,])
#logical matrix, genes which are both present in Female_up and filter_DF1
logical_vector <-(rownames(Female_Up) %in% (fold_2_genes))
#count number of TRUES which (I think, hope, presume) = number of genes FBG FC >2
sum(logical_vector, na.rm = TRUE) # best way to count TRUE values

######################
#FC >4
#####################

#extract row names (genes) from filter_DF1 with FC > 4
fold_4_genes <- row.names(filter_df1[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >4,])
#logical matrix, genes which are both present in Female_up and filter_DF1
logical_vector <-(rownames(Female_Up) %in% (fold_4_genes))
#count number of TRUES which (I think, hope, presume) = number of genes FBG FC >4
sum(logical_vector, na.rm = TRUE) # best way to count TRUE values

######################
#FC >6
#####################

#extract row names (genes) from filter_DF1 with FC > 4
fold_6_genes <- row.names(filter_df1[filter_df1$padj<.05 & abs(filter_df1$log2FoldChange) >6,])
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
sanity_check <- row.names(filter_df3[filter_df3$padj<.05 & abs(filter_df3$log2FoldChange) <=2,])

#####this is fucking correct!!!! I think....
