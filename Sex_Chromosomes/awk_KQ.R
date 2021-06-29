#awk calculations to get candidate scaffolds KQ method

######## MALES ##########

#importing data
KQ_table <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/genomeX.size.100kb.sorted.cov", header = FALSE)

#adding new awk variable which is (V4/V6)*1000
KQ_table$awk_value <- KQ_table[,"V4"]/KQ_table[,"V6"]*1000

#filtering
KQ_table$awk_value >= 250
KQ_sexLink <- KQ_table[KQ_table$awk_value >= 250,]


######## FEMALES ##########

#importing data
KQ_table_F <- read.table("/Users/goncaloleiria/Desktop/Sex_biased_GeneExpression/Sexes_in_the_Rockpools/genome_F.size.100kb.sorted.cov", header = FALSE)

#adding new awk variable which is (V4/V6)*1000
KQ_table_F$awk_value <- KQ_table_F[,"V4"]/KQ_table_F[,"V6"]*1000

#filtering
KQ_table_F$awk_value >= 250
KQ_sexLink_F <- KQ_table_F[KQ_table_F$awk_value >= 250,]
